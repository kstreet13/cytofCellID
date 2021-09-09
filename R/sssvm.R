

# source('R/utils.R')
# fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
# sce <- prepData(fname)



runSSSVM <- function(sce, datatype = c('Bagwell','Crompton')){
    datatype <- match.arg(datatype)
    out <- rep('cell', ncol(sce))
    
    # setup
    if(datatype == 'Crompton'){
        bead_channels <- c(19, 21, 30, 32, 54)
        tech <- t(assay(sce,'exprs')[bead_channels,])
        tech_bead_channels <- 1:ncol(tech)
        tech <- cbind(tech, t(assay(sce,'exprs')[rownames(sce) %in% c('DNA1','DNA2','Viability','VeriCells'),])) # 'VeriCells' ??
        tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))
    }else{
        bead_channels <- grep('Bead', rownames(sce))
        tech <- t(assay(sce,'exprs')[bead_channels,])
        tech_bead_channels <- 1:ncol(tech)
        tech <- cbind(tech, t(assay(sce,'exprs')[c('DNA1','DNA2','Live_Dead'),]))
        tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))
    }
    
    # we need fast, approximate KNN
    s.tech <- scale(tech)
    require(BiocNeighbors)
    k.dex <- buildKmknn(s.tech) 
    knn <- findKNN(s.tech, 10, BNINDEX=k.dex)
    
    # STEP 1: MARK GDP ZEROS
    out[int_colData(sce)$Center == 0 |
            int_colData(sce)$Offset == 0 |
            int_colData(sce)$Width == 0 |
            int_colData(sce)$Residual == 0] <- 'GDPzero'
    
    
    
    # approximate the label
    # pick events near the boundary
    # SVM on those events (+ some random ones)
    # optional: re-pick events, re-fit SVM
    
    #####################
    ### STEP 2: BEADS ###
    #####################
    X <- s.tech[which(out == 'cell'), ]
    
    # approximate the label
    {
        isBeadMat <- sapply(tech_bead_channels, function(col){
            x <- X[,col]
            g <- find_groups(x[x>0])
            mns <- by(x[x>0], g, mean)
            beadgp <- which.max(mns)
            out <- rep(FALSE, nrow(X))
            if(length(levels(g)) > 1){
                out[x>0][g == beadgp] <- TRUE
            }
            return(out)
        })
        init <- rowMeans(isBeadMat) > 0.5
    } # init
    
    # pick "true" events (all "neighbors of neighbors" share the same label)
    {
        mismatch <- matrix(init[knn$index[which(out=='cell'),]], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        neighbormismatch <- matrix(rowMeans(mismatch)[knn$index[which(out=='cell'),]], ncol = 10)
        neighbormismatch[is.na(neighbormismatch)] <- 0
        #ind.label <- which(rowSums(neighbormismatch) == 0)
        ind.unlab <- which(rowSums(neighbormismatch) != 0)
    } # init.unlab
    
    # semi-supervised SVM to classify all events
    {
        library(ssc)
        library(e1071)
        ytrain <- as.integer(init)
        ytrain[ind.unlab] <- NA
        ssvm <- selfTraining(x=X, y = ytrain, x.inst = TRUE, learner = svm,
                             learner.pars =list(kernel ="linear", scale = FALSE, probability = TRUE),
                             pred = function(m, x){
                                 attr(predict(m, x, probability = TRUE), "probabilities")
                             })
        pred <- predict(ssvm, X)
    } # pred
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(as.numeric(as.character(pred)))))] <- 'bead'
    
    
    
    ########################
    ### STEP 3: DOUBLETS ###
    ########################
    X <- s.tech[which(out == 'cell'), ]
    
    # approximate the label
    {
        doubletscore <- X[,'Residual'] + 
            X[,'Event_length'] +
            X[,'DNA2'] +
            X[,'DNA1'] -
            X[,'Offset'] -
            .5 * X[,'Width']
        g <- find_groups(doubletscore)
        doubletclus <- max(as.numeric(levels(g)))
        init <- g == doubletclus
    } # init
    
    # pick "true" events (all "neighbors of neighbors" share the same label)
    {
        mismatch <- matrix(init[knn$index[which(out=='cell'),]], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        neighbormismatch <- matrix(rowMeans(mismatch)[knn$index[which(out=='cell'),]], ncol = 10)
        neighbormismatch[is.na(neighbormismatch)] <- 0
        #ind.label <- which(rowSums(neighbormismatch) == 0)
        ind.unlab <- which(rowSums(neighbormismatch) != 0)
    } # init.unlab
    
    # semi-supervised SVM to classify all events
    {
        library(ssc)
        library(e1071)
        ytrain <- as.integer(init)
        ytrain[ind.unlab] <- NA
        ssvm <- selfTraining(x=X, y = ytrain, x.inst = TRUE, learner = svm,
                             learner.pars =list(kernel ="linear", scale = FALSE, probability = TRUE),
                             pred = function(m, x){
                                 attr(predict(m, x, probability = TRUE), "probabilities")
                             })
        pred <- predict(ssvm, X)
    } # pred
    
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(as.numeric(as.character(pred)))))] <- 'doublet'
    
    
    ######################
    ### STEP 4: DEBRIS ###
    ######################
    X <- s.tech[which(out == 'cell'), ]
    
    # approximate the label
    {
        debriscore <- X[,'Residual'] - 
            X[,'Event_length'] -
            X[,'DNA2'] -
            X[,'DNA1'] -
            X[,'Center'] -
            .5 * X[,'Width']
        require(mclust)
        g <- Mclust(debriscore, G=2)$classification
        # check for "1 wide, 1 narrow peak" issue
        r1 <- range(debriscore[which(g==1)])
        r2 <- range(debriscore[which(g==2)])
        m1 <- mean(debriscore[which(g==1)])
        m2 <- mean(debriscore[which(g==2)])
        if(r1[1] < r2[1] & r1[2] > r2[2]){ # group 1 is wide
            if(m1 > m2){
                g[which(g==1 & debriscore < m2)] <- 2
            }else{
                g[which(g==1 & debriscore > m2)] <- 2
            }
        }
        if(r2[1] < r1[1] & r2[2] > r1[2]){ # group 2 is wide
            if(m2 > m1){
                g[which(g==2 & debriscore < m1)] <- 1
            }else{
                g[which(g==2 & debriscore > m1)] <- 1
            }
        }
        
        debrisclus <- which.max(c(m1,m2))
        init <- g == debrisclus
    } # init
    
    # pick "true" events (all "neighbors of neighbors" share the same label)
    {
        mismatch <- matrix(init[knn$index[which(out=='cell'),]], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        neighbormismatch <- matrix(rowMeans(mismatch)[knn$index[which(out=='cell'),]], ncol = 10)
        neighbormismatch[is.na(neighbormismatch)] <- 0
        #ind.label <- which(rowSums(neighbormismatch) == 0)
        ind.unlab <- which(rowSums(neighbormismatch) != 0)
    } # init.unlab
    
    # semi-supervised SVM to classify all events
    {
        library(ssc)
        library(e1071)
        ytrain <- as.integer(init)
        ytrain[ind.unlab] <- NA
        ssvm <- selfTraining(x=X, y = ytrain, x.inst = TRUE, learner = svm,
                             learner.pars =list(kernel ="linear", scale = FALSE, probability = TRUE),
                             pred = function(m, x){
                                 attr(predict(m, x, probability = TRUE), "probabilities")
                             })
        pred <- predict(ssvm, X)
    } # pred
    
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(as.numeric(as.character(pred)))))] <- 'debris'
    
    out <- factor(out, levels = c("cell","debris","doublet","bead","GDPzero"))
    
    return(out)
}



# out3 <- runSSSVM(sce)



# hist(tech[which(out=='cell'),'Offset'])
# hist(tech[which(out=='cell'),'Width'])
# hist(tech[which(out=='cell'),'Center'])
# hist(tech[which(out=='cell'),'Residual'])
# hist(tech[which(out=='cell'),'DNA1'], breaks = 50)
# 
# 
# ind <- sample(which(out == 'cell'), 5000)
# pairs(tech[ind, c('DNA1','Width','Event_length','Residual')], col = rgb(0,0,0,.1))
# 
# 
# 
# ind <- sample(which(out == 'cell'), 5000)
# pairs(tech[,c('Center','Event_length','DNA1')][ind,], col = rgb(0,0,0,.1))
# 
# plot(density(tech[which(out=='cell'),'DNA1']))
