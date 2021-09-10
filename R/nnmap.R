

# source('R/utils.R')
# fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
# sce <- prepData(fname)



runNNSVM <- function(sce, secondguess = FALSE, datatype = c('Bagwell','Crompton')){
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
    unclassified.ind <- which(out == 'cell')
    X <- s.tech[unclassified.ind, ]
    
    # approximate the label (init)
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
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- rowMeans(isBeadMat) > 0.5
    }
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 1)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
    }
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = s.tech[ind,], y = factor(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, s.tech[unclassified.ind, ])
        if(secondguess){
            # re-pick
            {
                init <- rep(FALSE, ncol(sce))
                init[unclassified.ind] <- as.logical(pred)
                mismatch <- matrix(init[knn$index], ncol = 10)
                mismatch[init, ] <- !mismatch[init, ]
                mismatch[is.na(mismatch)] <- FALSE
                poss.ind <- which(rowSums(mismatch) >= 2)
                # ~ equal representation
                poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
                ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
            }
            svmfit2 <- svm(x = s.tech[ind,], y = factor(init[ind]),
                           kernel = "linear", scale = FALSE)
            pred <- predict(svmfit2, s.tech[unclassified.ind, ])
        }
    }
    
    # update labels
    out[unclassified.ind][as.logical(pred)] <- 'bead'
    
    
    
    ########################
    ### STEP 3: DOUBLETS ###
    ########################
    unclassified.ind <- which(out == 'cell')
    X <- s.tech[unclassified.ind, ]
    
    # approximate the label (init)
    {
        doubletscore <- X[,'Residual'] + 
            X[,'Event_length'] +
            X[,'DNA2'] +
            X[,'DNA1'] -
            X[,'Offset'] -
            .5 * X[,'Width']
        g <- find_groups(doubletscore)
        doubletclus <- max(as.numeric(levels(g)))
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- g == doubletclus
    } # init
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 2)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
    } # ind
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = s.tech[ind,], y = factor(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, s.tech[unclassified.ind, ])
        if(secondguess){
            # re-pick
            {
                init <- rep(FALSE, ncol(sce))
                init[unclassified.ind] <- as.logical(pred)
                mismatch <- matrix(init[knn$index], ncol = 10)
                mismatch[init, ] <- !mismatch[init, ]
                mismatch[is.na(mismatch)] <- FALSE
                poss.ind <- which(rowSums(mismatch) >= 2)
                # ~ equal representation
                poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
                ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
            }
            svmfit2 <- svm(x = s.tech[ind,], y = factor(init[ind]),
                           kernel = "linear", scale = FALSE)
            pred <- predict(svmfit2, s.tech[unclassified.ind, ])
        }
    } # pred
    
    # update labels
    out[unclassified.ind][as.logical(pred)] <- 'doublet'
    
    
    ######################
    ### STEP 4: DEBRIS ###
    ######################
    unclassified.ind <- which(out == 'cell')
    X <- s.tech[unclassified.ind, ]
    
    # approximate the label (init)
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
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- g == debrisclus
    } # init
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 2)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
    } # ind
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = s.tech[ind,], y = factor(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, s.tech[unclassified.ind, ])
        if(secondguess){
            # re-pick
            {
                init <- rep(FALSE, ncol(sce))
                init[unclassified.ind] <- as.logical(pred)
                mismatch <- matrix(init[knn$index], ncol = 10)
                mismatch[init, ] <- !mismatch[init, ]
                mismatch[is.na(mismatch)] <- FALSE
                poss.ind <- which(rowSums(mismatch) >= 2)
                # ~ equal representation
                poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
                ind <- sample(poss.ind, min(c(4000, length(poss.ind))), prob = poss.wt)
            }
            svmfit2 <- svm(x = s.tech[ind,], y = factor(init[ind]),
                           kernel = "linear", scale = FALSE)
            pred <- predict(svmfit2, s.tech[unclassified.ind, ])
        }
    } # pred
    
    # update labels
    out[unclassified.ind][as.logical(pred)] <- 'debris'
    
    out <- factor(out, levels = c("cell","debris","doublet","bead","GDPzero"))
    
    return(out)
}




# out1 <- runNNSVM(sce)
# out2 <- runNNSVM(sce, secondguess = TRUE)

#barplot(table(out))






# rd <- uwot::umap(X[,c('Residual','Event_length','DNA1','DNA2','Center','Width')])
# tmpsce <- SingleCellExperiment(assays = list(exprs = assay(sce,'exprs')[,which(out == 'cell')]),
#                                reducedDims = list(umap = rd))
# tmpsce <- make_hexbin(tmpsce, dimension_reduction = 'umap', nbins = 40)
# 
# plot_hexbin_density(tmpsce)
# tmpsce$init <- factor(init)
# plot_hexbin_meta(tmpsce, col="init", action="prop")

#plot(rd, col=alpha(c('blue','red')[1+init],alpha=.5),asp=1, cex=.5)

