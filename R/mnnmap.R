

# source('R/utils.R')
# fname <- '../ctc/data/Crompton/U2OS_NT-1_GFP_Spike_NORMAL_01_1.fcs'
# sce <- prepData(fname)



runMNNSVM <- function(sce, datatype = c('Bagwell','Crompton')){
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
    knn <- findKNN(s.tech, 20, BNINDEX=k.dex)
    mnn <- knn$index
    for(i in 1:ncol(sce)){
        x <- knn$index[i,]
        keep <- rowMeans(knn$index[x,] == i) > 0
        mnn[i, !keep] <- NA
    }
    
    
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
        init <- rowMeans(isBeadMat) > 0.5
    }
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[mnn[which(out=='cell'),]], ncol = 20)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 1)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, 4000, prob = poss.wt)
    }
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = X[ind,], y = as.numeric(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, X)
    }
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(pred)))] <- 'bead'
    
    
    
    ########################
    ### STEP 3: DOUBLETS ###
    ########################
    X <- s.tech[which(out == 'cell'), ]
    
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
        init <- g == doubletclus
    } # init
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[mnn[which(out=='cell'),]], ncol = 20)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 2)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, 4000, prob = poss.wt)
    } # ind
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = X[ind,], y = as.numeric(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, X)
    } # pred
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(pred)))] <- 'doublet'
    
    
    ######################
    ### STEP 4: DEBRIS ###
    ######################
    X <- s.tech[which(out == 'cell'), ]
    
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
        init <- g == debrisclus
    } # init
    
    # pick events near the boundary
    {
        mismatch <- matrix(init[mnn[which(out=='cell'),]], ncol = 20)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        poss.ind <- which(rowSums(mismatch) >= 2)
        # ~ equal representation
        poss.wt <- (1000/table(init[poss.ind]))[1+init[poss.ind]]
        ind <- sample(poss.ind, 4000, prob = poss.wt)
    } # ind
    
    
    # SVM on those events
    {
        require(e1071)
        svmfit <- svm(x = X[ind,], y = as.numeric(init[ind]),
                      kernel = "linear", scale = FALSE)
        pred <- predict(svmfit, X)
    } # pred
    
    # update labels
    out[which(out == 'cell')][which(as.logical(round(pred)))] <- 'debris'
    
    out <- factor(out, levels = c("cell","debris","doublet","bead","GDPzero"))
    
    return(out)
}




# out <- runMNNSVM(sce)

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

