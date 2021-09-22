

# source('R/utils.R')
# fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
# sce <- prepData(fname)



runNNprop <- function(sce, datatype = c('Bagwell','Crompton')){
    datatype <- match.arg(datatype)
    out <- rep(NA, ncol(sce))
    
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
    
    
    
    ### GDP ZEROS ###
    #################
    sure.gdpz <- rep('cell', ncol(sce))
    sure.gdpz[int_colData(sce)$Center == 0 |
                  int_colData(sce)$Offset == 0 |
                  int_colData(sce)$Width == 0 |
                  int_colData(sce)$Residual == 0] <- 'GDPzero'
    # label GDP zeros
    out[which(sure.gdpz == 'GDPzero')] <- 'GDPzero'

    
    # approximate the label
    # pick events we're relatively sure of
    # let nearest neighbors vote
    # (if 9/10 are "bead", it's a bead, etc.)

    ### BEADS ###
    #############
    unclassified.ind <- which(is.na(out))
    X <- s.tech[unclassified.ind, ]
    
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
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- rowMeans(isBeadMat) > 0.5
    } # init
    
    # pick events we're "sure of"
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        sure.ind <- which(rowSums(mismatch) == 0)
        sure.bead <- rep(NA, ncol(sce))
        sure.bead[sure.ind] <- c('cell','bead')[1+init[sure.ind]]
    } # sure.bead
    
    # label beads we're sure of
    out[which(sure.bead == 'bead')] <- 'bead'
    
    ### DOUBLETS ###
    ################
    unclassified.ind <- which(is.na(out))
    X <- s.tech[unclassified.ind, ]
    
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
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- g == doubletclus
    } # init
    
    # pick events we're "sure of"
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        sure.ind <- which(rowSums(mismatch) == 0)
        sure.doub <- rep(NA, ncol(sce))
        sure.doub[sure.ind] <- c('cell','doublet')[1+init[sure.ind]]
    } # sure.doub
    
    # label doublets we're sure of
    out[which(sure.doub == 'doublet')] <- 'doublet'

    ### DEBRIS ###
    ##############
    unclassified.ind <- which(is.na(out))
    X <- s.tech[unclassified.ind, ]
    
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
        init <- rep(FALSE, ncol(sce))
        init[unclassified.ind] <- g == debrisclus
    } # init
    
    # pick events we're "sure of"
    {
        mismatch <- matrix(init[knn$index], ncol = 10)
        mismatch[init, ] <- !mismatch[init, ]
        mismatch[is.na(mismatch)] <- FALSE
        sure.ind <- which(rowSums(mismatch) == 0)
        sure.debr <- rep(NA, ncol(sce))
        sure.debr[sure.ind] <- c('cell','debris')[1+init[sure.ind]]
    } # sure.debr
    
    # label debris we're sure of
    out[which(sure.debr == 'debris')] <- 'debris'
    
    
    ### CELLS ###
    #############
    # label cells we're sure of
    out[which(sure.gdpz == 'cell' &
              sure.bead == 'cell' &
              sure.doub == 'cell' &
              sure.debr == 'cell')] <- 'cell'
    
    
    ########################
    ### PROPOGATE LABELS ###
    ########################
    # slowly decrease the number of neighbors that need to agree
    for(cutoff in 9:1){
        print(cutoff)
        out.new <- out
        update_ind <- 1
        # let the labels fully propagate at this cutoff value
        while(length(update_ind) > 0){
            out <- out.new
            unlab <- which(is.na(out))
            # get neighbors of unlabeled events
            neighbor_labels <- matrix(out[knn$index[unlab,]], ncol = 10)
            inagreement <- apply(neighbor_labels, 1, function(x){
                suppressWarnings(max(table(x)))
            })
            # propogate labels where more than {cutoff} neighbors agree on a label
            update_ind <- which(inagreement >= cutoff)
            out.new[unlab[update_ind]] <- sapply(update_ind, function(ui){
                names(sort(table(neighbor_labels[ui,]), decreasing = TRUE))[1]
            })
            print(c(length(update_ind), sum(is.na(out.new))))
        }
    }
    
    out <- factor(out, levels = c("cell","debris","doublet","bead","GDPzero"))
    
    return(out)
}




# out1 <- runNNSVM(sce, datatype = 'Bagwell')
# out2 <- runNNprop(sce, datatype = 'Bagwell')
