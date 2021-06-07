

# things we want to identify

# maybe: zeros for GDPs    # easy
# beads                    # 
# bead+cell doublets       # 
# debris
# doublets

require(CATALYST)
fname <- 'data/Bagwell/REP_1_deid.fcs'
sce <- prepData(fname)
out <- rep('cell', ncol(sce))

# setup
bead_channels <- grep('Bead', rownames(sce))
tech <- t(assay(sce,'exprs')[bead_channels,])
tech <- cbind(tech, t(assay(sce,'exprs')[c('DNA1','DNA2','Live_Dead'),]))
tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))


# STEP 1: MARK GDP ZEROS
out[int_colData(sce)$Center == 0 |
        int_colData(sce)$Offset == 0 |
        int_colData(sce)$Width == 0 |
        int_colData(sce)$Residual == 0] <- 'GDPzero'


# STEP 2: FIND BEADS & BEAD+CELL DOUBLETS
# initialize with peak finder, then SVM
find_groups <- function(x){
    d <- density(x, adjust = 2)
    # initial peaks:
    #  (1) higher than it's 10 neighbors
    #  (2) at leat 1% of the maximum peak height
    peakTF <- sapply(1:length(d$y), function(i){
        neighbors <- c((i-5):(i-1), (i+1):(i+5))
        neighbors <- neighbors[neighbors %in% 1:length(d$y)]
        
        return(all(d$y[i] > d$y[neighbors]) &
                   d$y[i] > .01*max(d$y))
    })
    peaks <- which(peakTF)
    # find valleys (lowpoints between peaks)
    if(sum(peakTF) > 1){
        valleys <- sapply(1:(sum(peakTF)-1), function(i){
            p1 <- which(peakTF)[i]
            p2 <- which(peakTF)[i+1]
            btwn <- p1:p2
            return(btwn[which.min(d$y[btwn])])
        })
    }
    # peak filtering
    # remove peak if not >120% of both nearby valleys
    for(ii in 1:length(peaks)){
        ht <- d$y[peaks[ii]]
        if(ii == 1){
            lv <- 1
        }else{
            lv <- valleys[ii-1]
        }
        if(ii == length(peaks)){
            rv <- length(d$x)
        }else{
            rv <- valleys[ii]
        }
        lv <- d$y[lv]
        rv <- d$y[rv]
        if(ht < 1.2*lv | ht < 1.2*rv){
            peakTF[peaks[ii]] <- FALSE
            print(paste('removing',ii))
        }
    }
    peaks <- which(peakTF)
    # re-find valleys
    if(sum(peakTF) > 1){
        valleys <- sapply(1:(sum(peakTF)-1), function(i){
            p1 <- which(peakTF)[i]
            p2 <- which(peakTF)[i+1]
            btwn <- p1:p2
            return(btwn[which.min(d$y[btwn])])
        })
        out <- cut(x, breaks = c(-Inf, d$x[valleys], Inf))
    }else{
        out <- factor(rep(1, length(x)))
    }
    levels(out) <- 1:length(peaks)
    return(out)
}

isBeadMat <- sapply(bead_channels, function(ch){
    x <- assay(sce,'exprs')[ch,]
    g <- find_groups(x[x>0])
    mns <- by(x[x>0], g, mean)
    beadgp <- which.max(mns)
    out <- rep(FALSE, ncol(sce))
    if(length(levels(g)) > 1){
        out[x>0][g == beadgp] <- TRUE
    }
    return(out)
})
include <- rowMeans(isBeadMat) > 0.01
# add random events, in case initial classification is perfect
include[sample(which(!include), 1000)] <- TRUE
init <- rowMeans(isBeadMat)[include] > 0.5

#pairs(t(assay(sce,'exprs')[c(grep('Bead', rownames(sce)),43:44), include]), col = colorby(init, alpha=.5))

require(e1071)
svmfit <- svm(x = t(assay(sce,'exprs')[bead_channels, include]), y = as.numeric(init),
              kernel = "linear", scale = FALSE)
pred <- predict(svmfit, t(assay(sce,'exprs')[bead_channels, include]))

#pairs(t(assay(sce,'exprs')[c(grep('Bead', rownames(sce)),43:44), include]), col = colorby(round(pred), alpha=.5))

out[include][pred > 0.5] <- 'bead'
out[int_colData(sce)$Center == 0 |
        int_colData(sce)$Offset == 0 |
        int_colData(sce)$Width == 0 |
        int_colData(sce)$Residual == 0] <- 'GDPzero'



# STEP 2A: SEPARATE BEADS FROM BEAD+CELL DOUBLETS

# pairs(tech[which(out=='bead'),c('DNA1','DNA2','Event_length','Center','Offset','Width','Residual')], col = rgb(0,0,0,.3))
# subset = possible bead+cell doublets
subset <- which(out == 'bead' & 
                    assay(sce,'exprs')['DNA1',] > 0 &
                    assay(sce,'exprs')['DNA1',] > 0)

require(mclust)
mc <- Mclust(tech[subset, c('DNA1','DNA2','Event_length','Center','Offset','Width','Residual')],
             G = 2)
doubletclus <- which.max(sapply(1:2,function(clID){
    mean(assay(sce,'exprs')[c('DNA1','DNA2'), subset[mc$classification == clID]])
}))
out[subset[mc$classification == doubletclus]] <- 'beadCell'




