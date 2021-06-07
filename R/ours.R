
# things we want to identify

# maybe: zeros for GDPs    # easy
# beads                    # easy
# bead+cell doublets       # fairly easy
# debris                   # harder
# doublets                 # medium, but not sure if we want to remove them all

source('R/utils.R')
fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
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



# STEP ?: IDENTIFY DOUBLETS (CELL + CELL)

#       Center: low+high bimodal   
#       Offset: LOW   
#        Width: low-ish    
#     Residual: HIGH 
# Event_length: HIGH 
#          DNA: HIGH


hist(tech[which(out=='cell'),'Offset'])
hist(tech[which(out=='cell'),'Width'])
hist(tech[which(out=='cell'),'Center'])
hist(tech[which(out=='cell'),'Residual'])
hist(tech[which(out=='cell'),'DNA1'], breaks = 50)

foo <- tech[,'Width'] / tech[,'Event_length']

ind <- sample(which(out == 'cell'), 5000)

pairs(cbind(foo, tech[,c('Width','Event_length','DNA1')])[ind,], col = rgb(0,0,0,.1))


# STEP ?: IDENTIFY DEBRIS

#       Center: LOW
#       Offset: highly variable
#        Width: low-ish
#     Residual: high-ish 
# Event_length: LOW
#          DNA: LOW