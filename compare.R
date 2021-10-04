
fname <- '../ctc/data/Crompton/U2OS_NT-1_GFP_Spike_NORMAL_01_1.fcs'
datatype <- 'Crompton' # 'Bagwell'

base <- basename(fname)
base <- gsub('.fcs','', base)

source('R/utils.R')
source('R/catalyst.R')
source('R/bagwell.R')
source('R/nnmap.R')
source('R/mnnmap.R')
source('R/sssvm.R')
source('R/nnprop.R')
sce <- prepData(fname)


# setup

if(datatype == 'Crompton'){
    bead_channels <- c(19, 30, 32, 44, 54) # "dvs" from normCytof
    tech <- t(assay(sce,'exprs')[bead_channels,])
    tech_bead_channels <- 1:ncol(tech)
    tech <- cbind(tech, t(assay(sce,'exprs')[rownames(sce) %in% c('DNA1','DNA2','Viability'),])) # 'VeriCells' ??
    tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))
}else{
    bead_channels <- grep('Bead', rownames(sce))
    tech <- t(assay(sce,'exprs')[bead_channels,])
    tech_bead_channels <- 1:ncol(tech)
    tech <- cbind(tech, t(assay(sce,'exprs')[c('DNA1','DNA2','Live_Dead'),]))
    tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))
}
umap <- uwot::umap(tech)


t1 <- system.time(catalyst <- runCATALYST(sce))
t2 <- system.time(bagwell <- runBagwell(sce, datatype = datatype))
t3 <- system.time(nnsvm <- runNNSVM(sce, datatype = datatype))
#t4 <- system.time(nnsvm2 <- runNNSVM(sce, secondguess = TRUE, datatype = datatype))
t4 <- system.time(mnnsvm <- runMNNSVM(sce, datatype = datatype))
t5 <- system.time(nnprop <- runNNprop(sce, datatype = datatype))
#t5 <- system.time(sssvm <- runSSSVM(sce, datatype = datatype))

res <- list(umap = umap, times = list(catalyst = t1, bagwell = t2, nnsvm = t3, mnnsvm = t4, nnprop = t5),
            labels = data.frame(catalyst, bagwell, nnsvm, mnnsvm, nnprop))

saveRDS(res, file=paste0('../ctc/RESULTS_',base,'.rds'))

