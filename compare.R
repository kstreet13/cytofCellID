
fname <- '../ctc/data/20210427 GFP staining protocol determination U20S NT-1GFP & EW8 GFP w ctrls - Normalized/EW8_01_1 (1).fcs'
datatype <- 'Crompton' # 'Bagwell'

base <- basename(fname)
base <- gsub('.fcs','', base)

source('R/utils.R')
source('R/catalyst.R')
source('R/bagwell.R')
source('R/nnmap.R')
source('R/sssvm.R')
sce <- prepData(fname)


# setup

if(datatype == 'Crompton'){
    bead_channels <- c(19, 21, 30, 32, 54)
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


t1 <- system.time(catalyst <- runCATALYST(sce, datatype = datatype))
t2 <- system.time(bagwell <- runBagwell(sce, datatype = datatype))
t3 <- system.time(nnsvm <- runNNSVM(sce, datatype = datatype))
t4 <- system.time(nnsvm2 <- runNNSVM(sce, secondguess = TRUE, datatype = datatype))
t5 <- system.time(sssvm <- runSSSVM(sce, datatype = datatype))

res <- list(umap = umap, times = list(catalyst = t1, bagwell = t2, nnsvm = t3, nnsvm2 = t4, sssvm = t5),
            labels = data.frame(catalyst, bagwell, nnsvm, nnsvm2, sssvm))

saveRDS(res, file=paste0('../ctc/RESULTS_',base,'.rds'))

