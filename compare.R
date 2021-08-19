
fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'

base <- basename(fname)
base <- gsub('.fcs','', base)

source('R/utils.R')
source('R/catalyst.R')
source('R/bagwell.R')
source('R/nnmap.R')
source('R/sssvm.R')
sce <- prepData(fname)


# setup
bead_channels <- grep('Bead', rownames(sce))
tech <- t(assay(sce,'exprs')[bead_channels,])
tech <- cbind(tech, t(assay(sce,'exprs')[c('DNA1','DNA2','Live_Dead'),]))
tech <- cbind(tech, log1p(as.matrix(int_colData(sce)[,c('Event_length','Center','Offset','Width','Residual')])))
umap <- uwot::umap(tech)

t1 <- system.time(catalyst <- runCATALYST(sce))
t2 <- system.time(bagwell <- runBagwell(sce))
t3 <- system.time(nnsvm <- runNNSVM(sce))
t4 <- system.time(nnsvm2 <- runNNSVM(sce, secondguess = TRUE))
t5 <- system.time(sssvm <- runSSSVM(sce))

res <- list(umap = umap, times = list(catalyst = t1, bagwell = t2, nnsvm = t3, nnsvm2 = t4, sssvm = t5),
            labels = data.frame(catalyst, bagwell, nnsvm, nnsvm2, sssvm))

saveRDS(res, file=paste0('~/RESULTS_',base,'.rds'))

