# do all the dimensionality reductions

filenames <- system('ls ~/Desktop/Bagwell_data', intern=TRUE)

for(fname in filenames){
    print(fname)
    fname <- paste0('~/Desktop/Bagwell_data/', fname)
    datatype <- 'Bagwell' #'Crompton'
    
    base <- basename(fname)
    base <- gsub('.fcs','', base)
    require(CATALYST)
    sce <- prepData(fname)
    
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
    s.tech <- scale(tech)
    
    require(BiocSingular)
    pca <- BiocSingular::runPCA(s.tech, rank = 3)
    require(uwot)
    set.seed(1)
    umap1 <- uwot::umap(s.tech)
    set.seed(1)
    umap2 <- uwot::umap(s.tech, n_neighbors = 30, metric = "cosine")
    require(Rtsne)
    set.seed(1)
    tsne <- Rtsne::Rtsne(s.tech, check_duplicates = FALSE)$Y
    require(snifter)
    set.seed(1)
    fitsne <- fitsne(s.tech)[,]
    
    out <- list(pca = pca$x,
                umap1 = umap1,
                umap2 = umap2,
                tsne = tsne,
                fitsne = fitsne,
                meta = list(pca.sdev = pca$sdev))
    
    saveRDS(out, file = paste0('~/Desktop/dimRed/DR_',base,'.rds'))
    
}

