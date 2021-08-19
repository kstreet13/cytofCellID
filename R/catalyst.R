
# require(CATALYST)
# fname <- '../ctc/data/Bagwell/REP_1_deid.fcs'
# sce <- prepData(fname)

runCATALYST <- function(sce){
    sce <- normCytof(sce, beads = 'dvs', remove_beads = FALSE)
    
    out <- rep('cell', ncol(sce$data))
    out[colData(sce$data)$remove] <- 'remove'
    out[colData(sce$data)$is_bead] <- 'bead'
    
    return(out)
}

