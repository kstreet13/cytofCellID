
require(CATALYST)
fname <- 'data/Bagwell/REP_1_deid.fcs'
sce <- prepData(fname)
sce <- normCytof(sce, beads = 'dvs')

