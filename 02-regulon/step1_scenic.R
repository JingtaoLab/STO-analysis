rm(list = ls()) 
library(SCENIC)

## Load data
loomPath <- system.file(package="SCENIC", "scenic.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
exprMat[1:4,1:4]
head(cellInfo)
table(cellInfo$CellType)

dim(exprMat)
write.csv(exprMat,"scenic.auc.csv",quote=False)

