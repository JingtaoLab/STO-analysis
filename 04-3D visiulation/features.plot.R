setwd("E:/project/CS7/bin50_new/figures/figureS2/")
load("E:/project/CS7/bin50_new/cluster_RDS/CS7.rds")
cs7_2$final=Idents(cs7_2)

genes=c( "SOX2", "SOX17",'BMP4', "GABRP", "TBXT","NANOS3","SHH","IGF2")
exp=FetchData(cs7_2,genes)
exp$barcode=rownames(exp)

exp1=reshape2::melt(exp)
exp2=subset(exp1,value>0)

data=cs7_2@meta.data
data1=merge(data,exp2,by="barcode",all=T)

library(RColorBrewer)
library(plotly)
library(geomorph)
library(Seurat)
library(RColorBrewer)
library(plotly)
library(geomorph)
library(Seurat)



library(plotly)
SOX17=subset(data1,variable=='SOX17')
p=plot_ly(x=SOX17$newz_paste*10,y=SOX17$newy_paste,z=SOX17$newx_paste,
            #size = SOX17$size,sizes = c(1,200),
            
            color = as.numeric(SOX17$value),###400
            colors = colorRamp(c("lightgrey", "red")),
            marker = list(#color="red",
              size=3,
              line = list(width=0,opacity = 0.1)),
            type = "scatter3d", mode = "markers")
