#R

####拆分slice
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(ggplot2)
flist <- list.files(path = ".",pattern = "*h5ad$")
for (i in flist){
  Convert(i, "h5seurat",
          overwrite = TRUE,assay = "RNA")
}


###slice1
slice1 = LoadH5Seurat("CP1-1.bin50_seurat.h5seurat")
sliceA=subset(slice1,x>5000&x<11000&y>5500&y<10000) ####according the distribution of stereo-seq
sliceB=subset(slice1,x>17000&x<21000&y>6000&y<13000)
sliceC=subset(slice1,x>17000&x<19000&y>15900&y<18000)
sliceD=subset(slice1,x>9000&x<12000&y>15000&y<22000)

plot2=function(obj){
  data=data.frame(x=obj$newx,y=obj$newy,cluster=obj$spatial_leiden)
  p=ggplot(data,aes(x=x,y=y,color=cluster))+geom_point(size=0.5)+theme_bw()+
    scale_x_continuous(breaks = seq(0,30000,2000))+scale_y_continuous(breaks = seq(0,25000,2000))
  return(p)
}

sliceA$newy=sliceA$y-min(sliceA$y)
sliceA$newx=sliceA$x-min(sliceA$x)
sliceA$sample="S1"
plot2(sliceA)

sliceB$newy=sliceB$y-min(sliceB$y)
sliceB$newx=sliceB$x-min(sliceB$x)
sliceB$sample="S2"
plot2(sliceB)

sliceC$newy=abs(max(sliceC$y)-sliceC$y)
sliceC$newx=abs(max(sliceC$x)-sliceC$x)
sliceC$sample="S3"
plot2(sliceC)

sliceD$newy=abs(max(sliceD$y)-sliceD$y)
sliceD$newx=abs(max(sliceD$x)-sliceD$x)
sliceD$sample="S4"
plot2(sliceD)

slice1=merge(sliceA,y=list(sliceB,sliceC,sliceD),
                   add.cell.ids = c("S1","S2","S3","S4"))

######merge
obj=merge(x=slice1,y=list(slice2,slice3,slice4,slice5,slice6,slice7,slice8,slice9,slice10,slice11,slice12),add.cell.ids = c("slice1","slice2","slice3","sice4","slice5","slice6","slice7","slice8","slice9","slice10","slice11","slice12"))

####
###convert obj to loom for cluster
obj=CreateSeuratObject(obj@assays$RNA@counts,meta.data = obj@meta.data)
main.loom <- as.loom(obj, filename = "obj.loom", 
                     verbose = T)
