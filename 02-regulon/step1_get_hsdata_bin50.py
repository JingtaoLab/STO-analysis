import os,sys,re
import pandas as pd
import scanpy as sc
import numpy as np

in_h5ad=sys.argv[1]
sample_name=sys.argv[2]
outdir=sys.argv[3]
organ=sys.argv[4]

adata=sc.read_h5ad(in_h5ad)
adata=adata[adata.obs['annotation'].isin([organ])]


adata.X[np.isnan(adata.X)]=0



exp=pd.DataFrame(adata.X)


exp.index=adata.obs.index.tolist()

gene_name=[i.replace("\"","") for i in adata.var.index.tolist()]
exp.columns=gene_name

#physical coords
x_coords=adata.obs['x'].tolist() #um
y_coords=adata.obs['y'].tolist()  #um

z_coords=adata.obs['z'].tolist() #um

coords=pd.DataFrame({'physical_x':x_coords,'physical_y':y_coords,'physical_z':z_coords})
coords.index=adata.obs.index.tolist()

exp.T.to_csv(outdir+"/"+sample_name+".exp.csv")
coords.to_csv(outdir+"/"+sample_name+".phy.coords.csv")
