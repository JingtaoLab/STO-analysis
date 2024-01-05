import os,sys,re
import pandas as pd
import scanpy as sc
import numpy as np

in_h5ad=sys.argv[1]
sample_name=sys.argv[2]
outdir=sys.argv[3]

adata=sc.read_h5ad(in_h5ad)


adata.layers['raw_counts'][np.isnan(adata.layers['raw_counts'])]=0


exp=pd.DataFrame(adata.layers['raw_counts'])
exp.index=adata.obs.index.tolist()

gene_name=[i.replace("\"","") for i in adata.var.index.tolist()]
exp.columns=gene_name

#physical coords
###The z was from PASTE

x_coords=adata.obs['newx'].tolist() #um
y_coords=adata.obs['newy'].tolist()  #um
adata.obs['new_z']=adata.obs['new_z']*0.7 #
z_coords=adata.obs['new_z'].tolist() #um

coords=pd.DataFrame({'physical_x':x_coords,'physical_y':y_coords,'physical_z':z_coords})
coords.index=adata.obs.index.tolist()

exp.T.to_csv(outdir+"/"+sample_name+".exp.csv")
coords.to_csv(outdir+"/"+sample_name+".phy.coords.csv")
