import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scanpy as sc
import paste as pst

# Load Slices
data_dir = './data/' # change this path to the data you wish to analyze

adata=sc.read_h5ad('cs8.input.h5ad')
adata.obsm['spatial']=np.array(adata.obs[["newx", "newy"]])

adata1=adata[adata.obs['sample'] == 'S1'].copy()
adata2=adata[adata.obs['sample'] == 'S2'].copy()
adata3=adata[adata.obs['sample'] == 'S2'].copy()
adata4=adata[adata.obs['sample'] == 'S2'].copy()

pi1= pst.pairwise_align(S1, S2)
pi2= pst.pairwise_align(S2, S3)
pi3= pst.pairwise_align(S3, S4)

pis = [pi1, pi2, pi3]
slices = [adata1, adata2, adata3, adata4]
new_slices = pst.stack_slices_pairwise(slices, pis)

slice1,slice2,slice3,slice4=new_slices
names = new_slices.obs_names.tolist()
coor=new_slices.obsm['spatial']



