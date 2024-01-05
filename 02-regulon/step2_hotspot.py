import os,sys,re
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpb
mpb.use("Agg")
import seaborn as sns

auc_file=sys.argv[1]###scenic.auc.csv
pos_file=sys.argv[2]
outdir=sys.argv[3]
project_name=sys.argv[4]


pos = pd.read_csv(pos_file, index_col=0,sep=",")
counts = pd.read_csv(auc_file, index_col=0,sep=",")

counts = counts.loc[:, pos.index]
barcodes = pos.index.values

num_umi = counts.sum(axis=0)

gene_counts = (counts > 0).sum(axis=1)
valid_genes = gene_counts >= 50
counts = counts.loc[valid_genes]

hs = hotspot.Hotspot(counts, model='bernoulli', latent=pos, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=300,
)
hs_results = hs.compute_autocorrelations(jobs=16)
# Select the genes with significant spatial autocorrelation
hs_genes = hs_results.index[hs_results.FDR < 0.05]

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=20)
modules = hs.create_modules(
    min_gene_threshold=30, core_only=False, fdr_threshold=0.05
)
modules.value_counts()
hs.plot_local_correlations(yticklabels = True)

plt.savefig(outdir+"/"+project_name+"_"+"module_corr.pdf")

import pickle

f=open(outdir+"/"+project_name+'.hs','wb')
pickle.dump(hs,f)
