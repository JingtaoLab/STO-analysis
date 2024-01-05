
adata = sc.read_loom('obj.loom', sparse=True, cleanup=False, dtype='float32')
adata.obsm['spatial']=np.array(adata.obs[["newx", "newy"]])
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, percent_top=None, qc_vars=["mt"], log1p=False, inplace=True)

adata
##filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) ###2000 HVG,n_top_genes=5000
#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,n_top_genes=4500)
sc.pp.regress_out(adata,['nCount_RNA'])

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_pcs=20)##35,
sq.gr.spatial_neighbors(adata)
###SSC
conn = adata.obsp['connectivities']
conn.data[conn.data > 0] = 1
adj = conn + adata.obsp['spatial_connectivities']
adj.data[adj.data  > 0] = 1
sc.tl.umap(adata)
sc.tl.leiden(adata, adjacency=adj, resolution = 1.3, key_added = "clusters" )

sc.pl.spatial(adata, color='clusters', spot_size=50, legend_fontsize='x-small',save='cluster.pdf')

adata.write("data.h5ad")
