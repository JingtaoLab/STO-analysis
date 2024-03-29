####extract gef from mapping redults
#!bin/python
import warnings

warnings.filterwarnings('ignore')
import stereo as st
###bin50
# read the gef file information
data_path = '/tissuecut/CP1-1.tissue.gef'
st.io.read_gef_info(data_path)

data = st.io.read_gef(file_path=data_path, bin_size=50)
data.tl.cal_qc()
#Filtering
data.tl.filter_cells(min_gene=20, min_n_genes_by_counts=3, pct_counts_mt=10, inplace=True)
data.tl.filter_genes(min_cell=2, inplace=True)
#Normalization
data.tl.raw_checkpoint()
data.tl.raw.exp_matrix
data.tl.normalize_total(target_sum=10000)
data.tl.log1p()

##
#HVG
data.tl.highly_variable_genes(min_mean=0.0125, max_mean=3,min_disp=0.5,
                              n_top_genes=2000, res_key='highly_variable_genes')
data.tl.scale(max_value=10, zero_center=False)
#pca
data.tl.pca(use_highly_genes=False, n_pcs=30, res_key='pca')
data.tl.neighbors(pca_res_key='pca', n_pcs=30, res_key='neighbors')
data.tl.spatial_neighbors(neighbors_res_key='neighbors',res_key='spatial_neighbors')
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')

data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')
data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')
data.tl.louvain(neighbors_res_key='neighbors', res_key='louvain')
data.tl.phenograph(phenograph_k=30, pca_res_key='pca', res_key='phenograph')

adata = st.io.stereo_to_anndata(data,flavor='seurat',output='CP1-1.bin50_seurat.h5ad')
