import scanpy as sc

def calculate_qc_metrics(adata):
    sc.pp.calculate_qc_metrics(adata, inplace=True)

def filter_and_normalize(adata, area_mask, count_mask):
    filtered_adata = adata[area_mask & count_mask, :]
    sc.pp.calculate_qc_metrics(filtered_adata, inplace=True)
    sc.pp.normalize_total(filtered_adata, inplace=True)
    sc.pp.log1p(filtered_adata)
    return filtered_adata

def perform_clustering(adata, resolution=0.35):
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="clusters")
    return adata

