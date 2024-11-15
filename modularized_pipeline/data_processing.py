import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
from anndata import AnnData

def create_geodataframe(polys):
    geometries = [Polygon([(y, x) for x, y in zip(poly[0], poly[1])]) for poly in polys]
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i + 1}" for i, _ in enumerate(gdf.index)]
    return gdf

def merge_tissue_positions(adata, tissue_positions_file):
    print(tissue_position_file)
    # Read the tissue positions file into a DataFrame
    #tissue_positions = pd.read_parquet(tissue_positions_file, engine = 'pyarrow')
    # Set the 'barcode' column as the index
    #tissue_positions = tissue_positions.set_index('barcode')
    # Add a new column 'index' with the index value
    #tissue_positions['index'] = tissue_positions.index
    # Merge the tissue positions with adata.obs (the metadata)
    #adata.obs = pd.merge(adata.obs, tissue_positions, left_index=True, right_index=True)
    return adata


def filter_spatial_overlap(gdf, adata):
    # Perform spatial joins and remove overlapping nuclei logic here
    pass

def sum_gene_counts(adata, group_by_key):
    # Your logic for summing gene counts here
    pass

