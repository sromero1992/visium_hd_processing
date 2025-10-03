import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, Point
from scipy import sparse
from anndata import AnnData

from shapely.geometry import Polygon
import geopandas as gpd

def create_geodataframe(polys):
    """
    Creates a GeoDataFrame by performing the (X, Y) -> (Y, X) coordinate swap 
    on the polygon vertices, matching the explicit swap in the 10x pipeline code.
    
    It assumes polys[0] contains all X-coordinates and polys[1] contains all Y-coordinates.
    """
    
    # This line implements the swap: [(y, x) for x, y in zip(poly[0], poly[1])]
    # Where:
    # poly[0] is X-coordinates (assigned to 'x' in the inner loop)
    # poly[1] is Y-coordinates (assigned to 'y' in the inner loop)
    # The output tuple is (y, x) -> resulting in (Y, X) order for the polygon vertices.
    geometries = [Polygon([(y, x) for x, y in zip(poly[0], poly[1])]) for poly in polys]
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i + 1}" for i, _ in enumerate(gdf.index)]
    return gdf

def merge_tissue_positions(adata, tissue_positions_file):
    # Read the tissue positions file into a DataFrame
    df_tissue_positions = pd.read_parquet(tissue_positions_file, engine='pyarrow')
    # Set 'barcode' as the index in df_tissue_positions
    df_tissue_positions = df_tissue_positions.set_index('barcode')
    # Create an index in the dataframe to check joins
    df_tissue_positions['index'] = df_tissue_positions.index
    # Adding the tissue positions to the meta data
    adata.obs =  pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True)
    return adata, df_tissue_positions

def create_geodf_coords(df_tissue_positions):
    # Create a GeoDataFrame from the DataFrame of coordinates
    geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'], df_tissue_positions['pxl_row_in_fullres'])]
    gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)
    return gdf_coordinates


def filter_spatial_overlap(adata, gdf, gdf_coordinates):
    # Perform a spatial join to check which coordinates are in a cell nucleus
    result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how='left', predicate='within')
    
    # Identify nuclei associated barcodes and find barcodes that are in more than one nucleus
    result_spatial_join['is_within_polygon'] = ~result_spatial_join['index_right'].isna()
    barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
    result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlaping_polygons)
    
    # Remove barcodes in overlapping nuclei
    barcodes_in_one_polygon = result_spatial_join[result_spatial_join['is_within_polygon'] & result_spatial_join['is_not_in_an_polygon_overlap']]
    
    # The AnnData object is filtered to only contain the barcodes that are in non-overlapping polygon regions
    filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
    filtered_adata = adata[filtered_obs_mask,:]
    
    # Add the results of the point spatial join to the Anndata object
    filtered_adata.obs =  pd.merge(
            filtered_adata.obs, 
            barcodes_in_one_polygon[['index','geometry','id','is_within_polygon','is_not_in_an_polygon_overlap']], 
            left_index=True, right_index=True)
    
    return filtered_adata


def filter_spatial_overlap_assign(adata: AnnData, gdf: gpd.GeoDataFrame, gdf_coordinates: gpd.GeoDataFrame) -> AnnData:
    """
    Filters AnnData object to retain only bins whose center falls within a nucleus polygon.
    In case of overlap (bin center is in multiple nuclei), the bin is assigned to the 
    first matching nucleus (less strict than removing all overlapping bins).

    Args:
        adata: The AnnData object containing the binned gene counts.
        gdf: GeoDataFrame of nucleus polygons (must be in the same coordinate system as gdf_coordinates).
        gdf_coordinates: GeoDataFrame of bin center points (barcodes).

    Returns:
        A new AnnData object filtered to include only assigned bins.
    """
    
    # 1. Perform a spatial join: Find which bin centers fall within a nucleus
    # The result has multiple rows per barcode if it overlaps multiple polygons.
    # The 'index' column of the result is the bin barcode.
    # The 'index_right' column is the index/ID of the nucleus polygon.
    result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how='left', predicate='within')
    
    # 2. Identify and filter for bins inside ANY nucleus polygon
    # index_right is NaN if the bin center is not within any polygon
    barcodes_in_any_polygon = result_spatial_join[~result_spatial_join['index_right'].isna()]
    
    # 3. Resolve overlaps: Keep only the FIRST match for each bin barcode.
    # This assigns the bin to a single nucleus, even if it overlapped multiple.
    barcodes_assigned_to_one = barcodes_in_any_polygon.drop_duplicates(subset=['index'], keep='first')
    
    # 4. Filter the AnnData object
    # The 'index' column here holds the filtered barcode list.
    filtered_obs_mask = adata.obs_names.isin(barcodes_assigned_to_one['index'])
    filtered_adata = adata[filtered_obs_mask, :].copy() # .copy() avoids AnnData view warnings
    
    # 5. Add the spatial join results (nucleus ID and geometry) to the AnnData object
    # Rename 'index_right' to something meaningful like 'nucleus_match_id'
    barcodes_to_merge = (
        barcodes_assigned_to_one[['index', 'geometry', 'id', 'index_right']]
        .rename(columns={'index_right': 'nucleus_match_id', 'id': 'nucleus_id'})
        .set_index('index') # Set the barcode as the index for merging
    )
    
    # Merge the nucleus assignment information back into filtered_adata.obs
    filtered_adata.obs = pd.merge(
        filtered_adata.obs, 
        barcodes_to_merge, 
        left_index=True, 
        right_index=True,
        how='left' # Use left merge to keep all filtered adata obs
    )
    
    return filtered_adata

def sum_gene_counts_binned(filtered_adata):
    # Group the data by unique nucleous IDs
    groupby_object = filtered_adata.obs.groupby(['id'], observed=True)
    
    # Extract the gene expression counts from the AnnData object
    counts = filtered_adata.X
    
    # Obtain the number of unique nuclei and the number of genes in the expression data
    N_groups = groupby_object.ngroups
    N_genes = counts.shape[1]
    
    # Initialize a sparse matrix to store the summed gene counts for each nucleus
    summed_counts = sparse.lil_matrix((N_groups, N_genes))
    
    # Lists to store the IDs of polygons and the current row index
    polygon_id = []
    row = 0
    
    # Iterate over each unique polygon to calculate the sum of gene counts.
    for polygons, idx_ in groupby_object.indices.items():
        summed_counts[row] = counts[idx_].sum(0)
        row += 1
        polygon_id.append(polygons)
    
    # Create and AnnData object from the summed count matrix
    summed_counts = summed_counts.tocsr()
    grouped_filtered_adata = AnnData(X=summed_counts,obs=pd.DataFrame(polygon_id,columns=['id'],index=polygon_id),var=filtered_adata.var)
    
    return grouped_filtered_adata


