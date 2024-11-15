import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from shapely.geometry import Polygon
from plotting_utils import crop_image, filter_geodataframe_by_bbox, plot_cropped_image

# General image plotting functions
def plot_mask_and_save_image(title, gdf, img, cmap, output_name=None, bbox=None):
    """Plots and saves a mask image with optional bounding box."""
    cropped_img = crop_image(img, bbox)
    filtered_gdf = filter_geodataframe_by_bbox(gdf, bbox)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plot_cropped_image(axes[0], cropped_img, title)
    
    filtered_gdf.plot(cmap=cmap, ax=axes[1])
    axes[1].axis('off')

    if output_name:
        plt.savefig(output_name, bbox_inches='tight')
    else:
        plt.show()

def merge_gene_with_geodataframe(gdf, gene, adata):
    """Merges gene expression data with the GeoDataFrame."""
    gene_expression = adata[:, gene].to_df()
    gene_expression['id'] = gene_expression.index
    return gdf.merge(gene_expression, on='id')

def plot_gene_and_save_image(title, gdf, gene, img, adata, bbox=None, output_name=None):
    """Plots and saves a gene expression image with optional bounding box."""
    cropped_img = crop_image(img, bbox)
    merged_gdf = merge_gene_with_geodataframe(gdf, gene, adata)
    filtered_gdf = filter_geodataframe_by_bbox(merged_gdf, bbox)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plot_cropped_image(axes[0], cropped_img, title)
    
    filtered_gdf.plot(column=gene, cmap='inferno', legend=True, ax=axes[1])
    axes[1].set_title(gene)
    axes[1].axis('off')

    if output_name:
        plt.savefig(output_name, bbox_inches='tight')
    else:
        plt.show()

def create_custom_colormap(color_list, num_categories):
    """Creates a custom colormap with the specified color list."""
    if len(color_list) >= num_categories:
        return ListedColormap(color_list[:num_categories], name='custom_cmap')
    return ListedColormap(plt.cm.tab20.colors[:num_categories], name='custom_tab20_cmap')

def merge_clusters_with_geodataframe(gdf, color_by_obs, adata):
    """Merges clustering data with the GeoDataFrame."""
    return gdf.merge(adata.obs[color_by_obs].astype('category'), left_on='id', right_index=True)

def plot_clusters_and_save_image(title, gdf, img, adata, bbox=None, color_by_obs=None, output_name=None, color_list=None):
    """Plots and saves clustered nuclei with optional bounding box."""
    cropped_img = crop_image(img, bbox)
    merged_gdf = merge_clusters_with_geodataframe(gdf, color_by_obs, adata)
    filtered_gdf = filter_geodataframe_by_bbox(merged_gdf, bbox)
    
    num_categories = len(adata.obs[color_by_obs].astype('category').cat.categories)
    custom_cmap = create_custom_colormap(color_list or [], num_categories)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plot_cropped_image(axes[0], cropped_img, title)
    
    filtered_gdf.plot(column=color_by_obs, cmap=custom_cmap, ax=axes[1], legend=True)
    axes[1].axis('off')

    if output_name:
        plt.savefig(output_name, bbox_inches='tight')
    else:
        plt.show()

def plot_nuclei_area(gdf, area_cut_off):
    """Plots histograms of nuclei areas with a cutoff for filtering."""
    fig, axs = plt.subplots(1, 2, figsize=(15, 4))
    axs[0].hist(gdf['area'], bins=50, edgecolor='black')
    axs[0].set_title('Nuclei Area')

    axs[1].hist(gdf[gdf['area'] < area_cut_off]['area'], bins=50, edgecolor='black')
    axs[1].set_title(f'Nuclei Area Filtered: {area_cut_off}')

    plt.tight_layout()
    plt.show()

def total_umi(adata_, cut_off):
    """Plots the distribution of total UMI counts."""
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    axs[0].boxplot(adata_.obs["total_counts"], vert=False, patch_artist=True, boxprops=dict(facecolor='skyblue'))
    axs[0].set_title('Total Counts')

    axs[1].boxplot(adata_.obs["total_counts"][adata_.obs["total_counts"] > cut_off], vert=False,
                   patch_artist=True, boxprops=dict(facecolor='skyblue'))
    axs[1].set_title(f'Total Counts > {cut_off}')

    for ax in axs:
        ax.get_yaxis().set_visible(False)

    plt.tight_layout()
    plt.show()

