import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from shapely.geometry import Polygon
from visium_hd_functions.plotting_utils import crop_image, filter_geodataframe_by_bbox, plot_cropped_image

# General image plotting functions
def plot_mask_and_save_image(title, gdf, img, cmap, output_name=None, bbox=None, dpi=300):
    """
    Plots and saves a mask image with optional bounding box, ensuring high-quality output.

    :param title: Title of the plot.
    :param gdf: GeodataFrame to plot.
    :param img: Image array to display.
    :param cmap: Colormap for the geodataframe plot.
    :param output_name: File name to save the output image (optional).
    :param bbox: Bounding box to crop the image and geodataframe (optional).
    :param dpi: Dots per inch for high-quality saving.
    """
    cropped_img = crop_image(img, bbox)
    filtered_gdf = filter_geodataframe_by_bbox(gdf, bbox)

    # Create subplots with a larger figure size for better resolution
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))  # Adjust figsize as needed
    plot_cropped_image(axes[0], cropped_img, title)

    filtered_gdf.plot(cmap=cmap, ax=axes[1])
    axes[1].axis('off')

    if output_name:
        # Save the figure with higher DPI for better quality
        plt.savefig(output_name, bbox_inches='tight', dpi=dpi)
    else:
        # Show the plot with better interactive quality
        plt.show()

    # Close the figure to free memory
    plt.close(fig)



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
    color_list=["#7f0000","#808000","#483d8b","#008000","#bc8f8f","#008b8b","#4682b4","#000080","#d2691e","#9acd32","#8fbc8f","#800080","#b03060","#ff4500","#ffa500","#ffff00","#00ff00","#8a2be2","#00ff7f","#dc143c","#00ffff","#0000ff","#ff00ff","#1e90ff","#f0e68c","#90ee90","#add8e6","#ff1493","#7b68ee","#ee82ee"]
    if bbox is not None:
        cropped_img = img[bbox[1]:bbox[3], bbox[0]:bbox[2]]
    else:
        cropped_img = img

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].imshow(cropped_img, cmap='gray', origin='lower')
    axes[0].set_title(title)
    axes[0].axis('off')

    if bbox is not None:
        bbox_polygon = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])

    unique_values = adata.obs[color_by_obs].astype('category').cat.categories
    num_categories = len(unique_values)

    if color_list is not None and len(color_list) >= num_categories:
        custom_cmap = ListedColormap(color_list[:num_categories], name='custom_cmap')
    else:
        # Use default tab20 colors if color_list is insufficient
        tab20_colors = plt.cm.tab20.colors[:num_categories]
        custom_cmap = ListedColormap(tab20_colors, name='custom_tab20_cmap')

    merged_gdf = gdf.merge(adata.obs[color_by_obs].astype('category'), left_on='id', right_index=True)

    if bbox is not None:
        intersects_bbox = merged_gdf['geometry'].intersects(bbox_polygon)
        filtered_gdf = merged_gdf[intersects_bbox]
    else:
        filtered_gdf = merged_gdf

    # Plot the filtered polygons on the second axis
    plot = filtered_gdf.plot(column=color_by_obs, cmap=custom_cmap, ax=axes[1], legend=True)
    axes[1].set_title(color_by_obs)
    legend = axes[1].get_legend()
    legend.set_bbox_to_anchor((1.05, 1))
    axes[1].axis('off')

    # Move legend outside the plot
    plot.get_legend().set_bbox_to_anchor((1.25, 1))

    if output_name is not None:
        plt.savefig(output_name, bbox_inches='tight')
    else:
        plt.show()


def plot_nuclei_area(gdf, area_cut_off, output_name=None):
    """Plots histograms of nuclei areas with a cutoff for filtering."""
    fig, axs = plt.subplots(1, 2, figsize=(15, 4))
    axs[0].hist(gdf['area'], bins=50, edgecolor='black')
    axs[0].set_title('Nuclei Area')

    axs[1].hist(gdf[gdf['area'] < area_cut_off]['area'], bins=50, edgecolor='black')
    axs[1].set_title(f'Nuclei Area Filtered: {area_cut_off}')

    plt.tight_layout()

    if output_name:
        # Ensure the directory exists
        directory = os.path.dirname(output_name)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        plt.savefig(output_name, bbox_inches='tight')
        print(f"Plot saved as: {output_name}")
    else:
        plt.show()

    plt.close(fig)


def total_umi(adata_, cut_off, output_name=None):
    """Plots the distribution of total UMI counts with an option to save the plot."""
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    
    # Boxplot for all total counts
    axs[0].boxplot(adata_.obs["total_counts"], vert=False, patch_artist=True, boxprops=dict(facecolor='skyblue'))
    axs[0].set_title('Total Counts')

    # Boxplot for total counts greater than the cutoff
    axs[1].boxplot(adata_.obs["total_counts"][adata_.obs["total_counts"] > cut_off], vert=False,
                   patch_artist=True, boxprops=dict(facecolor='skyblue'))
    axs[1].set_title(f'Total Counts > {cut_off}')

    # Hide y-axis
    for ax in axs:
        ax.get_yaxis().set_visible(False)

    plt.tight_layout()

    # Save or show the plot
    if output_name:
        # Ensure the directory exists
        directory = os.path.dirname(output_name)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        plt.savefig(output_name, bbox_inches='tight')
        print(f"Plot saved as: {output_name}")
    else:
        plt.show()

    # Free resources
    plt.close(fig)


