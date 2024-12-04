# Import necessary modules
import scanpy as sc
from visium_hd_functions.image_processing import load_image, normalize_image, load_stardist_model, predict_nuclei, enhance_contrast
from visium_hd_functions.data_processing import create_geodataframe, merge_tissue_positions, filter_spatial_overlap, sum_gene_counts_binned, create_geodf_coords
from visium_hd_functions.plotting import plot_mask_and_save_image, plot_gene_and_save_image, plot_clusters_and_save_image, plot_nuclei_area, total_umi
from visium_hd_functions.clustering import calculate_qc_metrics, filter_and_normalize, perform_clustering
from matplotlib.colors import ListedColormap

# Set up configurations
#Save results
dir_base = '/mnt/SCDC/Optimus/selim_working_dir/visium_hd_mouse_colon/'
samp_name = 'AfGC_513_'

# Open files
dir_files = '/mnt/SCDC/Bumblebee/2024_AfGC_AdipoRon_Visium_HD_Training_Data/'
dir_sub = 'High-Res_Stiched_Images/btf_files/'
img_file = dir_files + dir_sub + 'AfGC#513.btf'

dir_sub2 = 'mapped-data/AfGC-513/outs/binned_outputs/square_002um/' 
dir_hd5 = dir_files + dir_sub2 
h5_file = dir_hd5  + 'filtered_feature_bc_matrix.h5'

dir_tissue = dir_hd5 + 'spatial/'
tissue_position_file = dir_tissue+'tissue_positions.parquet'

# Load and process image
print("Loading and processing image...")
img = load_image(img_file)
img = normalize_image(img, min_percentile = 5, max_percentile = 95)

# Ensure the normalized image is within the range [0, 1]
img = img.clip(0, 1)

# Enhance image contrast
print("Enhancing constrast and brightness...")
img = enhance_contrast(img, output_dir=dir_base, img_name=samp_name, clip_limit=0.02, brightness_factor=2.5)

print("Loading pre-trained StarDist model...")
# Try to download and read from default directory, if fails, use pre-trained model in repository
model_dir = 'default'
model = load_stardist_model(model_dir)

print("Predicting nuclei...") # This is crucial
#labels, polys = predict_nuclei(model, img, block_size = 4096, prob_thresh = 0.01, nms_thresh = 0.001)
labels, polys = predict_nuclei(model, img, block_size = 4096, prob_thresh = 0.005, nms_thresh = 0.0005, min_overlap = 64, context = 128)

# Create GeoDataFrame using Polygon geometries
print("Creating GeoDataFrame...")
gdf = create_geodataframe(polys['coord'])

# Display summary of the GeoDataFrame
gdf.head()

# Plot the nuclei segmentation (ROI)
# bbox=(x min,y min,x max,y max)
#bbox=(5662,1708,6766,2974) # This ROI is blurry
bbox=(8604,4776,9708,5796) 

# Define a single color cmap
cmap=ListedColormap(['grey'])

# Create Plot
plot_mask_and_save_image(title = "Region of Interest 1", gdf = gdf, bbox = bbox,
                         cmap = cmap,img = img, output_name = dir_base+samp_name+"image_mask.ROI1.tif")


# Process AnnData from Visium HD data
print("Reading AnnData...")
print(h5_file)
adata = sc.read_10x_h5(h5_file)

# Display summary of AnnData object
print(adata)


# Adding tissue positions to the meta data + df with tissue positions and barcodes
print("Merging tissue positions...")
adata, df_tissue_positions = merge_tissue_positions(adata, tissue_position_file)

# Create a GeoDataFrame from the DataFrame of coordinates 
gdf_coordinates = create_geodf_coords(df_tissue_positions)

# Spatial join check of coordinates that are in a cell nucleus
filtered_adata = filter_spatial_overlap(adata, gdf, gdf_coordinates)

# Add counts in bins
grouped_filtered_adata = sum_gene_counts_binned(filtered_adata)

# Store the area of each nucleus in the GeoDataframe
gdf['area'] = gdf['geometry'].area

# Perform quality control and normalization
print("Calculating QC metrics...")
calculate_qc_metrics(grouped_filtered_adata)

# Plot the nuclei area distribution before and after filtering
plot_nuclei_area(gdf=gdf,area_cut_off=500, output_name = dir_base+samp_name+"nuclei_area.png")

# Plot total UMI distribution
total_umi(grouped_filtered_adata, 100,  output_name = dir_base+samp_name+"total_umi.png")


# Create a mask based on the 'id' column for values present in 'gdf' with 'area' less than 500
mask_area = grouped_filtered_adata.obs['id'].isin(gdf[gdf['area'] < 500].id)
# Create a mask based on the 'total_counts' column for values greater than 100
mask_count = grouped_filtered_adata.obs['total_counts'] > 100

# Apply both masks to the original AnnData to create a new filtered AnnData object
grouped_filtered_adata = filter_and_normalize(grouped_filtered_adata, mask_area, mask_count)

# Perform clustering
print("Performing clustering...")
grouped_filtered_adata = perform_clustering(grouped_filtered_adata)

# Plot and save the clustering results
plot_clusters_and_save_image(title="Region of interest 1", gdf=gdf, img=img, 
                             adata=grouped_filtered_adata, bbox=bbox, color_by_obs='clusters', 
                             output_name=dir_base+"image_clustering.ROI1.tiff")


print("Plotting gene expression...")
# Plot Lyz1 gene expression
plot_gene_and_save_image(title="Region of interest 1", gdf=gdf, gene='Epcam', img=img, adata=grouped_filtered_adata, bbox= bbox,output_name=dir_base+"image_Epcam.ROI1.tiff")
plot_gene_and_save_image(title="Region of interest 1", gdf=gdf, gene='Muc2', img=img, adata=grouped_filtered_adata, bbox= bbox,output_name=dir_base+"image_Muc2.ROI1.tiff")


# Save results
print("Saving processed data...")
adata.write_h5ad(dir_base + 'processed_data.h5ad')
print("Pipeline complete!")

print("Pipeline completed successfully.")

