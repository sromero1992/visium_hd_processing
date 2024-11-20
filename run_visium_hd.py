# Import necessary modules
import scanpy as sc
from visium_hd_functions.image_processing import load_image, normalize_image, load_stardist_model, predict_nuclei
from visium_hd_functions.data_processing import create_geodataframe, merge_tissue_positions, filter_spatial_overlap, sum_gene_counts_binned 
from visium_hd_functions.plotting import plot_mask_and_save_image, plot_gene_and_save_image, plot_clusters_and_save_image
from matplotlib.colors import ListedColormap

# Set up configurations
#Save results
dir_base = '/mnt/SCDC/Optimus/selim_working_dir/visium_hd_mouse_colon/'

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

print("Loading pre-trained StarDist model...")
# Try to download and read from default directory, if fails, use pre-trained model in repository
model_dir = 'default'
model = load_stardist_model(model_dir)

print("Predicting nuclei...") # This is crucial
labels, polys = predict_nuclei(model, img, block_size = 4096, prob_thresh = 0.01, nms_thresh = 0.001)

# Create GeoDataFrame using Polygon geometries
print("Creating GeoDataFrame...")
gdf = create_geodataframe(polys['coord'])

# Display summary of the GeoDataFrame
gdf.head()

# Plot the nuclei segmentation (ROI)
# bbox=(x min,y min,x max,y max)
bbox=(1840,3296,2672,4160)

# Define a single color cmap
cmap=ListedColormap(['grey'])

# Create Plot
plot_mask_and_save_image(title = "Region of Interest 1", gdf = gdf, bbox = bbox,
                         cmap = cmap,img = img, output_name = dir_base+"image_mask.ROI1.tif")
