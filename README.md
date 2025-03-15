# visium_hd_processing
Package installation
-CPU:
    - Windows: 
      conda env create -f stardist-cpu-env-windows.yml (Tested)
    - Linux: 
      conda env create -f stardist-cpu-env-linux.yml   (Tested)
-GPU:
    - Windows: 
      conda env create -f stardist-gpu-env-windows.yml (Tested)
    - Linux: 
      conda env create -f stardist-gpu-env-linux.yml   (Working on this section)

Activate environment: 
conda activate stardist-cpu-env

To have global access, in the currect directory of the git repository (visium_hd_processing) install: 
pip install -e path_to/visium_hd_functions

Example:
pip install -e .

Install kernel for linux jupyter: 
python -m ipykernel install --user --name=stardist-cpu-env --display-name "Python (stardist-cpu-env)"

Project Directory Hierarchy
visium_hd_processing/
├── README.md
├── check_gpu_tf.py              # Script to check TensorFlow GPU availability
├── pipeline_visium_hd.ipynb     # Main pipeline notebook for Visium HD processing
├── stardist-cpu-env-linux.yml   # Conda environment file for Linux (CPU)
├── stardist-cpu-env-windows.yml # Conda environment file for Windows (CPU)
├── stardist-gpu-env-linux.yml   # Conda environment file for Linux (GPU)
├── stardist-gpu-env-windows.yml # Conda environment file for Windows (GPU)
├── visium_hd_functions/         # Core functions package
│   ├── __init__.py              # Makes this directory a package
│   ├── clustering.py            # Clustering-related functions
│   │   └── perform_clustering(adata, resolution=0.35)
│   ├── data_processing.py       # Functions for preprocessing and QC
│   │   ├── calculate_qc_metrics(adata)
│   │   └── filter_and_normalize(adata, area_mask, count_mask)
│   ├── image_processing.py      # Image processing functions
│   │   ├── load_image(file_path)
│   │   ├── normalize_image(img, min_percentile=5, max_percentile=95)
│   │   ├── load_stardist_model(model_dir='default')
│   │   └── predict_nuclei(model, img, block_size=4096, prob_thresh=0.01, nms_thresh=0.001, min_overlap=128, context=128, normalizer=None, n_tiles=(4, 4, 1))
│   ├── plotting.py              # Plotting functions
│   │   ├── plot_mask_and_save_image(title, gdf, img, cmap, output_name=None, bbox=None)
│   │   ├── plot_gene_and_save_image(title, gdf, gene, img, adata, bbox=None, output_name=None)
│   │   ├── plot_clusters_and_save_image(title, gdf, img, adata, bbox=None, color_by_obs=None, output_name=None, color_list=None)
│   │   ├── plot_nuclei_area(gdf, area_cut_off)
│   │   └── plot_total_umi(adata_, cut_off)
│   ├── plotting_utils.py        # Helper functions for plotting
│   │   ├── crop_image(img, bbox)
│   │   ├── filter_geodataframe_by_bbox(gdf, bbox)
│   │   └── plot_cropped_image(ax, img, title)
│   └── setup.py                 # setup.py for the package
├── website_pipeline/            # Example notebooks or tutorials
│   └── segmentation_nuclei_tutorial2.ipynb # Example tutorial for segmentation
└── pretrained/                  # Pre-trained models
    └── 2D_versatile_he/         # Example pre-trained model directory

