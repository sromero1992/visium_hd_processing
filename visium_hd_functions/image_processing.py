#from tifffile import imread
import os
from skimage.io import imread
from csbdeep.utils import normalize
from stardist.models import StarDist2D

# Function to load an image
def load_image(file_path):
    """
    Load an image from a given file path.
    :param file_path: Path to the image file.
    :return: Loaded image.
    """
    return imread(file_path)

# Function to normalize an image
def normalize_image(img, min_percentile=5, max_percentile=95):
    """
    Normalize an image based on percentile values.
    :param img: Input image to normalize.
    :param min_percentile: Minimum percentile for normalization (default: 5%).
    :param max_percentile: Maximum percentile for normalization (default: 95%).
    :return: Normalized image.
    """
    return normalize(img, min_percentile, max_percentile)

# Function to load the StarDist model from a specific directory
def load_stardist_model(model_dir='default'):
    """
    Load the StarDist2D model from a custom path.
    :param model_dir: Path to the pretrained model directory or 'default' to download.
    :return: Loaded StarDist2D model.
    """
    # Try loading the model from the default pre-trained location
    if model_dir == 'default':
        try:
            print("Trying to load pre-trained model from the default source...")
            model = StarDist2D.from_pretrained('2D_versatile_he')
            print("Pre-trained model loaded from default location.")
        except Exception as e:
            print(f"Error downloading from the default location: {e}")
            print(f"Please, try to use pre-trained model included in the visium_hd_processing repo")
            model = None
            #print(f"Falling back to local model at {model_dir}.")
            #if not os.path.exists(model_dir):
            #    raise ValueError(f"Model directory {model_dir} not found.")
            #model = StarDist2D.from_pretrained(model_dir)
    else:
        # If a specific directory is provided, try to load from there
        if not os.path.exists(model_dir):
            raise ValueError(f"Model directory {model_dir} not found.")
        print(f"Loading pre-trained model from provided directory: {model_dir}")
        model = StarDist2D(None, name='2D_versatile_he',basedir=model_dir)

    return model

# Function to predict nuclei using a trained model
def predict_nuclei(model, img, block_size=4096, prob_thresh=0.01, nms_thresh=0.001, 
                   min_overlap=128, context=128, normalizer=None, n_tiles=(4, 4, 1)):
    """
    Predict nuclei instances using a StarDist2D model.
    :param model: Trained StarDist2D model.
    :param img: Input image to predict on.
    :param block_size: Size of the image blocks to process (default: 4096).
    :param prob_thresh: Probability threshold for predicting nuclei (default: 0.01).
    :param nms_thresh: Non-Maximum Suppression threshold (default: 0.001).
    :param min_overlap: Minimum overlap for combining predictions (default: 128).
    :param context: Context padding for tile predictions (default: 128).
    :param normalizer: Normalizer function (default: None).
    :param n_tiles: Number of tiles in X, Y, and Z for processing large images (default: (4, 4, 1)).
    :return: Predicted nuclei instances.
    """
    return model.predict_instances_big(img, axes='YXC', block_size=block_size, prob_thresh=prob_thresh,
                                       nms_thresh=nms_thresh, min_overlap=min_overlap, context=context, 
                                       normalizer=normalizer, n_tiles=n_tiles)

