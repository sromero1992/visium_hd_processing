from tifffile import imread
from csbdeep.utils import normalize
from stardist.models import StarDist2D

def load_image(file_path):
    return imread(file_path)

def normalize_image(img, min_percentile=5, max_percentile=95):
    return normalize(img, min_percentile, max_percentile)

def load_stardist_model():
    return StarDist2D.from_pretrained('2D_versatile_he')

def predict_nuclei(model, img, block_size=4096, prob_thresh=0.01, nms_thresh=0.001, min_overlap=128, context=128, normalizer=None, n_tiles=(4, 4, 1)):
    return model.predict_instances_big(img, axes='YXC', block_size=block_size, prob_thresh=prob_thresh,
                                        nms_thresh=nms_thresh, min_overlap=min_overlap, context=context,normalizer=None, n_tiles=n_tiles)

