import matplotlib.pyplot as plt
from shapely.geometry import Polygon

def crop_image(img, bbox):
    """Crops the image to the given bounding box."""
    if bbox:
        return img[bbox[1]:bbox[3], bbox[0]:bbox[2]]
    return img

def filter_geodataframe_by_bbox(gdf, bbox):
    """Filters a GeoDataFrame to only include geometries within the bounding box."""
    if bbox:
        bbox_polygon = Polygon([(bbox[0], bbox[1]), (bbox[2], bbox[1]), (bbox[2], bbox[3]), (bbox[0], bbox[3])])
        return gdf[gdf['geometry'].intersects(bbox_polygon)]
    return gdf

def plot_cropped_image(ax, img, title):
    """Plots a cropped image."""
    ax.imshow(img, cmap='gray', origin='lower')
    ax.set_title(title)
    ax.axis('off')

