```
Python code for an improved bathtub model considering hydrological connectivity and attenuation. 
Date: 2025/02/11

Wang, Y., Ye, Y., Nicholls R., et al. (2025) Development policy affects coastal risk in China more than sea-level rise. [Unpublished manuscript].
```

from osgeo import gdal
from queue import Queue
from scipy.ndimage import convolve
import numpy as np


def read_img(filename, is_convert_nan=True, is_verbose=False):
    dataset = gdal.Open(filename)
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_proj = dataset.GetProjection()
    im_geotrans = dataset.GetGeoTransform()
    no_data = dataset.GetRasterBand(1).GetNoDataValue()
    im_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, im_width, im_height)
    if is_convert_nan:
        im_data[im_data==no_data] = np.nan
    if is_verbose:
        return im_data, im_proj, im_geotrans
    else:
        return im_data


def write_img(file_path, im_proj, im_geotrans, im_data, dtype=None, nodata=None):
    if dtype is None:
        if 'int8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32
    else:
        datatype = dtype

    if len(im_data.shape) == 2:
        im_bands, (im_height, im_width) = 1, im_data.shape
    else:
        im_bands, im_height, im_width = im_data.shape 

    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(file_path, im_width, im_height, im_bands, datatype, options=['COMPRESS=LZW'])
    dataset.SetGeoTransform(im_geotrans)
    dataset.SetProjection(im_proj)

    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)
        if nodata:
            dataset.GetRasterBand(1).SetNoDataValue(nodata)
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i+1).WriteArray(im_data[i])
            if nodata:
                dataset.GetRasterBand(i+1).SetNoDataValue(nodata)
    del dataset


def inf2nan(x):
    x[np.isinf(x)] = np.nan
    return x


def nan2neginf(x):
    x[np.isnan(x)] = -np.inf
    return x


def convolve_sealand_edge(mask):
    window = [[-1,-1,-1],
              [-1, 8,-1],
              [-1,-1,-1]]
    mask = mask - 1
    edges = np.where(convolve(mask, window, mode='constant') > 1)
    return np.column_stack(edges).tolist()


def get_land_mask(dem, slr, mask):
    g1 = (mask == 1)
    g2 = (dem <= np.nanmax(slr))
    return g1 & g2


def initialize_queue(all_mask):
    ini_list = convolve_sealand_edge(all_mask)
    q = Queue()
    for idx in ini_list:
        q.put(idx)
    return q


def fast_attenuated_bathtub(dem, slr, all_mask, atte_factor=0.02):
    land_mask = get_land_mask(dem, slr, all_mask)
    update_slr = slr.copy()
    update_slr = nan2neginf(update_slr)
    queue = initialize_queue(all_mask)
    while not queue.empty():
        x, y = queue.get()
        x8 = np.array([x] * 8)
        y8 = np.array([y] * 8)
        dx = np.array([-1, 1, 0, 0, -1, -1, 1, 1])
        dy = np.array([0, 0, -1, 1, -1, 1, -1, 1])
        nx, ny = x8 + dx, y8 + dy
        valid_indices = np.logical_and(np.logical_and(nx >= 0, nx < dem.shape[0]), 
                                       np.logical_and(ny >= 0, ny < dem.shape[1]))
        nx, ny = nx[valid_indices], ny[valid_indices]
        land_mask_nxny, dem_nxny, update_slr_x_y = land_mask[nx, ny], dem[nx, ny], update_slr[x, y]
        new_slr = update_slr_x_y - atte_factor
        update_indices = np.logical_and(land_mask_nxny, dem_nxny < update_slr_x_y - atte_factor)
        update_indices = np.logical_and(update_indices, new_slr > update_slr[nx, ny])
        update_slr[nx[update_indices], ny[update_indices]] = new_slr
        [queue.put(i) for i in list(zip(nx[update_indices], ny[update_indices])) if i is not None]
    
    flood_depth = np.where(all_mask == 1, update_slr, np.nan)
    flood_depth = inf2nan(flood_depth) - dem
    flood_map = flood_depth.copy()
    flood_map[~np.isnan(flood_map)] = 1
    flood_map[np.isnan(flood_map)] = 255
    flood_map = flood_map.astype(int)
    return flood_depth, flood_map