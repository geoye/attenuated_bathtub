Code for an improved geometric inundation model considering hydrological connectivity and attenuation (Wang et al., 2025).

## Introduction

The standard ‘bathtub’ model is a passive and static method. It computes flood depths by taking the difference between the water level and ground elevation. It assumes that the water surface is uniform, just like water rising in a bathtub. To prevent depressions which are not connected to the ocean being incorrectly considered as flooded, bathtub models that consider hydrological connectivity are widely used to constrain flood extent.

Nonetheless, the flood depth and extent can still be overestimated as the relevant hydrodynamic processes that constrain flow during a flood event are ignored, such as bed friction, surface roughness, and water flow direction. In this way, the attenuation factor was introduced to improve the understanding on hydraulic characteristics of coastal flooding. Attenuation refers to the reduction of inundation depth compared to a bathtub, and the attenuation factor is defined as the reduction of water level (in m) per unit of distance traveled (per pixel). The attenuation factor was specified as a uniform value across all pixels.

In this model, ocean water propagates from the current ocean and inundates the cells whose elevation is below a specific value subtracting the attenuation from extreme sea level rise, and keeps iterating until no cells can be inundated. The deterministic eight (D8) neighborhood was used to constraint connectivity. In each iteration, the algorithm computes the inundation depth from the ocean-land border, and updates the ocean-land border through finding their neighboring cells. When multiple cells share the same target cell to be inundated, its inundation depth will be determined as the highest value among all possible candidates. This process will keep iterating until no more cells could be flooded.

## Dependencies
Dependencies:`numpy`, `gdal`, `scipy`
## Example usage
```python
from atte_bathtub import fast_atte_bathtub, nan2neginf, read_img, write_img
from osgeo import gdal
import os
import numpy as np
import time

if __name__ == '__main__':
    for atte_factor in [0, 0.01, 0.1]:
        out_root = "./output"
        if not os.path.exists(out_root):
            os.makedirs(out_root)
        out_path = f"{out_root}/flood_depth_esl_high_end_2030_{str(atte_factor).replace('.', 'p')}.tif"
        if not os.path.exists(out_path):
            mask = read_img("./input/sea_land_mask.tif")
            dem, proj, geotrans = read_img("./input/dem_sub2030.tif", is_verbose=True)
            slr_data = read_img("./input/ts_slr_high_end_2030.tif")
            slr_data = nan2neginf(slr_data)
            fd = fast_atte_bathtub(dem, slr_data, mask, atte_factor=atte_factor)
            write_img(out_path, proj, geotrans, fd)
            print(f"Finish: {out_path}")
            del fd, dem, proj, geotrans
```

Where:

- `dem` is digital elevation model data of the study area.
- `slr_data` is the interpolated sea level rise data **on the land-ocean border** for the coastal zone.
- `mask` is a binary raster indicating whether a pixel is seawater or land.
- `atte_factor` is the attenuation factor (unit: m/pixel)
- The output file in the `out_path` represents the inundated water depth (in meter).

## License
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
