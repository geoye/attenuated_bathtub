Code for an improved bathtub model considering hydrological connectivity and attenuation (Wang et al., 2025).

## Introduction

The standard ‘bathtub’ model is a passive and static method. It computes flood depths by taking the difference between the water level and ground elevation. It assumes that the water surface is uniform, just like water rising in a bathtub. To prevent depressions which are not connected to the ocean being incorrectly considered as flooded, bathtub models that consider hydrological connectivity are widely used to constrain flood extent.

Nonetheless, the flood depth and extent can still be overestimated as the relevant hydrodynamic processes that constrain flow during a flood event are ignored, such as bed friction, surface roughness, and water flow direction. In this way, the attenuation factor was introduced to improve the understanding on hydraulic characteristics of coastal flooding. Attenuation refers to the reduction of inundation depth compared to a bathtub, and the attenuation factor is defined as the reduction of water level (in m) per unit of distance traveled (per pixel). The attenuation factor was specified as a uniform value across all pixels.

In this model, ocean water propagates from the current ocean and inundates the cells whose elevation is below a specific value subtracting the attenuation from extreme sea level rise, and keeps iterating until no cells can be inundated. The deterministic eight (D8) neighborhood was used to constraint connectivity. In each iteration, the algorithm computes the inundation depth from the ocean-land border, and updates the ocean-land border through finding their neighboring cells. When multiple cells share the same target cell to be inundated, its inundation depth will be determined as the highest value among all possible candidates. This process will keep iterating until no more cells could be flooded.

## Dependencies
Dependencies:`numpy`, `gdal`, `scipy`
## Example usage
```python
from attenuated_bathtub import fast_attenuated_bathtub, read_img, write_img
import os


if __name__ == '__main__':
    for atte_factor in [0, 0.01, 0.1]:
        fd_path = f"./output/flood_depth_esl_g3_ssp585_2030_{str(atte_factor).replace('.', 'p')}.tif"
	if not os.path.exists(fd_path):
	    mask = read_img("./sea_land_mask_g3.tif")
            dem, proj, geotrans = read_img("./dem_sub2030_g3.tif", is_verbose=True)
            slr_data = read_img("./esl_g3_ssp585_2030.tif")
            fd, fm = fast_attenuated_bathtub(dem, slr_data, mask, atte_factor=atte_factor)
            write_img(fd_path, proj, geotrans, fd)
            # `fm` represents the binary flooding output. 
            # One can also write it to disk using `write_img(fm_path, proj, geotrans, fm)`
```

Where:

- `dem` is digital elevation model data of the research area.
- `slr_data` is the interpolated sea level rise data for the coastal zone. **Note that only pixels on the land-ocean border will be used to model SLR!**
- `land_sea_mask` is a binary raster indicating whether a pixel is seawater or land.
- `atte_factor` is the attenuation factor (unit: m/pixel)
- The output file `fd_path` represents the inundated water depth (in meter).

## License
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)

## To cite the article

> Wang, Y., Ye, Y., Nicholls R., et al. (2025) Development policy affects coastal risk in China more than sea-level rise. [Unpublished manuscript].
