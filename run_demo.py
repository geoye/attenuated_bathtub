```
Python code for an improved geometric inundation model considering hydrological connectivity and attenuation. 
Date: 2025/06/06
```

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