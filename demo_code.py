#!/usr/bin/env python
# -*- coding:utf-8 -*-

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