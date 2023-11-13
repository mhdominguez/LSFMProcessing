# BDV H5toTIF
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# this script converts an entire H5 dataset to TIF, excluding dataset.xml
# usage example 1, 8-bit processing on target directory: python3 dataset_folder_export_all_h5_to_TIFs_scikit-image.py /path/to/target/dataset 8
# usage example 2, 16-bit processing on current working directory: python3 dataset_folder_export_all_h5_to_TIFs_scikit-image.py 16


import sys
import os
import h5py
import numpy as np
#from PIL import Image
#from tifffile import imsave
from skimage.io._plugins import tifffile_plugin as tp
import re

# handle in arguments
outbits = "16" #default
if len(sys.argv) > 1:
	if os.path.exists(os.path.dirname(sys.argv[1])):
		os.chdir(sys.argv[1])
	elif sys.argv[1].startswith("8") or ( len(sys.argv) > 2 and sys.argv[2].startswith("8") ):
		outbits = "8"

for root, dirs, files in os.walk("."):
	for filename in files:
		if re.search("\.h5",filename): #( filename.endswith("00.h5") or filename.endswith("01.h5") ):
			file_to_modify = filename
			print( "Exporting TIFs from file", file_to_modify)
			
			f = h5py.File(file_to_modify, "r")
			
			group_names = f.keys()

			for this_group_name in list(group_names):
				if not this_group_name.startswith('t0'):
					continue  # Skip to the next iteration if the group name does not start with 't0'
				print( "Working on group:", this_group_name)
				subgroup_names = f["/"+this_group_name+"/"].keys()
				for this_name in list(subgroup_names):
					print( " subgroup:", this_name)

					tif_filename = this_group_name + "_" + this_name + ".tif"

					# Check if the TIFF file already exists
					if os.path.exists(tif_filename):
						print(f"Skipping, TIFF file already exists: {tif_filename}")
						continue  # Skip to the next iteration

					image = np.array(f["/"+this_group_name+"/"+this_name+"/0/cells"]).astype("uint"+outbits)
					tp.imsave(tif_filename, image, compression='zlib', imagej=True)

			f.close()
