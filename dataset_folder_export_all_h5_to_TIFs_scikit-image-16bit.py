# BDV H5toTIF, 16bit
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# this script converts an entire 16-bit H5 dataset to TIF, excluding dataset.xml
# usage: python3 dataset_folder_export_all_h5_to_TIFs_scikit-image-16bit.py



import sys
import os
import h5py
import numpy as np
#from PIL import Image
#from tifffile import imsave
from skimage.io._plugins import tifffile_plugin as tp
import re

for root, dirs, files in os.walk("."):
	for filename in files:
		if re.search("\d\.h5",filename): #( filename.endswith("00.h5") or filename.endswith("01.h5") ):
			#print(filename)
			#print "Please note that renaming rules are coded in this python script; please modify script to change rules for renaming if desired!"
			file_to_modify = filename
			print( "Exporting TIFs from file", file_to_modify)
			
			f = h5py.File(file_to_modify, "r")
			
			group_names = f.keys()
			last_group_name = list(group_names)[len(group_names)-1]
			print( "Working on group:", last_group_name)
			last_group = f["/"+last_group_name+"/"]
			
			subgroup_names = last_group.keys()
			#print "Subgroups: ", subgroup_names
			
			for this_name in list(subgroup_names):
				print( " subgroup:", this_name)
				image = np.array(f["/"+last_group_name+"/"+this_name+"/0/cells"]).astype("uint16")
				tp.imsave(last_group_name + "_" + this_name + ".tif", image)
				#fi.write_multipage(image, last_group_name + "_" + this_name + ".tif")
				
			f.close()
