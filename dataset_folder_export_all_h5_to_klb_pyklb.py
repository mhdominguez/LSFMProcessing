# BDV H5toKLB
# 2019-2023 Martin H. Dominguez
# Gladstone Institutes

# this script converts an entire H5 dataset to KLB, excluding dataset.xml
# usage example 1, 8-bit processing on target directory: python3 dataset_folder_export_all_h5_to_klb_pyklb.py /path/to/target/dataset 8
# usage example 2, 16-bit processing on current working directory: python3 dataset_folder_export_all_h5_to_klb_pyklb.py 16


import sys
import os
import h5py
import numpy as np
import pyklb #requires pyklb (https://github.com/bhoeckendorf/pyklb)
import re

# handle in arguments
outbits = "16" #default
if len(sys.argv) > 1:
	if os.path.exists(os.path.dirname(sys.argv[1])):
		os.chdir(sys.argv[1])
		if len(sys.argv) > 2 and sys.argv[2].startswith("8"):
			outbits = "8"
	elif sys.argv[1].startswith("8"):
		outbits = "8"

for root, dirs, files in os.walk("."):
	for filename in files:
		if re.search("\.h5",filename): #( filename.endswith("00.h5") or filename.endswith("01.h5") ):
			file_to_modify = filename
			print( "Exporting KLBs from file", file_to_modify)
			
			f = h5py.File(file_to_modify, "r")
			
			group_names = f.keys()

			for this_group_name in list(group_names):
				if not this_group_name.startswith('t0'):
					continue  # Skip to the next iteration if the group name does not start with 't0'
				print( "Working on group:", this_group_name)
				subgroup_names = f["/"+this_group_name+"/"].keys()
				for this_name in list(subgroup_names):
					print( " subgroup:", this_name)

					klb_filename = this_group_name + "_" + this_name + ".klb"

					# Check if the KLB file already exists
					if os.path.exists(klb_filename):
						print(f"Skipping, KLB file already exists: {klb_filename}")
						continue  # Skip to the next iteration

					pyklb.writefull(np.array(f["/"+this_group_name+"/"+this_name+"/0/cells"]).astype("uint"+outbits),klb_filename)

			f.close()
