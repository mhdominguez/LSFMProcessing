# BDV H5toKLB
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# this script converts an entire H5 dataset to KLB, excluding dataset.xml (needs to be done manually with BigDataViewer KLB import)
# usage example 1, 8-bit processing on target directory: python3 dataset_folder_export_all_h5_to_klb_pyklb.py /path/to/target/dataset 8
# usage example 2, 16-bit processing on current working directory: python3 dataset_folder_export_all_h5_to_klb_pyklb.py 16

import sys
import os
import h5py
import numpy as np
import pyklb #requires pyklb (https://github.com/bhoeckendorf/pyklb); remember to install Cython and run setup.py as root, then pip install as root; get libklb.so and libklb_static.a from /tgmm-paper/build/keller-lab-block-filetype/src and copy to /usr/lib
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
		if re.search("\d\.h5",filename): #( filename.endswith("00.h5") or filename.endswith("01.h5") ):
		 	#print(filename)
			#print "Please note that renaming rules are coded in this python script; please modify script to change rules for renaming if desired!"
			file_to_modify = filename
			print( "Exporting KLBs from file", file_to_modify)
			
			f = h5py.File(file_to_modify, "r")
			
			group_names = f.keys()
			last_group_name = list(group_names)[len(group_names)-1]
			print( "Working on group:", last_group_name)
			last_group = f["/"+last_group_name+"/"]
			
			subgroup_names = last_group.keys()
			#print "Subgroups: ", subgroup_names
			
			for this_name in list(subgroup_names):
				print( " subgroup:", this_name)
				pyklb.writefull(np.array(f["/"+last_group_name+"/"+this_name+"/0/cells"]).astype("uint"+outbits),last_group_name + "_" + this_name + ".klb")
				
			f.close()
