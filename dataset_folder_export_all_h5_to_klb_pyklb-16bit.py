# BDV H5toKLB, 8bit
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# this script converts an entire 8-bit H5 dataset to KLB, excluding dataset.xml (needs to be done manually with BigDataViewer KLB import)
# usage: python3 dataset_folder_export_all_h5_to_klb_pyklb-8bit.py

import sys
import os
import h5py
import numpy as np
import pyklb #requires pyklb (https://github.com/bhoeckendorf/pyklb); remember to install Cython and run setup.py as root, then pip install as root; get libklb.so and libklb_static.a from /tgmm-paper/build/keller-lab-block-filetype/src and copy to /usr/lib
import re

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
				pyklb.writefull(np.array(f["/"+last_group_name+"/"+this_name+"/0/cells"]).astype("uint16"),last_group_name + "_" + this_name + ".klb")
				
			f.close()
