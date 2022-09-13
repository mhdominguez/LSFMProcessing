import sys
import os
import h5py
import numpy as np
import pyklb #requires pyklb (https://github.com/bhoeckendorf/pyklb); remember to install Cython and run setup.py as root, then pip install as root; get libklb.so and libklb_static.a from /tgmm-paper/build/keller-lab-block-filetype/src and copy to /usr/lib
from skimage.io._plugins import tifffile_plugin as tp
import re

for root, dirs, files in os.walk("."):
	for filename in files:
		if re.search("\d\.klb",filename): #( filename.endswith("00.h5") or filename.endswith("01.h5") ):
			filename_root, file_extension = os.path.splitext(filename)
			
			file_to_modify = filename
			print( "Exporting TIFs from file", file_to_modify)
			

			image = pyklb.readfull(filename)	
			tp.imsave(filename_root + ".tif", image, compress=6)
