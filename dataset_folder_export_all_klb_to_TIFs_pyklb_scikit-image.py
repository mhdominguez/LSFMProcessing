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
			print( "Exporting TIF from file", file_to_modify)
			
			
			image = pyklb.readfull(filename)	
			shape = image.shape
			print(  shape)
			tp.imsave(filename_root + ".tif", image.reshape(shape[0],1,shape[1],shape[2]), compress=6, resolution=(1./0.380490284561, 1./0.380490284561), imagej=True, metadata={'spacing': 1.52196113824,'unit': 'microns'} )
