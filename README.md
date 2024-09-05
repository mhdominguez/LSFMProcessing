# LSFMProcessing
A collection of methods for deconvolving, filtering, and rendering LSFM images. This collection is currently designed to take Zeiss ZEN .czi, Lavision UM2 .ome.tiff, or Luxendo .h5/.json input files, though it could be generalized to accept any [proprietary] input format that BioFormats allows.
<br><br>
## Main Fiji Methods (LSFM Processing Macros.ijm)
### Installation / Requirements:
* Fiji, using ImageJ 1.53f or greater
* BigStitcher, Bio-Formats, SiMView packages installed (checked) in Fiji Update
* [PSF generator](http://bigwww.epfl.ch/algorithms/psfgenerator/) 18.12.2017
  - download Fiji plugin and place PSF_Generator.jar in `Fiji.app/plugins` folder
* [Parallel Spectral Deconvolution](https://sourceforge.net/projects/spectraldeconv/files/spectraldeconv/1.12/parallel_spectral_deconvolution-1.12-bin.zip/download) 1.12
  - download plugin binary zip file and copy all *.jar files to `Fiji.app/plugins` folder
  - delete file jars/jtransforms-2.4.jar before using deconvolution
* [Parallel Iterative Deconvolution (if using iterative deconvolution)](https://sourceforge.net/projects/iterativedeconv/files/iterativedeconv/1.12/parallel_iterative_deconvolution-1.12-bin.zip/download) 1.12
  - download plugin binary zip file and copy all *.jar files to `Fiji.app/plugins` folder
  - delete file jars/jtransforms-2.4.jar before using deconvolution
* h5py and scikit-image installed on system for h5/tif functions
  - on Ubuntu, can use `sudo pip3 install h5py scikit-image` at console to install
* pyklb installed on system for klb functions and TGMM
  - on Ubuntu, can use console to install...
```
git clone https://github.com/bhoeckendorf/pyklb
cd pyklb
sudo pip3 install cython==0.29.35
python3 setup.py build
sudo python3 setup.py install
```

Use `Plugins->Macros->Install` to add `LSFM Processing Macros.ijm` to the macrons menu in Fiji.<br>
Then, use `Plugins->Macros->0. Change LSFM Processing Settings` to adjust user settings, including filter parameters and deconvolution block size.<br>
<br><br>
## BigStitcher dataset copy view transformations (dataset_folder_copy_view_transformations.pl)
 this script opens the BigStitcher dataset.xml in the working directory, and for each timepoint assuming two view setups per angle, finds the view setup with the most number of transforms for each angle (i.e. the angle used to register that view), and copies them to the other view so that view is now registered.<br>
 usage: `perl dataset_folder_copy_view_transformations.pl`
<br><br>
## BigDataViewer dataset format conversion scripts
 convert an entire H5 dataset to TIF, excluding dataset.xml:<br>
  usage: `python3 dataset_folder_export_all_h5_to_tif_scikit-image.py [optional path/to/dataset] 16`, 16 can be replaced with 8 for 8-bit<br><br>
 convert an entire H5 dataset to KLB, excluding dataset.xml:<br>
  usage: `python3 dataset_folder_export_all_h5_to_klb_pyklb.py [optional path/to/dataset] 16`, 16 can be replaced with 8 for 8-bit<br><br>
