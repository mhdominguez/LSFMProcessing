# LSFMProcessing
A collection of methods for deconvolving, filtering, and rendering LSFM images. This collection is currently designed to take Zeiss ZEN .czi or Lavision UM2 .ome.tiff input files, though it could be generalized to accept any [proprietary] input format that BioFormats allows.
<br><br>
## Main Fiji Methods (LSFM Processing Macros.ijm)
### Installation / Requirements:
* Fiji, using ImageJ 1.53f or greater
* BigStitcher, Bio-Formats, SiMView packages installed (checked) in Fiji Update
* [PSF generator](http://bigwww.epfl.ch/algorithms/psfgenerator/) 18.12.2017
  - download Fiji plugin and place PSF_Generator.jar in plugins folder
* [Parallel Spectral Deconvolution ](https://sites.google.com/site/piotrwendykier/software/deconvolution/parallelspectraldeconvolution) 1.12
  - download plugin binary zip file and copy all *.jar files to plugins folder
  - delete file jars/jtransforms-2.4.jar before using deconvolution
* pyklb and h5py installed on system for h5/klb functions
  - on Ubuntu, can use `sudo pip3 install git+https://github.com/bhoeckendorf/pyklb.git@skbuild h5py` at console to install

Use `Plugins->Macros->Install` to add `LSFM Processing Macros.ijm` to the macrons menu in Fiji.<br>
Then, use `Plugins->Macros->0. Change LSFM Processing Settings` to adjust user settings, including filter parameters and deconvolution block size.<br>
<br><br>
## BigStitcher dataset copy view transformations (dataset_folder_copy_view_transformations.pl)
 this script opens the BigStitcher dataset.xml in the working directory, and for each timepoint assuming two view setups per angle, finds the view setup with the most number of transforms for each angle (i.e. the angle used to register that view), and copies them to the other view so that view is now registered.<br>
 usage: `perl dataset_folder_copy_view_transformations.pl`
<br><br>
## BigDataViewer dataset format conversion scripts
 convert an entire 8-bit H5 dataset to TIF, excluding dataset.xml:<br>
  usage: `python3 dataset_folder_export_all_h5_to_TIFs_scikit-image.py [optional path/to/dataset] 16`, 16 can be replaced with 8 for 8-bit<br><br>
 convert an entire 8-bit H5 dataset to KLB, excluding dataset.xml (needs to be done manually with BigDataViewer KLB import):<br>
  usage: `python3 dataset_folder_export_all_h5_to_klb_pyklb.py [optional path/to/dataset] 8`, 8 can be replaced with 16 for 16-bit<br><br>
