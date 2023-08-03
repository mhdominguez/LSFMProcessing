// Zeiss Light Sheet Processing Macro Scripts
// 2019-2023 Martin H. Dominguez
// Gladstone Institutes

// requires: 
//     Fiji, using ImageJ 1.53f or greater
//     BigStitcher, Bio-Formats, SiMView packages installed (checked) in Fiji Update
//     PSF generator 18.12.2017 (http://bigwww.epfl.ch/algorithms/psfgenerator/)
//       - download Fiji plugin and place PSF_Generator.jar in plugins folder
//     Parallel Spectral Deconvolution 1.12 (https://sites.google.com/site/piotrwendykier/software/deconvolution/parallelspectraldeconvolution)
//       - download plugin binary zip file and copy all *.jar files to plugins folder
//       - delete file jars/jtransforms-2.4.jar before using deconvolution
//     Parallel Iterative Deconvolution 1.12 (if using iterative deconvolution, https://sites.google.com/site/piotrwendykier/software/deconvolution/paralleliterativedeconvolution)
//       - download plugin binary zip file and copy all *.jar files to plugins folder
//       - delete file jars/jtransforms-2.4.jar before using deconvolution
//     pyklb and h5py installed separately for h5/klb functions
//       - on Ubuntu, can use "sudo pip3 install git+https://github.com/bhoeckendorf/pyklb.git@skbuild h5py" at console to install

function show_instructions () {
	Dialog.create("LSFM processing install info/instructions...");
	//Dialog.setInsets(0,0,0);
	Dialog.addMessage("LSFMProcessing requires:",16,"Black");
    Dialog.addMessage("   Fiji, using ImageJ 1.53f or greater",14,"Black");
    Dialog.addMessage("   BigStitcher, Bio-Formats, SiMView packages installed (checked) in Fiji Update",14,"Black");
    Dialog.addMessage("   PSF generator 18.12.2017",14,"Black");
    Dialog.addMessage("       - download Fiji plugin (http://bigwww.epfl.ch/algorithms/psfgenerator/)",12,"Black");
	Dialog.addMessage("       - place PSF_Generator.jar in Fiji plugins folder",12,"Black");
    Dialog.addMessage("   Parallel Spectral Deconvolution 1.12",14,"Black");
    Dialog.addMessage("       - download plugin binary zip file (https://sites.google.com/site/piotrwendykier/software/deconvolution/parallelspectraldeconvolution)",12,"Black");
    Dialog.addMessage("       - unzip and copy all *.jar files to Fiji plugins folder",12,"Black");
    Dialog.addMessage("       - delete file jtransforms-2.4.jar in Fiji jars folder before using deconvolution",12,"Black");
    Dialog.addMessage("   Parallel Iterative Deconvolution 1.12 (if using iterative deconvolution)",14,"Black");
    Dialog.addMessage("       - download plugin binary zip file (https://sites.google.com/site/piotrwendykier/software/deconvolution/paralleliterativedeconvolution)",12,"Black");
    Dialog.addMessage("       - unzip and copy all *.jar files to Fiji plugins folder",12,"Black");
    Dialog.addMessage("       - delete file jtransforms-2.4.jar in Fiji jars folder before using deconvolution",12,"Black");
    Dialog.addMessage("   pyklb and h5py installed on system, for h5/klb functions",14,"Black");
    Dialog.addMessage("       - on Ubuntu, can use \"sudo pip3 install git+https://github.com/bhoeckendorf/pyklb.git@skbuild h5py\" at console to install",12,"Black");
	Dialog.show();
}

var file_sep = File.separator();
var max_slice_depth = 480; //heap allocation in MB should be at least 200 times max slice depth
var rolling_ball_radius = 200;
var fixed_rolling_ball_radius = 500;
var fixed_blob_removal = false;
var fixed_precipitate_removal = true;
var convert_to_8_bit = true;
var generate_mips_when_filtering = true;
var deconvolution_regression_parameter = 50;
var deconvolution_iterations = 0;
var deconvolution_subtract_camera_noise = true;
var deconvolve_yes = true;
var detection_NA_penalty = 10;
var path_dataset_folder_export_all_h5_to_klb_pyklb = "dataset_folder_export_all_h5_to_klb_pyklb.py";

var default_tissue_refractive_index = 1.4;
var default_immersion_refractive_index = 1.333;

function removeLeadingWhitespace(str) {
    while (startsWith(str, " ")) {
        str = substring(str, 1);  // Remove the first character
    }
    return str;
}

function listFilesRecursive(directory) {
	list = getFileList(directory);
	fileList = newArray(list.length);

	for (i = 0; i < list.length; i++) {
		filePath = directory + file_sep + list[i];
		
		if (File.isDirectory(filePath)) {
			// Recursively call the function for subdirectories
			subdirectoryFiles = listFilesRecursive(filePath);
			
			// Concatenate the subdirectory files with the current directory files
			fileList = concat(fileList, subdirectoryFiles);
		} else {
			// Add the file path to the array
			fileList[i] = filePath;
		}
	}
	
	// Remove null entries from the array (if any)
	fileList = Array.trim(fileList);
	
	// Return the array of file names
	return fileList;
}


function LS_PSF_generator ( NA_objective, NA_lightsheet, RI_immersion, RI_sample, lambda_lightsheet, lambda_detection, res_XY, res_Z, Z_radius_lightsheet, WD_return, CS_return ) {
	output_file = "#PSFGenerator\nPSF-shortname=GL\n"; // start us off
	output_file = output_file + "ResLateral=" + d2s(res_XY,2) + "\nResAxial=" + d2s(res_Z,2) + "\n"; //voxel units
    
	//penalize NA on objective to oversize PSF for better handling of axial spherical aberration -- this is a user-specified parameter
	if ( detection_NA_penalty > 0 &&  detection_NA_penalty < 75 ) {
		NA_objective -= ( NA_objective * detection_NA_penalty / 100 );
	}
	
    //if no lambda_detection, use lambda_lightsheet
    if ( lambda_detection == NaN || lambda_detection <= 0 ) {
		lambda_detection = lambda_lightsheet;
    }
    
    //if no lambda_lightsheet, use lambda_detection
    if ( lambda_lightsheet == NaN || lambda_lightsheet <= 0 ) {
		lambda_lightsheet = lambda_detection;
    }    
	
	//now, calculate smallest PSF size
	//first, estimate PSF Z dimension as at least 4 waists; the beam waist is related to lightsheet NA and wavelength: w0 = 2n*(lambda_lightsheet)/pi*NA
	
	//start with w(0), beam waist at center
	beam_waist = RI_immersion * lambda_lightsheet / ( NA_lightsheet * PI );
	//print( "  ..beam waist at center (nm): " + d2s(beam_waist,2) + " due to " + d2s(RI_immersion,3) + " " + d2s(NA_lightsheet,3) + " " + d2s(lambda_lightsheet,2) + " " + d2s(PI,4) + " " + d2s(Z_radius_lightsheet,4) + " " + d2s(res_XY,4) + "\n" );

	//Gaussian light sheet calculations:
	//https://en.wikipedia.org/wiki/Rayleigh_length
	//https://www.photometrics.com/wp-content/uploads/2019/10/Light-Sheet-Microscopy-App-Note.pdf
	//https://core.ac.uk/download/pdf/157810694.pdf
	//https://opg.optica.org/aop/fulltext.cfm?uri=aop-10-1-111&id=381035
	
	//Rayleigh distance is not really used post-hoc, but may be helpful for manually setting up light sheet during imaging
	//Zr = ( PI * beam_waist ) / lambda_lightsheet;
	//print( "  ..Rayleigh radius (nm): " + d2s(Zr,2) );

	// determine light sheet thickness as an X-invariant single value, so choose thickness w(z) at some reasonable X value based on Z_radius_lightsheet and res_XY
	//now calculate with w(z), the beam waist at distance Z_radius_lightsheet (in pixels) from the center of the focus, where lightsheet is narrowest at width beam_waist
	beam_waist *= beam_waist; // square beam waist for this
	beam_width_z_squared = lambda_lightsheet * Z_radius_lightsheet * res_XY / ( PI * beam_waist );
	beam_width_z_squared *= beam_width_z_squared;
	beam_width_z_squared = beam_waist * ( 1 + beam_width_z_squared ); // this is actually w(z) squared, but we will use this instead of w(z)	
	//print( "w(z): " + d2s(sqrt(beam_width_z_squared),4) + " due to " + d2s(res_XY,2) + " and "  + d2s(Z_radius_lightsheet,1) + " and " + d2s(res_Z,4) + "\n" );
	
	if ( beam_waist > beam_width_z_squared ) { //take as lightsheet thickness the greater of width at Rayleigh distance, or calculated distance based on image width input into this function (Z_radius_lightsheet)
		//print ( "using beam_waists squared = " + d2s(beam_waist,8) + ", beam_width_z_squared = " + d2s(beam_width_z_squared,8)  );
		beam_width_z_squared = beam_waist;
	} else {
		//print ( "using beam_width_z_squared squared = " + d2s(beam_width_z_squared,8) + ", beam_waists = " + d2s(beam_waist,8) );
	}
	
	//estimate the depth of the PSF in pixels, smallest to get the job done
	Z_dim = round( 3 * sqrt(beam_width_z_squared) / res_Z );
	//print( "  NAobj = " + d2s(NA_objective,4) + ", resZ = " + d2s(res_Z,4) + ", RIimm = " + RI_immersion );
	//print ( "Zdim: " + d2s(Z_dim,4) + " based on " + d2s(beam_width_z_squared,6) + "\n" );
	Z_dim = return_odd_pixel_dimension( Z_dim );
	if ( isNaN(Z_dim) ) {
		Z_dim = 31;
	}
	
	//next estimate smallest XY size possible to get the job done
	XY_dim = return_odd_pixel_dimension( 6500 / (NA_objective * res_XY * RI_immersion) ); //30 = k * / 0.9(NA) * 250(XY) * 1(RIi)
	if ( isNaN(XY_dim) ) {
		XY_dim = 31 ;
	}
	
	output_file = output_file + "NY=" + d2s(XY_dim,0) + "\nNX=" + d2s(XY_dim,0) + "\nNZ=" + d2s(Z_dim,0) + "\nType=16-bits\n"; //dimensions and color depth of output PSF
	
	//finish working on the file -- this allows PSF generator to create detection PSF
	output_file = output_file + "NA=" + d2s(NA_objective,6) + "\nLUT=Grays\nLambda=" + d2s(lambda_detection,4) + "\nScale=Linear\n";
	output_file = output_file + "psf-GL-NI=" + d2s(RI_immersion,4) + "\npsf-GL-NS=" + d2s(RI_sample,4) + "\npsf-GL-accuracy=Good\npsf-GL-ZPos=" + d2s(CS_return,1) + "\npsf-GL-TI=" + d2s(WD_return,1) + "\n";	
	
	//now, work on illumination PSF
	rand_num = d2s(floor(random() * 1000000 ),0);
	LS_ill_PSF = "Ill-PSF-" + rand_num;
	//name_newimage = "Ill-PSF-" + rand_num;
	//midpoint_Z_dim = floor(Z_dim/2) + (WD_return/res_Z); // what is the center of the beam?
	midpoint_Z_dim = floor(Z_dim/2);
	newImage("Ill-PSF-" + rand_num, "16-bit black", 1, 1, Z_dim );
	
	//create gaussian gradient to finish illumination PSF
	selectWindow( LS_ill_PSF );
	for (a=0; a<Z_dim; a++ ) {
		//exp(-2r^2)/w0^2
		setZCoordinate(a);
		pos_Z = ( a-midpoint_Z_dim ) * res_Z;
		setPixel( 0, 0, 255 * exp( -2 * pos_Z * pos_Z  / beam_width_z_squared ) );
		//print( "setting pixel at " + a + " with value " + d2s(getPixel(0,0),0) + " from midpoint_Z_dim " + d2s(midpoint_Z_dim,3) + " and pos_Z " + d2s(pos_Z,3) + "\n" );
	}	
	run("Size...", "width=" + XY_dim + " height=" + XY_dim + " depth=" + Z_dim + " constrain average interpolation=Bilinear");
	
	//now write config file for  -- this template will allow PSF Generator to create detection PSF
	filehandle = File.open( "Det-PSF-" + rand_num + ".txt" );
	print( filehandle, output_file );
	File.close( filehandle);
	
	//run PSF generator
	//setBatchMode(false);
	//run_script = "importClass(Packages.PSFGenerator);\nnew_imp = PSFGenerator.computeImagePlus( \"" + "Det-PSF-" + rand_num + ".txt" + "\" );\nnew_imp.show();\nnew_imp.getID();\n";
	run_script = "importClass(Packages.PSFGenerator);\nPSFGenerator.computeImagePlus( \"" + "Det-PSF-" + rand_num + ".txt" + "\" ).show();\n";
	//eval_result = eval( "js", run_script ); 
	eval( "js", run_script );  //throw away result
	//print ( "Eval result: " + eval_result );
	
	LS_det_PSF = getInfo("window.title");
	//imageCalculator("multiply stack create 32-bit", parseInt(eval_result), LS_ill_PSF );
	imageCalculator("multiply stack create 32-bit", LS_det_PSF, LS_ill_PSF );
	//selectWindow(LS_det_PSF);
	//exit();
	name_newimage = getInfo("window.title");
  	run("Properties...", "unit=micron pixel_width=" + d2s(res_XY/1000,10) + " pixel_height=" + d2s(res_XY/1000,10) + " voxel_depth=" + d2s(res_Z/1000,10) );
	//exit();
	//selectWindow( parseInt(eval_result) ); close();
	selectWindow(LS_det_PSF); close(); 
	selectWindow(LS_ill_PSF); close();
	delete_result = File.delete( "Det-PSF-" + rand_num + ".txt" ); //throw away result
	
	return name_newimage;
}

function Zeiss_WDforObjective(Olightsheet) {
	//returns working distance and coverslip thickness (i.e. depth to sample) of the objective
	return_array = newArray(3); // WD, coverslip, nominal RIimmersion
	Array.fill(return_array, 0 );
	if (startsWith( Olightsheet, "W " ) ) {
		//probably water dipping objective
		if (matches(Olightsheet, ".*20\(x\|X\).*" ) ) {
			return_array[0] = 2400;	
		} else if (matches(Olightsheet, ".*40\(x\|X\).*" ) ) {
			return_array[0] = 2500;
		} else if (matches(Olightsheet, ".*5\(x\|X\).*" ) ) {
			return_array[0] = 5100;			
		} else if (matches(Olightsheet, ".*10\(x\|X\).*" ) ) {
			return_array[0] = 3700;			
		} else if (matches(Olightsheet, ".*63\(x\|X\).*" ) ) {
			return_array[0] = 2100;
		}
		return_array[2] = 1.333;
	} else if (startsWith(Olightsheet, "Clr " ) || startsWith(Olightsheet, "LSFM cl" )) {
		//probably clearing objective
		if (matches(Olightsheet, ".*20\(x\|X\).*" ) ) {
			return_array[0] = 5600;	
		}	
		return_array[2] = 1.45;
	} else if (startsWith(Olightsheet, "EC " )) {
		//probably clearing objective
		if (matches(Olightsheet, ".*5\(x\|X\).*" ) ) {
			return_array[0] = 18500;
			return_array[1] = 170;
			return_array[2] = 1.333;
		}
	} else {
		print ("The working distance and coverslip are not defined for " + Olightsheet + " \n");
	}
	
	return return_array;
}

function UM2_WDforObjective(Olightsheet) {
	//returns working distance and coverslip thickness (i.e. depth to sample) of the objective
	return_array = newArray(3); // WD, coverslip (always 0), NA
	Array.fill(return_array, 0 );
	
	if (startsWith( Olightsheet, "mi plan " ) ) { // LaVision UM2
		//probably water dipping objective
		if (matches(Olightsheet, ".*1\.1\(x\|X\).*" ) ) {
			if (matches(Olightsheet, ".*dc57.*" ) ) {
				return_array[0] = 17600;
			} else if (matches(Olightsheet, ".*dc40.*" ) ) {
				return_array[0] = 16000;
			}
			return_array[2] = 0.1;	
		} else if (matches(Olightsheet, ".*4\(x\|X\).*" ) ) {
			if (matches(Olightsheet, ".*dc57.*" ) ) {
				return_array[0] = 16000;
			} else if (matches(Olightsheet, ".*dc49.*" ) ) {
				return_array[0] = 16000;
			} else if (matches(Olightsheet, ".*dc33.*" ) ) {
				return_array[0] = 16000;
			}
			return_array[2] = 0.35;	
		} else if (matches(Olightsheet, ".*12\(x\|X\).*" ) ) {
			if (matches(Olightsheet, ".*dc57.*" ) ) {
				return_array[0] = 10000;
			} else if (matches(Olightsheet, ".*dc49.*" ) ) {
				return_array[0] = 10900;
			} else if (matches(Olightsheet, ".*dc33.*" ) ) {
				return_array[0] = 8500;
			}	
			return_array[2] = 0.53;	
		}
	} else {
		print ("The working distance and coverslip are not defined for " + Olightsheet + " \n");
	}
	//print ("LaVision: " + d2s(return_array[0],0) + " " + d2s(return_array[1],0) );
	return return_array;
}

function return_parse_xml_line_for_value ( inString ) {
	//print( "checking for value in " + inString );
	inString = toLowerCase(inString);
	
	start_index = indexOf(inString, "value=\"") +7;
	//print( "   ...at " + d2s(start_index+6,0) );
	
	if ( start_index < 0 ) {
		return ""
	}
	subline = substring(inString, start_index );
	//print( "   ... " + subline );
	stop_index = indexOf(subline, "\"" );
	
	if ( stop_index < 0 ) {
		return ""
	}
	
	subline = substring(subline, 0, stop_index );
	//print( "   ... " + subline );
	
	return subline ;
}

function return_odd_pixel_dimension( dimension ) {

	//print( "return_odd_pixel_dimension: " + d2s(dimension,2) );
	if ( dimension <= 0 ) {
		print ("Problem with negative or zero dimension in return_odd_pixel_dimension\n");
		return NaN;
	//} else if ( dimension < 4 ) {
	//	return 3;
	} else if ( dimension < 8 ) {
		return 7;
	} else if ( dimension < 16 ) {
		return 15;
	} else if ( dimension < 24 ) {
		return 23;		
	} else if ( dimension < 32 ) {
		return 31;
	} else if ( dimension < 40 ) {
		return 39;		
	} else if ( dimension < 48 ) {
		return 47;		
	} else {
		return 63; //don't go larger than this
	}
}
	
function main_ijm_deconvolve_large_stack(folder_in_batch,directory,outputDirectory) {	
	
	processList = newArray(0);
	unique_PSF_filenames = newArray(0);
	unique_PSF_parameters = newArray(0);	
	run("Bio-Formats Macro Extensions");
	type = "";
	
	if ( folder_in_batch ) {
		if ( directory == "" || !File.exists(directory) ) {
			//Ask user to choose the input and output directories
			directory = getDirectory("Choose input directory");
		}
		fileList = getFileList(directory);

		// identify type / format of dataset
		for (i=0; i<fileList.length; i++) {
			if (endsWith(fileList[i], ".czi")) { // Zeiss Z.1 CZI format
				type = ".czi";
				break;
			} else if (endsWith(fileList[i], "Z0000.ome.tif")) { // OME format, possibly UM2
				type = "Z0000.ome.tif";
				break;
			}
		}
		for (i=0; i<fileList.length; i++) {
			if (endsWith(fileList[i], type )) {
				//now, remove .czi from end of filename
				
				filepaddedname = substring(fileList[i], 0, indexOf(fileList[i], type) );
				
				/*name_ext = split(fileList[i],".");
				filepaddedname = "";
				for (n=0; n<name_ext.length-1; n++) {
					filepaddedname += name_ext[n];
				}*/
				
				//print( "    " + name_ext[0] + " " + name_ext[0] + " " + name_ext[0] + " " + name_ext[0] + " " + name_ext[0] + " " +);
	
				//get ready to rename on save
				if ( endsWith(filepaddedname,".") ) { //additional trailing periods should be removed
					while( endsWith(filepaddedname,".") ) {
						filepaddedname = substring(filepaddedname,0,filepaddedname.length);					
					}
				}
				
				processList = Array.concat( processList, filepaddedname + "///" + fileList[i] ); //new name will consist of padded old name so we can sort correctly by timepoint
			}
		}
	} else {
		path = File.openDialog("Select a File");
  		directory = File.getParent(path);
  		name = File.getName(path);		
		//name_ext = split(name,".");

		if (endsWith(name, ".czi")) { // Zeiss Z.1 CZI format
			type = ".czi";
		} else if (endsWith(name, "Z0000.ome.tif")) { // OME format, possibly UM2
			type = "Z0000.ome.tif";
		}
		
		//now, remove .czi from end of filename
		/*filepaddedname = "";
		for (n=0; n<name_ext.length-1; n++) {
			filepaddedname += name_ext[n];
		}*/
		
		filepaddedname = substring(name, 0, indexOf(name, type) );
		
		//get ready to rename on save
		if ( endsWith(filepaddedname,".") ) { //additional trailing periods should be removed
			while( endsWith(filepaddedname,".") ) {
				filepaddedname = substring(filepaddedname,0,filepaddedname.length);					
			}
		}
		processList = Array.concat( processList, filepaddedname + "///" + name ); //new name will consist of padded old name so we can sort correctly by timepoint		
	}
	
	if ( outputDirectory == "" || !File.exists(outputDirectory) ) {
		outputDirectory = getDirectory("Choose output directory");
	}	
	
	Array.sort(processList); //sort files by their intended name, not true name
	//RIlightsheet = 0;
	//Olightsheet = "";
	max_channels = 0;
    max_time = 0;
	this_time = 0;
	
	//TODO:uncomment
	setBatchMode(true);

	for (i=0; i<processList.length; i++) {
		name_ext = split(processList[i],"(///)");
		file = directory + file_sep + name_ext[1];
		print( "Preparing " + directory + file_sep + name_ext[1] + " for deconvolution..." );
		
		Ext.setId(file);
		Ext.getSeriesCount(nPositions);

		RIlightsheet = 0;
		Olightsheet = "";
		WD_CS_return = newArray(3);
		Array.fill(WD_CS_return,0);

		//get view angles, again metadata permeates all series and these values can be extracted with any series open
		string_angles = newArray(nPositions);
		angles = newArray(nPositions);

		Array.fill(string_angles,NaN); // default angle is NaN
		Array.fill(angles,NaN); // default angle is NaN
		
		if ( type == ".czi" ) {
			//get refractive index
			//although this will retrieve metadata only from active series, refractive index is not different for each channel or view but instead is the same for everything
			Ext.getMetadataValue("Information|Image|RefractiveIndex #1", RIlightsheet);
			if ( RIlightsheet == 0 ) {
				Ext.getMetadataValue("Information|Image|RefractiveIndex", RIlightsheet);
			}
			//print ( "Metadata: \n" + metadata_value  + "\n" ); exit();

			//get objective name
			Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|Objective #1", Olightsheet);
			if ( Olightsheet == "" || Olightsheet == 0 ) {
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|Objective", Olightsheet);
			}
			
			WD_CS_return = Zeiss_WDforObjective( Olightsheet );
			if ( WD_CS_return[2] > RIlightsheet ) {
				RIlightsheet = WD_CS_return[2];
			}

			//grab initial angles from file
			if ( nPositions < 1) {
				print ( "Fewer than one view/position present in file " + directory + file_sep + name_ext[1] + "!\n" );
				continue;
			} else {
				//fill angles array
				for(a=0; a<nPositions; a++) {
					Ext.getMetadataValue("Information|Image|V|View|Offset #" + d2s((a+1),0), string_angles[a]);
				}

				//fill actual angles array with numbers
				for(a=0; a<nPositions; a++) {
					angles[a] = parseFloat( string_angles[a]);
				}
			}

			//modify angles to try to put front at zero degrees
			if ( nPositions == 1 ) {
				//set only available angle to zero degrees
				angles[0] = 0;
			} else { // 2 or more angles

				//first, figure out whether we should use angles from -180 to +180, or 0 to 360 -- to do this, if the difference between the largest and smallest angles is greater than 180, use -180 to 180 rather than 0 to 360
				//fill sort angles array
				sort_angles = Array.copy(angles);
				Array.sort(sort_angles);

				if ( sort_angles[sort_angles.length-1] - sort_angles[0] > 180 ) {
					//modify angles to use -180 to 180
					for ( a=0; a<angles.length; a++ ) {
						if ( angles[a] > 180 ) {
							angles[a] -= 360;
						}
					}

					//redo sort_angles array now since angles have been modified
					sort_angles = Array.copy(angles);
					Array.sort(sort_angles);
				}

				if( nPositions == 2 ) {
					abs_diff = sort_angles[1] - sort_angles[0];

					if ( abs_diff < 180.2 && abs_diff > 179.8 ) { // views are opposing, so assume one is front and one is back
						//set view 1 to angle 0 and view 2 to angle 180, then be done
						angles[0] = 0;
						angles[1] = 180;
					} else if (angles[1]>=angles[0])  { // views are not opposing and second view has greater angle than first view, so assume that front is actually midpoint of the two views -- and set that to 0 degrees
						angles[0] = 360 - (abs_diff /2);
						angles[1] = abs_diff/2;
					} else { // views are not opposing and first view has greater angle than second view, so assume that front is actually midpoint of the two views -- and set that to 0 degrees
						angles[1] = 360 - (abs_diff /2);
						angles[0] = abs_diff/2;
					}
				} else { // 3 or more positions, so zero the middle view if odd or view right before middle view (assume views done left-front-right-back or right-front-left-back)
					//zero the second view
					middle_view = floor((angles.length / 2) - 0.5);
	                		diff_angle = angles[middle_view];
					for ( a=0; a<angles.length; a++ ) {
						angles[a] -= diff_angle;

						if ( angles[a] < -180 ) {
							angles[a] += 360;
						} else if ( angles[a] >= 360 ) {
							angles[a] -= 360;
						}
					}
				}

				//okay, now correct angles for 0 to 360, regardless of which scheme we used above
				for ( a=0; a<angles.length; a++ ) {
					if ( angles[a] >= 360 ) {
						angles[a] -= 360;
					} else if ( angles[a] < 0 ) {
						angles[a] += 360;
					}
				}
			}
		}

		if ( nPositions == 1 ) {
			angles[0] = 0; // first angle is 0
		}
		
		//get lightsheet data for each channel
		num_channels = 0;
	    num_time = 0;
		max_time = 0;
		max_width = 0;
		for(a=0; a<nPositions; a++) {
			dim_width = 0;
				
			Ext.setSeries(a);
			Ext.getSizeC(num_channels);
	        Ext.getSizeT(num_time);
			Ext.getSizeX(dim_width);
			
			if ( num_channels > max_channels ) {
				max_channels = num_channels;
			}
			if ( num_time > max_time ) {
				max_time = num_time;
			}
			if ( dim_width > max_width ) {
				max_width = dim_width;
			}
		}
		
		NAlightsheets = newArray(max_channels);
		NAdetections = newArray(max_channels);
		WLlightsheets = newArray(max_channels);
		WLdetections = newArray(max_channels);
		Rayleigh_dist_calc_waist = max_width/5; // default waist is 20% of total x dimension of image -- really should consider this per channel per position (not just per channel), but most likely positions will all have same dim_width

		Array.fill(NAlightsheets,0); // default NA is 0
		Array.fill(NAdetections,0); // default NA is 0
		Array.fill(WLlightsheets,0); // default WL is 0
		Array.fill(WLdetections,0); // default WL is 0

		if ( type == ".czi" ) {
			for(a=0; a<max_channels; a++) {
				Ext.getMetadataValue("Information|Image|Channel|NALightSheet #" + d2s((a+1),0), NAlightsheets[a]);
				Ext.getMetadataValue("Information|Image|Channel|NADetection #" + d2s((a+1),0), NAdetections[a]);
				Ext.getMetadataValue("Information|Image|Channel|IlluminationWavelength|SinglePeak #" + d2s((a+1),0), WLlightsheets[a]);
				Ext.getMetadataValue("Information|Image|Channel|DetectionWavelength|SinglePeak #" + d2s((a+1),0), WLdetections[a]);
				
				print( "  ..examining " + name_ext[0] + " ...channel " + d2s(a,0) + " parameters: NAill " + d2s(NAlightsheets[a],6) + ", NAdet " + d2s(NAdetections[a],6) + ", WLill " + d2s(WLlightsheets[a],2) + ", WLdet " + d2s(WLdetections[a],2) );
			}
		}

		for (t=0; t<max_time; t++ ) {		
			for(a=0; a<nPositions; a++) {
				if ( nPositions == 1 ) {
					run("Bio-Formats", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT " + " t_begin=" + d2s(t+1,0) + " t_end=" + d2s(t+1,0) + " t_step=1" );	
				} else {
					run("Bio-Formats", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT series_"+ d2s(a+1,0) + " t_begin_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_end_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_step_" + d2s(a+1,0) + "=1" );
				}
				
				//Get data from this particular stack if UM2
				if ( type == "Z0000.ome.tif" ) {
					open( file );
					infoData=getMetadata("Info");
					close();
					/*
			 		* fname="UltraII LRI" nTy="3" nId="524292" Value="1.550000"
			 		* fname="UltraII ObjectiveName" nTy="5" nId="524292" Value="MI Plan 1.1x DC57"
			 		* fname="UltraII ExBeamWaist" nTy="1" nId="262148" Value="3.300000"/>
			 		* label="UltraII Sheet thickness FWHM micron" fname="UltraII SheetThickness" nTy="3" nId="524292" Value="3.889772"/>
			 		* label="UltraII Na Value" fname="UltraII NAValue" nTy="3" nId="524292" Value="0.156093"/>
			 		* fname="UltraII Wavelength0" nTy="2" nId="524292" Value="561"/>
				 	* fname="UltraII EmWavelength0" nTy="2" nId="524292" Value="620"/>
				 	* <ExcitationWL ExcitationWL="488"/>
			 		* <EmissionWL EmissionWL="525"/>
			 		* <Pixels ID="Pixels:25492249-EF27-4C5C-9DC5-D4343416BCD4" DimensionOrder="XYZTC" PixelType="uint16" BigEndian="false" SizeX="1772" SizeY="1742" SizeZ="1117" SizeT="1" SizeC="1" PhysicalSizeX="3.559694" PhysicalSizeY="3.559693" PhysicalSizeZ="1.700000">
				 	*/
					infoString = split(infoData, "\n");
					beam_waist = "0";
					for ( l=0; l<infoString.length; l++ ) {
						if (indexOf(infoString[l], "fname=\"UltraII LRI\"") >= 0) {
							RIlightsheet = return_parse_xml_line_for_value(infoString[l]);
						} else if (indexOf(infoString[l], "fname=\"UltraII ObjectiveName\"") >= 0) {
							Olightsheet = return_parse_xml_line_for_value(infoString[l]);
						} else if (indexOf(infoString[l], "fname=\"UltraII CurWLExcitation\"") >= 0) {
							WLlightsheets[0] = return_parse_xml_line_for_value(infoString[l]);
						} else if (indexOf(infoString[l], "fname=\"UltraII CurWLEmission\"") >= 0) {
							WLdetections[0] = return_parse_xml_line_for_value(infoString[l]);							
						} else if (indexOf(infoString[l], "fname=\"UltraII ExBeamWaist\"") >= 0) {
							beam_waist = return_parse_xml_line_for_value(infoString[l]);
						} 
					}
		
					WD_CS_return = UM2_WDforObjective( Olightsheet );
					NAdetections[0] = d2s(WD_CS_return[2],1);

					//calculate lightsheet NA					
					// beam_waist = RI_immersion * lambda_lightsheet / ( NA_lightsheet * PI );
					//NA = (1.55) * .488 / ( PI * 3.889772 )
					//	( PI * 3.889772 ) = 12.22						
					NAlightsheets[0] = d2s( (parseFloat(RIlightsheet) * parseFloat(WLlightsheets[0])/1000) / ( PI * parseFloat(beam_waist)),10);
					
				}
				
				//Get name of opened stack	
				master_title = getTitle();
				rand_num = d2s( floor(random() * 1000000 ), 0 );
				rename( "Image" + rand_num );
				title = getTitle();
				getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
				getVoxelSize(vox_width, vox_height, vox_depth, vox_unit);
				voxel_definition_output = "unit=" + vox_unit + " pixel_width=" + d2s(vox_width,10) + " pixel_height=" + d2s(vox_height,10) + " voxel_depth=" + d2s(vox_depth,10); //will use later when we output the file
				
				//make all units nanometers
				if (vox_unit=="nm" || vox_unit=="nanometers" || vox_unit=="nanometer") {
					//do nothing
				} else if (vox_unit==getInfo("micrometer.abbreviation") || vox_unit=="um" || vox_unit=="microns" || vox_unit=="micron") {
					vox_width * = 1000;
					vox_height * = 1000;
					vox_depth * = 1000;
				} else if (vox_unit=="mm" || vox_unit=="millimeters" || vox_unit=="millimeter"){
					vox_width * = 1000000;
					vox_height * = 1000000;
					vox_depth * = 1000000;
				} else if (vox_unit=="cm" || vox_unit=="centimeters" || vox_unit=="centimeter"){
					vox_width * = 10000000;
					vox_height * = 10000000;
					vox_depth * = 10000000;
				} else if (vox_unit=="m" || vox_unit=="meters" || vox_unit=="meter"){
					vox_width * = 1000000000;
					vox_height * = 1000000000;
					vox_depth * = 1000000000;
				} else {
					print ( "Cannot interpret voxel unit " + vox_unit + " for " + file + ", series " + d2s(a,0) + "\n" );
					continue;
				}
				
				//check height and width scales
				if ( vox_width > 1.001 * vox_height || vox_height > 1.001 * vox_width ) {
					print ( "Unexpected XY scale difference: " + d2s(vox_height,5) + " vs. " + d2s(vox_width,5) + " for " + file + ", series " + d2s(a,0) + "\n" );
					continue;
				}

				print( "  ..examining " + name_ext[0] + " ...series " + d2s(a,0) + " attributes: Xdim (nm): " + d2s(dim_width*vox_width,2) + ", Ydim (nm): " + d2s(dim_height*vox_height,2) + ", Zdim (nm): " + d2s(dim_slices*vox_depth,2) );
		
				//Split channels and record names of each new image stack
				channelList = newArray(0);
				if ( dim_channels > 1 ) {
					run("Split Channels");
					for (b=1; b<=dim_channels; b++ ) {
						channelList = Array.concat( channelList, "C" + IJ.pad(b,1) + "-" + title ); 
					}
				} else if ( dim_channels == 1 ) {
					channelList = Array.concat( channelList, title );
				} else {
					continue; //skip this series altogether	
				}
				
				if ( type == "Z0000.ome.tif" ) {
					for (b=0; b<channelList.length; b++ ) {
						if ( b > 0 ) {
							WLlightsheets[b] = WLlightsheets[0];
							WLdetections[b] = WLdetections[0];
							NAdetections[b] = NAdetections[0];
							NAlightsheets[b] = NAlightsheets[0];
						}
						print( "  ..examining " + name_ext[0] + " ...channel " + d2s(b,0) + " parameters: NAill " + d2s(NAlightsheets[b],6) + ", NAdet " + d2s(NAdetections[b],6) + ", WLill " + d2s(WLlightsheets[b],2) + ", WLdet " + d2s(WLdetections[b],2) );
					}
				}
		
				fileoutname = newArray(channelList.length);
				
				psf_id = newArray(max_channels);
				psf_depth = newArray(max_channels); // Z-depth for PSF
				Array.fill(psf_id,NaN); // default is no fill
				Array.fill(psf_depth,NaN); // default is no fill
				
				//lambdas = newArray(max_channels);
				//Array.fill(lambdas,0); //no regularization by default
				
				//figure out file out names
				for (b=0; b<channelList.length; b++ ) {
					fileoutname[b] = "LSFM__" + replace(name_ext[0],' ','_') + "__C" + IJ.pad(b,1) + "__A" + IJ.pad(d2s(round(angles[a]),0),3);
				}
				
				//do PSFs separately since they require batchmode to be off
				for (b=0; b<channelList.length; b++ ) {
					if ( NAlightsheets[b] > 0 && WLlightsheets[b] > 0 && RIlightsheet > 0 && !(isNaN(angles[a]) ) ) {
						
						this_PSF_identifier = NAdetections[b] + "//" + NAlightsheets[b] + "//" + RIlightsheet + "//" + WLlightsheets[b] + "//" + WLdetections[b] + "//" + d2s(vox_width,8) + "//" + d2s(vox_depth,8) + "//" + d2s(Rayleigh_dist_calc_waist,8) + "//" + Olightsheet;
						
						//see if we have done this PSF before
						is_match = -1;
						for (m=0;m<unique_PSF_parameters.length; m++ ) {
							if ( this_PSF_identifier == unique_PSF_parameters[m] ) {
								is_match = m;
								break;
							}
						}
						
						if ( is_match >= 0 ) {
							//copy file and be done; we've done it before
							//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif is match with " + unique_PSF_filenames[is_match] + "...\n" + this_PSF_identifier );
							//File.copy( unique_PSF_filenames[is_match], outputDirectory + fileoutname[b] + "_psf.tif" );
							open( unique_PSF_filenames[is_match] );
							
							getDimensions(_, _, _, psf_depth[b], _);
							psf_id[b] = getImageID();
							
						} else { //new PSF, need to generate
							//generate and save lightsheet  PSF
							//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif is no match...\n" + this_PSF_identifier );
							
							//TODO:uncomment below two
							setBatchMode("exit and display"); //needed for PSF generation
							setBatchMode(false);

							tissue_RI = default_tissue_refractive_index;
							if ( RIlightsheet > default_tissue_refractive_index ) {
								tissue_RI = RIlightsheet;
							}
							PSF_name = LS_PSF_generator ( parseFloat(NAdetections[b]), parseFloat(NAlightsheets[b]), parseFloat(RIlightsheet), tissue_RI, parseFloat(WLlightsheets[b]), parseFloat(WLdetections[b]), vox_width, vox_depth, Rayleigh_dist_calc_waist, WD_CS_return[0], WD_CS_return[1] ); //estimate lightsheet radius from w(0) as being 1/4th the total width of the acquired image -- this should cut it about midway
							selectWindow( PSF_name );
							
							
							run("Enhance Contrast...", "saturated=0 process_all use");
							run("16-bit");
							
							//print( "Delete result " + outputDirectory + fileoutname[b] + "_psf.tif: " + delete_result ); 
							getDimensions(_, _, _, psf_depth[b], _);
							psf_id[b] = getImageID();
							
							delete_result = File.delete( outputDirectory + file_sep + fileoutname[b] + "_psf.tif" ); //throw away result
							saveAs("Tiff", outputDirectory + file_sep + fileoutname[b] + "_psf.tif" );
							//close();
							//exit();
							
							//TODO:uncomment
							setBatchMode(true);
							
							//now store in arrays for future use
							unique_PSF_filenames = Array.concat( unique_PSF_filenames, outputDirectory + file_sep + fileoutname[b] + "_psf.tif" );
							unique_PSF_parameters = Array.concat( unique_PSF_parameters, this_PSF_identifier );
						}

					} else {
						//this_PSF_identifier = NAdetections[b] + "//" + NAlightsheets[b] + "//" + RIlightsheet + "//" + WLlightsheets[b] + "//" + WLdetections[b] + "//" + d2s(vox_width,8) + "//" + d2s(vox_depth,8) + "//" + d2s(dim_width/4,8) + "//" + Olightsheet;
						//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif cannot be saved...\n" + this_PSF_identifier );
						psf_id[b] = NaN;
						psf_depth[b] = NaN;
						
					}
				}				
				
				//now process actual images				
				//remove blank images -- from all channels together to prevent incorrect registrations later
				selectWindow(channelList[0]); //use channel 0 as primary channel
				setSlice(dim_slices);
				getStatistics( area, mean );
				if ( mean < 10 ) { //we potentially have a
					//print( "Mean is : " + d2s(mean,4) + " on slice " + d2s(dim_slices,0) );
					
					final_slice = dim_slices-1;
					
					for (ss=dim_slices-10; ss>0; ss-=10 ) { //decrement 10 to find the first slice with data
						setSlice(ss);
						getStatistics( area, mean );
						//print( " ...checking mean: " + d2s(mean,4) + " on slice " + d2s(ss,0) );
						if ( mean > 10 ) {
							//okay pinpoint exact end of stack now
							for (tt=ss+9; tt>=ss; tt-- ) { //decrement individual slices to find the first slice with data
								setSlice(tt);
								getStatistics( area, mean );
								//print( "    ...checking mean: " + d2s(mean,4) + " on slice " + d2s(tt,0) );
								if ( mean > 10 ) { //okay we found the end
									final_slice = tt;
									break;
								}
							}
							break;
						}
					}
					
					for (b=0; b<channelList.length; b++ ) {
						selectWindow(channelList[b]);
						img_id = getImageID();
						if (isNaN(psf_id[b]) || isNaN(img_id) ) {
							continue; //don't process this image
						}						
						run("Slice Remover", "first=" + d2s(final_slice+1,0) + " last=" + d2s(dim_slices,0) + " increment=1");
					}
				}
					
				//preprocess images, first by trying to remove noisy in the background
				for (b=0; b<channelList.length; b++ ) {
					selectWindow(channelList[b]);
					img_id = getImageID();
					if (isNaN(psf_id[b]) || isNaN(img_id) ) {
						continue; //don't process this image
					}					
					
					if ( deconvolution_subtract_camera_noise ) {
						//calculate minimum projection and subtract from stack (attempt to remove noise)
						run("Z Project...", "projection=[Min Intensity]");
						average_img = getTitle();
						imageCalculator("Subtract stack",channelList[b],average_img);
						selectWindow(average_img);
						close();
					}
					
					//save image to free heap space for deconvolution
					if ( channelList.length > 2 && b > 0 ) {
						selectWindow(channelList[b]);
						saveAs("Tiff", outputDirectory + file_sep + channelList[b] + ".tif" );
						close();
					}
				}
				
				
				//do deconvolution here
				for (b=0; b<channelList.length; b++ ) {
					if ( channelList.length > 2 && b > 0 ) {
						open( outputDirectory + file_sep + channelList[b] + ".tif" );
					} else {
						selectWindow(channelList[b]);
					}
					img_id = getImageID();
					//print( "Processing image " + d2s(img_id,0) + " with title " + channelList[b] );
					if (isNaN(psf_id[b]) || isNaN(img_id) ) {
						continue; //don't process this image
					}
					//do deconvolution here
					
					//try to estimate regularization parameter based on image histogram -- this should be revised
					Stack.getStatistics( area, mean, min, max, std );
					
					//arbitrarily set regularization parameter, based on 16-bit images input
					//lambda = 75 / ( std + 100); // 10/std was the original, but produced ringing artifacts -- therefore try 50/std , but in truth we should do two iterations stepping down from lambda 1000/std then 100/std
					lambda = deconvolution_regression_parameter / ( std + deconvolution_regression_parameter);
					
					//decide if stack is too large for single run with deconvolution plugin (plugin does not use block-wise)
					getDimensions(dim_width, dim_height, _, dim_slices, _);
					if ( dim_slices > max_slice_depth ) {  // if ( dim_slices * dim_width * dim_height > 1769472000 ) {
						//chunk image into 480-slice components (including PSF depth padding), then deconvolve each one and recompose
						num_blocks = floor((dim_slices/(max_slice_depth-(3*psf_depth[b])))+1); //internal blocks: one PSF padding and half PSF on each side
						
						if ( num_blocks < 2 ) { //case num_blocks == 1 -- shouldn't ever get here
							//32-bit output (with single precision calculation)
							run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n";
							if ( deconvolution_iterations > 1 ) {
								run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
							}
							//print (run_script);
							//call("edu.emory.mathcs.restoretools.spectral.ParallelSpectralDeconvolution3D.deconvolve", pathToBlurredImage, pathToPsf, pathToDeblurredImage, method, stencil, resizing, output, precision, threshold, regParam, nOfThreads, showPadded);
							if ( deconvolve_yes ) {
								eval( "js", run_script );  //throw away result
							} else {
								run("Duplicate...", "duplicate");
							}
						} else {
							img_ids = newArray(num_blocks); //create array with all image IDs of blocks
							Array.fill(img_ids,NaN);
							//print( "DEBUG: Dim slices (" +d2s(dim_slices,0)+ "), block depth (" + d2s(block_depth,0) + "), PSF depth (" + d2s(psf_depth[b],0) + ") for channel "+d2s(b,0) +", max slice depth (" + d2s(max_slice_depth,0) +")." );
							block_depth = floor((dim_slices/num_blocks)+0.5); //if odd, round up
							if ( block_depth <= psf_depth[b] ) { //declare a problem if block depth is not bigger than two PSFs
								print( "Block depth (" + d2s(block_depth,0) + ") is not greater than PSF depth (" + d2s(psf_depth[b],0) + ") for channel "+d2s(b,0) +", cannot deconvolve " +master_title+"!" );
								return;
							}
							
							//fill img_ids array with parent stack
							rename( channelList[b] + "-block0" ); //these renaming lines aren't necessary, deconvolution is by imageID, not window title
							img_ids[0] = img_id; //last and parent block maintains image title, will make substacks from this image (concatenate will look like 1,2,3,...,n,0)
							
							//take off first substack, block_depth slices (including single PSF pad)
							run("Make Substack...", "slices=1-" + d2s(block_depth+psf_depth[b],0) ); 
							rename( channelList[b] + "-block1" );
							img_ids[1] = getImageID();
							
							//go back to parent and delete non-padding portion of first block
							//selectWindow(channelList[b] + "-block0");
							selectImage(img_ids[0]);
							run("Slice Remover", "first=1 last=" + d2s(block_depth-psf_depth[b],0) + " increment=1");

							//deconvolve first child
							run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_ids[1], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n"; //print (run_script);								
							if ( deconvolution_iterations > 1 ) {
								run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_ids[1], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
							}
							if ( deconvolve_yes ) {
								eval( "js", run_script );  //throw away result
							} else {
								run("Duplicate...", "duplicate");
							}
							run("Slice Remover", "first=" + d2s(block_depth+1,0) + " last=" + d2s(block_depth+psf_depth[b],0) + " increment=1");//delete padding
							rand_num = d2s(floor(random() * 1000000 ),0);
							rename("block1-" + rand_num);
							cat_text = "  image1=block1-"+rand_num+" ";
							
							selectImage(img_ids[1]); //close undeconvolved image
							close();
								
							for(ii=2; ii<num_blocks; ii++) { //first and last block already taken care of
								//take off intermediate substack, block_depth slices (including two PSF pad)
								selectImage(img_ids[0]);
								run("Make Substack...", "slices=1-" + d2s(block_depth+(2*psf_depth[b]),0) ); //take off an internal substack, block_depth slices (plus double PSF pad)
								rename( channelList[b] + "-block" + d2s(ii,0) );
								img_ids[ii] = getImageID();
								
								//go back to parent and delete non-padding portion of intermediate block
								//selectWindow(channelList[b] + "-block0");
								selectImage(img_ids[0]);
								run("Slice Remover", "first=1 last=" + d2s(block_depth,0) + " increment=1");
									
								//deconvolve intermediate child
								run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_ids[ii], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n"; //print (run_script);								
								if ( deconvolution_iterations > 1 ) {
									run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_ids[ii], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
								}
								if ( deconvolve_yes ) {
									eval( "js", run_script );  //throw away result
								} else {
									run("Duplicate...", "duplicate");
								}
								run("Slice Remover", "first=1 last=" + d2s(psf_depth[b],0) + " increment=1");//delete padding
								run("Slice Remover", "first=" + d2s(block_depth+1,0) + " last=" + d2s(block_depth+psf_depth[b],0) + " increment=1");//delete padding
								rand_num = d2s(floor(random() * 1000000 ),0);
								rename("block"+d2s(ii,0)+"-" + rand_num);
								//cat_text = "  image1=block"+d2s(ii,0)+"-"+rand_num+" ";
								
								//concatenate intermediate children
								run("Concatenate...", cat_text + " image2=block"+d2s(ii,0)+"-"+rand_num+" image3=[-- None --]" );
								rename("block"+d2s(ii,0)+"-" + rand_num); //rename the concatenated stack
								cat_text = "  image1=block"+d2s(ii,0)+"-"+rand_num+" "; //refresh cat_text with concatenated stack
								
								//clean up
								selectImage(img_ids[ii]);
								close();
								call("java.lang.System.gc");
							}
								
							//deconvolve parent -- note this script DOES NOT deconvolve in place
							run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_ids[0], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n"; //print (run_script);								
							if ( deconvolution_iterations > 1 ) {
								run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_ids[0], 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
							}
							if ( deconvolve_yes ) {
								eval( "js", run_script );  //throw away result
							} else {
								run("Duplicate...", "duplicate");
							}
							run("Slice Remover", "first=1 last=" + d2s(psf_depth[b],0) + " increment=1");//delete padding
							rand_num = d2s(floor(random() * 1000000 ),0);
							rename("block0-" + rand_num);
							//cat_text_end = "image"+d2s(num_blocks,0)+"=block0-"+rand_num+" image"+d2s(num_blocks+1,0)+"=[-- None --]";
							//image1=block1-"+rand_num+" "; //
							cat_text = cat_text + " image2=block0-"+rand_num+" image3=[-- None --]";

							//put final stack back together
							run("Concatenate...", "  " + cat_text);

							img_id = img_ids[0];
							
							//clean up residual mess
							/*
							for(ii=1; ii<num_blocks; ii++) {
								if (isOpen(img_ids[ii]) ) {
									selectImage(img_ids[ii]);
									close();
								}
							}*/
						}
					} else if ( dim_slices < 2 * psf_depth[b] )  { //declare a problem if block depth is not bigger than two PSFs
						print( "Image depth (" + d2s(dim_slices,0) + ") is not greater than PSF depth (" + d2s(psf_depth[b],0) + ") for channel "+d2s(b,0) +", cannot deconvolve " +master_title+"!" );
						return;
					} else { //normal flow through this else statement -- images not greater than 480 Z-depth
						//32-bit output (with single precision calculation)
						run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n";
						if ( deconvolution_iterations > 1 ) {
							run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
						}
						//print (run_script);
						if ( deconvolve_yes ) {
							eval( "js", run_script );  //throw away result
						} else {
							run("Duplicate...", "duplicate");
						}
					}
					//run("Enhance Contrast...", "saturated=0 normalize process_all use");
					setMinAndMax(0, 65535); //clip the upper bounds of signal (lower bounds already clipped for nonnegativity)
					run("16-bit");
					run("Properties...", voxel_definition_output );
					
					//16-bit output (with double, not single precision calculation)
					//run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.DoubleReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new DoubleReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.DoubleStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.SHORT, false, 0.01, 0).deconvolve().show();\n";
					//print (run_script);
					//eval( "js", run_script );  //throw away result
					
					//okay write TIF file
					delete_result = File.delete( outputDirectory + file_sep + fileoutname[b] + ".tif" ); //throw away result
					
					saveAs("Tiff", outputDirectory + file_sep + fileoutname[b] + ".tif" );
					close();
					print("  ..saved "+fileoutname[b]+".tif");
					
					//close original image
					//if ( isOpen(img_id) ) {
						selectImage(img_id);
						close();
					//}
					//if ( isOpen(psf_id[b]) ) {
						selectImage(psf_id[b]);
						close();
					//}
					//delete image block on disk
					if ( channelList.length > 2 && b > 0 ) {
						delete_result = File.delete( outputDirectory + file_sep + channelList[b] + ".tif" );
					} 
					
					//write to logfile
					File.append( file + ": angle " + string_angles[a] + ", time " + d2s(t+1,0) + " of " + d2s(max_time,0) + ", channel " + d2s(b,0) + " >> " + outputDirectory + file_sep + fileoutname[b] + ".tif", outputDirectory + file_sep + "CZI_to_TIF_with_deconvolution.txt" );
					//close("*");
					call("java.lang.System.gc");
				}
			}
			//Close original concatenated stack
			//close("C*");   
			this_time++; //increment time counter
		}
		//close();
		Ext.close();
		call("java.lang.System.gc");
		print( "  ..finished saving deconvolved images to " + outputDirectory ); 
	}
}

function getValueByKeyJson(key, jsonText) {
	index = indexOf(jsonText, '"' + key + '":');
	if (index != -1) {
		startIndex = index + lengthOf(key) + 3;
		endIndex = -1;
		str = replace(String.trim(substring(jsonText, startIndex)), "\n", "");
		startIndex = 0;
		if (startsWith(str, '"')) {
			startIndex = 1;
			endIndex = indexOf(str, '"', startIndex);
		} else if (startsWith(str, "'")) {
			startIndex = 1;
			endIndex = indexOf(str, "'", startIndex);
		} else if (startsWith(str, "{")) {
			startIndex = 1;
			endIndex = indexOf(str, '}', startIndex);
		} else if (startsWith(str, "[")) {
			startIndex = 1;
			endIndex = indexOf(str, ']', startIndex);
		} else {
			endIndex = indexOf(str, ",", startIndex);
			if (endIndex == -1) {
				endIndex = indexOf(str, "}", startIndex);
			}
		}

		/*if (startsWith(str, '"') && endsWith(str, '"')) {
 			str = substring(str, 1, lengthOf(str) - 1);
		}*/
		//return str;
		if ( endIndex < 0 ) {
			return replace(str, " ", "");
		} else {
			return replace(substring(str, startIndex, endIndex), " ", "");
		}
	}
	return "";
}

function main_ijm_deconvolve_series(directory,outputDirectory,do_twice) {	

	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);
	
	if ( outputDirectory == "" || !File.exists(outputDirectory) ) {
		outputDirectory = getDirectory("Choose output directory");
	}

	//Count the maximum number of positions and slices in dataset
	run("Bio-Formats Macro Extensions");

	//print ( "directory: " + directory );
    add_this_file = false;
	processList = newArray(0);
	unique_PSF_filenames = newArray(0);
	unique_PSF_parameters = newArray(0);
	
	//determine if Z.1 or MuVi SPIM
	type = "";
	if (File.isDirectory(directory + "raw") && File.exists(directory + "bdv.h5")) {
		type = "json/h5";
	} 
	for (i=0; i<fileList.length; i++) {
		if (endsWith(fileList[i], ".czi")) { // Zeiss Z.1 CZI format
			type = ".czi";
			break;
		}
	}
	//print ("type = " + type );
	if ( type == ".czi" ) {
		for (i=0; i<fileList.length; i++) {
			if (endsWith(fileList[i], ".czi")) {
				
				//now, remove .czi from end of filename
				filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
				//print( fileList[i] + " -> " + filepaddedname );
				
				//filepaddedname = name_ext[0];
				//print( "File " + d2s(i,0) + " ..." + filepaddedname + " " + fileList[i] );
				//need to check for presence of \((\d+)\) at end of filename, and if not present, add "(00000)"
	
				if ( endsWith(filepaddedname,".") ) { //additional trailing periods should be removed
					while( endsWith(filepaddedname,".") ) {
						filepaddedname = substring(filepaddedname,0,filepaddedname.length);					
					}
				}
							
				if ( endsWith(filepaddedname,")") && matches(filepaddedname, ".*\\(\\d+\\)$") ) {
					//make sure the digits are padded for name sort
					//print ( "Matches...\n" );
	                
	                //this is a file that is a later file in a time sequence, so we don't need to add it
	                add_this_file = false;
	                
					filepaddedname = replace(filepaddedname, "\\(\(\\d\)\\)", "(" + "0000$1" + ")"); //replace single digits with padding
					filepaddedname = replace(filepaddedname, "\\(\(\\d\\d\)\\)", "(" + "000$1" + ")" ); //replace double digits with paddi
					filepaddedname = replace(filepaddedname, "\\(\(\\d\\d\\d\)\\)", "(" + "00$1" + ")" ); //replace triple digits with padding
					filepaddedname = replace(filepaddedname, "\\(\(\\d\\d\\d\\d\)\\)", "(" + "0$1" + ")" ); //replace quadruple digits with padding
		
				} else {
					//add that stuff to end
					filepaddedname = filepaddedname + "(00000)";
	                
	                //this is a principle file, so use it
	                add_this_file = true;
				}
				
				filepaddedname = replace(filepaddedname, "\(\\D\)\(\\d\)\(\\D\)", "$1" + "0000" + "$2$3" ); //replace single digits with padding
				filepaddedname = replace(filepaddedname, "\(\\D\)\(\\d\\d\)\(\\D\)", "$1" + "000" + "$2$3" ); //replace double digits with paddi
				filepaddedname = replace(filepaddedname, "\(\\D\)\(\\d\\d\\d\)\(\\D\)", "$1" + "00" + "$2$3" ); //replace triple digits with padding
				filepaddedname = replace(filepaddedname, "\(\\D\)\(\\d\\d\\d\\d\)\(\\D\)", "$1" + "0" + "$2$3" ); //replace quadruple digits with padding
	            
	            if ( add_this_file ) {
	                processList = Array.concat( processList, filepaddedname + "///" + fileList[i] ); //new name will consist of padded old name so we can sort correctly by timepoint
	            } 
			}
		}

	} else if ( type == "json/h5" ) {
		//print( "getting :" + directory + "raw" );
		folderList = getFileList(directory + "raw");
		fileList = newArray(0);
		
		for (i=0; i<folderList.length; i++) {
			//print( " checking1: " + directory + "raw" + file_sep + folderList[i] );
			if ( File.isDirectory(directory + "raw" + file_sep + folderList[i]) ) {
				if ( endsWith(folderList[i],"/") || endsWith(folderList[i],file_sep) ) {
					folderList[i] = substring(folderList[i],0,lengthOf(folderList[i])-1);
				} 
				thisFileList = getFileList(directory + "raw" + file_sep + folderList[i]);
				
				for (ii=0; ii<thisFileList.length; ii++) {
					fileList = Array.concat( fileList, "raw" + file_sep + folderList[i] + file_sep + thisFileList[ii] );
					//print( " go " + "raw" + file_sep + folderList[i] );
				}
			}
		}
		
		for (i=0; i<fileList.length; i++) {
			if (endsWith(fileList[i], ".h5")) {

				//now, remove .h5 from end of filename
				filepaddedname = substring(fileList[i], 0, fileList[i].length-3);
				
				//remove .lux-XXX extension from end of filename
				dotIndex = lastIndexOf(filepaddedname, ".");
				if (dotIndex >= 0) {
					filepaddedname = substring(filepaddedname, 0, dotIndex);
				}

				view_setups = newArray(0);
				//confirm json is present
				//print ( "  check " + directory + filepaddedname + ".json" );
				
	            if ( File.exists(directory + filepaddedname + ".json") ) {
					//open json to extract timepoint and view angles
					input_json = File.openAsString( directory + filepaddedname + ".json");

					// Extract the desired information
					cam_left_to_right =  getValueByKeyJson("cam_left_to_right", input_json);
					time_point =  parseInt( getValueByKeyJson("time_point", input_json) );
					channel =  parseInt( getValueByKeyJson("channel", input_json) );
					camera =  getValueByKeyJson("camera", input_json);
					start_deg =  parseFloat( getValueByKeyJson("start_deg", input_json) );
					cam_left_to_right_list = split(cam_left_to_right, ",");
					cam_left_to_right = cam_left_to_right_list[0];

					filepaddedname = "LSFM_T" + IJ.pad(d2s(time_point,0),4) + "_C" + IJ.pad(d2s(channel,0),1);
					this_angle = round(start_deg);
					if ( this_angle > 360 ) {
						this_angle -= 360;
					} else if ( this_angle < 0 ) {
						this_angle += 360;
					}
					this_view = -1;
					if ( camera == "left" || camera == "right"  ) {
						this_view = 0; // left follows far-to-near acquisition, right is opposite -- will use this later to reverse PSF for right camera
						if ( camera == "right" ) {
							//this_angle -= 180;
							this_view = 1;
						}

						/*if ( this_angle > 360 ) {
							this_angle -= 360;
						} else if ( this_angle < 0 ) {
							this_angle += 360;
						}
						*/

						//now, segment viewsetup into AXXX_VX and then create view setup with that identity
						filepaddedname = filepaddedname + "_A" + IJ.pad(d2s(this_angle,0),3) + "_V" + IJ.pad(d2s(this_view,0),1);
						//print( "filepaddedname: " + filepaddedname );
						processList = Array.concat( processList, filepaddedname + "///" + fileList[i] + "///" + d2s(time_point,0) + "///" + d2s(channel,0) + "///" + d2s(this_angle,0) + "///" + d2s(this_view,0) ); //new name will consist of padded old name so we can sort correctly by timepoint
						//print ("   adding: " + filepaddedname + "///" + fileList[i] );
					}


	            }
			}
		}

	}
	
	Array.sort(processList); //sort files by their intended name, not true name
	//Array.print(processList);
	
	//NAlightsheet = 0;
	//WLlightsheet = 0;
	nPositions = 1;
	string_angles = newArray(nPositions);
	Array.fill(string_angles,NaN); 
	
	angles = newArray(1);
	angles[0] = 0;
	RIlightsheet = default_immersion_refractive_index;
	max_channels = 0;
    max_time = 0;
	this_time = 0;
	Olightsheet = "";
	//exit();
	
	setBatchMode(true);
	input_json = "";
	/*	for (i=0; i<processList.length; i++) {
			print ( "process list: " + processList[i] );
		}*/

	for (i=0; i<processList.length; i++) {
	//for (i=0; i<1; i++) {
		
		name_ext = split(processList[i],"(///)");
		//file = directory + name_ext[1];
		file = directory + name_ext[1];
		if ( type == ".czi" ) {
			
			Ext.setId(file);
			Ext.getSeriesCount(nPositions);
		
			//get refractive index
			//although this will retrieve metadata only from active series, refractive index is not different for each channel or view but instead is the same for everything
			RIlightsheet = 0;
			Ext.getMetadataValue("Information|Image|RefractiveIndex #1", RIlightsheet);
			if ( RIlightsheet == 0 ) {
				Ext.getMetadataValue("Information|Image|RefractiveIndex", RIlightsheet);
			}
			//print ( "Metadata: \n" + metadata_value  + "\n" ); exit();

			//get view angles, again metadata permeates all series and these values can be extracted with any series open
			string_angles = newArray(nPositions);
			angles = newArray(nPositions);
			Array.fill(string_angles,NaN); // default angle is NaN
			Array.fill(angles,NaN); // default angle is NaN
			
			//get objective name
			Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|Objective #1", Olightsheet);
			if ( Olightsheet == "" || Olightsheet == 0 ) {
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|Objective", Olightsheet);
			}
			//grab initial angles from file
			if ( nPositions < 1) {
				print ( "Fewer than one view/position present in file " + directory + name_ext[1] + "!\n" );
				continue;
			} else {
				//fill angles array
				for(a=0; a<nPositions; a++) {
					Ext.getMetadataValue("Information|Image|V|View|Offset #" + d2s((a+1),0), string_angles[a]);
				}

				//fill actual angles array with numbers
				for(a=0; a<nPositions; a++) {
					angles[a] = parseFloat( string_angles[a]);
				}
			}
		} else if ( type == "json/h5" ) {
			//filepaddedname = substring(name_ext[0], 0, name_ext[0].length-4);
			filepaddedname = name_ext[1];
			//print( "tried to open: " + name_ext[1] );
			dotIndex = lastIndexOf(filepaddedname, ".");
			if (dotIndex >= 0) {
				filepaddedname = substring(name_ext[1], 0, dotIndex);
			}
			dotIndex = lastIndexOf(filepaddedname, ".");
			if (dotIndex >= 0) {
				filepaddedname = substring(name_ext[1], 0, dotIndex);
			}
			
			
			//print( "tried to open: " + directory + filepaddedname + ".json" );
			if ( type == "json/h5" && File.exists(directory + filepaddedname + ".json") ) {
				//open json to extract timepoint and view angles
				input_json = File.openAsString( directory + filepaddedname + ".json");
			}
			
			//print( input_json); exit(0);
		
			// Extract the desired information
			refractive_index =  parseFloat(getValueByKeyJson("refractive_index", input_json));
			if ( isNaN(refractive_index) || refractive_index < 1 ) {
				refractive_index = default_immersion_refractive_index;
			}
			angles[0] =  parseInt(name_ext[4]); //round(parseFloat( getValueByKeyJson("start_deg", input_json) ));
			string_angles[0] = name_ext[4] + "_" + name_ext[5];

		}

		/*
		//modify angles to try to put front at zero degrees
		if ( nPositions == 1 ) {
			//set only available angle to zero degrees
			if ( type == ".czi" ) {
				angles[0] = 0;
			} else if () {


			}
		} else { // 2 or more angles
		
			//first, figure out whether we should use angles from -180 to +180, or 0 to 360 -- to do this, if the difference between the largest and smallest angles is greater than 180, use -180 to 180 rather than 0 to 360
			//fill sort angles array
			sort_angles = Array.copy(angles);
			Array.sort(sort_angles);
			
			if ( sort_angles[sort_angles.length-1] - sort_angles[0] > 180 ) {
				//modify angles to use -180 to 180
				for ( a=0; a<angles.length; a++ ) {
					if ( angles[a] > 180 ) {
						angles[a] -= 360;
					}
				}
				
				//redo sort_angles array now since angles have been modified
				sort_angles = Array.copy(angles);
				Array.sort(sort_angles);
			}
			
			if( nPositions == 2 ) {
				abs_diff = sort_angles[1] - sort_angles[0];
				
				if ( abs_diff < 180.2 && abs_diff > 179.8 ) { // views are opposing, so assume one is front and one is back
					//set view 1 to angle 0 and view 2 to angle 180, then be done
					angles[0] = 0;
					angles[1] = 180;
				} else if (angles[1]>=angles[0])  { // views are not opposing and second view has greater angle than first view, so assume that front is actually midpoint of the two views -- and set that to 0 degrees
					angles[0] = 360 - (abs_diff /2);
					angles[1] = abs_diff/2;
				} else { // views are not opposing and first view has greater angle than second view, so assume that front is actually midpoint of the two views -- and set that to 0 degrees
					angles[1] = 360 - (abs_diff /2);
					angles[0] = abs_diff/2;
				}
			} else { // 3 or more positions, so zero the middle view if odd or view right before middle view (assume views done left-front-right-back or right-front-left-back)
				//zero the second view 
				middle_view = floor((angles.length / 2) - 0.5);
                		diff_angle = angles[middle_view];
				for ( a=0; a<angles.length; a++ ) {
					angles[a] -= diff_angle;
					
					if ( angles[a] < -180 ) {
						angles[a] += 360;	
					} else if ( angles[a] >= 360 ) {
						angles[a] -= 360;	
					}
				}
				
			}
			
			//okay, now correct angles for 0 to 360, regardless of which scheme we used above
			for ( a=0; a<angles.length; a++ ) {
				if ( angles[a] >= 360 ) {
					angles[a] -= 360;
				} else if ( angles[a] < 0 ) {
					angles[a] += 360;
				}
			}
		}
		*/
		
		//get lightsheet data for each channel
		num_channels = 1;
        num_time = 0;
		max_time = 1;

		if ( type == ".czi" ) {
			for(a=0; a<nPositions; a++) {
				Ext.setSeries(a);
				Ext.getSizeC(num_channels);
	            Ext.getSizeT(num_time);

				if ( num_channels > max_channels ) {
					max_channels = num_channels;
				}
				if ( num_time > max_time ) {
					max_time = num_time;
				}
			}
		}

		NAlightsheets = newArray(max_channels);
		NAdetections = newArray(max_channels);
		WLlightsheets = newArray(max_channels);
		WLdetections = newArray(max_channels);

		Array.fill(NAlightsheets,0); // default NA is 0
		Array.fill(NAdetections,0); // default NA is 0
		Array.fill(WLlightsheets,0); // default WL is 0
		Array.fill(WLdetections,0); // default WL is 0

		objective_wd = 0;
		objective_immersion = "";
		if ( type == ".czi" ) {
			for(a=0; a<max_channels; a++) {
				Ext.getMetadataValue("Information|Image|Channel|NALightSheet #" + d2s((a+1),0), NAlightsheets[a]);
				Ext.getMetadataValue("Information|Image|Channel|NADetection #" + d2s((a+1),0), NAdetections[a]);			
				Ext.getMetadataValue("Information|Image|Channel|IlluminationWavelength|SinglePeak #" + d2s((a+1),0), WLlightsheets[a]);
				Ext.getMetadataValue("Information|Image|Channel|DetectionWavelength|SinglePeak #" + d2s((a+1),0), WLdetections[a]);		
				
				print ("RIlightsheet: " + RIlightsheet );
				print ("NAlightsheets["+d2s(a,0)+"]: " + NAlightsheets[a] );
			}
		} else if ( type == "json/h5" ) {
			// detection objective information -- cue json_obj text to beginning of selected objective
			objective =  getValueByKeyJson("objective", input_json);
			index = indexOf(input_json, '"objectives":');
			json_obj = substring(input_json, index + lengthOf("objectives") + 3);
			index = indexOf(json_obj, '"name": ' + '"' + objective + '"' );
			json_obj = substring(json_obj, index + lengthOf(objective) + 10);

			// now gather objective data
			NAdetections[0] = parseFloat(getValueByKeyJson("na", json_obj));
			objective_wd = 1000 * parseFloat( replace( getValueByKeyJson("wd_mm", json_obj), "[^0-9.]", "") );
			objective_immersion = getValueByKeyJson("immersion", json_obj);

			if ( objective_immersion == "water" ) {
				RIlightsheet = 1.333;
			} else if ( objective_immersion == "air" ) {
				RIlightsheet = 1;
			} else if ( objective_immersion == "glycerol" ) {
				RIlightsheet = 1.47;
			}

			// laser information
			index = indexOf(input_json, '"lasers":');
			//json_laser = substring(input_json, index + lengthOf("lasers") + 3);
			
			//startIndex = indexOf(jsonString, "\"lasers\"");
			endIndex = indexOf(input_json, "]", index) + 1;
			lasersString = substring(input_json, index, endIndex);

			// Split the "lasers" string into individual laser objects
			laserList = split(lasersString, "{");

			wavelength_sum = 0.0;
			wavelength_number = 0;
			
			
			// Loop through each laser object
			for (ii = 1; ii < lengthOf(laserList); ii++) {
			    laserObject = laserList[ii];
			
			    // Extract the "name" and "on" values from the laser object
			    nameIndex = indexOf(laserObject, "\"name\"");
			    onIndex = indexOf(laserObject, "\"on\"");
			    onValue = substring(laserObject, onIndex + 6, indexOf(laserObject, ",", onIndex + 7));
			
			    // Parse the "wavelength_nm" if the laser is on
			    if (onValue == "true") {
			        wavelengthIndex = indexOf(laserObject, "\"wavelength_nm\"");
			        wavelengthValue = parseFloat( substring(laserObject, wavelengthIndex + 17, indexOf(laserObject, ",", wavelengthIndex + 19)) );
			
			        // Print the laser name and wavelength
			        //print("Wavelength: " + wavelengthValue);

					wavelength_sum += wavelengthValue;
					wavelength_number++;
			    }
			}
			
			// Average overall laser wavelength, ignoring intensities, etc
			excitation = wavelength_sum / wavelength_number;

			// default emission is 30nm + excitation
			emission = excitation + 30;
			
			// filter wheel information
			index = indexOf(input_json, '"filterwheels":');
			json_fw = substring(input_json, index + lengthOf("filterwheels") + 3);
			index = indexOf(json_fw, '"name": ' + '"' + 'fw_' + objective + '"' );
			json_fw = substring(json_fw, index + lengthOf(objective) + 10);

			// gather emission data
			this_filter = getValueByKeyJson("selection", json_fw);
			if ( startsWith(this_filter, "BP")  ) {
				this_filter = substring(this_filter, 2);
				wavelengths = split(this_filter,"-");
				emission0 = parseFloat(replace(wavelengths[0], "[^0-9.]", ""));
				emission1 = parseFloat(replace(wavelengths[1], "[^0-9.]", ""));
				
				if ( isNaN(emission0) || emission0 < 0 ) {
					if ( !(isNaN(emission1)) && emission1 > 0 ) {
						emission = emission1;
					}
				} else if ( isNaN(emission1) || emission1 < 0 ) {
					if ( !(isNaN(emission0)) && emission0 > 0 ) {
						emission = emission0;
					}
				} else {
					emission = ( emission0 + emission1 ) / 2;
				}

			} else if ( startsWith(this_filter, "LP")  ) {
					emission = parseFloat( replace( this_filter, "[^0-9.]", "") );
					emission = emission + 30;
			} else if ( startsWith(this_filter, "SP")  ) {
					emission = parseFloat( replace( this_filter, "[^0-9.]", "") );
					emission = emission - 30;
			} else {
				emission = parseFloat( replace( this_filter, "[^0-9.]", "") );
			}
			
			if ( isNaN(excitation) || excitation < 0 ) {
				if ( !(isNaN(emission)) && emission > 0 ) {
					excitation = emission - 30;
				}
			}
		
			//print ( "emission: " + emission );
			
			WLlightsheets[0] = d2s(excitation,3);
			WLdetections[0] = d2s(emission,3);
			
			// lightsheet thickness is located in beamexpander section -- this is w0 or the radius at beam waist
			index = indexOf(input_json, '"beamexpander":');
			json_beamexpander = substring(input_json, index + lengthOf("beamexpander") + 3);
			index = indexOf(json_beamexpander, '"name": "beamexp"' );
			json_beamexpander = substring(json_beamexpander, index + lengthOf("beamexp") + 10);
			this_waist = getValueByKeyJson("selection", json_beamexpander);
			//print( "json_beamexpander : " + this_waist ); 
			this_waist = replace( this_waist, "&micro;", "u");
			this_waist_multiplier = 1.0;
			
			if ( indexOf(this_waist,"um") >= 0 ) {
				//microns -- do nothing
			} else if ( indexOf(this_waist,"mm") >= 0 ) {
				this_waist_multiplier = 1000.0;
			} else if ( indexOf(this_waist,"cm") >= 0 ) {
				this_waist_multiplier = 10000.0;
			} else if ( indexOf(this_waist,"nm") >= 0 ) {
				this_waist_multiplier = 0.001;
			}
			
			// calculate lightsheet NA					
			// beam_waist = RI_immersion * lambda_lightsheet / ( NA_lightsheet * PI );
			// NA = (1.55) * .488 / ( PI * 3.889772 )
			//	( PI * 3.889772 ) = 12.22						
			//NAlightsheets[0] = d2s(RIlightsheet * WLlightsheets[0]/1000 / ( PI * beam_waist),10);
			NAlightsheets[0] = d2s(RIlightsheet * (parseFloat(WLlightsheets[0])/1000) / ( PI * parseFloat( replace(this_waist, "[^0-9.]", ""))  * this_waist_multiplier ),10);
			//print ("RIlightsheet: " + RIlightsheet );
			//print ("NAlightsheets[0]: " + NAlightsheets[0] + " based on " + d2s(parseFloat( replace(this_waist, "[^0-9.]", "")),6) );
		}
		
		
		
		for (t=0; t<max_time; t++ ) {		
			for(a=0; a<nPositions; a++) {
				//print ("About to open " + "Bio-Formats", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT series_"+ d2s(a+1,0) + " t_begin_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_end_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_step_" + d2s(a+1,0) + "=1" );
				//exit();                
				//run("Bio-Formats", "open=[/media/martin/Dominguez Data/Onset of F6-nGFP/2019-01-24 embryo B -- keep, 2v/1-24 embryo B 00.czi] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT series_2 t_begin_2=3 t_end_2=3 t_step_2=1");
				if ( type == ".czi" ) {
					if ( nPositions == 1 ) {
						run("Bio-Formats", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT " + " t_begin=" + d2s(t+1,0) + " t_end=" + d2s(t+1,0) + " t_step=1" );	
					} else {
						run("Bio-Formats", "open=[" + file + "] color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT series_"+ d2s(a+1,0) + " t_begin_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_end_" + d2s(a+1,0) + "=" + d2s(t+1,0) + " t_step_" + d2s(a+1,0) + "=1" );
					}
				} else if ( type == "json/h5" ) {
					//run("N5", "n5=["+file+"/Data]");
					file_js = replace(file, "\\", "/");
					js = 'importClass(Packages.loci.plugins.in.ImporterOptions);importClass(Packages.org.janelia.saalfeldlab.n5.ij.N5Importer);var filePath = '+ '"' + file_js + '/Data' + '"' + ';var importer = new N5Importer();importer.process(filePath, false );';
					//print( "JS=" + js );
					result = eval("js",js);
					//print( "eval'd" );
				}
				
				//Get name of opened stack	
				rand_num = floor(random() * 1000000 );
				rename( "Image" + rand_num );
				title = getTitle();
				getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
				
				getVoxelSize(vox_width, vox_height, vox_depth, vox_unit);
				if ( type == "json/h5" ) {
					voxel_size_um = getValueByKeyJson("voxel_size_um", input_json);
					vox_width = parseFloat( getValueByKeyJson("width", voxel_size_um) );
					vox_height =  parseFloat( getValueByKeyJson("height", voxel_size_um) );
					vox_depth =  parseFloat( getValueByKeyJson("depth", voxel_size_um) );
					vox_unit = "um";
					//print ( "vox_json: " + voxel_height_um + " " + voxel_depth_um );
				}
				
				voxel_definition_output = "unit=" + vox_unit + " pixel_width=" + d2s(vox_width,10) + " pixel_height=" + d2s(vox_height,10) + " voxel_depth=" + d2s(vox_depth,10); //will use later when we output the file
				//print ( " vox -- " + voxel_definition_output );
				//make all units nanometers
				if (vox_unit=="nm" || vox_unit=="nanometers" || vox_unit=="nanometer") {
					//do nothing
				} else if (vox_unit==getInfo("micrometer.abbreviation") || vox_unit=="um" || vox_unit=="microns" || vox_unit=="micron") {
					vox_width * = 1000;
					vox_height * = 1000;
					vox_depth * = 1000;
				} else if (vox_unit=="mm" || vox_unit=="millimeters" || vox_unit=="millimeter"){
					vox_width * = 1000000;
					vox_height * = 1000000;
					vox_depth * = 1000000;
				} else if (vox_unit=="cm" || vox_unit=="centimeters" || vox_unit=="centimeter"){
					vox_width * = 10000000;
					vox_height * = 10000000;
					vox_depth * = 10000000;
				} else if (vox_unit=="m" || vox_unit=="meters" || vox_unit=="meter"){
					vox_width * = 1000000000;
					vox_height * = 1000000000;
					vox_depth * = 1000000000;
				} else {
					print ( "Cannot interpret voxel unit " + vox_unit + " for " + file + ", series " + d2s(a,0) + "\n" );
					continue;
				}
				
				//check height and width scales
				if ( vox_width > 1.001 * vox_height || vox_height > 1.001 * vox_width ) {
					print ( "Unexpected XY scale difference: " + d2s(vox_height,5) + " vs. " + d2s(vox_width,5) + " for " + file + ", series " + d2s(a,0) + "\n" );
					continue;
				}
		
				//Split channels and record names of each new image stack
				channelList = newArray(0);
				if ( dim_channels > 1 ) {
					run("Split Channels");
					for (b=1; b<=dim_channels; b++ ) {
						channelList = Array.concat( channelList, "C" + IJ.pad(b,1) + "-" + title ); 
					}
				} else if ( dim_channels == 1 ) {
					channelList = Array.concat( channelList, title );
				} else {
					continue; //skip this series altogether	
				}
		
				fileoutname = newArray(channelList.length);
				viewpsfname = newArray(channelList.length);
				psf_id = newArray(max_channels);
				Array.fill(psf_id,NaN); // default is no fill
				
				//lambdas = newArray(max_channels);
				//Array.fill(lambdas,0); //no regularization by default
				
				//figure out file out names
				
				if ( type == ".czi" ) {
					for (b=0; b<channelList.length; b++ ) {
						fileoutname[b] = "LSFM_T" + IJ.pad(this_time,4) + "_C" + IJ.pad(b,1) + "_A" + IJ.pad(d2s(round(angles[a]),0),3);
						viewpsfname[b] = 0;
					}
				} else if ( type == "json/h5" ) {
					for (b=0; b<channelList.length; b++ ) {
						fileoutname[b] = "LSFM_T" + IJ.pad(name_ext[2],4) + "_C" + IJ.pad(name_ext[3],1) + "_A" + IJ.pad(d2s(round(angles[a]),0),3) + "_V" + IJ.pad(b,1) + IJ.pad(name_ext[5],1);
						viewpsfname[b] = name_ext[5];
					}
				}
				
				//do PSFs separately since they require batchmode to be off
				for (b=0; b<channelList.length; b++ ) {
					//print( "c"+b + ": " + d2s(NAlightsheets[b],4) + " " +  d2s(WLlightsheets[b],4) + " " + d2s(RIlightsheet,4) + " " + d2s(angles[a],4) );
					if ( NAlightsheets[b] > 0 && WLlightsheets[b] > 0 && RIlightsheet > 0 && !(isNaN(angles[a]) ) ) {
						
						this_PSF_identifier = NAdetections[b] + "//" + NAlightsheets[b] + "//" + RIlightsheet + "//" + WLlightsheets[b] + "//" + WLdetections[b] + "//" + d2s(vox_width,8) + "//" + d2s(vox_depth,8) + "//" + d2s(dim_width/5,8) + "//" + Olightsheet + "//" + viewpsfname[b];
						
						//see if we have done this PSF before
						is_match = -1;
						for (m=0;m<unique_PSF_parameters.length; m++ ) {
							if ( this_PSF_identifier == unique_PSF_parameters[m] ) {
								is_match = m;
								break;
							}
						}
						
						if ( is_match >= 0 ) {
							//copy file and be done; we've done it before
							//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif is match with " + unique_PSF_filenames[is_match] + "...\n" + this_PSF_identifier );
							//File.copy( unique_PSF_filenames[is_match], outputDirectory + fileoutname[b] + "_psf.tif" );
							open( unique_PSF_filenames[is_match] );
							psf_id[b] = getImageID();
							
						} else { //new PSF, need to generate
							//generate and save lightsheet  PSF
							//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif is no match...\n" + this_PSF_identifier );
							
							setBatchMode("exit and display"); //needed for PSF generation
							setBatchMode(false);

							WD_CS_return = newArray(3); // WD, coverslip, nominal RIimmersion

							if ( type == ".czi" ) {
								WD_CS_return = Zeiss_WDforObjective( Olightsheet );
							} else if ( type == "json/h5" ) {
								WD_CS_return[0] = objective_wd;
								WD_CS_return[1] = 0;
								
								if ( indexOf( objective_immersion, "air") >= 0 || indexOf( objective_immersion, "dry") >= 0 ) {
									WD_CS_return[1] = 170;
								}
							}
							if ( WD_CS_return[2] > RIlightsheet ) {
								RIlightsheet = WD_CS_return[2];
							}			
							tissue_RI = default_tissue_refractive_index;
							if ( RIlightsheet > tissue_RI ) {
								tissue_RI = RIlightsheet;
							}
							
							PSF_name = LS_PSF_generator ( parseFloat(NAdetections[b]), parseFloat(NAlightsheets[b]), parseFloat(RIlightsheet), tissue_RI, parseFloat(WLlightsheets[b]), parseFloat(WLdetections[b]), vox_width, vox_depth, dim_width/5, WD_CS_return[0], WD_CS_return[1] ); //estimate lightsheet radius from w(0) as being 1/4th the total width of the acquired image -- this should cut it about midway
							selectWindow( PSF_name );
							
							// reverse PSF for right-camera views since they have far-field late in the stack rather than early
							if ( type == "json/h5" && parseInt(viewpsfname[b]) > 0 ) {
								run("Reverse");
							}
							
							run("Enhance Contrast...", "saturated=0 process_all use");
							run("16-bit");
							
							//print( "Delete result " + outputDirectory + fileoutname[b] + "_psf.tif: " + delete_result ); 
							psf_id[b] = getImageID();
							
							delete_result = File.delete( outputDirectory + file_sep + fileoutname[b] + "_psf.tif" ); //throw away result
							saveAs("Tiff", outputDirectory + file_sep + fileoutname[b] + "_psf.tif" );
							//close();
							//exit();
							setBatchMode(true);
							
							//now store in arrays for future use
							unique_PSF_filenames = Array.concat( unique_PSF_filenames, outputDirectory + file_sep + fileoutname[b] + "_psf.tif" );
							unique_PSF_parameters = Array.concat( unique_PSF_parameters, this_PSF_identifier );
						}

					} else {
						//this_PSF_identifier = NAdetections[b] + "//" + NAlightsheets[b] + "//" + RIlightsheet + "//" + WLlightsheets[b] + "//" + WLdetections[b] + "//" + d2s(vox_width,8) + "//" + d2s(vox_depth,8) + "//" + d2s(dim_width/4,8) + "//" + Olightsheet;
						//print ( "File " + outputDirectory + fileoutname[b] + "_psf.tif cannot be saved...\n" + this_PSF_identifier );
						psf_id[b] = NaN;
						
					}
				}				
				
				//remove blank images -- from all channels together to prevent incorrect registrations later
				selectWindow(channelList[0]); //use channel 0 as primary channel
				setSlice(dim_slices);
				getStatistics( area, mean );
				if ( mean < 10 ) { //we potentially have a
					//print( "Mean is : " + d2s(mean,4) + " on slice " + d2s(dim_slices,0) );
					
					final_slice = dim_slices-1;
					
					for (ss=dim_slices-10; ss>0; ss-=10 ) { //decrement 10 to find the first slice with data
						setSlice(ss);
						getStatistics( area, mean );
						//print( " ...checking mean: " + d2s(mean,4) + " on slice " + d2s(ss,0) );
						if ( mean > 10 ) {
							//okay pinpoint exact end of stack now
							for (tt=ss+9; tt>=ss; tt-- ) { //decrement individual slices to find the first slice with data
								setSlice(tt);
								getStatistics( area, mean );
								//print( "    ...checking mean: " + d2s(mean,4) + " on slice " + d2s(tt,0) );
								if ( mean > 10 ) { //okay we found the end
									final_slice = tt;
									break;
								}
							}
							break;
						}
					}
					
					for (b=0; b<channelList.length; b++ ) {
						selectWindow(channelList[b]);
						img_id = getImageID();
						if (isNaN(psf_id[b]) || isNaN(img_id) ) {
							continue; //don't process this image
						}						
						run("Slice Remover", "first=" + d2s(final_slice+1,0) + " last=" + d2s(dim_slices,0) + " increment=1");
					}
				}
				
				if ( deconvolution_subtract_camera_noise ) {
					//preprocess images, first by trying to remove noisy in the background
					for (b=0; b<channelList.length; b++ ) {
						selectWindow(channelList[b]);
						img_id = getImageID();
						if (isNaN(psf_id[b]) || isNaN(img_id) ) {
							continue; //don't process this image
						}					
						
						//calculate minimum projection and subtract from stack (attempt to remove noise)
						run("Z Project...", "projection=[Min Intensity]");
						average_img = getTitle();
						imageCalculator("Subtract stack",channelList[b],average_img);
						selectWindow(average_img);
						close();
					}
				}
				
				
				//do deconvolution here
				for (b=0; b<channelList.length; b++ ) {
					selectWindow(channelList[b]);
					img_id = getImageID();
					//print( "Processing image " + d2s(img_id,0) + " with title " + channelList[b] );
					if (isNaN(psf_id[b]) || isNaN(img_id) ) {
						continue; //don't process this image
					}
					//do deconvolution here
					
					//try to estimate regularization parameter based on image histogram -- this should be revised
					Stack.getStatistics( area, mean, min, max, std );
					
					//
					//
					//  Deconvolution Round #1
					//
					//arbitrarily set regularization parameter, based on 16-bit images input
					if ( do_twice ) {
						//two sequential (iterative) rounds of deconvolution
						lambda = 4 * deconvolution_regression_parameter / ( std + (4.6*deconvolution_regression_parameter) ); // 10/std was the original, but produced ringing artifacts -- therefore try 50/std , but in truth we should do two iterations stepping down from lambda 1000/std then 100/std
						
					} else {
						//lambda = 75 / ( std + 100); // 10/std was the original, but produced ringing artifacts -- therefore try 50/std , but in truth we should do two iterations stepping down from lambda 1000/std then 100/std
						lambda = deconvolution_regression_parameter / ( std + deconvolution_regression_parameter);
					}

					//run deconvolution
					//32-bit output (with single precision calculation)
					run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n";
					if ( deconvolution_iterations > 1 ) {
						run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
					}
					//print (run_script);
					//print ( "1_orig ID " + img_id + ", next_img_ID " + getImageID() );
					if ( deconvolve_yes ) {
						eval( "js", run_script );  //throw away result
					} else {
						run("Duplicate...", "duplicate");
					}
					
					if ( do_twice ) {
						next_img_id = getImageID(); //get new image ID
						
						//close original image
						selectImage(img_id);
						close();
						img_id = next_img_id;
						selectImage(img_id); //go back to deconvolved image
						
						
					}
					
					//adjust properties then save
					setMinAndMax(0, 65535); //clip the upper bounds of signal (lower bounds already clipped for nonnegativity)
					run("16-bit");
					run("Properties...", voxel_definition_output );

					if ( do_twice ) {	
						//
						//
						//  Deconvolution Round #2
						//
						//arbitrarily set regularization parameter, based on 16-bit images input
						lambda = 4 * deconvolution_regression_parameter / ( std + (4.6*deconvolution_regression_parameter) ); // 10/std was the original, but produced ringing artifacts -- therefore try 50/std , but in truth we should do two iterations stepping down from lambda 1000/std then 100/std
	
						//run deconvolution
						run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.gtik.FloatReflexiveGeneralizedTikhonov3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.spectral.SpectralEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar frgt = new FloatReflexiveGeneralizedTikhonov3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), SpectralEnums.FloatStencil3DType.LAPLACIAN.stencil, SpectralEnums.ResizingType.NONE, Enums.OutputType.FLOAT, false, " + d2s(lambda,8) + ", 0).deconvolve().show();\n";
						if ( deconvolution_iterations > 1 ) {
							run_script = "importClass(WindowManager);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver3D);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions);\nimportClass(Packages.edu.emory.mathcs.restoretools.iterative.IterativeEnums);\nimportClass(Packages.edu.emory.mathcs.restoretools.Enums);\nvar options = new WPLOptions(); options.setNormalize(true);\nvar frgt = new WPLFloatIterativeDeconvolver3D(WindowManager.getImage(" + d2s(img_id, 0) + "), WindowManager.getImage(" + d2s(psf_id[b], 0) + "), IterativeEnums.BoundaryType.REFLEXIVE, IterativeEnums.ResizingType.AUTO, Enums.OutputType.FLOAT, " + d2s(deconvolution_iterations,0) +", false, options ).deconvolve().show();\n";
						}
						//print (run_script);
						if ( deconvolve_yes ) {
							eval( "js", run_script );  //throw away result
						} else {
							run("Duplicate...", "duplicate");
						}
						next_img_id = getImageID(); //get new image ID
						
						//close original image
						selectImage(img_id);
						close();
						img_id = next_img_id;
						selectImage(img_id); //go back to deconvolved image
						
						//adjust properties then save
						setMinAndMax(0, 65535); //clip the upper bounds of signal (lower bounds already clipped for nonnegativity)
						run("16-bit");
						run("Properties...", voxel_definition_output );
						delete_result = File.delete( outputDirectory + file_sep + fileoutname[b] + ".tif" ); //throw away result
						saveAs("Tiff", outputDirectory + file_sep + fileoutname[b] + ".tif" );
	
						//close image
						close();
					
					} else {
						//one round of deconvolution
						
						//okay write TIF file
						delete_result = File.delete( outputDirectory + file_sep + fileoutname[b] + ".tif" ); //throw away result
						
						saveAs("Tiff", outputDirectory + file_sep + fileoutname[b] + ".tif" );
						close();
						
						//close original image
						selectImage(img_id);
						close();
					}
					
					//close PSF
					selectImage(psf_id[b]);
					close();
					
					//write to logfile
					File.append( file + ": angle " + string_angles[a] + ", time " + d2s(t+1,0) + " of " + d2s(max_time,0) + ", channel " + d2s(b,0) + " >> " + outputDirectory + file_sep + fileoutname[b] + ".tif", outputDirectory + file_sep + "CZI_to_TIF_with_deconvolution.txt" );
				}

				//print ( XXXX );
			}
			//Close original concatenated stack
			close("C*");   
			this_time++; //increment time counter
		}
		//close();
		Ext.close();
		call("java.lang.System.gc");
		
		//print( "file: " + file +", " + d2s(i+1,0) + " of " + processList.length );
	}
}

function main_series_tif_to_anaglyphs(directory,outputDirectory) {
	default_degree_separation = 5; //this could be user-changable, but 5 is usually a good number
	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);

	//Count the maximum number of positions and slices in dataset
	run("Bio-Formats Macro Extensions");

	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
	
	//process_view_list = newArray(0);
	
	for (i=0; i<fileList.length; i++) {
		if ((startsWith(fileList[i], "LSFM_T") || startsWith(fileList[i], "Fused_T")) && endsWith(fileList[i], ".tif")) {

			//now, remove .tif from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
			
			//ignore PSF files
			if ( endsWith(filepaddedname, "_psf") ) {
				continue;
			}
			
			//now, segment filename into TXXXX_CX_AXXX and then create view setup with "CX_AXXX"
			name_ext = split(filepaddedname,"_");
			if ( name_ext.length > 3 ) {
				this_view = name_ext[2] + "_" + name_ext[3];
			} else if ( name_ext.length > 2 ) {
				this_view = name_ext[2];
			} else {
				this_view = "V00";
			}
			
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				view_setups = Array.concat( view_setups, this_view );
				
			}
			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] );
		}
	}
	
	Array.sort(view_setups); //sort files by their intended name, not true name
	File.makeDirectory(directory + file_sep + "Anaglyphs");
	
	max_channels = 0;
    max_time = 0;
	view_sizeZ = newArray(view_setups.length);
	Array.fill(view_sizeZ,0);
	
	//okay, first things first, need to figure out how big each stack should be
	setBatchMode(true);
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		outfileList = newArray(0);
		max_sizeZ = 0;
		
		for (i=0; i<process_file_list.length; i++) {
		//for (i=0; i<1; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			file = directory + file_sep + name_ext[1];
			processList = Array.concat( processList, file );
			outfileList = Array.concat( outfileList, directory + file_sep + "Anaglyphs" + file_sep + name_ext[1] );
			
			Ext.setId(file);
			Ext.getSizeZ(sizeZ);
			//print( "   file: " + name_ext[1] + " sizeZ = " + d2s(sizeZ,0) );
			if ( sizeZ > max_sizeZ ) {
				max_sizeZ = sizeZ;
			}
			Ext.close();
			
		}
		
		//determine the angle "degree separation" automatically based on Z-position of mean signal in the middle of the stack
		degree_separation = default_degree_separation; 
		run("TIFF Virtual Stack...", "open=[" + processList[floor(processList.length/2)] + "]"); //open an image in the middle of time series
		OrigStack = getTitle();
		//run("Reslice [/]...", "output=2.000 start=Top rotate avoid"); //get ready to find Z center of mass, first reslice
		getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
		run("Size...", "width=" + round(dim_width/2) + " height=" + round(dim_height/2) + " depth=" + round(dim_slices/2) + " constrain average interpolation=Bilinear");
		run("Reslice [/]...", "start=Left rotate avoid");
		run("Z Project...", "projection=[Max Intensity]");
		run("Enhance Contrast...", "saturated=0.001 normalize");
		run("8-bit");
		ctr_of_mass_z = NaN;
		run("Set Measurements...", "area center redirect=None decimal=3");
		run("Measure"); //measure center of mass
		row_num = nResults() - 1;
		ctr_of_mass_z = getResult("XM", row_num);
		getVoxelSize(vox_width, vox_height, vox_depth, vox_unit); //get voxel and image scaling factors
		getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
		if ( isNaN(ctr_of_mass_z) ) {
			ctr_of_mass_z = 0;
		} else {
			ctr_of_mass_z /= vox_width;
			ctr_of_mass_z /= dim_width;
			
			//now, scale things up, where 0->45deg, 1->2deg (y=45 -43x)
			if ( ctr_of_mass_z > 1 || ctr_of_mass_z < 0 ) {
				//do nothing; leave as default
			} else {
				degree_separation = 45 - (43 * ctr_of_mass_z);
			}
		}		
		
		projection_options = "projection=[Brightest Point] axis=Y-Axis";
		projection_options += " initial=-" + d2s(degree_separation/2,2);
		projection_options += " total=" + d2s(degree_separation,2);
		projection_options += " rotation=" + d2s(degree_separation,2);
		projection_options += " interpolate";		
		
		print( "View setup: " + view_setups[m] + " sizeZ = " + d2s(max_sizeZ,0) + " center of mass Z: " + d2s(ctr_of_mass_z,4) + " degree separation:" + d2s(degree_separation,2) );
		close("*");
		//exit();
		call("java.lang.System.gc");
		
		for (i=0; i<processList.length; i++) {
			//run("TIFF Virtual Stack...", "open=[" + processList[i] + "]");
			open(processList[i]);
			//OrigStack = getTitle();
			
			bits = bitDepth();
			
			if ( bitDepth == 8 ) {
				//do nothing, already 8-bit
			} else {
				OrigStack = getTitle();
				getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
				Stack.getStatistics( area, mean, min, max, std );
				new_max = mean + ((10*std) + (2500*log(std)))/2;
				new_min = mean-(0.5*std);			
				
				new_max = mean+(10*std);
				new_min = mean-(0.5*std);
				if ( new_min < min ) {
					new_min = min;	
				}
				if ( new_max > max ) {
					new_max = max;
				}
				setMinAndMax(new_min,new_max);
				run("Apply LUT", "stack");			
				run("Duplicate...", "duplicate");
				NewStack = getTitle();
				
				selectWindow(OrigStack);
				close();
				selectWindow(NewStack);
				run("Gamma...", "value=0.80 stack");
				run("Enhance Contrast...", "saturated=0.001 normalize process_all use");			
				run("8-bit");
			}
			OrigStack = getTitle();
			run( "3D Project...", projection_options );
			NewStack_12 = getTitle();
			run("Make Substack...", "  slices=2-2");
			NewStack_3 = getTitle();
			run("Concatenate...", "open image1=["+NewStack_12+"] image2=["+NewStack_3+"] image3=[-- None --]");
			run("Make Composite", "display=Composite");
			run("Stack to RGB");
			
			saveAs("Tiff", outfileList[i] );
			close(); 
			
			if ( isOpen(OrigStack) ) {
				close(OrigStack);
			}
			if ( isOpen(NewStack_12) ) {
				close(NewStack_12);
			}
			if ( isOpen(NewStack_3) ) {
				close(NewStack_3);
			}
			
			Ext.close();
			call("java.lang.System.gc");
		}
	}
}

function main_time_series_tif_to_mip(directory,outputDirectory) {
	
	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);
	
	//Count the maximum number of positions and slices in dataset
	run("Bio-Formats Macro Extensions");


	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
		
	for (i=0; i<fileList.length; i++) {
		if (startsWith(fileList[i], "LSFM_") && endsWith(fileList[i], ".tif")) {

			//now, remove .tif from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
			
			//ignore PSF files
			if ( endsWith(filepaddedname, "_psf") ) {
				continue;
			}
			
			//now, segment filename into TXXXX_CX_AXXX and then create view setup with "CX_AXXX"
			name_ext = split(filepaddedname,"_"); //default for time series runs
			if (startsWith(fileList[i], "LSFM__")) { //output for large stack deconvolution runs contains __
				name_ext = split(filepaddedname,"(__)");
			}
			this_view = name_ext[2] + "_" + name_ext[3];
			
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				view_setups = Array.concat( view_setups, this_view );
				
			}
			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] );
		}
	}
	
	Array.sort(view_setups); //sort files by their intended name, not true name
	
	File.makeDirectory(directory + file_sep + "MIPs");
	
	max_channels = 0;
    max_time = 0;
	view_sizeZ = newArray(view_setups.length);
	Array.fill(view_sizeZ,0);
	
	//okay, first things first, need to figure out how big each stack should be
	setBatchMode(true);
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		outfileList = newArray(0);
		max_sizeZ = 0;
		
		for (i=0; i<process_file_list.length; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			file = directory + file_sep + name_ext[1];
			processList = Array.concat( processList, file );
			outfileList = Array.concat( outfileList, directory + file_sep + "MIPs" + file_sep + name_ext[1] );
			
			Ext.setId(file);
			Ext.getSizeZ(sizeZ);
			if ( sizeZ > max_sizeZ ) {
				max_sizeZ = sizeZ;	
			}
			Ext.close();
		}
		print( "View setup: " + view_setups[m] + " sizeZ = " + d2s(max_sizeZ,0) );
		//exit();
		
		Ext.close();
		call("java.lang.System.gc");
		
		for (i=0; i<processList.length; i++) {
			//run("TIFF Virtual Stack...", "open=[" + processList[i] + "]");
			open(processList[i]);
			OrigStack = getTitle();
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("Tiff", outfileList[i] );
			close(); 
			if ( isOpen(OrigStack) ) {
				close(OrigStack);
			}
			call("java.lang.System.gc");
		}
	}
}


function main_klb_to_mip(directory) {
	
	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);

	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
	
	//process_view_list = newArray(0);
	Array.sort(fileList); //very critical step, can be done here or done to array processList below
	
	for (i=0; i<fileList.length; i++) {
		if (startsWith(fileList[i], "t0") && endsWith(fileList[i], ".klb")) {

			//now, remove .klb from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);

			//now, segment filename into TXXXX_SX and then create view setup with "SX"
			name_ext = split(filepaddedname,"_");
			this_view = name_ext[1];
			
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				view_setups = Array.concat( view_setups, this_view );
				
			}

			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] + "///" + filepaddedname );
		}
	}
	
	Array.sort(view_setups); //sort files by their intended name, not true name

	max_channels = 0;
    max_time = 0;
	
	//set up default parameters, and ask user about options for what we are going to do
	do_full_zmip = true;
	do_half_zmips = true;
	intellegent_halfway_finder = false;
	oblique_proj_mean_cleanup = true;
	do_oblique_proj_1 = true;
	oblique_proj_angle_1 = -45;
	do_oblique_proj_2 = true;
	oblique_proj_angle_2 = 45;
	do_oblique_proj_3 = true;
	oblique_proj_angle_3 = -45;
	do_anaglyph = true;
	do_anaglyph_oblique_proj_1 = true;
	do_anaglyph_oblique_proj_2 = true;
	do_anaglyph_oblique_proj_2 = true;
	default_degree_separation = 5;
	Dialog.create("KLB processing settings...");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("KLB z-MIP settings...",16,"Black");
	Dialog.addCheckbox("Do full z-MIPs", do_full_zmip);
	Dialog.addCheckbox("Do two halves z-MIPs", do_half_zmips);
	Dialog.addCheckbox("Determine z-halfway from image features", intellegent_halfway_finder);
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("KLB oblique MIP settings...",16,"Black");
	Dialog.addCheckbox("Subtract mean projection for cleanup", oblique_proj_mean_cleanup);
	Dialog.addCheckbox("Do oblique MIP #1 (y-axis)", do_oblique_proj_1);	
	Dialog.addNumber("Oblique MIP #1 angle:", oblique_proj_angle_1);
	Dialog.addCheckbox("Do oblique MIP #2 (y-axis)", do_oblique_proj_2);	
	Dialog.addNumber("Oblique MIP #2 angle:", oblique_proj_angle_2);
	Dialog.addCheckbox("Do oblique MIP #3 (x-axis)", do_oblique_proj_3);	
	Dialog.addNumber("Oblique MIP #3 angle:", oblique_proj_angle_3);	
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("KLB anaglyph settings...",16,"Black");
	Dialog.addCheckbox("Do full z-anaglyph", do_anaglyph);	
	Dialog.addCheckbox("Do oblique anaglyph #1 (as above)", do_anaglyph_oblique_proj_1);	
	Dialog.addCheckbox("Do oblique anaglyph #2 (as above)", do_anaglyph_oblique_proj_2);	
	Dialog.addCheckbox("Do oblique anaglyph #3 (as above)", do_anaglyph_oblique_proj_2);		
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("KLB custom partial MIP settings...",16,"Black");
    Dialog.addString("Partial z-MIP #1 (slice range):", "###-###", 7);
    Dialog.addString("Partial z-MIP #2 (slice range):", "###-###", 7);

	//Dialog.addNumber("Degrees separation:", default_degree_separation);	
	Dialog.show();
	do_full_zmip = Dialog.getCheckbox();
	do_half_zmips = Dialog.getCheckbox();
	intellegent_halfway_finder = Dialog.getCheckbox();
	oblique_proj_mean_cleanup = Dialog.getCheckbox();
	do_oblique_proj_1 = Dialog.getCheckbox();
	oblique_proj_angle_1 = Dialog.getNumber();
	do_oblique_proj_2 = Dialog.getCheckbox();
	oblique_proj_angle_2 = Dialog.getNumber();
	do_oblique_proj_3 = Dialog.getCheckbox();
	oblique_proj_angle_3 = Dialog.getNumber();	
	do_anaglyph = Dialog.getCheckbox();
	do_anaglyph_oblique_proj_1 = Dialog.getCheckbox();
	do_anaglyph_oblique_proj_2 = Dialog.getCheckbox();
	do_anaglyph_oblique_proj_3 = Dialog.getCheckbox();	
	do_zmip1 = Dialog.getString();
	do_zmip2 = Dialog.getString();

	zmip1_range = split(do_zmip1, "-");
	zmip2_range = split(do_zmip2, "-");

	zmip1_range_int = newArray(2);
	zmip2_range_int = newArray(2);
	
	zmip1_range_int[0] = parseInt(zmip1_range[0]);
	zmip1_range_int[1] = parseInt(zmip1_range[1]);
	zmip2_range_int[0] = parseInt(zmip2_range[0]);
	zmip2_range_int[1] = parseInt(zmip2_range[1]);

	//explore dataset.xml files to get unit information
	X_voxel_micron = 0;
	Y_voxel_micron = 0;
	Z_voxel_micron = 0;
	voxel_units = "";
	input_MVL_lines = newArray(100);
	//print( "looking for... " + directory+ "dataset_klb.xml" );
	if ( File.exists(directory+ "dataset.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb_s0.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb_s0.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb_s1.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb_s1.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	}
	unit_line = -1;
	size_line = -1;
	for ( l=0; l<input_MVL_lines.length; l++ ) {
		if (indexOf(input_MVL_lines[l], "<size>") >= 0) {
			size_line = l;
		} else if (indexOf(input_MVL_lines[l], "<unit>") >= 0) {
			unit_line = l;
			size_line = -1;
		}
		if ( unit_line >= 0 && size_line >= 0 ) {
			break;
		}
	}	
	if ( unit_line >= 0 && size_line >= 0 ) {
		//print("got lines");
		start_index = indexOf(input_MVL_lines[size_line], "<size>");
		subline = substring(input_MVL_lines[size_line], start_index + 6 );
		stop_index = indexOf(subline, "</size>" );
		size_line_data = substring(subline, 0, stop_index );
		sizes_parsed = split( size_line_data, " " );
		X_voxel_micron = parseFloat(sizes_parsed[0]);
		Y_voxel_micron = parseFloat(sizes_parsed[1]);
		Z_voxel_micron = parseFloat(sizes_parsed[2]);
		
		start_index = indexOf(input_MVL_lines[unit_line], "<unit>");
		subline = substring(input_MVL_lines[unit_line], start_index + 6 );
		stop_index = indexOf(subline, "</unit>" );
		voxel_units = substring(subline, 0, stop_index );
	}
	
	File.makeDirectory(directory + file_sep + "MIPs" );
	
	//iterate over channels
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		
		for (i=0; i<process_file_list.length; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			processList = Array.concat( processList, name_ext[2] );
		}
		call("java.lang.System.gc");
		
		
		//determine the angle "degree separation" automatically based on Z-position of mean signal in the middle of the stack
		degree_separation = default_degree_separation; 
		projection_options = "projection=[Brightest Point]";
		intelligent_mid_slice = 0; // used for half z-MIPs		
		if ( do_anaglyph || do_anaglyph_oblique_proj_1 || do_anaglyph_oblique_proj_2 || do_anaglyph_oblique_proj_3 || intellegent_halfway_finder ) {
			setBatchMode(false); //KLB does not work with batch mode
			
			//open an image in the middle of time series
			run("KLB...", "open=[" + directory + processList[floor(processList.length/2)] + ".klb]");
			VirtStack = getImageID();//("window.title");
			
			setBatchMode(true);
			selectImage(VirtStack);
			run("Duplicate...", "duplicate");
			OrigStack = getTitle();//("window.title");
			selectImage(VirtStack);
			close();
			selectImage(OrigStack);
			
			//run("Reslice [/]...", "output=2.000 start=Top rotate avoid"); //get ready to find Z center of mass, first reslice
			getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
			run("Size...", "width=" + round(dim_width/2) + " height=" + round(dim_height/2) + " depth=" + round(dim_slices/2) + " constrain average interpolation=Bilinear");
			run("Reslice [/]...", "start=Left rotate avoid");
			VirtStack = getImageID();
			
			run("Z Project...", "projection=[Max Intensity]");
			NewProj = getImageID();
			
			run("Enhance Contrast...", "saturated=0.001 normalize");
			
			run("8-bit");
			ctr_of_mass_z = NaN;
			run("Set Measurements...", "area center redirect=None decimal=3");
			run("Measure"); //measure center of mass
			row_num = nResults() - 1;
			ctr_of_mass_z = getResult("XM", row_num);
			getVoxelSize(vox_width, vox_height, vox_depth, vox_unit); //get voxel and image scaling factors
			getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
			if ( isNaN(ctr_of_mass_z) ) {
				ctr_of_mass_z = 0;
			} else {
				ctr_of_mass_z /= dim_width;
				if ( intellegent_halfway_finder ) {
					intelligent_mid_slice = floor((dim_slices*ctr_of_mass_z) + 0.5);
				}					
				//ctr_of_mass_z /= vox_width;
				
				//now, scale things up, where 0->45deg, 1->2deg (y=45 -43x)
				if ( ctr_of_mass_z > 1 || ctr_of_mass_z < 0 ) {
					//do nothing; leave as default
				} else {
					degree_separation = 7 - (6.2 * ctr_of_mass_z);
				}
			}		
			projection_options += " total=" + d2s(degree_separation,2);
			projection_options += " rotation=" + d2s(degree_separation,2);
			//projection_options += " interpolate";		
			
			//print( "View setup: " + view_setups[m] + " center of mass Z: " + d2s(ctr_of_mass_z,4) + " degree separation:" + d2s(degree_separation,2) );
			
			if (isOpen(OrigStack)) {
				close(OrigStack);
			}
			if (isOpen(VirtStack)) {
				selectImage(VirtStack);
				close();	
			}
			if (isOpen(NewProj)) {
			selectImage(NewProj);
			close();	
			}
			if (isOpen("Results")) {
				selectWindow("Results"); 
				run("Close" );
			}
			//exit();
		}
		
		call("java.lang.System.gc");		
		
		for (i=0; i<processList.length; i++) {
			setBatchMode(false); //KLB does not work with batch mode
			
			run("KLB...", "open=[" + directory + processList[i] + ".klb]");
			VirtStack = getImageID();//("window.title");
			
			setBatchMode(true);
			selectImage(VirtStack);
			run("Duplicate...", "duplicate");
			OrigStack = getImageID();//("window.title");
			selectImage(VirtStack);
			close();
			
			
			//set up units
			selectImage(OrigStack);
			if ( !(voxel_units == "") ) {
				Stack.setXUnit(voxel_units);
				run("Properties...", "pixel_width="+d2s(X_voxel_micron,8)+" pixel_height="+d2s(Y_voxel_micron,8)+" voxel_depth="+d2s(Z_voxel_micron,8) );
			}
			
			
			//full z-MIP
			if ( do_full_zmip ) {
				run("Z Project...", "projection=[Max Intensity]");
				saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_complete" );
				close();
				selectImage(OrigStack);
			}

			//custom z-MIP
			if ( !(isNaN(zmip1_range_int[0])) && !(isNaN(zmip1_range_int[1])) && zmip1_range_int[0] > 0 && zmip1_range_int[1] > zmip1_range_int[0] ) {
				run("Z Project...", "start=" + d2s(zmip1_range_int[0],0) + " stop=" + d2s(zmip1_range_int[1],0) + " projection=[Max Intensity]" );
				saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_zmip" + zmip1_range_int[0] + "to" + zmip1_range_int[1] );
				close();
				selectImage(OrigStack);
			}
			if ( !(isNaN(zmip2_range_int[0])) && !(isNaN(zmip2_range_int[1])) && zmip2_range_int[0] > 0 && zmip2_range_int[1] > zmip2_range_int[0] ) {
				run("Z Project...", "start=" + d2s(zmip2_range_int[0],0) + " stop=" + d2s(zmip2_range_int[1],0) + " projection=[Max Intensity]" );
				saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_zmip" + zmip2_range_int[0] + "to" + zmip2_range_int[1] );
				close();
				selectImage(OrigStack);
			}

			//partial z-MIPs
			if ( do_half_zmips ) {
				//get per-slice pixel values
				getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
				mid_slice = floor((dim_slices+0.5)/2); //default
				/*
				if ( i == 0 ) {
					last_mid_slice = mid_slice;
				}
				
				if ( intellegent_halfway_finder ) {
					first_slice = 1;
					last_slice = dim_slices;
					last_mean = 0;
					for (ss=1; ss<=dim_slices; ss++ ) { //decrement slices to fill array with max pixel values for each slice
						setSlice(ss);
						getStatistics( area, mean );
						if ( mean > 1 && mean > last_mean ) {
							first_slice = ss;
							break;
						} else {
							last_mean = mean;
						}
					}
					for (ss=dim_slices; ss>first_slice; ss-- ) { //decrement slices to fill array with max pixel values for each slice
						setSlice(ss);
						getStatistics( area, mean );
						if ( mean > 1 && mean > last_mean ) {
							last_slice = ss;
							break;
						} else {
							last_mean = mean;
						}
					}
					
					if ( first_slice +2 < last_slice ) {
						mid_slice = first_slice + floor(((last_slice-first_slice)+0.5)/2)
					}
					
					//check to make sure we aren't wildly swinging around the image from timepoint to timepoint, and restrict midway movement to no more than 5% of the total stack thickness
					if ( abs(mid_slice - last_mid_slice) > dim_slices/15 ) {
						if ( mid_slice > last_mid_slice ) { //current midway is greater than last, add 5%
							mid_slice = last_mid_slice + floor(dim_slices/20);
						} else { //current midway is greater than last, subtract 5%
							mid_slice = last_mid_slice - floor(dim_slices/20);
						}
					}
					
					//these should not happen
					if ( mid_slice >= dim_slices ) {
						mid_slice = dim_slices -1;
					} else if ( mid_slice < 2 ) {
						mid_slice = 2;	
					}
					
					//store for next timepoint
					last_mid_slice = mid_slice;
				}
				*/
				run("Z Project...", "start=1 stop=" + d2s(mid_slice,0) + " projection=[Max Intensity]");
				saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_front_half" );
				close(); 
				
				selectImage(OrigStack);
				run("Z Project...", "stop=" + d2s(dim_slices,0) + " start=" + d2s(mid_slice,0) + " projection=[Max Intensity]");
				saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_back_half" );
				close(); 
			
				selectImage(OrigStack);
			}
			
			//oblique MIPs
			getVoxelSize(vox_width, vox_height, vox_depth, vox_unit);
			if ( do_oblique_proj_1 || do_oblique_proj_2 || do_oblique_proj_3 ) {
				for ( ii=1; ii<=3; ii++ ) {
					my_proj_angle = 0;
					my_axis = "Y";
					if ( ii==1 && do_oblique_proj_1 ) {
						my_proj_angle = oblique_proj_angle_1;
					} else if ( ii==2 && do_oblique_proj_2 ) {
						my_proj_angle = oblique_proj_angle_2;
					} else if ( ii==3 && do_oblique_proj_3 ) {
						my_proj_angle = oblique_proj_angle_3;
						my_axis = "X";						
					}
					if ( my_proj_angle == 0 ) {
						continue;
					}
					
					run("3D Project...", "projection=[Brightest Point] axis="+my_axis+"-Axis slice="+d2s(vox_depth,4)+" initial="+d2s(my_proj_angle,4)+" total=0 rotation=10 lower=1 upper=255 opacity=0 surface=100 interior=50");
					if (oblique_proj_mean_cleanup) {
						rename("Image" + d2s( floor(random() * 1000000 ), 0 ));
						Brightest = getTitle();
						
						selectImage(OrigStack);
						run("3D Project...", "projection=[Mean Value] axis="+my_axis+"-Axis slice="+d2s(vox_depth,4)+" initial="+d2s(my_proj_angle,4)+" total=0 rotation=10 lower=1 upper=255 opacity=0 surface=100 interior=50"); 
						rename("Image" + d2s( floor(random() * 1000000 ), 0 )); 
						Mean = getTitle();
						
						if ( Mean == Brightest ) { //is this collision-proof?
							rename( "MeanImage" );
							Mean = getTitle();
						}
						
						imageCalculator("Subtract create",Brightest,Mean);
						saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_oblique" + d2s(ii,0) );
						close();
						
						selectImage(Brightest);
						saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_ob"+d2s(ii,0)+"_max" );
						close();
						
						selectImage(Mean);
						saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_ob"+d2s(ii,0)+"_mean" );
						close();					
					} else {
						saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_oblique" + d2s(ii,0) );
						close(); 
					}
				}
				selectImage(OrigStack);
			}
			
			if ( do_anaglyph || do_anaglyph_oblique_proj_1 || do_anaglyph_oblique_proj_2 || do_anaglyph_oblique_proj_3 ) {
				for ( ii=0; ii<=3; ii++ ) {
					my_proj_angle = 0;
					my_axis = "Y";
					if ( ii == 0 ) {
						if (!( do_anaglyph )) {
							continue;
						}
					} else if ( ii==1 && do_anaglyph_oblique_proj_1 ) {
						my_proj_angle = oblique_proj_angle_1;
					} else if ( ii==2 && do_anaglyph_oblique_proj_2 ) {
						my_proj_angle = oblique_proj_angle_2;
					} else if ( ii==3 && do_anaglyph_oblique_proj_3 ) {
						my_proj_angle = oblique_proj_angle_3;
						my_axis = "X";	
					}
					if ( my_proj_angle == 0 && ii > 0) {
						continue;
					}
					
					run( "3D Project...", projection_options + " axis=" + my_axis + "-Axis slice="+d2s(vox_depth,4) + " initial=" + d2s(my_proj_angle-(degree_separation/2),2) );
					NewStack_12 = getTitle();
					run("Make Substack...", "  slices=2-2");
					NewStack_3 = getTitle();
					run("Concatenate...", "open image1=["+NewStack_12+"] image2=["+NewStack_3+"] image3=[-- None --]");
					run("Make Composite", "display=Composite");
					NewProj = getImageID();
					run("Stack to RGB");
				
					saveAs("Tiff", directory + file_sep + "MIPs" + file_sep + processList[i] + "_anaglyph" + d2s(ii,0) );
					close(); 
					if (isOpen(NewProj)) {
						selectImage(NewProj);
						close();	
					}
					//setBatchMode("exit and display");
					//exit();
				}
			}
			
			
			close();
			call("java.lang.System.gc");
		}
	}
	setBatchMode(false);
}

function main_tifseries_16bit_to_8bit(folder_in_batch,directory,outputDirectory) {	

	fileList = newArray(0);
	name_ext = newArray(0);
	
	if ( folder_in_batch ) {
		if ( directory == "" || !File.exists(directory) ) {
			//Ask user to choose the input and output directories
			directory = getDirectory("Choose input directory");
		}
		fileList = getFileList(directory);

	} else {
		path = File.openDialog("Select a File");
  		directory = File.getParent(path) + file_sep;
  		name = File.getName(path);
		
		//now, remove .tif from end of filename
		/*name_ext = split(name,".");
		filepaddedname = "";
		for (n=0; n<name_ext.length-1; n++) {
			filepaddedname += name_ext[n];
		}*/
		filepaddedname = substring(name, 0, name.length-4);

		//fileList = Array.concat(fileList,name);
		fileList_preview = getFileList(directory);
		name_ext_preview = split(filepaddedname,"(__)");
		
		for (i=0; i<fileList_preview.length; i++) {
			if (startsWith(fileList_preview[i], "LSFM__") && (endsWith(fileList_preview[i], ".tif")||endsWith(fileList_preview[i], ".zip"))) {
				/*name_ext = split(fileList_preview[i],".");
				filepaddedname = "";
				for (n=0; n<name_ext.length-1; n++) {
					filepaddedname += name_ext[n];
				}*/
				filepaddedname = substring(fileList_preview[i], 0, fileList_preview[i].length-4);
				
				name_ext = split(filepaddedname,"(__)");
				//print( "looking into adding " + filepaddedname + " to jobs list...");
				//print( "  " + name_ext[0] + " == " + name_ext_preview[0]  + " && " +  name_ext[1]  + " == " +  name_ext_preview[1] );
				if (name_ext[0] == name_ext_preview[0] && name_ext[1] == name_ext_preview[1] ) {
					fileList = Array.concat(fileList,fileList_preview[i]);
					//print( "  ..added");
				}
			}
		}
	}
	
	File.makeDirectory(directory + "Filtered");
	File.makeDirectory(directory + "Tmp");
	if ( generate_mips_when_filtering ) {
		File.makeDirectory(directory + "Filtered" + file_sep + "MIPs" );
	}
	this_view = "";
	
	//Count the maximum number of positions and slices in dataset
	run("Bio-Formats Macro Extensions");

	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
	
	//print( "about to iterate on filelist, " + fileList.length );
	
	for (i=0; i<fileList.length; i++) {
		if (startsWith(fileList[i], "LSFM_") && (endsWith(fileList[i], ".tif")||endsWith(fileList[i], ".zip"))) {

			//now, remove .tif from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
			
			//ignore PSF files
			if ( endsWith(filepaddedname, "_psf") ) {
				continue;
			}
			
			//now, segment filename into TXXXX_CX_AXXX and then create view setup with "CX_AXXX"
			if (startsWith(fileList[i], "LSFM__")) { //output for large stack deconvolution runs contains __
				name_ext = split(filepaddedname,"(__)");
				this_view = name_ext[2]; //group by angles only
			} else {
				name_ext = split(filepaddedname,"_"); //default for time series runs
				this_view = name_ext[2] + "_" + name_ext[3];
			}
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				//print( "adding view setup: " + this_view );
				view_setups = Array.concat( view_setups, this_view );
				
			}
			//print( "adding file " + filepaddedname + " as part of view setup " + this_view );
			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] + "///" + filepaddedname );
		}
	}
	
	Array.sort(view_setups); //sort files by their intended name, not true name
	
	max_channels = 0;
    max_time = 0;
	view_sizeZ = newArray(view_setups.length);
	Array.fill(view_sizeZ,0);
	
	//okay, first things first, need to figure out how big each stack should be
	setBatchMode(true);
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		processList_mean = newArray(0);
		processList_std = newArray(0);
		processList_max = newArray(0);
		max_sizeZ = 0;
		
		for (i=0; i<process_file_list.length; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			//file = directory + name_ext[1];
			processList = Array.concat( processList, name_ext[2] );

			//run("TIFF Virtual Stack...", "open=[" + directory + name_ext[1] + "]");
			open(directory + name_ext[1]);
			
			//take care of dimension issue for black image padding below
			getDimensions(dim_width, dim_height, dim_channels, sizeZ, dim_frames);		
			if ( sizeZ > max_sizeZ ) {
				max_sizeZ = sizeZ;	
			}			

			//subject average image intensity across the stack
			orig_stack = getTitle();
			if (startsWith(name_ext[2], "LSFM__")) { //output for large stack deconvolution runs contains __, and don't using sliding parabola here
				if ( fixed_rolling_ball_radius > 10 ) {
					run("Subtract Background...", "rolling="+d2s(fixed_rolling_ball_radius,0)+" stack");				
				}
			} else if ( rolling_ball_radius > 10 ) {
				run("Subtract Background...", "rolling="+d2s(rolling_ball_radius,0)+" sliding disable stack");
			}			
			
			if ( deconvolution_subtract_camera_noise ) {
				if (startsWith(name_ext[2], "LSFM__")) { 
					run("Z Project...", "projection=[Min Intensity]");
					run("Gamma...", "value=0.9");
				} else {
					run("Z Project...", "projection=[Average Intensity]");
					run("Gamma...", "value=0.85");
				}
				
				average_img = getTitle();
				imageCalculator("Subtract stack", orig_stack,average_img);
				selectWindow(average_img);
				close();
			}
			
			//remove bright blobs	
			if ( fixed_precipitate_removal || fixed_blob_removal ) {
				bright_blob_remover_this_image(fixed_precipitate_removal,fixed_blob_removal);
			}
			
			//okay, now grab final data for compiling for modification of images the second time around
			Stack.getStatistics( area, mean, min, max, std );
			processList_mean = Array.concat( processList_mean, mean );
			processList_std = Array.concat( processList_std, std );
			processList_max = Array.concat( processList_max, max );
			//print( "item " + i + ", mean: " + mean + ", stdev: " + std );
			
			saveAs("Tiff", directory + "Tmp" + file_sep + name_ext[2] + ".tif" );
			close();
			call("java.lang.System.gc");			
			
		}
		
		Array.getStatistics(processList_mean, min, max, grand_mean, grand_mean_std);
		total_sum_squares = 0;
		combine_stdev = 0;
		for (i=0; i<processList.length; i++) {
			total_sum_squares = processList_mean[i] - grand_mean;
			combine_stdev += (processList_std[i] * processList_std[i]) + (total_sum_squares * total_sum_squares);
		}
		combine_stdev = sqrt(combine_stdev/processList.length);
		Array.getStatistics(processList_max, min, max, max_mean, max_std);
		
		new_min = 10 + grand_mean-(0.5*combine_stdev);		
		new_max = max_mean + max_std;
		if (max < new_max) {	
			new_max = max;
		}

		print( "Processing view setup: " + view_setups[m] + " sizeZ = " + d2s(max_sizeZ,0) + "..." );

		for (i=0; i<processList.length; i++) {
			open(directory + "Tmp" + file_sep + processList[i] + ".tif" );
			getVoxelSize(vox_width, vox_height, vox_depth, vox_unit);
			//run("Duplicate...", "duplicate");
			orig_stack = getTitle();

			//run("Gamma...", "value=0.75 stack");
			//run("Enhance Contrast...", "saturated=0.0001 normalize process_all use");
			//go ahead and modify histogram here, convert to 8-bit
			Stack.getStatistics( area, mean, min, max, std );
			if ( new_min > min ) { //if desired minimim value is actually less than this particular frame's minimum, don't add background intensity just keep this frame's minimum
				print( "   ..For " + processList[i] + ": adjusting timepoint minimum " + d2s(min,0) + " to series minimum " + d2s(new_min,0) + "." );
				min = new_min;
			} else {
				print( "   ..For " + processList[i] + ": keeping timepoint minimum " + d2s(min,0) + ". Not using series minimum " + d2s(new_min,0) + "." );
			}
			if ( new_max < max ) { //if desired maximum value is less than this frame's maximum, please clip off more intense pixels to max cutoff
				print( "   ..For " + processList[i] + ": adjusting timepoint maximum " + d2s(max,0) + " to series maximum " + d2s(new_max,0) + "." );
				max = new_max;
			} else {
				print( "   ..For " + processList[i] + ": keeping timepoint maximum " + d2s(max,0) + ". Not using series maximum " + d2s(new_max,0) + "." );
			}
			//setPixel(0,0,max);
			//setPixel(1,0,min);			
			setMinAndMax(min,max);
			run("Apply LUT", "stack");
			//run("Enhance Contrast...", "process_all use");
			run("Gamma...", "value=0.75 stack");
			
			if ( convert_to_8_bit ) {
				run("8-bit");
			} else {
				run("16-bit");
			}
			
			getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);

			//make all stacks the same z-depth in this view setup, but only if doing a time series
			if ( !(startsWith(processList[i], "LSFM__")) && folder_in_batch && dim_slices < max_sizeZ ) {
				if ( convert_to_8_bit ) {
					newImage("__AddStack", "8-bit black", dim_width, dim_height, max_sizeZ-dim_slices);
				} else {
					newImage("__AddStack", "16-bit black", dim_width, dim_height, max_sizeZ-dim_slices);
				}
				
				run("Concatenate...", "  title=__ConcatStack image1=__AddStack image2=" + orig_stack + " image3=[-- None --]");
			}
			
			setVoxelSize(vox_width, vox_height, vox_depth, vox_unit);
			saveAs("Tiff", directory + "Filtered" + file_sep + processList[i] + ".tif" );
			delete_result = File.delete( directory + "Tmp" + file_sep + processList[i] + ".tif" ); //throw away result
			
			//Close concatenated stack
			if ( generate_mips_when_filtering ) {
				run("Z Project...", "projection=[Max Intensity]");
				saveAs("Tiff", directory + "Filtered" + file_sep + "MIPs" + file_sep + processList[i] + ".tif" );
				close(); 
			}			
			close();
			
			Ext.close();
			call("java.lang.System.gc");
		}
	}
}

function main_multichannel_to_channels_and_merge () {
	//allow user to select image
	//showMessageWithCancel("Select an image and click OK or cancel to stop processing.");
	//imageId = getImageId();
	//selectImage(imageId);
	imageTitle = getTitle();
	
	getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames); //XY and channels will be important
	
	if ( dim_slices > 1 || dim_frames > 1 ) {
	    print( "The selected image, " + imageTitle + " has too many slices("+d2s(dim_slices,0)+") or frames("+d2s(dim_frames,0)+"), reduce dimensionality to one each." );    
	    exit();
	}
	
	if ( dim_channels < 2 ) {
	    print( "The selected image, " + imageTitle + " has too few channels("+d2s(dim_channels,0)+") and therefore does not require this processing." );    
	    exit();
	}
	
	//now, prompt for user input on channel name and pseudocolor
	channel_titles = newArray(dim_channels);
	channel_colors = newArray(dim_channels);
	//overlay_items = Overlay.size();
	Dialog.create("Enter channel parameters...");
	for (b=1; b<=dim_channels; b++ ) {
	    def_color = "Red";
		def_title = "RFP";
	    if ( b == 2) { def_color = "Green"; def_title = "GFP"; }
	    else if ( b == 3) { def_color = "Blue"; def_title = "DAPI"; }
	    else if ( b == 4) { def_color = "Grays"; def_title = "DAPI"; }    
	    Dialog.addChoice("Channel "+d2s(b,0)+" LUT:", newArray("Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "Grays"),def_color);
	    Dialog.addToSameRow();
	    Dialog.addString("Title:", def_title);
	}
	Dialog.addString("Merge Title:", "Merge");
	//Dialog.addCheckbox("Combine images vertically", false);
	layout_array = newArray("Vertical", "Box", "Horizontal");
	Dialog.addChoice("Montage layout", layout_array, "Box" );
	Dialog.addCheckbox("Auto-level", false);
	Dialog.addNumber("Font size:", floor(dim_height/12));
	Dialog.addNumber("White spacing:", floor(dim_height/60));
	scalebar_array = newArray("None", "Upper Right", "Lower Right", "Upper Left", "Lower Left" );
	Dialog.addChoice("Scale bar in Merge image", scalebar_array, "None");
	Dialog.show();
	for (b=0; b<dim_channels; b++ ) {
	    channel_colors[b] = Dialog.getChoice();
	    channel_titles[b] = Dialog.getString();
	}
	merge_title = Dialog.getString();
	//combine_vertical = Dialog.getCheckbox();
	combine_layout = Dialog.getChoice();
	auto_level = Dialog.getCheckbox();
	chosen_font_size = Dialog.getNumber();
	chosen_spacing = Dialog.getNumber();
	add_scale_bar = Dialog.getChoice();
	
	
	
	//duplicate image upfront
	is_overlay = false;
	if ( Overlay.size > 0 ) {
		Overlay.copy;
		is_overlay = true;
	}
	run("Duplicate...", "duplicate");
	//working_imageId = getImageId();
	working_imageTitle = getTitle();
	//print( "working image title: " + working_imageTitle );
	
	//Split channels and record names of each new image stack
	channelList = newArray(0);
	if ( dim_channels > 1 ) {
		run("Split Channels");
		for (b=1; b<=dim_channels; b++ ) {
			channelList = Array.concat( channelList, "C" + IJ.pad(b,1) + "-" + working_imageTitle ); 
			//print( "new channel: " + channelList[b-1] );
		}
	} else if ( dim_channels == 1 ) {
			channelList = Array.concat( channelList, working_imageTitle );
	} else {
	    print( "Encountered an error splitting channels" );
	    close();
		exit();
	}
	
	//now process each channel image
	for (b=0; b<dim_channels; b++ ) {
		selectWindow(channelList[b]);
		run(channel_colors[b]);
		//print( "top " + d2s(b,0) );
	
		if ( auto_level ) {
			setMinAndMax(0, 65535); //clip the upper bounds of signal (lower bounds already clipped for nonnegativity)
			run("Apply LUT", "stack");
			run("Gamma...", "value=0.80 stack");	
			run("Enhance Contrast...", "saturated=0.001 normalize");
			run("8-bit");
		}
	
		//duplicate to create image for final merge
		run("Duplicate...", "duplicate");
		rename("ch" + d2s(b,0) + "_temp" );
	
		//go back to original image and continue processing individual channel
		selectWindow(channelList[b]);
	
		//draw channel title
		setFont("SansSerif", chosen_font_size, " antialiased");
		if (channel_colors[b]=="Grays") {
			setColor("#dddddd");
		} else {
			setColor(channel_colors[b]);
		}
	    //setBackgroundColor("Black");
	    setBackgroundColor(0, 0, 0);
		drawString(channel_titles[b], floor(chosen_font_size / 3), floor(chosen_font_size * 1.3333) );	
		//print( "midA " + d2s(b,0) );
		//merge to create this channel's image as an RGB instead of single color
		if (channel_colors[b]=="Grays") {
			run("Merge Channels...", "c1=["+channelList[b]+"] c2=["+channelList[b]+"] c3=["+channelList[b]+"]");
			//print ("Proceed Gray");
			rename(channelList[b]);
		} else {
			//create extra black channel for merge to create RGB
			depth = bitDepth();
			if ( depth == 32 ) {
				newImage("blank_temp", "32-bit black", dim_width, dim_height, 1);
			} else if ( depth == 16 ) {
				newImage("blank_temp", "16-bit black", dim_width, dim_height, 1);
			} else if ( depth == 8 ) {
				newImage("blank_temp", "8-bit black", dim_width, dim_height, 1);
			}
			
			//merge to create channel image as an RGB
			run("Merge Channels...", "c1=["+channelList[b]+"] c2=blank_temp create");
			//print( "midC " + d2s(b,0) );
			rename("merge_temp");
			run("Stack to RGB");
			rename(channelList[b]);
			close("merge_temp");		
		}
		//print( "midB " + d2s(b,0) );
	
		
		//print( "bottom " + d2s(b,0) );
	
		if ( is_overlay ) {
			rename("merge_temp");
			Overlay.paste;
			Overlay.flatten;
			rename(channelList[b]);
			close("merge_temp");		
		}
	
		//now pad the image
		setBackgroundColor(255,255,255);
		if ( combine_layout == "Vertical" ) {
			run("Canvas Size...", "width="+d2s(dim_width,0)+" height="+d2s(dim_height+chosen_spacing,0)+" position=Top-Left");
		} else if ( combine_layout == "Horizontal" ) {
			run("Canvas Size...", "width="+d2s(dim_width+chosen_spacing,0)+" height="+d2s(dim_height,0)+" position=Top-Left");
		} else { //box is default
			run("Canvas Size...", "width="+d2s(dim_width+chosen_spacing,0)+" height="+d2s(dim_height+chosen_spacing,0)+" position=Top-Left");
		}
		
		//
		//makeText(channel_titles[b], 16, 16);
		//run("Add Selection...", "stroke="+channel_colors[b]+" fill=#00000000 new");
		//run("Draw", "slice");
		//run("Select None");	
	}
	
	//overlay
	
	//now create final merge
	merge_text_list = "c1=[ch0_temp] c2=[ch1_temp]";
	for (b=2; b<dim_channels; b++ ) {
		merge_text_list = merge_text_list + " c"+d2s(b+1,0)+"=[ch" +d2s(b,0)+ "_temp]";
	}
	run("Merge Channels...", merge_text_list + " create");
	
	//if ( is_overlay ) {
	//	Overlay.copy;
	//}
	rename("final_temp");
	//exit();
	//selectWindow(imageTitle);
	run("Stack to RGB");
	if ( is_overlay ) {
		rename("stack_temp");
		Overlay.paste;
		Overlay.flatten;
		rename("merge_temp");
		close("stack_temp");
	} else {
		rename("merge_temp");
	}
	close("final_temp");
	
	//run("Duplicate...", "duplicate");
	//working_imageTitle = getTitle();
	setFont("SansSerif", chosen_font_size, " antialiased");
	setColor("White");
	drawString(merge_title, floor(chosen_font_size / 3), floor(chosen_font_size * 1.3333) );
	if (add_scale_bar != "None") {
		run("Scale Bar...", "width=100 height="+d2s(chosen_spacing/2+4,0)+" font="+d2s(floor(chosen_font_size/2.8+18),0)+" color=White background=None location=["+add_scale_bar+"]");
	}
	if ( combine_layout == "Box" ) { //box is default
			run("Canvas Size...", "width="+d2s(dim_width+chosen_spacing,0)+" height="+d2s(dim_height+chosen_spacing,0)+" position=Top-Left");
		}
	
	
	
	//now combine all
	if ( combine_layout == "Vertical" ) {
		run("Combine...", "stack1=["+channelList[0]+"] stack2=["+channelList[1]+"] combine");
		rename("combine_temp");
		for (b=2; b<dim_channels; b++ ) {
			run("Combine...", "stack1=combine_temp stack2=["+channelList[b]+"] combine");
			rename("combine_temp");
		}
		run("Combine...", "stack1=combine_temp stack2=merge_temp combine");
	} else if ( combine_layout == "Horizontal" ) {
		run("Combine...", "stack1=["+channelList[0]+"] stack2=["+channelList[1] + "]");
		rename("combine_temp");
		for (b=2; b<dim_channels; b++ ) {
			run("Combine...", "stack1=combine_temp stack2=["+channelList[b] + "]");
			rename("combine_temp");	
		}
		run("Combine...", "stack1=combine_temp stack2=merge_temp");
	} else if ( dim_channels < 5 ) { //box is default; can be 2x2 or 3x3
		//first make temporary rows then combine those
		for (b=0; b<dim_channels; b+=2 ) {
			if (b+1<dim_channels) { //enough for two images to combine
				run("Combine...", "stack1=["+channelList[b]+"] stack2=["+channelList[b+1] + "]");
				//run("Combine...", "stack1=combine_temp stack2=["+channelList[b] + "]");
				rename("combine_temp_" + d2s(b/2,0));	
			} else {
				selectWindow(channelList[b]);
				rename("combine_temp_" + d2s(b/2,0));
			}
		}
		//add merge to end of last combined row
		if( dim_channels == 2 ) {
			selectWindow("merge_temp");
			rename("combine_temp_" + d2s(floor(dim_channels/2),0));	
	
			//now combine rows
			run("Combine...", "stack1=[combine_temp_0] stack2=[combine_temp_1] combine");
		} else {
			run("Combine...", "stack1=[combine_temp_" + d2s(floor((dim_channels-1)/2),0)+ "] stack2=merge_temp");
			rename("combine_temp_" + d2s(floor((dim_channels-1)/2),0));
	
			//now combine rows
			for (b=2; b<dim_channels; b+=2 ) {
				run("Combine...", "stack1=[combine_temp_0] stack2=[combine_temp_"+d2s(b/2,0) + "] combine");
				//run("Combine...", "stack1=combine_temp stack2=["+channelList[b] + "]");
				rename("combine_temp_0");
			}		
		}
		getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
		run("Canvas Size...", "width="+d2s(dim_width-chosen_spacing,0)+" height="+d2s(dim_height-chosen_spacing,0)+" position=Top-Left");
	} else { //if ( dim_channels < 10 ) {
		print( "Not capable of box montage with so many channels.  Use horizontal or vertical montage, then manually crop and re-combine as desired." );
	}
}


function bright_blob_remover_this_image ( precipitate, blob ) {
	
	//remove small precipitates
	Stack.getStatistics( area, mean, min, max, std );
	
	if ( precipitate ) {
		run("Remove Outliers...", "radius=4 threshold="+d2s(mean+std,0)+" which=Bright stack");
	}
	
	if ( !(blob) ) {
		return;
	}
	
	//duplicate stack and get basic data
	orig_id = getTitle();
	run("Duplicate...", "duplicate");
	background_id = getTitle();
	getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
	
	
	//get max Z projection and find cutoff value for top 0.X% of pixels, where X is dependent on series maximum value
	//denominator = -280*log(mean)+2280;
	if ( mean < 100 ) {
		print( "Mean does not meet cutoff for bright blob removal for image " + orig_id );
		close(background_id);
		return;
	}
	denominator = ( -0.144 * max ) + 10015;
	percent_px = floor( dim_width * dim_height / denominator );
	
	run("Z Project...", "projection=[Max Intensity]");
	getHistogram(value,count,256);
	
	close();
	total = 0;
	cutoff = 65535;
	for ( i=255; i>=0; i-- ) {
		total += count[i];
		if ( total > percent_px ) {
			cutoff = value[i];
			break;	
		}
	}
	
	//get per-slice pixel values
	max_px_values = newArray(dim_slices);
	mean_px_values = newArray(dim_slices);
	std_px_values = newArray(dim_slices);
	for (ss=1; ss<=dim_slices; ss++ ) { //decrement slices to fill array with max pixel values for each slice
		setSlice(ss);
		getStatistics( area, mean_px_values[ss-1], min, max_px_values[ss-1], std_px_values[ss-1] );
	}
	Array.getStatistics(mean_px_values, minimum, maximum, average_mean_px_values);
	Array.getStatistics(std_px_values, minimum, maximum, average_std_px_values);
	
	//print( "Average values: avg: " + d2s(average_mean_px_values,0) + ", std:" + d2s(average_std_px_values,0) );
	average_mean_px_values /= average_std_px_values;
	//print( "Average values: avg: " + d2s(average_mean_px_values,4) + ", std:" + d2s(average_std_px_values,0) );
	//average_mean_px_values = 1;
	//iterate over image slices and set up outlier pixels using min/max and apply LUT
	for (ss=0; ss<dim_slices; ss++ ) { //decrement slices to fill array with max pixel values for each slice
		setSlice(ss+1);
		//print( "slice " + d2s(ss+1,0) + ", avg " + d2s(mean_px_values[ss],0) + ", max " + d2s(max_px_values[ss],0) + ", std " + d2s(std_px_values[ss],0) );
		//this_cutoff = floor(cutoff*(average_mean_px_values*std_px_values[ss]/mean_px_values[ss]));
		//print( "slice " + d2s(ss+1,0) + ", range " + d2s(this_cutoff,0) + " to " + d2s(max_px_values[ss],0) );
		setMinAndMax(floor(cutoff*(average_mean_px_values*std_px_values[ss]/mean_px_values[ss])), max_px_values[ss]);
		run("Apply LUT", "slice");
	}
	//print(d2s(cutoff,0)); exit();
	//morpology processing to capture large area around blob that may not have triggered by threshold value
	run("Enhance Contrast...", "saturated=25.0 normalize process use");
	run("Minimum...", "radius=1 stack");
	run("Square", "stack");
	run("Gaussian Blur 3D...", "x=0 y=0 z=2");
	run("Square", "stack");
	run("Maximum...", "radius=16 stack");
	
	//now subtract
	imageCalculator("Subtract stack", orig_id,background_id);
	close(background_id);
	selectWindow(orig_id);
}


var default_microns_per_z_plane = 1.52196113824;
function main_stacks_to_partial_collapse (directory,collapse_Z_size_micron,overlap_redundancy) {
	file_sep = File.separator();
	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);
	
	microns_per_Z_unit = 10; //default
	overlap_factor = 2; //default	
	
	if ( directory == "" || !File.exists(directory) ) {
		print( "Cannot find directory " + directory + "!" );
		return;	
	}
	
	if ( collapse_Z_size_micron == "" || overlap_redundancy == "" ) {
		//prompt for collapse Z size and overlap factor, default is 10um and 1
		Dialog.create("Set Z plane size and overlap redundancy level...");
		Dialog.addNumber("Microns per Z-plane if not specified in metadata",default_microns_per_z_plane);
		Dialog.addChoice("Desired Z-plane depth in micron:", newArray("32","200","100","50","24", "20", "16", "10", "8", "6", "2"));
		Dialog.addChoice("Overlap fold redundancy:", newArray("1", "0", "2", "3", "4", "5", "8"));
		Dialog.show();
		default_microns_per_z_plane = Dialog.getNumber();
		collapse_Z_size_micron = Dialog.getChoice();
		overlap_redundancy = Dialog.getChoice();
	}
	
	microns_per_Z_unit = parseInt(collapse_Z_size_micron);
	overlap_factor = parseInt(overlap_redundancy) + 1;

	//PSFgen_path = File.openDialog("Select PSFGenerator.jar ...");
	
	//Count the maximum number of positions and slices in dataset
	run("Bio-Formats Macro Extensions");
	
	//newPosition = 0;
	//newSlice = 0;
	//maxPosition = 0;
	//maxSlice = 0;
    //add_this_file = false;

	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
	
	//process_view_list = newArray(0);
	
	for (i=0; i<fileList.length; i++) {
		if ((startsWith(fileList[i], "LSFM_") || startsWith(fileList[i], "t000")) && endsWith(fileList[i], ".tif")) {

			//now, remove .tif from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
			
			//ignore PSF files
			if ( endsWith(filepaddedname, "_psf") ) {
				continue;
			}
			
			//now, segment filename into TXXXX_CX_AXXX and then create view setup with "CX_AXXX"
			/*name_ext = split(filepaddedname,"_");
			if ( name_ext.length > 3 ) {
				this_view = name_ext[2] + "_" + name_ext[3];
			} else if ( name_ext.length > 2 ) {
				this_view = name_ext[2];
			} else {
				this_view = "V00";
			}*/
			name_ext = split(filepaddedname,"_"); //default for time series runs
			if (startsWith(fileList[i], "LSFM__")) { //output for large stack deconvolution runs contains __
				name_ext = split(filepaddedname,"(__)");
				this_view = name_ext[2] + "_" + name_ext[3];
			} else if (startsWith(fileList[i], "t000")) {
				this_view = name_ext[1];
			} else {
				this_view = name_ext[2] + "_" + name_ext[3];
			}
			
			
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				view_setups = Array.concat( view_setups, this_view );
				
			}
			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] );
		}
	}
	
	useKLB = false;
	if ( process_file_list.length < 1 ) {
		view_setups = newArray(0);
		useKLB = true;
		for (i=0; i<fileList.length; i++) {
			if (startsWith(fileList[i], "t0") && endsWith(fileList[i], ".klb")) {
	
				//now, remove .tif from end of filename
				/*name_ext = split(fileList[i],".");
				filepaddedname = "";
				for (n=0; n<name_ext.length-1; n++) {
					filepaddedname += name_ext[n];
				}*/
				filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
	
				//now, segment filename into TXXXX_SX and then create view setup with "SX"
				name_ext = split(filepaddedname,"_");
				this_view = name_ext[1];
				
				//see if we have seen this setup before
				is_match = false;
				for (m=0;m<view_setups.length; m++ ) {
					if ( this_view == view_setups[m] ) {
						is_match = true;
						break;
					}
				}
				if ( is_match ) {
					//do nothing	
				} else {
					//add unique view setup
					view_setups = Array.concat( view_setups, this_view );
					
				}
	
				process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] ); //+ "///" + filepaddedname );
			}
		}	
	}
	
	Array.sort(view_setups); //sort files by their intended name, not true name
	//Array.print(view_setups);
	//Array.print(process_file_list);
	
	File.makeDirectory(directory + file_sep + "Partial_Z_Collapsed");
	
	max_channels = 0;
    max_time = 0;
	//this_time = 0;
	//exit();
	view_sizeZ = newArray(view_setups.length);
	Array.fill(view_sizeZ,0);
	
	//okay, first things first, need to figure out how big each stack should be
	setBatchMode(true);
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		outfileList = newArray(0);
		max_sizeZ = 0;
		max_scaleZ = 0;
		
		for (i=0; i<process_file_list.length; i++) {
		//for (i=0; i<1; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			file = directory + file_sep + name_ext[1];
			processList = Array.concat( processList, file );
			outfileList = Array.concat( outfileList, directory + file_sep + "Partial_Z_Collapsed" + file_sep + name_ext[1] );
			
			Ext.setId(file);
			Ext.getSizeZ(sizeZ);
			Ext.getPixelsPhysicalSizeZ(scaleZ); // spacing between Z sections in microns, or NaN if the spacing is not stored in the original file.
			Ext.close();
			
			//print( "   file: " + name_ext[1] + " sizeZ = " + d2s(sizeZ,0) );
			if ( sizeZ > max_sizeZ ) {
				max_sizeZ = sizeZ;	
			}
			if ( isNaN(scaleZ) ) {
				scaleZ = default_microns_per_z_plane;
			} else if ( scaleZ > max_scaleZ ) {
				max_scaleZ = scaleZ;	
			}			
			
		}
		if (max_scaleZ == 0 ) {
			max_scaleZ = default_microns_per_z_plane; //just use 1 micron as default size if size isn't known
		}
		
		//exit();
		
		Ext.close();
		call("java.lang.System.gc");
		
		//okay, now figure out how many original Z-stacks get collapsed into a new Z-stack, and what offset to use on consecutive collapsed images
		images_per_frame = microns_per_Z_unit / max_scaleZ; //total micros of image depth divided by size desired per frame
		
		if (isNaN(images_per_frame)) {
			images_per_frame = 1; //default
			offset_between_frames = 1; //default
		} else {
			//now, figure out redundancy factor
			if ( overlap_factor>1 ) {
				offset_between_frames = round(images_per_frame / overlap_factor);
			} else {			
				offset_between_frames = images_per_frame; //default for no redundancy
			}
		}
		
				
		
		print( "View setup: " + view_setups[m] + " sizeZ = " + d2s(max_sizeZ,0)  + " scaleZ = " + d2s(max_scaleZ,5) + " images per frame = " + d2s(images_per_frame,0) + " frame offset = " + d2s(offset_between_frames,0) );
		
		for (i=0; i<processList.length; i++) {
			if ( images_per_frame < 2 ) {
				//just copy TIFs
				File.copy(processList[i], outfileList[i]);
				
			} else {
				if ( useKLB ) {
					setBatchMode(false); //KLB does not work with batch mode
					run("KLB...", "open=[" + processList[i] + "]");
					VirtStack = getImageID();//("window.title");
					setBatchMode(true);
					selectImage(VirtStack);
					run("Duplicate...", "duplicate");
					OrigStack = getTitle();
					selectImage(VirtStack);
					close();
				} else {
					setBatchMode(false); //KLB does not work with batch mode
					open(processList[i]);
					OrigStack = getTitle();
				}
				//OrigStack = getInfo("window.title");
				
				here = 1; //starting Z coordinate
				img_number = 1; //increments for image tites
				running_list = "";
				
				while( here < max_sizeZ ) {
					//print( "here1");
					selectWindow(OrigStack);
					//print( "here3");
					run("Z Project...", "start=" + d2s(here,0) + " stop=" + d2s(here+images_per_frame-1,0) + " projection=[Max Intensity]"); //if no start given, start=1; if no stop given, stop=Z_dim
					newname = "Z" + d2s(here,0) + "_" + OrigStack;
					rename( newname );
					running_list = running_list + " image" + d2s(img_number,0) + "=" + newname;
					here += offset_between_frames;
					img_number++;
				}
				
				run("Concatenate...", " " + running_list + " " + "image" + d2s(img_number,0) + "=[-- None --]");
				
				
				//print( "About to save... " + outfileList[i] );
				
				saveAs("Tiff", outfileList[i] );
				close(); 
				
				if ( isOpen(OrigStack) ) {
					selectImage(OrigStack);
					close();
				}

				close("*");

			}
			call("java.lang.System.gc");
			
		}
	}
}



function main_extract_defined_slices(directory) {
	
	if ( directory == "" || !File.exists(directory) ) {
		//Ask user to choose the input and output directories
		directory = getDirectory("Choose input directory");
	}
	fileList = getFileList(directory);

	//declare expandable arrays
	process_file_list = newArray(0);
	view_setups = newArray(0);
	
	//process_view_list = newArray(0);
	Array.sort(fileList); //very critical step, can be done here or done to array processList below
	
	for (i=0; i<fileList.length; i++) {
		if ((startsWith(fileList[i], "LSFM_") || startsWith(fileList[i], "t000")) && endsWith(fileList[i], ".tif")) {

			//now, remove .tif from end of filename
			/*name_ext = split(fileList[i],".");
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}*/
			filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
			
			//ignore PSF files
			if ( endsWith(filepaddedname, "_psf") ) {
				continue;
			}
			
			//now, segment filename into TXXXX_CX_AXXX and then create view setup with "CX_AXXX"
			name_ext = split(filepaddedname,"_"); //default for time series runs
			if (startsWith(fileList[i], "LSFM__")) { //output for large stack deconvolution runs contains __
				name_ext = split(filepaddedname,"(__)");
				this_view = name_ext[2] + "_" + name_ext[3];
			} else if (startsWith(fileList[i], "t000")) {
				this_view = name_ext[1];
			} else {
				this_view = name_ext[2] + "_" + name_ext[3];
			}
			
			//see if we have seen this setup before
			is_match = false;
			for (m=0;m<view_setups.length; m++ ) {
				if ( this_view == view_setups[m] ) {
					is_match = true;
					break;
				}
			}
			if ( is_match ) {
				//do nothing	
			} else {
				//add unique view setup
				view_setups = Array.concat( view_setups, this_view );
				
			}

			process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] + "///" + filepaddedname );
		}
	}
	useKLB = false;
	if ( process_file_list.length < 1 ) {
		view_setups = newArray(0);
		useKLB = true;
		for (i=0; i<fileList.length; i++) {
			if (startsWith(fileList[i], "t0") && endsWith(fileList[i], ".klb")) {
	
				//now, remove .tif from end of filename
				/*name_ext = split(fileList[i],".");
				filepaddedname = "";
				for (n=0; n<name_ext.length-1; n++) {
					filepaddedname += name_ext[n];
				}*/
				filepaddedname = substring(fileList[i], 0, fileList[i].length-4);
	
				//now, segment filename into TXXXX_SX and then create view setup with "SX"
				name_ext = split(filepaddedname,"_");
				this_view = name_ext[1];
				
				//see if we have seen this setup before
				is_match = false;
				for (m=0;m<view_setups.length; m++ ) {
					if ( this_view == view_setups[m] ) {
						is_match = true;
						break;
					}
				}
				if ( is_match ) {
					//do nothing	
				} else {
					//add unique view setup
					view_setups = Array.concat( view_setups, this_view );
					
				}
	
				process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] + "///" + filepaddedname );
			}
		}	
	}
	Array.sort(view_setups); //sort files by their intended name, not true name

	max_channels = 0;
    max_time = 0;
	
	startFrame = newArray( view_setups.length );
	stopFrame = newArray( view_setups.length );
	
	//set up default parameters, and ask user about options for what we are going to do
	Dialog.create("Choose the slices you want to extract for this series...");
	for ( m=0; m<view_setups.length; m++ ) {
		Dialog.addMessage("View setup: " + view_setups[m],16,"Black");
		Dialog.addNumber("Start frame (leave as 0 to ignore this view setup: ", 0);
		Dialog.addNumber("Stop frame (leave as 0 to only collect start frame: ", 0);
		//Dialog.addCheckbox("Do oblique MIP #3 (x-axis)", do_oblique_proj_3);	
	}
	Dialog.show();
	for ( m=0; m<view_setups.length; m++ ) {
		startFrame[m] = Dialog.getNumber();
		stopFrame[m] = Dialog.getNumber();
	}

	//explore dataset.xml files to get unit information
	X_voxel_micron = 0;
	Y_voxel_micron = 0;
	Z_voxel_micron = 0;
	voxel_units = "";
	input_MVL_lines = newArray(100);
	//print( "looking for... " + directory+ "dataset_klb.xml" );
	if ( File.exists(directory+ "dataset.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb_s0.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb_s0.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	} else if ( File.exists(directory+ "dataset_klb_s1.xml") ) { 
		input_MVL = File.openAsString(directory+ "dataset_klb_s1.xml");
		input_MVL_lines = split(input_MVL, "\n");
		//print( "got XML");
	}
	unit_line = -1;
	size_line = -1;
	for ( l=0; l<input_MVL_lines.length; l++ ) {
		if (indexOf(input_MVL_lines[l], "<size>") >= 0) {
			size_line = l;
		} else if (indexOf(input_MVL_lines[l], "<unit>") >= 0) {
			unit_line = l;
			size_line = -1;
		}
		if ( unit_line >= 0 && size_line >= 0 ) {
			break;
		}
	}	
	if ( unit_line >= 0 && size_line >= 0 ) {
		//print("got lines");
		start_index = indexOf(input_MVL_lines[size_line], "<size>");
		subline = substring(input_MVL_lines[size_line], start_index + 6 );
		stop_index = indexOf(subline, "</size>" );
		size_line_data = substring(subline, 0, stop_index );
		sizes_parsed = split( size_line_data, " " );
		X_voxel_micron = parseFloat(sizes_parsed[0]);
		Y_voxel_micron = parseFloat(sizes_parsed[1]);
		Z_voxel_micron = parseFloat(sizes_parsed[2]);
		
		start_index = indexOf(input_MVL_lines[unit_line], "<unit>");
		subline = substring(input_MVL_lines[unit_line], start_index + 6 );
		stop_index = indexOf(subline, "</unit>" );
		voxel_units = substring(subline, 0, stop_index );
	}
	
	File.makeDirectory(directory + file_sep + "ExtractedFrames" );
	
	//iterate over channels
	for ( m=0; m<view_setups.length; m++ ) {
		
		//initialize variables for processing this view setup
		processList = newArray(0);
		if ( startFrame[m] < 1 ) {
			continue;
		}
		
		for (i=0; i<process_file_list.length; i++) {
			name_ext = split(process_file_list[i],"(///)");
			if ( name_ext[0] != view_setups[m] ) {
				continue;
			}
			
			processList = Array.concat( processList, name_ext[2] );
		}
		
		for (i=0; i<processList.length; i++) {
			setBatchMode(false); //KLB does not work with batch mode
			
				if ( useKLB ) {
					setBatchMode(false); //KLB does not work with batch mode
					run("KLB...", "open=[" + directory + processList[i] + ".klb]");
					VirtStack = getImageID();//("window.title");
					setBatchMode(true);
					selectImage(VirtStack);
					run("Duplicate...", "duplicate");
					OrigStack = getTitle();
					selectImage(VirtStack);
					close();
				} else {
					setBatchMode(false); //KLB does not work with batch mode
					open(processList[i]);
					OrigStack = getTitle();
				}
			
			
			//set up units
			selectImage(OrigStack);
			if ( !(voxel_units == "") ) {
				Stack.setXUnit(voxel_units);
				run("Properties...", "pixel_width="+d2s(X_voxel_micron,8)+" pixel_height="+d2s(Y_voxel_micron,8)+" voxel_depth="+d2s(Z_voxel_micron,8) );
			}
			
			if ( stopFrame[m] < 1 ) {
				//single frame
				run("Duplicate...", "duplicate range=" + d2s(startFrame[m],0));
				saveAs("Tiff", directory + file_sep + "ExtractedFrames" + file_sep + processList[i] + "_f" + d2s(startFrame[m],0) );
				selectImage(OrigStack);	
			} else {
				//incomplete z-MIP
				run("Z Project...", "stop=" + d2s(stopFrame[m],0) + " start=" + d2s(startFrame[m],0) + " projection=[Max Intensity]");
				saveAs("Tiff", directory + file_sep + "ExtractedFrames" + file_sep + processList[i] + "_f" + d2s(startFrame[m],0) + "-" + d2s(stopFrame[m],0) + "_zmax" );
				close();  // close projection
				selectImage(OrigStack);	
				saveAs("Tiff", directory + file_sep + "ExtractedFrames" + file_sep + processList[i] + "_f" + d2s(startFrame[m],0) + "-" + d2s(stopFrame[m],0) );
			}
			close(); // close origina image
			call("java.lang.System.gc");
		}
	}
	setBatchMode(false);
}

function main_mip_to_avi () {
	directory = getDirectory("Choose input directory");
	fileList = getFileList(directory);
	//declare expandable arrays
	//process_file_list = newArray(0);
	view_setups = newArray(0);
	chan_setups = newArray(0);
	
	//process_view_list = newArray(0);
	Array.sort(fileList); //very critical step, can be done here or done to array processList below
	min_time = NaN;
	max_time = -1;
	min_chan = NaN;
	max_chan = 0;
	
	for (i=0; i<fileList.length; i++) {
		if (startsWith(fileList[i], "t0") && endsWith(fileList[i], ".tif")) {

			//now, remove .tif from end of filename
			name_ext = split(fileList[i],".");
			if (  endsWith(fileList[i], ".c.tif") ) {
				name_ext = split(fileList[i],"(\\.c\\.)");
			}
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}

			//now, segment filename into TXXXX_SX and then create view setup with "SX"
			name_ext = split(filepaddedname,"_");
			this_view = name_ext[2];
			for (n=3; n<name_ext.length; n++ ) {
				this_view = this_view + "_" + name_ext[n];
			}			
			this_chan = name_ext[1];
			
			//capture max_time for dialog box
			this_time = parseInt( substring( name_ext[0], 1 ) );
			if ( this_time > max_time ) {
				max_time = this_time;
			}
			if ( isNaN(min_time) || this_time < min_time ) {
				min_time = this_time;	
			}
			
			//at time min, get all the possible views and channels
			if ( this_time == min_time ) {
				//see if we have seen this setup before
				is_match = false;
				for (m=0;m<view_setups.length; m++ ) {
					if ( this_view == view_setups[m] ) {
						is_match = true;
						break;
					}
				}
				if ( is_match ) {
					//do nothing	
				} else {
					//add unique view setup
					view_setups = Array.concat( view_setups, this_view );
				}
				is_match = false;
				for (m=0;m<chan_setups.length; m++ ) {
					if ( this_chan == chan_setups[m] ) {
						is_match = true;
						break;
					}
				}
				if ( is_match ) {
					//do nothing	
				} else {
					//add unique view setup
					chan_setups = Array.concat( chan_setups, this_chan );
				}
			}

			//process_file_list = Array.concat( process_file_list, this_view + "///" + fileList[i] + "///" + filepaddedname );
		}
	}
	
	Array.sort(view_setups);
	Array.sort(chan_setups);
	
	gamma_filter = 0.8; //1 is no gamma adjustment
	contrast_adjust = 0.01; //0 is no contrast adjustment
	clip_lower_end = 50; // percent of mean histogram value for the stack, default is 50%
	save_AVI = true;
	merge_chan = true;
	AVI_framerate = 12;
	
	view_setups_do = newArray(view_setups.length);
	Array.fill( view_setups_do, true );
	chan_setups_do = newArray(chan_setups.length);
	Array.fill( chan_setups_do, true );
	
	Dialog.create("MIPs (folder in batch) to AVI settings...");
	
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("MIP to AVI settings...",16,"Black");

	Dialog.addNumber("Gamma filter (1 = filter off):", gamma_filter, 2, 5, "");
	Dialog.addNumber("Contrast clip % (0 = filter off):", contrast_adjust, 3,6, "");
	Dialog.addNumber("Clip below (% of histogram mean; 0 = off):", clip_lower_end, 0, 3, "");
	Dialog.addCheckbox("Merge channels (except anaglyphs):", merge_chan );
	Dialog.addCheckbox("Save AVI (remains open if unchecked)", save_AVI );
	Dialog.addNumber("AVI framerate:", AVI_framerate );

	Dialog.addMessage("\n\n");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("\n\nChoose timepoints...",16,"Black");
	Dialog.addNumber("Start:", min_time);
	Dialog.addNumber("End:", max_time);
	
	Dialog.addMessage("\n\n");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("\n\nChoose channels...",16,"Black");
	for (i=0; i<chan_setups_do.length; i++) {
		Dialog.addCheckbox(chan_setups[i], chan_setups_do[i] );
	}
	Dialog.addMessage("\n\n");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("\n\nChoose views...",16,"Black");
	for (i=0; i<view_setups_do.length; i++) {
		Dialog.addCheckbox(view_setups[i], view_setups_do[i] );
	}
	Dialog.show();
	
	gamma_filter = Dialog.getNumber();
	contrast_adjust = Dialog.getNumber();
	clip_lower_end = Dialog.getNumber();
	merge_chan = Dialog.getCheckbox();	
	save_AVI = Dialog.getCheckbox();	
	AVI_framerate = Dialog.getNumber();

	min_time = Dialog.getNumber();
	max_time = Dialog.getNumber();
	
	for (i=0; i<chan_setups_do.length; i++) {
		chan_setups_do[i] = Dialog.getCheckbox();
	}
	for (i=0; i<view_setups_do.length; i++) {
		view_setups_do[i] = Dialog.getCheckbox();
	}
	
	//getDir("luts") 
	chan_setups_lut = newArray(chan_setups.length);
	Array.fill( chan_setups_lut, "" );
	channels_to_merge = 0;
	for (i=0; i<chan_setups_do.length; i++) {
		if ( chan_setups_do[i] ) {
			chan_setups_lut[i] = File.openDialog("Select a LUT file for " + chan_setups[i] + "..." );
			channels_to_merge++;
		}
	}
	
	File.makeDirectory(directory + file_sep + "AVIs");
	last_chan = chan_setups.length -1;
	for (ii=last_chan; i>=0; ii--) {
		if ( chan_setups_do[ii] ) {
			last_chan = ii;
			break;
		}
	}

	for (i=0; i<view_setups.length; i++) {
		if ( view_setups_do[i] ) {
			//print( "View setup: " + view_setups[i] + ", last chan: " + last_chan );
			for (ii=0; ii<=last_chan; ii++) {
				//print( " testing " + d2s(ii,0) + ", " + d2s(chan_setups_do[ii],0) );
				if ( chan_setups_do[ii] ) {
					//print( "  Chan setup: " + chan_setups[ii] );
					type_string = "";
					is_anaglyph = false;
					if ( startsWith(view_setups[i], "anag" )  ) {
						//print( "   doing RGB" ) ;
						type_string = " type=RGB";
						is_anaglyph = true;
					}
					run("Image Sequence...", "dir=["+directory+"]"+type_string+" filter="+chan_setups[ii]+"_"+view_setups[i]+" start="+d2s(min_time,0)+" count="+d2s(max_time-min_time,0)+" sort");
					rand_num = d2s( floor(random() * 1000000 ), 0 );
					rename( "Movie_" + rand_num );
					run("Animation Options...", "speed="+d2s(AVI_framerate,0) );
					
					//Stack.getStatistics( area, mean, min, max, std );
					
					run("Gamma...", "value="+d2s(gamma_filter,2)+" stack");
					if ( is_anaglyph ) {
						run("Enhance Contrast...", "saturated="+d2s(contrast_adjust,3)+" process_all use");	
					} else {
						run("Enhance Contrast...", "saturated="+d2s(contrast_adjust,3)+" normalize process_all use");			
					}
					
					//getMinAndMax(min, max);
					Stack.getStatistics( area, mean, min, max, std );
					setMinAndMax(floor(clip_lower_end * mean / 100), max);
					run("Apply LUT", "stack");
					
					if ( !is_anaglyph ) {
						run("8-bit");
						run("LUT... ", "open=["+chan_setups_lut[ii]+"]");
					}
					
					if ( is_anaglyph || merge_chan==false || channels_to_merge < 2 ) {
						run("AVI... ", "compression=PNG frame=12 save=["+directory+file_sep + "AVIs"+file_sep+chan_setups[ii]+"_"+view_setups[i]+"_"+d2s(min_time,0)+"to"+d2s(max_time,0)+".avi]");
						close();
					} else {
						rename( "Movie_" + d2s(ii,0) );
						
						if ( ii == last_chan )  { //last channel for this setup
							cur_chan = 1;
							merge_text_list = "c" + d2s(cur_chan,0) + "=Movie_" + d2s(ii,0);
							cur_chan++;
							for ( jj=ii-1; jj>=0; jj-- ) {
								if (chan_setups_do[jj]) {
									merge_text_list = merge_text_list + " c" + d2s(cur_chan,0) + "=Movie_" + d2s(jj,0);
									cur_chan++;
								}
							}
							
							run("Merge Channels...", merge_text_list + " create");
							run("AVI... ", "compression=PNG frame=12 save=["+directory+file_sep + "AVIs"+file_sep+view_setups[i]+"_"+d2s(min_time,0)+"to"+d2s(max_time,0)+".avi]");
							close();
						}
						
					}
				}
			}			
		}
	}
	//run("Image Sequence...", "dir=[/media/martin/12TiB RAID0/Fused Live Datasets/2019-11-07 E6.50 F6nGFP+LmChMesp1 Dataset/Side/MIPs/] filter=s01_back count=151 sort use");
	//close();
	//run("Image Sequence...", "dir=[/media/martin/12TiB RAID0/Fused Live Datasets/2019-11-07 E6.50 F6nGFP+LmChMesp1 Dataset/Side/MIPs/] filter=s01_complete start=3 count=151 sort use");
	//run("Image Sequence...", "dir=[/media/martin/12TiB RAID0/Fused Live Datasets/2019-11-07 E6.50 F6nGFP+LmChMesp1 Dataset/Side/MIPs/] type=RGB filter=s01_anaglpyh1 count=151 sort use");
}

//--------------
//   MAIN
//--------------


//get files for processing -- for command line interface run of this file
arglist = getArgument();
if (arglist == "" ) {
	//do nothing if this is being run as a macro within IJ GUI
} else {
	list = split(arglist,"###");
	
	//batch processing from command line will favor deconvolution of large stacks rather than whole series, because whole series are huge and time-consuming and don't really need to be run in batch (each series is already a batch)
	main_ijm_deconvolve_large_stack(list[0],list[1]);
	//main_ijm_deconvolve_series(list[0],list[1],false);
}
macro "Show LSFMProcessing installation instructions" {
	show_instructions();	
}
macro "  " {
	// menu spacer
}
macro "0. Change LSFM processing settings..." {
	Dialog.create("LSFM processing settings...");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("LSFM .czi deconvolution settings...",16,"Black");
	Dialog.addCheckbox("Deconvolve stacks before writing in .tif format", deconvolve_yes );
	Dialog.addNumber("Max depth (slices) for large stack blockwise deconvolution:", max_slice_depth);
	//Dialog.addMessage("  heap allocation in MB should be 200 times max slice depth");
	Dialog.addNumber("# iterations (0 non-iterative OR >1 iterative):", deconvolution_iterations);
	Dialog.addNumber("Regression coefficient multiplier (non-iterative only):", deconvolution_regression_parameter);
	Dialog.addCheckbox("Attempt to subtract camera noise from stack", deconvolution_subtract_camera_noise );
	Dialog.addNumber("Detection NA penalty (% of NA, improves aberration handling):", detection_NA_penalty );
	Dialog.addNumber("Default tissue RI (when metadata does not specify):", default_tissue_refractive_index );
	Dialog.addNumber("Default immersion RI (when metadata does not specify):", default_immersion_refractive_index );
	Dialog.addMessage("\n\n");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("\n\nLSFM .tif filtering settings...",16,"Black");
	Dialog.addCheckbox("Generate z-MIPs for previewing results", generate_mips_when_filtering);
	Dialog.addCheckbox("Minimal-loss 8-bit conversion (single image and time series)", convert_to_8_bit);
	Dialog.addNumber("Element radius (px) for background subtraction (time series):", rolling_ball_radius);
	Dialog.addNumber("Element radius (px) for background subtraction (single):", fixed_rolling_ball_radius);
	
	Dialog.addCheckbox("Bright precipitate removal (single image only)", fixed_precipitate_removal);	
	Dialog.addCheckbox("Bright blob removal (single image only)", fixed_blob_removal);	
	Dialog.show();

	deconvolve_yes = Dialog.getCheckbox();
	max_slice_depth = Dialog.getNumber();
	deconvolution_iterations = floor( Dialog.getNumber() );
	deconvolution_regression_parameter = Dialog.getNumber();
	deconvolution_subtract_camera_noise = Dialog.getCheckbox();
	detection_NA_penalty = Dialog.getNumber();
	default_tissue_refractive_index = Dialog.getNumber();
	default_immersion_refractive_index = Dialog.getNumber();
	generate_mips_when_filtering = Dialog.getCheckbox();
	convert_to_8_bit = Dialog.getCheckbox();
	rolling_ball_radius = Dialog.getNumber();
	fixed_rolling_ball_radius = Dialog.getNumber();
	fixed_precipitate_removal = Dialog.getCheckbox();	
	fixed_blob_removal = Dialog.getCheckbox();	
}
macro "1. Deconvolve Z.1 or UM2 acquisitions with large stacks (folder in batch)." {
	main_ijm_deconvolve_large_stack(true,"","");
	setBatchMode("exit and display");
	print("DONE: Deconvolve Z.1 or UM2 acquisitions with large stacks (folder in batch).");
}
macro "1. Deconvolve single Z.1 or UM2 acquisitions with large stacks (single file)..." {
	main_ijm_deconvolve_large_stack(false,"","");
	setBatchMode("exit and display");
	print("DONE: Deconvolve single Z.1 or UM2 acquisitions with large stacks (single file).");
}
macro "1. Deconvolve Z.1 or MuVi time series files (folder in batch)..." {
	main_ijm_deconvolve_series("","",false);
	setBatchMode("exit and display");
	print("DONE: Deconvolve Z.1 or MuVi time series .czi files (folder in batch).");
}
macro "1. Deconvolve x2 (very slow) Z.1 or MuVi time series files (folder in batch)..." {
	main_ijm_deconvolve_series("","",true);
	setBatchMode("exit and display");
	print("DONE: Deconvolve x2 (very slow) Z.1 or MuVi time series .czi files (folder in batch).");
}
macro "2. Filter and unify z-depth of LSFM .tif files (best for time series; folder in batch)..." {
	main_tifseries_16bit_to_8bit(true,"","");
	setBatchMode("exit and display");
	print("DONE: Convert LSFM .tif files with large stacks (best for time series; folder in batch).");
}
macro "2. Filter LSFM .tif files (single file)..." {
	main_tifseries_16bit_to_8bit(false,"","");
	setBatchMode("exit and display");
	print("DONE: Convert and filter LSFM .tif files to 8-bit (single file).");
}
macro "3. BigStitcher (PLEASE split fused h5/xml images to 1 timepoint and 1 setup per partition)..." {
	run("BigStitcher");
}
macro "3. Convert LSFM .tif time series to z-MIP anaglyphs (folder in batch)..." {
	main_series_tif_to_anaglyphs("","");
	setBatchMode("exit and display");
	print("DONE: Convert LSFM .tif time series to z-MIP anaglyphs (folder in batch).");
}
macro "3. Convert LSFM .tif files to z-MIPs (folder in batch)..." {
	main_time_series_tif_to_mip("","");
	setBatchMode("exit and display");
	print("DONE: Convert LSFM .tif files to z-MIPs (folder in batch).");
}
macro "3. Convert .tif (folder in batch) to AVIs..." {
	directory = getDirectory("Choose input directory");
	fileList = getFileList(directory);
	fr = getNumber("Framerate (fps):", 12);
	File.makeDirectory(directory + file_sep + "AVIs");
	setBatchMode(true);
	for (i=0; i<fileList.length; i++) {
		if ( endsWith(fileList[i], ".tif") ) {
			//now, remove .tif from end of filename
			name_ext = split(fileList[i],".");
			if (  endsWith(fileList[i], ".c.tif") ) {
				name_ext = split(fileList[i],"(\\.c\\.)");
			}
			filepaddedname = "";
			for (n=0; n<name_ext.length-1; n++) {
				filepaddedname += name_ext[n];
			}
			
			open( directory + file_sep + fileList[i] );
			run("AVI... ", "compression=PNG frame="+d2s(fr,0)+" save=["+directory+file_sep + "AVIs"+file_sep+filepaddedname+".avi]");
			close();
		}
	}
	setBatchMode("exit and display");
	print("DONE: Convert .tif (folder in batch) to AVIs.");
}
macro "4. Convert fused .h5 to .klb (folder in batch, requires h5/xml split to 1 timepoint and 1 setup per partition)..." {
	directory = getDirectory("Choose h5/xml input directory");
	
	if ( !File.exists(path_dataset_folder_export_all_h5_to_klb_pyklb) ) {
		path_dataset_folder_export_all_h5_to_klb_pyklb = File.openDialog("Please locate python script: dataset_folder_export_all_h5_to_klb_pyklb.py");
	}
	
	Dialog.create("Choose bitdepth...");
	Dialog.addChoice("Bits per pixel:",newArray("8","16"), "16")
	Dialog.show();
	
	outbits = "16";
	outbits = Dialog.getChoice();
	print( "calling: python3 " + path_dataset_folder_export_all_h5_to_klb_pyklb + " " + directory + " " + outbits );
	//exec("python3 " + path_dataset_folder_export_all_h5_to_klb_pyklb); //,directory,outbits);	
	exec("python3", path_dataset_folder_export_all_h5_to_klb_pyklb, directory, outbits);
	print("DONE: Convert fused .h5 to .klb (folder in batch).");
}
macro "5. Create .klb BigDataViewer dataset.xml file..." {
	run("Open KLB");
	//print("DONE: Create .klb BigDataViewer dataset.xml file.");
}
macro "5. Convert t0XXXX .klb files to anaglyphs or MIPs (folder in batch)..." {
	main_klb_to_mip("");
	setBatchMode("exit and display");
	print("DONE: Convert t0XXXX .klb files to anaglyphs or MIPs (folder in batch).");
}
macro "5. Convert LSFM .tif or t0XXXX .klb files to partially collapsed Z-stacks..." {
	main_stacks_to_partial_collapse("","","");
	setBatchMode("exit and display");
	print("DONE: Convert LSFM .tif or t0XXXX .klb files to partially collapsed Z-stacks.");
}
macro "6. Extract defined slice(s) from LSFM .tif or t0XXXX .klb files in time series..." {
	main_extract_defined_slices ("");
	setBatchMode("exit and display");
	print("DONE: Extract defined slice(s) from LSFM .tif or t0XXXX .klb files in time series.");
}
macro "7. Process t0XXXX MIPs (folder in batch) to AVIs..." {
	main_mip_to_avi();	
	setBatchMode("exit and display");
	print("DONE: Process t0XXXX MIPs (folder in batch) to AVIs.");
}
macro "  " {
	// menu spacer
}
macro "Macro: Bright blob/precipitate remover (selected image ONLY)..." {
	if ( fixed_precipitate_removal || fixed_blob_removal ) {
	bright_blob_remover_this_image(fixed_precipitate_removal,fixed_blob_removal);
	}
}
macro "Macro: Enhance local contrast (stack, selected image ONLY)..." {
	blocksize = 127;
	histogram_bins = 256;
	maximum_slope = 2.0;
	mask = "*None*";
	fast = true;
	process_as_composite = true;
	
	Dialog.create("CLAHE processing settings...");
	Dialog.addNumber("Blocksize:", blocksize);
	Dialog.addNumber("Histogram bins:", histogram_bins);
	Dialog.addNumber("Maximum slope:", maximum_slope);
	Dialog.addCheckbox("Fast", fast );
	Dialog.show();
	
	blocksize = Dialog.getNumber();
	histogram_bins = Dialog.getNumber();
	maximum_slope = Dialog.getNumber();
	fast = Dialog.getCheckbox();	
	
	getDimensions( width, height, channels, slices, frames );
	isComposite = channels > 1;
	parameters =
	  "blocksize=" + blocksize +
	  " histogram=" + histogram_bins +
	  " maximum=" + maximum_slope +
	  " mask=" + mask;
	if ( fast )
	  parameters += " fast_(less_accurate)";
	if ( isComposite && process_as_composite ) {
	  parameters += " process_as_composite";
	  channels = 1;
	}
	  
	for ( f=1; f<=frames; f++ ) {
	  Stack.setFrame( f );
	  for ( s=1; s<=slices; s++ ) {
	    Stack.setSlice( s );
	    for ( c=1; c<=channels; c++ ) {
	      Stack.setChannel( c );
	      run( "Enhance Local Contrast (CLAHE)", parameters );
	    }
	  }
	}
	print("DONE: Enhance local contrast (stack, selected image ONLY).");
}
macro "Macro: Multichannel single slice image to montage (selected image ONLY)..." {
	main_multichannel_to_channels_and_merge();
}
macro "Macro: Montage/tiles to stack (selected image ONLY)..." {
	getDimensions(dim_width, dim_height, dim_channels, dim_slices, dim_frames);
	master_title = getTitle();
	rand_num = d2s( floor(random() * 1000000 ), 0 );
		
	n_wide = 2;
	n_tall = 2;
	border_internal = 12;
	border_external = 0;
	//fast = true;
	//process_as_composite = true;
	
	Dialog.create("Montage/Tiles to stack...");
	Dialog.addNumber("Horizontal panels:", n_wide);
	Dialog.addNumber("Vertical panels:", n_tall);
	Dialog.addNumber("Internal dividers width (px):", border_internal);
	Dialog.addNumber("External border width (px)", border_external );
	Dialog.show();
	
	n_wide = Dialog.getNumber();
	n_tall = Dialog.getNumber();
	border_internal = Dialog.getNumber();
	border_external = Dialog.getNumber();
	
	tile_width = ( dim_width - ( 2 * border_external ) - ( (n_wide-1) * border_internal ) ) / n_wide;
	tile_height = ( dim_height - ( 2 * border_external ) - ( (n_tall-1) * border_internal ) ) / n_tall;
	
	if ( tile_height != floor(tile_height) || tile_width != floor(tile_width) ) {
		print("Indivisible image dimensions for tile specifications!  Adjust number of tiles or border width(s), and try again.");
	} else if ( n_wide * n_tall < 2 || n_wide < 1 || n_tall < 1 ) {
		print("Not enough tiles specified!  Adjust number of tiles or border width(s), and try again.");
	} else {
		//proceed by counting down horizontal axis first, then vertical
		pos_x = border_external;
		pos_y = border_external;
		img_num = 0;
		for ( iii = 0; iii < n_tall; iii++ ) {
			for ( ii = 0; ii < n_wide; ii++ ) {
				selectWindow(master_title);
				run("Select None");
				pos_x = border_external + ( ii * ( tile_width + border_internal ) );
				pos_y = border_external + ( iii * ( tile_height + border_internal ) );
				//run("Specify...", "width="+d2s(tile_width,0)+" height="+d2s(tile_height,0)+" x="+d2s(pos_x,0)+" y="+d2s(pos_y,0)+" slice=1");
				makeRectangle(pos_x, pos_y, tile_width, tile_height);
				run("Duplicate...", " ");
				rename( "Image" + rand_num + "_" + d2s(img_num,0) );
				if ( img_num > 0 ) {
					run("Concatenate...", "  image1=Image" + rand_num + "_" + d2s(img_num-1,0) +" image2=Image" + rand_num + "_" + d2s(img_num,0) + " image3=[-- None --]" );
					rename( "Image" + rand_num + "_" + d2s(img_num,0) ); //rename the concatenated stack
				}
				img_num++;
			}
		}
		selectWindow("Image" + rand_num + "_" + d2s(img_num-1,0));
		rename( master_title + "_as_stack" );
		//selectWindow(master_title);
		//close();
	}
}
macro "Macro: Set Display Range..." {
	// Sets the display range of the active image.
	getMinAndMax(min, max);
	min = getNumber("Min:", min);
	max = getNumber("Max:", max);
	setMinAndMax(min, max);
}
