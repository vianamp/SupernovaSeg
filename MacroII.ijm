//
//          Matheus Viana - vianamp@gmail.com
//
// ===========================================================
//MainDir
//  |
//  |---- mito  (contains images of mito signal)
//  |      |
//  |      |------ file1.tif 
//  |              file2.tif
//  |              ...
//  |
//  |---- surface  (target of this routine and
//  |      |        contains images of surface signal + ROIs)
//  |      |
//  |      |------ file1.tif
//  |              file2.tif
//  |              ...
//  |              file1.zip
//  |              file2.zip
//  |              ...
//  |              
//  |---- MitoGraph  (this folder will be created)
//  |      |
//  |      |------ file1.centers
//  |              file2.centers
//  |              ...
// ===========================================================

// Should be the folder "surface" as described above

_RootFolder = getDirectory("Choose a Directory");

_SurfaceFolder = _RootFolder + "surface/";

_MitoFolder = _RootFolder + "mito/";

_dxy = 0.056;

_dz = 0.2;

// A folder called "MitoGraph" will be creasted

File.makeDirectory(_RootFolder+"/MitoGraph");

// List of files "surface"

_FileList = getFileList(_SurfaceFolder);

// Batch mode on

setBatchMode(true);

i = 0;
while (i < _FileList.length)  {

	// Only TIFF files

	if ( endsWith(_FileList[i],".tif") ) {

		// Image name

		_ImageName = split(_FileList[i],".");
		_ImageName = _ImageName[0];

		// Load image
		
		open(_SurfaceFolder+_ImageName+".tif");

		// Create coordinates file
		
		f = File.open(_RootFolder+"MitoGraph/"+_ImageName+".centers");

		// Write header

		print(f, "[RootFolder]");
		print(f, _RootFolder+"/MitoGraph/");
		print(f, "[MitoFolder]");
		print(f, _MitoFolder);
		print(f, "[SurfaceFolder]");
		print(f, _SurfaceFolder);
		print(f, "[Prefix]");
		print(f, _ImageName);
		print(f, "[SpacingXY]");
		print(f, _dxy);
		print(f, "[SpacingZ]");
		print(f, _dz);
		print(f, "[Centers]");

		// Load ROI file

		roiManager("Reset");
		roiManager("Open",_SurfaceFolder+_ImageName+".zip");

		// For each ROI
		
		for (roi = 0; roi < roiManager("count"); roi++) {
			
			roiManager("Select",roi);
			_z = getSliceNumber();
			getSelectionCoordinates(_x,_y);

			// Write the coordinates according to ROI type
			
			if (Roi.getType=="point") {
				
				print(f, d2s(roi,0) + ",1," + d2s(_x[0],0) + "," + d2s(_y[0],0) + "," + d2s(_z,0) + ",0,0");
				
			} else {
				
				if (_x.length == 3) {
					
					r1 = sqrt(pow(_x[1]-_x[0],2)+pow(_y[1]-_y[0],2));
					r2 = sqrt(pow(_x[2]-_x[0],2)+pow(_y[2]-_y[0],2));
					print(f, d2s(roi,0) + ",2," + d2s(_x[0],0) + "," + d2s(_y[0],0) + "," + d2s(_z,0) + "," + d2s(r1,3) + "," + d2s(r2,3));		
					
				} else {
					
					r1 = sqrt(pow(_x[1]-_x[0],2)+pow(_y[1]-_y[0],2));
					print(f, d2s(roi,0) + ",3," + d2s(_x[0],0) + "," + d2s(_y[0],0) + "," + d2s(_z,0) + "," + d2s(r1,3) + ",0");		
					
				}
			}
		}

		print(f, "[end]");
		File.close(f);
		
		// Image processing

		run("Select None");
		run("Median 3D...", "x=3 y=3 z=3");
		
		for (s = 0; s < nSlices; s++) {
			setSlice(s+1);
			run("Median...", "radius=2 slice");
			run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None*");
			run("Unsharp Mask...", "radius=10 mask=0.90");
		}

		// Save resulting image
		
		setSlice(0.5*nSlices);
		resetMinAndMax();
		saveAs("Tiff",_RootFolder+"MitoGraph/"+_ImageName+".tif");
		run("Close");
		
	}
	
	i++;
	
}

setBatchMode(false);
