parentout = 'C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/ExR_PAINT/ExR/';
parentin = 'V:/Jinyoung/2021.6.12 DNA-PAINT 3rd gel/';
regionfolders_in = newArray(
	//"0.05x PBS/",
	"1x PBS/"
	);

regionfolders_out = newArray(
	//"0.05x PBS/",
	"1x PBS/"
	);

wellfolders_in = newArray(
	//"well 1/",
	//"well 2/",
	//"well 3/",
	//"well 4/" //has ROI 1 subfolder only
	"well 5/" //has ROI1 subfolder only
	//"well 6/"
	);

roifolders_in = newArray(
	"ROI 1/"
	//"ROI 2/"
	);

for (j=0; j<regionfolders_in.length; j++){
	for (k = 0; k<wellfolders_in.length; k++){
		for (m = 0; m<roifolders_in.length; m++){
			run("Close All");
			//print(toString(k));
			image_folder = parentin + regionfolders_in[j] + wellfolders_in[k] + roifolders_in[m];
			print(image_folder);
			outfolder = parentout + regionfolders_out[j];
			File.makeDirectory(outfolder);
			list = getFileList(image_folder);									 //find all files in current data folder					
			extension = "nd2";							
		
			count = 0;
			for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
				if (endsWith(toLowerCase(list[i]), "." + extension)){
					//print(list[i]);
					if (indexOf(list[i], "40x") >= 0){
						count = count+1;
						name =  list[i];
			        	dotIndex = indexOf(name, ".");										//record position of dot that separates extension from filename	
			        	title = substring(name, 0, dotIndex); 								//record the abbreviated file name (excluding file extension; eg:'stack1')

			        	wellsplit = split(wellfolders_in[k],"/");
			        	wellname = wellsplit[0];

			        	roisplit = split(roifolders_in[m],"/");
			        	roiname = roisplit[0];
			        	
			        	newname = wellname + "_" + roiname;
			        	print(newname);
		
						impath = image_folder + list[i];
						//print(impath);
						run("Bio-Formats Importer","open=[" + impath + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
						run("Subtract Background...", "rolling=50 stack"); //subtract background
		
						saveAs("Tiff", outfolder + newname + ".tif");
						close();
					}
				}
			}
		}
	}
	run("Close All");
	run("Collect Garbage");
}
