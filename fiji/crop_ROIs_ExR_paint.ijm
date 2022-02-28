parent = 'C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/ExR_PAINT/unregistered/';
target = 'C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/ExR_PAINT/cropped_rois_doubles_unreg/';

ROIs = newArray("well1_roi1",
"well1_roi2",
"well2_roi1",
"well2_roi2",
"well4_roi1"
);

for (j=0; j<ROIs.length; j++){ //loop through each ROI

	//crop PAINT rois
	roiManager("Open", target + "ROIs/"+ ROIs[j] + "_PAINT_rois.zip"); //open the manually cropped fovs
						
    for (n=0; n<roiManager("count"); ++n) {
    	open(parent+ROIs[j]+"_PAINT.tif");
    	roiManager("Select", n);
    	run("Crop");
        saveAs("Tiff", target + ROIs[j] + "_PAINT_ROI"+n+".tif");
        close();
    }

    roiManager("Delete");

	if (isOpen("ROI Manager")) {
	 selectWindow("ROI Manager");
	 run("Close");
	  	}

  	//crop ExR
	roiManager("Open", target + "ROIs/"+ ROIs[j] + "_ExR_rois.zip"); //open the manually cropped fovs
						
    for (n=0; n<roiManager("count"); ++n) {
    	open(parent+ROIs[j]+"_ExR_reg.tif");
    	roiManager("Select", n);
    	run("Crop");
        saveAs("Tiff", target + ROIs[j] + "_ExR_ROI"+n+".tif");
        close();
    }

    roiManager("Delete");

	if (isOpen("ROI Manager")) {
	 selectWindow("ROI Manager");
	 run("Close");
	  	}

}
	


run("Close All");
run("Collect Garbage");