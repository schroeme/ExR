parent = '/Users/margaretschroeder/Dropbox (MIT)/2020.7.31_antigen_retrieval_crop/';
folders = newArray("AR_cav/","AR_homer/","AR_psd95/","noAR_cav/","noAR_homer/","noAR_psd95/");
parent_target = parent;

for (j=0; j<folders.length; j++) {
	run("Close All");
	findex = j;
	image_folder = parent + folders[findex];
	//print(image_folder);
	lines=split(folders[findex],"_");
	AR = lines[0];
	proteintemp=lines[1];
	lines2 = split(proteintemp,"/");
	protein=lines2[0];

	targetdir = parent_target + folders[findex] + "masks/";
	File.makeDirectory(targetdir);					  //create a folder called "masks" in the folder

	postCh = "C1";
	refCh = "C2";
	preCh = "C3";
	dapiCh = "C4";
	
	list = getFileList(image_folder);									 //find all files in current data folder					
	extension = "tif";							 //set file types to look at 'tif' stacks
	index=-1;													 //initialize index that will mark files with 
	 
	for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
		if (endsWith(toLowerCase(list[i]), "." + extension)){			//if its extension matches up with the _extensions of interest

			index=i;
			name =  list[index];												//record the full file name (including file extension; eg:'stack1.tif')
	        dotIndex = indexOf(name, ".t");										//record position of dot that separates extension from filename	
	        title = substring(name, 0, dotIndex); 								//record the abbreviated file name (excluding file extension; eg:'stack1')
	        //print(title);
	
			reffile=File.openAsString(parent + protein + "_ref_" + AR + ".txt");
			//print(parent + protein + "_ref_" + AR + ".txt");
			minref = parseInt(reffile)*7;
			//print(minref);

			prefile=File.openAsString(parent + protein + "_pre_" + AR + ".txt");
			//print(parent + protein + "_pre_" + AR + ".txt");
			minpre = parseInt(prefile)*7;
			//print(minpost);

			postfile=File.openAsString(parent + protein + "_post_" + AR + ".txt");
			minpost = parseInt(postfile)*7;
			//print(minpre);

			open(image_folder+list[i]);
			
			//run("Subtract Background...", "rolling=50 stack");

			//set different color channels to RGB, this way the channels can be separated again after registration
			Stack.setDisplayMode("color");
			Stack.setChannel(1);
			run("Red");
			Stack.setChannel(2);
			run("Green");
			Stack.setChannel(3);
			run("Blue");

			run("Split Channels"); //split the channels

			selectWindow(refCh + "-" + title + ".tif");
			saveAs("Tiff", targetdir + "Raw-refCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			run("Threshold...");
			setThreshold(minref, 65535);
			run("Convert to Mask", "method=Minimum background=Dark black");
			run("Auto Threshold", "method=Triangle white stack use_stack_histogram");
			saveAs("Tiff", targetdir + "Bin-refCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			close();

			selectWindow(preCh + "-" + title + ".tif");
			saveAs("Tiff", targetdir + "Raw-preCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			run("Threshold...");
			setThreshold(minpre, 65535);
			run("Convert to Mask", "method=Minimum background=Dark black");
			saveAs("Tiff", targetdir + "Bin-preCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			close();

			selectWindow(postCh + "-" + title + ".tif");
			saveAs("Tiff", targetdir + "Raw-postCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			run("Threshold...");
			setThreshold(minpost, 65535);
			run("Convert to Mask", "method=Minimum background=Dark black");
			saveAs("Tiff", targetdir + "Bin-postCh-" + protein + "_" + AR + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
			close();

		    run("Close All");
		    run("Collect Garbage");
		} 
	}


	run("Clear Results");
}