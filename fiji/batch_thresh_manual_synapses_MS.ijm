parent = 'C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/idExMDecrowdingSegs/';
sub = 'Decrowding crop image/images/';
layers = newArray("L1","L23","L4");
folders = newArray(
//					"A1_cacha_shank/",
					"A2_PSD95_homer/"
//					"A3_syngap_shank/",
//					"A4_homer_shank/",
//					"A5_rim_homer/",
//					"A6_bassoon_homer/",
//					"A7_shank_homer/",
//					"B1_cacha_shank/",
//					"B2_PSD95_homer/",
//					"B3_syngap_shank/",
//					"B4_homer_shank/",
//					"B5_rim_homer/",
//					"B6_bassoon_homer/",
//					"B7_shank_homer/",
//					"C1_cacha_shank/",
//					"C2_PSD95_homer/",
//					"C3_syngap_shank/",
//					"C4_homer_shank/",
//					"C5_rim_shank/",
//					"C6_bassoon_homer/",
//					"C7_shank_homer/"
					);
parent_target = 'C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/idExMDecrowdingSegs/';

for (j=0; j<folders.length; j++) {
	for (k=0; k<layers.length; k++) {
		run("Close All");
		findex = j;
		lindex = k;
		image_folder = parent + sub + folders[findex] + layers[lindex] + "/";
		lines=split(folders[findex],"_");
		sampleno = lines[0];
		samplecode = charCodeAt(sampleno,0);
		sample = fromCharCode(samplecode);
		teststain = lines[1];
		refstain_whole = lines[2];
		refstain_splits = split(refstain_whole,"/");
		refstain = refstain_splits[0];
		//print(sample);
		//print(teststain);
		//print(refstain);
	
		//File.makeDirectory(parent_target + folders[findex]);
		targetdir = parent_target + sub + folders[findex] + layers[lindex] + "/masks/";
		print(targetdir);
		File.makeDirectory(targetdir);					  //create a folder called "masks" in the folder
	
		//manually set these for now
		//for the 7/9/19 experiment, order is ref-pre-post (except for PSD95, which has order ref-post-pre)
		//for the 5/25/19 experiment, order is ref-post-pre
		//for rim 5/27, order is ref-post-pre

//		contains(seqname, 'A1_cacha_shank') | ...
//                contains(seqname, 'B1_cacha_shank') | ...
//                contains(seqname, 'A3_syngap_shank') | ...
//                contains(seqname, 'B3_syngap_shank') | ...
//                contains(seqname, 'A4_homer_shank') | ...
//                contains(seqname, 'B4_homer_shank')
                
		if ((indexOf(folders[findex],"A1_cacha_shank")>=0) ||
		(indexOf(folders[findex],"B1_cacha_shank")>=0) ||
		(indexOf(folders[findex],"A3_syngap_shank")>=0) ||
		(indexOf(folders[findex],"B3_syngap_shank")>=0) ||
		(indexOf(folders[findex],"A4_homer_shank")>=0) ||
		(indexOf(folders[findex],"B4_homer_shank")>=0))
		{
			postCh = "C3";
			refCh = "C1";
			preCh = "C2";
			dapiCh = "C4";
		} else {
			refCh = "C1";
			dapiCh = "C4";
			preCh = "C3";
			postCh = "C2";
		}
		
		list = getFileList(image_folder);									 //find all files in current data folder					
		extension = "tif";							 //set file types to look at 'tif' stacks
		index=-1;													 //initialize index that will mark files with 
		 
		for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
			if (endsWith(toLowerCase(list[i]), "." + extension) && (indexOf(list[i], "registered") < 0)){			//if its extension matches up with the _extensions of interest
	
				index=i;
				name =  list[index];												//record the full file name (including file extension; eg:'stack1.tif')
		        dotIndex = indexOf(name, ".t");										//record position of dot that separates extension from filename	
		        title = substring(name, 0, dotIndex); 								//record the abbreviated file name (excluding file extension; eg:'stack1')
		        //print(title);
		        //subtract background from each image in the stack
	
				//filesplits = split(title,"L");
				layer = layers[lindex];//"L" + filesplits[1];
				//print(layer);
				//print(parent_target + "thresholds/" + teststain + "_ref_" + layer + "_" + sample + ".txt");
				reffile=File.openAsString(parent_target + "thresholds/" + teststain + "_ref_" + layer + "_" + sample + ".txt");
				minref = parseInt(reffile)*7;
	
				//print(parent_target + teststain + "_post_" + layer + "_" + sample + ".txt");
				postfile=File.openAsString(parent_target + "thresholds/" + teststain + "_post_" + layer + "_" + sample + ".txt");
				minpost = parseInt(postfile)*7;
				//print(minpost);
	
				//print(parent_target + teststain + "_pre_" + layer + "_" + sample + ".txt");
				prefile=File.openAsString(parent_target + "thresholds/" + teststain + "_pre_" + layer + "_" + sample + ".txt");
				minpre = parseInt(prefile)*7;
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
				saveAs("Tiff", targetdir + "Raw-refCh-" + layer + "-" + sample + "_" + toString(index) + "_noBG.tif"); // save the channel we want to extract from separately
				run("Threshold...");
				setThreshold(minref, 65535);
				run("Convert to Mask", "method=Minimum background=Dark black");
				saveAs("Tiff", targetdir + "Bin-refCh-" + layer + "-" + sample + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
				close();
	
				selectWindow(preCh + "-" + title + ".tif");
				saveAs("Tiff", targetdir + "Raw-preCh-" + layer + "-" + sample + "_" + toString(index) + "_noBG.tif"); // save the channel we want to extract from separately
				run("Threshold...");
				setThreshold(minpre, 65535);
				run("Convert to Mask", "method=Minimum background=Dark black");
				saveAs("Tiff", targetdir + "Bin-preCh-" + layer + "-" + sample + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
				close();
	
				selectWindow(postCh + "-" + title + ".tif");
				saveAs("Tiff", targetdir + "Raw-postCh-" + layer + "-" + sample + "_" + toString(index) + "_noBG.tif"); // save the channel we want to extract from separately
				run("Threshold...");
				setThreshold(minpost, 65535);
				run("Convert to Mask", "method=Minimum background=Dark black");
				saveAs("Tiff", targetdir + "Bin-postCh-" + layer + "-" + sample + "_" + toString(index) + ".tif"); // save the channel we want to extract from separately
				close();
	
	
			    run("Close All");
			    run("Collect Garbage");
			} // end "if" loop selecting only files with the right extension
		} //end "for" loop going through all files in current folder 
	print(parseInt(postfile));
	}
	run("Clear Results");
}