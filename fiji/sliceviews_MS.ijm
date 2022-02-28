folder = 'D:/Margaret/ROIs/';
list = getFileList(folder);									 //find all files in current data folder					
extension = "tif";							 //set file types to look at 'tif' stacks												 //initialize index that will mark files with 
	 
	for (i=0; i<list.length; i++) {								 //go through all files in the current data folder
		if (endsWith(toLowerCase(list[i]), "." + extension)){			//if its extension matches up with the _extensions of interest
			run("Close All");
			open(folder+list[i]);
			name =  list[i];												//record the full file name (including file extension; eg:'stack1.tif')
	        dotIndex = indexOf(name, ".");										//record position of dot that separates extension from filename	
	        title = substring(name, 0, dotIndex); 								//record the abbreviated file name (excluding file extension; eg:'stack1')
	        
			// get number of slices in the stack
			getDimensions(width, height, channels, slices, frames);
			print(slices);
			
			//set position to that slice
			setSlice(round(slices/2));


			//adjust the gain automatically
			Stack.setDisplayMode("color");
			Stack.setChannel(1);
			resetMinAndMax();
//			run("Enhance Contrast", "saturated=0.35");
			Stack.setChannel(2);
			resetMinAndMax();
//			run("Enhance Contrast", "saturated=0.35");
			Stack.setDisplayMode("composite");

			//set scale to biological units for scale bar
			run("Properties...", "channels=2 slices=" + toString(slices) + " frames=1 unit=micron pixel_width=0.0114 pixel_height=0.0114 voxel_depth=0.02667");
			
			//create three orthogonal views and save each
			//run("Orthogonal Views");
			//windowlist = getList("window.titles");
//			print("Image windows:");
//		     for (i=0; i<windowlist.length; i++)
//		        print("   "+windowlist[i]);
//		  	}


			//reslice from left to right, varying X -> YZ view
			run("Reslice [/]...", "output=0.400 start=Left avoid");
			getDimensions(width, height, channels, slices, frames);
			n = slices;
			run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Upper Left] bold hide overlay label");
			//run("Make Substack...", "channels=1-2 slices=" + toString(round(n/2)));
			run("Duplicate...", "duplicate slices=" + toString(round(n/2)));
			//run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Lower Right] bold hide overlay");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_merge.jpeg");
			Stack.setActiveChannels("10");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_AB42.jpeg");
			Stack.setActiveChannels("01");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_Kv7.2.jpeg");
			close();

			//reslice from left to right, varying Y -> XZ view
			selectWindow(name);
			run("Reslice [/]...", "output=0.400 start=Top avoid");
			getDimensions(width, height, channels, slices, frames);
			n=slices;
			run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Upper Left] bold hide overlay label");
			run("Duplicate...", "duplicate slices=" + toString(round(n/2)));
			//run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Lower Right] bold hide overlay");
			//run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Lower Right] bold hide overlay");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_merge.jpeg");
			Stack.setActiveChannels("10");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_AB42.jpeg");
			Stack.setActiveChannels("01");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_Kv7.2.jpeg");
			close();

			selectWindow(name);
			getDimensions(width, height, channels, slices, frames);
			n=slices;
			run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Upper Left] bold hide overlay label");
			run("Duplicate...", "duplicate slices=" + toString(round(n/2)));
			//run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Lower Right] bold hide overlay");
			//run("Scale Bar...", "width=.1 height=2 font=14 color=White background=None location=[Lower Right] bold hide overlay");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_merge.jpeg");
			Stack.setActiveChannels("10");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_AB42.jpeg");
			Stack.setActiveChannels("01");
			saveAs("Jpeg", folder + "/orthogonal/" + title + "_YZ_Kv7.2.jpeg");
			close();

			run("Close All");
		}
	}