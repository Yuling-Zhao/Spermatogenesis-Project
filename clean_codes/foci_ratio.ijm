//This script is used to measure the intensity of the whole nuclei and the foci contained
//The output will be the ratio of the total foci intensity in the representative nuclei


//run("Set Measurements...", "area mean integrated redirect=None decimal=3");
input_path = getDirectory("");
fileList = getFileList(input_path); 

total_ratio = newArray(fileList.length);

for (f = 39; f < 108; f++) {
	//clean-up to prepare for analysis
	roiManager("reset");
	run("Close All");
	run("Clear Results");
	
	//open file
	open(input_path+fileList[f]);
	print(input_path+fileList[f]); //display file that is processed
    title = getTitle();
    
	//convert original image to 8-bit before measure
	run("8-bit");
	//duplicate the image without edit for intensity measurement
	run("Duplicate...", "title=copy_1");

	//process the original image to select the whole nuclear
	selectImage(title);
	run("Median...", "radius=8");
	run("Threshold...");
	waitForUser("Click OK when you are done");
	run("Convert to Mask");
	run("Analyze Particles...", "add");
	
	//select image without editing for whole nuclear intensity measure
	selectImage("copy_1");
	//select the nuclear for measure manually, click ok to proceed measure
	waitForUser("Click OK when you are done");
	roiManager("Measure");
	WN_IntDen = getResult("IntDen", 0);
	run("Clear Results");
	
	//duplicate the nuclear region for defining foci regions inside nuclear
	selectImage("copy_1");
	roiManager("Deselect");
	run("Duplicate...", "title=copy_2");
	roiManager("reset");
	//before looking for foci in the nuclear, duplicate again the nuclear for measuring foci intensity
	selectImage("copy_2");
	run("Duplicate...", "title=copy_3");
	
	//define foci with threshold
	selectImage("copy_2");
	run("Median...", "radius=1");
	run("Threshold...");
	waitForUser("Click OK when you are done");
	run("Convert to Mask");
	run("Analyze Particles...", "add");
	
	//measure foci intensity
	selectImage("copy_3");
	roiManager("Deselect");
	roiManager("Measure");
	total_speckle = 0;
	for(i=0; i<nResults; i++) {
    	total_speckle += getResult("IntDen", i);
    	}
	speckle_WN = total_speckle / WN_IntDen;
	total_ratio[f] = speckle_WN;
	print(speckle_WN);
	run("Clear Results");
	run("Close All");
}

Table.create("Total_Ratio");
Table.setColumn("Foci_Ratio", total_ratio);
