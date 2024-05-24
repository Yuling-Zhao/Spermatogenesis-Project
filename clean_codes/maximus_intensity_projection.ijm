//This script uses maximum intensity projection method to get the best focus slice from z stack of each time point;


baseDirectory = getDir("");
outputpath = "/Users/aktharlab/Desktop/experiment_records/PHF20_PHF20L1/PHF20/YLZ-119_live_imaging_analysis/NPC_D0/nontreated/2105/maxi_intensity/s09/region_1/";
for (t=112; t<=116; t++) {
	timePoint = d2s(t, 0);  // Convert number to string with 0 decimal places
    while (lengthOf(timePoint) < 3) {
        timePoint = "0" + timePoint;  // Pad with leading zeros
    }
	prefix = "NPC_multi_pos_GFP_s09t" + timePoint;
	File.openSequence(baseDirectory+prefix);
	makeRectangle(291, 425, 630, 630);
	run("Duplicate...", "duplicate");
    run("Z Project...", "projection=[Max Intensity]");
    selectImage("MAX_"+prefix+"-1");
    run("Brightness/Contrast...");
    waitForUser("Click OK when you are done");
    saveAs("Tiff", outputpath+"MAX_"+prefix+".tif");
	run("Close All");
}
