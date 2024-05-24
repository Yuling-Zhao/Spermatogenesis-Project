//This macron script is used for live imaging analysis: 48 hours with the interval of 15min, and 20 slices for each time point;
//This macron script separate the images from different time points into separate folders
//The next step will be using maximum intensity projection to pick the best focus slice from each time point

// Specify the base directory containing the files
baseDirectory = getDir("");
// Start the loop from 1 to 193
for (t=1; t<=116; t++) {
    // Format the time point into a three-digit string (e.g., 001, 002, ..., 193)
    timePoint = d2s(t, 0);  // Convert number to string with 0 decimal places
    while (lengthOf(timePoint) < 3) {
        timePoint = "0" + timePoint;  // Pad with leading zeros
    }
    
    // Create the prefix for the files
    prefix = "NPC_multi_pos_GFP_s09t" + timePoint;
    
    // Define the new directory path based on the prefix
    newDirPath = baseDirectory + prefix + "/";
    
    // Check if directory exists, if not, create it
    if (!File.exists(newDirPath)) {
        File.makeDirectory(newDirPath);
    }
    
    // List all files in the base directory
    list = getFileList(baseDirectory);
    
    // Move files with the specified prefix to the new directory
    for (i=0; i<list.length; i++) {
        if (startsWith(list[i], prefix)) {
            File.rename(baseDirectory + list[i], newDirPath + list[i]);
        }
    }
}
