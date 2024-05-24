//The script is used when there are multiple scenes in one experiment;
//The images of different scenes will be separated into different folders;

// Specify the base directory containing the files
baseDirectory = "/Users/aktharlab/Desktop/experiment_records/PHF20_PHF20L1/PHF20/YLZ-118 NPC transfection with H3K9me3 probe full-length for live imaging/image/2105/NPC_nondiff_multi_pos_GFP/";
// Start the loop from 1 to 10
for (s=1; s<=10; s++) {
    // Format the time point into a three-digit string (e.g., 001, 002, ..., 193)
    scene = d2s(s, 0);  // Convert number to string with 0 decimal places
    while (lengthOf(scene) < 2) {
        scene = "0" + scene;  // Pad with leading zeros
    }
    
    // Create the prefix for the files
    prefix = "NPC_multi_pos_GFP_s" + scene;
    
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
