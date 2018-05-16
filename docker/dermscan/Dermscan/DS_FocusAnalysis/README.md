Notes:

- The method currently works on 8bit per color channel (24 bit RGB images), and I have not yet tested them on 12 bit images. Therefore, please convert the images to 24bit RGB for now. The method works on the each spectral band separately(but loads them as white light images).

- The method looks for differences in focal depth images of the same band by making versions of blurred reference image (focal depth 0) and looking for differences with other focal depths. The final evaluation is visual by examining the difference images and leveled difference images.

- The command line to call the methods is DS_FocusAnalysis.exe <folder_path> . The calculated images will be stored in the "processed" subfolder of the folder path. Hence, make sure that "<folder_path>/processed" directory exists before running the app. Currently only loading and saving to PNG format is supported. 