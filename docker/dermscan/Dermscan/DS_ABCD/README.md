Notes on the current status of the methods:

- The method currently works on 8bit per color channel (24 bit RGB images), and I have not yet tested them on 12 bit images. Therefore, please convert the images to 24bit RGB for now. The method works on the white light images, so It would also be good that the images have corrected colors if needed.

- The asymmetry is calculated only on shape, the new version will probably have correlation of halves as a metric.

- The border contrast info is stores in the output arrays (see the main function of the project), you can add a line of code to print these out. Alternatively, I can add saving these values to a text file for external examination (e.g. using Matlab)

- The D feature is the diameter. A separate function segments the inner parts of the melanoma, it is a "work in progress" version, a bit slower and less accurate than the version I hope to get.

- The command line to call the methods is DS_ABCD.exe <file_name> . The calculated images will be stored in the folder of the original image. Currently only loading and saving to PNG format is supported. 