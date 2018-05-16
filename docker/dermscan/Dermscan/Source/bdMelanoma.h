/**************************************************************************
 Methods for ABCDE features of melanoma.

 Author: Danilo Babin
 File name: "bdMelanoma.h"
 ***************************************************************************/



#ifndef BD_MELANOMA_DEF
	#define BD_MELANOMA_DEF



#include "bdObject.h"
#include "bdBIP.h"
#include "bdGIP.h"
#include "bdCurveXY.h"




/// Color info and methods used for mapping colors (C) feature.

class bdColorCode
{
public:
    unsigned int m_red;
    unsigned int m_green;
    unsigned int m_blue;
    
    // this member is to be used as a helper (if need be).
    unsigned int m_intensity;
    
    bdColorCode(){ m_red = m_green = m_blue = 0; m_intensity = 0; };
    ~bdColorCode(){};
    
    void Set(unsigned int r, unsigned int g, unsigned int b, unsigned int intensity=0)
    {
        m_red = r;
        m_green = g;
        m_blue = b;
        m_intensity = intensity;
    };
    
    int operator==(bdColorCode &cc)
    {
        if(cc.m_red==m_red && cc.m_blue==m_blue && cc.m_green==m_green) return 1;
        return 0;
    };
    
    /// Distance from the given color
    double L2(unsigned int r, unsigned int g, unsigned int b)
    {
        return ( sqrt( (m_red-r)*(m_red-r) + (m_green-g)*(m_green-g) + (m_blue-b)*(m_blue-b)) );
    };
    
    double L1(unsigned int r, unsigned int g, unsigned int b)
    {
        double dr = (m_red-r);
        if(dr<0) dr = -dr;
        double dg = (m_green-g);
        if(dg<0) dg = -dg;
        double db = (m_blue-b);
        if(db<0) db = -db;

        return ( dr+dg+db );
    };
    
    double Lmax(unsigned int r, unsigned int g, unsigned int b)
    {
        double d = (m_red-r);
        if(d<0) d = -d;
        double dg = (m_green-g);
        if(dg<0) dg = -dg;
        if(d<dg) d=dg;
        double db = (m_blue-b);
        if(db<0) db = -db;
        if(d<db) d=db;
        
        return d;
    };
    
    double L2order(unsigned int r, unsigned int g, unsigned int b)
    {
        if((m_red>m_green && r<g) || ((m_red<m_green && r>g))) return 400;
        if((m_red>m_blue && r<b) || ((m_red<m_blue && r>b))) return 400;
        if((m_blue>m_green && b<g) || ((m_blue<m_green && b>g))) return 400;
        
        return ( sqrt( (m_red-r)*(m_red-r) + (m_green-g)*(m_green-g) + (m_blue-b)*(m_blue-b)) );
    };
};




///  Algorithms for calculation of ABCDE features of melanoma.

class bdMelanoma : public bdFunctionObject
{
public:
    
    //--- ASYMMETRY SCORES ---
    double m_asymmetry_soft_score_axis_1;
    double m_asymmetry_soft_score_axis_2;
    double m_asymmetry_hard_score_axis_1;
    double m_asymmetry_hard_score_axis_2;
    double m_asymmetry_soft_score;
    double m_asymmetry_hard_score;
    //------
    
    //--- BORDER SCORES ---
    double m_border_soft_score;
    double m_border_hard_score;
    //------
    
    //--- COLOR SCORES ---
    double m_color_brown_area;
    double m_color_brown_volume;
    double m_color_red_area;
    double m_color_red_volume;
    double m_color_blue_gray_area;
    double m_color_blue_gray_volume;

    double m_color_mapped_black_volume;
    double m_color_mapped_white_volume;
    double m_color_mapped_light_brown_volume;
    double m_color_mapped_dark_brown_volume;
    double m_color_mapped_red_volume;
    double m_color_mapped_blue_gray_volume;

    double m_color_soft_score;
    double m_color_hard_score;
    //------
    
    //--- DIFFERENTIAL STRUCTURES SCORES ---
    double m_DS_globules;
    double m_DS_dots;
    double m_DS_streaks;
    double m_DS_structureless_area;
    double m_DS_network;
    
    double m_DS_soft_score;
    double m_DS_hard_score;
    //-------------------------------------

    
    
    /// Checks if the pixel with coordinates r,c is masked in the 'mask' image. The input image is RGB and a pixel is considered masked if it is masked by at least 2 channels.
    int IsMasked(bdImage &mask, int r, int c);
    
    /// from the segmented melanoma image (which is an RGB image, where each channel is segmented) creates a mask image (single slice scalar pixel value image) using IsMasked() method.
    int CreateMaskImage(bdImage &segmented_image, bdImage &mask_image);
    
    /// Makes a projection of a 'position_to_project' on the vector. The output is the projected position and coefficient 'n' which defines where on the
    /// vector was the projection made (if the projection falls inside the vector, 'n' is in range [0,1], otherwise it falls outside. If 'n'=0.5, the
    /// projected position is exactly in the middle of the vector).
    int ProjectPositionToLine(double *line_pos1_src, double *line_pos2_src, double *position_to_project_src, double *output_projected_position_src, double &n);

    /// Fing in each channel histogram the peak value. Modify gray values in channels so that the histogram peaks fall at the same pixel gray value.
    int EqualizeBrightnessOfChannels(bdImage &input, bdImage &output);
    
    /// Performs dilation on (multiple) binary 2d imagesusing a circular structuring element (SE) of given squared radius.
    int DilateBinaryImageWithCircleSE(bdImage &input, bdImage &output, unsigned int squared_radius);
    
    /// Compute histogram of a 2d image slice (for given time t and slice s indexes) and within range [range_min, range_max].
    int Histogram(bdImage &input, bdArray<int> &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t=0, unsigned int s=0);
    
    /// Compute histogram of a 2d image slice (for given time t and slice s indexes) and within range [range_min, range_max].
    int Histogram(bdImage &input, bdCurveXY &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t=0, unsigned int s=0);
    
    /// Compute histogram of a masked part of the 2d image slice (for given time t and slice s indexes) and within range [range_min, range_max].
    int Histogram(bdImage &input, bdImage &mask, bdArray<int> &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t_input=0, unsigned int s_input=0, unsigned int t_mask=0, unsigned int s_mask=0);

    /// Compute histogram of a masked part of the 2d image slice (for given time t and slice s indexes) and within range [range_min, range_max].
    int Histogram(bdImage &input, bdImage &mask, bdCurveXY &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t_input=0, unsigned int s_input=0, unsigned int t_mask=0, unsigned int s_mask=0);
    
    /// Find the maxima in the histogram and compute for which area they are the dominant maximum - values are stored in the output maxima array (size of the histogram).
    int HistogramMaxima(bdArray<int> &input_histogram, bdArray<int> &output_maxima);
    
    /// For given input histogram make its envelope in range [start_index,end_index] with window size window_size.
    /// NOT TESTED.
    int EnvelopeOfHistogram(bdArray<int> &input_histogram, unsigned int start_index, unsigned int end_index, unsigned int window_size, bdArray<int> &output_envelope);

    /// Calculate the number of connected components and get the size of the component which is the closest to the input seed indexes.
    int ConnectedComponents_2D_WithSeedIndexes(bdImage &input_image, int threshold, unsigned int t, bdDiscreteCoordinates3D &seed_indexes, int &output_number_of_components, int &output_region_size);
    
    /// Each component of RGB image will be leveled to input level value (e.g. if level=25 all values above 25 will become 25).
    int LevelComponentsOfRGBImage(bdImage &inputRGB, unsigned int level_value, bdImage &output_leveled_RGB);
    
    /// For each RGB component separately perform mean (average) operator using the circular structuring element of given squared radius.
    /// The mean operator is performed only on voxels under mask.
    int MeanOperatorPerRGBComponent(bdImage &inputRGB, bdImage &mask, unsigned int squared_circular_SE_radius, bdImage &output_mean_RGB);
    
    /// Performs automatic normalization in range [lower_value, upper_value] for each color component of the image separately.
    int NormalizeComponentsOfRGBImage(bdImage &inputRGB, unsigned int lower_value, unsigned int upper_value, bdImage &output_normalized_RGB);

    /// Performs automatic normalization in range [lower_value, upper_value] for all color components together.
    int NormalizeWholeRGBRangeImage(bdImage &inputRGB, unsigned int lower_value, unsigned int upper_value, bdImage &output_normalized_RGB);
    
    /// For all voxels within the mask sum up voxel values over all components. Use this score to compare the level of information left in the subtracted images.
    double SumOfVoxelValuesForAllRGBComponents(bdImage &inputRGB, bdImage &mask);
    
    /// Downsample the input image by turning the window of size size_of_new_pixel*size_of_new_pixel into a single pixel in the output image.
    int Downsample(bdImage &input, bdImage &output, int size_of_new_pixel);
    
    /// Downsample the input image by averaging the pixel values in window of size size_of_new_pixel*size_of_new_pixel.
    int Downsample_Mean(bdImage &input, bdImage &output, int size_of_new_pixel);
    
    /// Downsample the input image by using max pixel value in window of size size_of_new_pixel*size_of_new_pixel.
    int Downsample_Max(bdImage &input, bdImage &output, int size_of_new_pixel);

    
    /// Using the downsampled image (obtained by bdMelanoma::Downsample() method) perform masking of the input image.
    int MaskWithDownsampledImage(bdImage &input, bdImage &downsampled_mask, bdImage &output, int size_of_downsampled_pixel);

    /// Fills in all the wholes in the segmentation result with values found in the original image.
    int FillInSegmentedMelanoma(bdImage &original_image, bdImage &segmented_image, bdImage &output);
    
    
	/// Segment melanoma in input 2d+Time image where each time instance is also image of different channel (images are channels of a multispectral image or just channels from a RGB image).
    /// We examine histogram to see which of the images is best suited for thresholding. The threshold is set similarily to Otsu method. The largest Conn.Comp. is extracted as the final segmentation.
    int SegmentMelanoma_Thresholding(bdImage &input, bdImage &output);
    
    /// Segment melanoma in input 2d+Time image where each time instance is also image of different channel (images are channels of a multispectral image or just channels from a RGB image).
    /// In the white light image select the channel with the highest variance. Find the most prominent (maximum) peak in the histogram. Search for the steepest slope before the maximum
    /// (highest difference between the 2 neighboring values in the histogram). Calculate where the line defined by 2 points of the steepest slope crosses the x axis. The crossing is the
    /// threshold value. The connected component that is the closest to the image center is the one that is extracted.
    int SegmentMelanoma_Thresholding2(bdImage &input, bdImage &output);

    /// In the white light image select the channel with the highest variance. Find the most prominent (maximum) peak in the histogram. Search for the steepest slope before the maximum
    /// (highest difference between the 2 neighboring values in the histogram). Calculate where the line defined by 2 points of the steepest slope crosses the x axis. The crossing is the
    /// threshold value. The connected component that is the closest to the image center is the one that is extracted. The extracted image is downsampled in such way that the new pixel
    /// size is larger than the radius of hairs (possibly found in the images). The downsampled image has non-zero value if the corresponding window in the input image had more than 2/3
    /// of non-zero pixels in it. The largest connected component in the downsampled image is extracted and dilated with 4-neighborhood. This result is used as a mask to the extracted
    /// image before downsampling.
    int SegmentMelanoma_Thresholding3(bdImage &input, bdImage &output);

    /// In the white light image select the channel with the highest variance. Find the most prominent (maximum) peak in the histogram. Search for the steepest slope before the maximum
    /// (highest difference between the 2 neighboring values in the histogram). Calculate where the line defined by 2 points of the steepest slope crosses the x axis. The crossing is the
    /// threshold value. The connected component that is the closest to the image center is the one that is extracted. The extracted image is downsampled in such way that the new pixel
    /// size is larger than the radius of hairs (possibly found in the images). The downsampled image has non-zero value if the corresponding window in the input image had more than 2/3
    /// of non-zero pixels in it. The largest connected component in the downsampled image is extracted and dilated with 4-neighborhood. This result is used as a mask to the extracted
    /// image before downsampling. The method calculates 2 segmentations from the 2 calculated threshold values (one is with the steepest rising slope to the highest peak in histogram,
    /// the other one is the steepest rising/falling slope between the highest peak and the second highest peak.)
    int SegmentMelanoma_Thresholding4(bdImage &input, bdImage &output, bdImage &output2);
    
    /// In the white light image select the channel with the highest variance. Find the most prominent (maximum) peak in the histogram. Search for the steepest slope before the maximum
    /// (highest difference between the 2 neighboring values in the histogram). Calculate where the line defined by 2 points of the steepest slope crosses the x axis. The crossing is the
    /// threshold value. The connected component that is the closest to the image center is the one that is extracted.
    /// If the (extracted) segmented region is too small, the other threshold value (steepest slope towards the nextmaximum peak in the histogram) is used.
    /// If the (extracted) segmented region touches the image border perform masking as described bellow.
    /// The extracted image is downsampled in such way that the new pixel
    /// size is larger than the radius of hairs (possibly found in the images). The downsampled image has non-zero value if the corresponding window in the input image had more than 2/3
    /// of non-zero pixels in it. The largest connected component in the downsampled image is extracted and dilated with 4-neighborhood. This result is used as a mask to the extracted
    /// image before downsampling. The method calculates 2 segmentations from the 2 calculated threshold values (one is with the steepest rising slope to the highest peak in histogram,
    /// the other one is the steepest rising/falling slope between the highest peak and the second highest peak.)
    int SegmentMelanoma_Thresholding5(bdImage &input, bdImage &output);
    

    /// In the white light image select the channel with the highest variance. Find the most prominent (maximum) peak in the histogram. Search for the steepest slope before the maximum
    /// (highest difference between the 2 neighboring values in the histogram). Calculate where the line defined by 2 points of the steepest slope crosses the x axis. The crossing is the
    /// threshold value. The connected component that is the closest to the image center is the one that is extracted.
    /// If the (extracted) segmented region is too small, the other threshold value (steepest slope towards the nextmaximum peak in the histogram) is used.
    /// If the (extracted) segmented region touches the image border perform masking as described bellow.
    /// The extracted image is downsampled in such way that the new pixel
    /// size is larger than the radius of hairs (possibly found in the images). The downsampled image has non-zero value if the corresponding window in the input image had more than 2/3
    /// of non-zero pixels in it. The largest connected component in the downsampled image is extracted and dilated with 4-neighborhood. This result is used as a mask to the extracted
    /// image before downsampling. The method calculates 2 segmentations from the 2 calculated threshold values (one is with the steepest rising slope to the highest peak in histogram,
    /// the other one is the steepest rising/falling slope between the highest peak and the second highest peak.)

    /// The method result is filled in (to surpass getting a porous segmentation). It can go wrong with very large lesions (that are larger than skin region).
    int SegmentMelanoma_Thresholding6(bdImage &input, bdImage &output);

    
    /// Based on version 5, just uses combined image of RGB channels to perform segmentation.
    int SegmentMelanoma_Thresholding7(bdImage &input, bdImage &output);
 
    
    /// Based on version 5, just uses the BLUE channel and discards zero values in the histogram.
    int SegmentMelanoma_Thresholding8(bdImage &input, bdImage &output);

    
    /// Based on version 5, just uses the B channel and the bending point (foot of the curve) between the peaks in the histogram (discards zero-valued entries in the histogram).
    int SegmentMelanoma_Thresholding9(bdImage &input, bdImage &output);


    
    /// Segments the inner structure of presegmented melanoma by increasing the threshold and finding 'white' isolated connected components.
    int SegmentMelanomaInnerStructure(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    /// Finds the local maxima values.
    int SegmentMelanomaInnerStructure2(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    /// Segments the brighter inner structures on a local level of certain size. Makes iterative thresholding and closing the thresholded structure to see which parts are left as local
    /// bright spots (the ones that are closed in the closed image but are not in the thresholded image)
    int SegmentMelanomaInnerStructure3(bdImage &original, bdImage &segmented_melanoma, bdImage &output);
  
    /// Segments the brighter inner structures on a local level for connected components of any size. Makes iterative thresholding and considers isolated connected components (the ones not
    /// touching lesion egde) as local bright spots.
    int SegmentMelanomaInnerStructure4(bdImage &original, bdImage &segmented_melanoma, bdImage &output);
    
    /// Segments the brighter inner structures on a local level for connected components of any size. Makes iterative thresholding and considers isolated connected components (the ones not
    /// touching lesion egde) as local bright spots. Improvement in speed compared to the previous version by checking if a region is isolated by examining its edge (grow edge, see if it's
    /// larger or smaller than the actual edge).
    int SegmentMelanomaInnerStructure5(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    int SegmentMelanomaInnerStructure6(bdImage &original, bdImage &segmented_melanoma, bdImage &output);
    
    
    /// Same algorithm as DARK v. 2 just with MAX operator when cobining the images.
    int SegmentMelanomaInnerStructure7(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    
    /// Same algorithm as SegmentStreaks with adjusted percentage and detection of network are by its density.
    int SegmentMelanomaInnerStructure8(bdImage &original, bdImage &segmented_melanoma, bdImage &output);


    
    /// Segments the inner DARK structure of presegmented melanoma by increasing the threshold and finding 'dark' isolated connected components.
    int SegmentMelanomaInnerDarkStructure(bdImage &original, bdImage &segmented_melanoma, bdImage &output);
    
    /// Segments the inner DARK structure of presegmented melanoma by increasing the threshold and finding 'dark' isolated connected components. Extracts a single combined gray-scale slice
    /// to perform segmentation on.
    int SegmentMelanomaInnerDarkStructure2(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    /// Segments the inner DARK structure of presegmented melanoma by using SegmentStructures Extracts a single gray-scale slice
    /// to perform segmentation on.
    int SegmentMelanomaInnerDarkStructure3(bdImage &original, bdImage &segmented_melanoma, bdImage &output);

    
    /// Segments the structures (both dark and white, depending on the setting of 'percentage' parameter).
    /// Pertorm segmentation on a single gray-scale slice given by index 't'.
    int SegmentStructures(bdImage &original, unsigned int t, bdImage &segmented_melanoma, unsigned int SE_squared_radius, unsigned int percentage, bdImage &output, unsigned short segmented_label=255);

    
    /// Segments the inner DARK structure as streaks using SegmentStructures method. Extracts a single gray-scale slice
    /// to perform segmentation on.
    int SegmentStreaks(bdImage &original, bdImage &segmented_melanoma, unsigned int percentage, bdImage &output);

    /// Segments the structureless area as the area of signifficant blur.
    int SegmentStructurelessArea(bdImage &original, bdImage &segmented_melanoma, unsigned int percentage, bdImage &output);

    /// Segments the structureless area as the area of of at least 10% of the whole in between segmented structures.
    int SegmentStructurelessArea2(bdImage &mask, bdImage &border_image, bdImage &dark_structures_image, bdImage &bright_structures_image, bdImage &streaks_image, bdImage &output);


    /// Calculates the A (Asymmetry) coefficient from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// For a given number of axes, asymmetry is ratio of the number of asymmetric pixels against the number of pixels in the union of overlapping halves of the melanoma for all axes.
    int AsymmetryOfMelanoma(bdImage &original, bdImage &mask, unsigned int number_of_symmetry_axes, bdArray<double> &output_min_symmetry_coeff, bdArray<double> &output_max_symmetry_coeff, bdArray<double> &output_mean_symmetry_coeff);
  
    /// Takes lists of colors and compares them to determine the asymmetry of colors in these regions. Each bdColorCode in these lists has m_intensity value equal to the number of pixels with that color in the region.
    /// The output text is the one used in the print out (use "structure" if you analyze structures and "color" if you analyze colors).
    /// 'max_diff' parameter is updated to a new value if the newly calculated value is higher than current max_diff.
    int AsymmetryOfColorsInRegionLists_Helper(bdList<bdColorCode> &list_of_colors_in_region1, bdList<bdColorCode> &list_of_colors_in_region2, double n_of_pixels_in_region1, double n_of_pixels_in_region2, double &max_diff, const char *text_for_output);
    
    /// Calculates the A (Asymmetry) coefficient from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// For a given number of axes, asymmetry is ratio of the number of asymmetric pixels against the number of pixels in the union of overlapping halves of the melanoma for all axes.
    int AsymmetryOfMelanoma2(bdImage &original, bdImage &mask, bdImage &mapped_chart_colors_image, bdImage &combined_structure_image, unsigned int number_of_symmetry_axes, bdImage &output);

    
    
    /// Calculates B (Border) coefficient from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// Follows the (largest found) border of the melanoma and in sequential way along the border calculates the contrast of averaged pixels within a circle of 'squared_radius' size.
    int BorderOfMelanoma(bdImage &original, bdImage &mask, int squared_radius, bdArray<int> &output_contrast_along_border);
    
    /// Calculates B (Border) coefficient from segmented (masked) melanoma image (must contain segmentation for each channel, so R,G,B).
    /// Follows the (largest found) border of the melanoma and in sequential way along the border calculates the minimum squared distance to the 2 other borders.
    int BorderOfMelanoma2(bdImage &segmented_melanoma, bdArray<int> &output_contrast_along_border, bdImage &output);
    
    
    /// Calculates B (Border) coefficient from segmented and mask melanoma image (must contain segmentation for each channel, so R,G,B; mask is just a single scalar slice).
    /// Follows the mask border of the melanoma and in sequential way along the border calculates the minimum squared distance to the 2 other borders.
    /// If the cumulative squared distance to RGB borders is smaller than 'threshold_squared_distance', the border is considered to be sharp.
    int BorderOfMelanoma3(bdImage &segmented_melanoma, bdImage &mask, bdArray<int> &output_contrast_along_border, unsigned int threshold_squared_distance, bdImage &output);


    
    /// Calculates C (Color) feature from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// The method uses l2 norm to map real colors to pre-defined set of colors. The best matching color is added to output image.
    int ColorsOfMelanoma(bdImage &original, bdImage &mask, bdImage &output_mapped_colors_image);
    
    /// Calculates C (Color) feature from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// Takes as input aslo color charts for the searched colors.
    /// The method uses l2 norm to map real colors to pre-defined set of colors. The best matching color is added to output image.
    int ColorsOfMelanoma2(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &output_mapped_colors_image);
    
    /// Calculates C (Color) feature from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// Takes as input color charts for the searched colors and color maps that indicate how the charts map to searched colors.
    /// The method uses l2 norm to map real colors to pre-defined set of colors. The best matching color is added to output image.
    int ColorsOfMelanoma3(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image);

    /// Calculates C (Color) feature from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// Takes as input color charts for the searched colors and color maps that indicate how the charts map to searched colors.
    /// The method defined the chart to search from by looking at relation between RGB values and only after uses l2 norm to map real colors to pre-defined set of colors. The best matching color is added to output image.
    int ColorsOfMelanoma4(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image);

    /// Calculates C (Color) feature from input original and masked melanoma image (masked image is the segmentation image of the melanoma).
    /// Takes as input color charts for the searched colors and color maps that indicate how the charts map to searched colors.
    /// The method defined the chart to search from by looking at relation between RGB values and only after uses l2 norm to map real colors to pre-defined set of colors. The best matching color is added to output image.
    int ColorsOfMelanoma5(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image);
    
    /// Calculate scores for the D feature based on segmented images.
    int Quantify_D_Feature(bdImage &mask, bdImage &segmented_melanoma, bdImage &dark_structures_image, bdImage &bright_structures_image, bdImage &streaks_image, bdImage &structureless_area_image);


    /// Calculates D (!!!Diameter!!!) feature from input border melanoma image (obtained from BorderOfMelanoma() method).
    /// The method uses measures minimum, maximum and mean diameter given in pixels (spacing of pixels is not taken into account).
    int DiameterOfMelanoma(bdImage &border_image, double *output_min_diameter_in_pixels, double *output_max_diameter_in_pixels, double *output_mean_diameter_in_pixels);
    
    /// Subtract the input images and record the absolute values in the output image. Considers both input images to be 2D images.
    int SubtractImages(bdImage &image1, bdImage &image2, bdImage &mask, bdImage &output_image);
    
    /// For visualization makes a white RGB image from segmented_melanoma as RED, melanoma_bright_structure as GREEN, melanoma_dark_structure as BLUE. Combinations are possible.
    int MakeCombinedSegmentedImage(bdImage &segmented_melanoma, bdImage &melanoma_bright_structure, bdImage &melanoma_dark_structure, bdImage &melanoma_dark_streaks, bdImage &output_combined_image);
    
    /// Extract the RED channel from input white image.
    int Extract_RedImage(bdImage &white_image, bdImage &output_image);
    
    /// Extract the GREEN channel from input white image.
    int Extract_GreenImage(bdImage &white_image, bdImage &output_image);

    /// Extract the BLUE channel from input white image.
    int Extract_BlueImage(bdImage &white_image, bdImage &output_image);
    
    /// Extract the max of GREEN and BLUE channels from input white image.
    int Extract_MaxOfGreenAndBlueImage(bdImage &white_image, bdImage &output_image);
    
    /// Extract the max of RED, GREEN and BLUE channels from input white image.
    int Extract_MaxOfRGBImage(bdImage &white_image, bdImage &output_image);

    /// for each pixel new value is |red-blue| + |blue-green| + |green-red|.
    int Extract_DifferenceOfRGB_Image(bdImage &white_image, bdImage &output_image);
    
    /// Results are printed to screen and to a text file with given file name.
    int PrintResults(const char *file_name);


};

#endif