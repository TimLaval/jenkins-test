/**************************************************************************
 Class of operations defined on binary images

 Author: Danilo Babin
 File name: "bdBIP.h"
***************************************************************************/


#if defined(BD_BIP_USE_SOURCE_CODE)
    #define BD_BIP_API
#else
    #if defined(_MSC_VER) //  Microsoft
        #ifdef BD_BIP_EXPORTS
            #define BD_BIP_API __declspec(dllexport)
        #else
            #define BD_BIP_API __declspec(dllimport)
        #endif
    #else // GCC
        #ifdef BD_BIP_EXPORTS
            #define BD_BIP_API __attribute__((visibility("default")))
        #else
            #define BD_BIP_API
        #endif
    #endif
#endif


#ifndef BD_BIP_DEF
	#define BD_BIP_DEF







#include "bdObject.h"
#include "bdImage.h"
#include "bdGeometry.h"
#include "bdVoxel.h"
#include "bdList.h"
#include "bdRegion3D.h"




//#pragma warning(disable: 4661)



// Class of operations on 3D image set. CAUTION: THIS DOES NOT MEAN THAT THEY PROCESS ONLY SPACE
// DIMENSIONS, THEY CAN ALSO PROCESS 2D TIME SEQUENCES IF THE INPUT IMAGE IS DEFINED IN THAT WAY!!!

class BD_BIP_API bdBIP : public bdFunctionObject
{
public:
   
    /// Like the center of mass, just every non-zero pixel has the same weight.
    int CenterOfObjectVolume(bdImage &image, unsigned int &output_center_s, unsigned int &output_center_r, unsigned int &output_center_c, unsigned int t=0);
	
	/// Area of the 3D object in the image.
	int AreaOfObject_3D(bdImage &image, unsigned int t=0);

	/// Area of the 2D object in the given slice of the image.
	int AreaOfObject_2D(bdImage &image, unsigned int slice_index, unsigned int t=0);

	/// Threshold the image with the given threshold value.
	int Threshold(bdImage &input_image, int threshold_value, bdImage &output_image);

	/// Get number of non-zero value voxels in the 26/27 neighborhood.
	int NumberOfNonZeroVoxelsIn26Neighborhood(bdImage &image, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVoxelsIn26Neighborhood(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVoxelsIn27Neighborhood(bdImage &image, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVoxelsIn27Neighborhood(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c);

    /// For given input seed indexes and image finds the closest non-zero voxel in the image and gives its indexes as output. If none found return fail 0, else success 1.
    int GetClosestNonZeroVoxelToSeedIndexes(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c, unsigned int &out_s, unsigned int &out_r, unsigned int &out_c);
    int GetClosestNonZeroVoxelToSeedIndexes(bdImage &image, unsigned int s, unsigned int r, unsigned int c, unsigned int &out_s, unsigned int &out_r, unsigned int &out_c);
    

	/// Extracts the largest found 26-connected component to output image.
	/// !!!!!NEEDS TO BE ADJUSTED TO 4D PROPERLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int ExtractLargest_26_ConnectedComponent(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t=0);

    //NOT IMPLEMENTED YET!!!
//    /// Extracts the CC(26-neighborhood) that is the largest and closest to the input seed. The proximity to seed acts as a coefficient as: distance_from seed/max_distance.
//    int Extract_CC26_LargestAndClosestToSeed(bdImage &input_image, int threshold, bdImage &output_image, unsigned int seed_s, unsigned int seed_r, unsigned int seed_c, double max_distance, unsigned int t=0);

	/// Mark every connected component in the image (CCs are made by connected voxels larger than 'threshold'.
	int ConnectedComponents_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t=0, int *output_number_of_components = NULL);

    /// Extract the 26-connected component from Seed points.
    int Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, bdList<bdDiscreteCoordinates3D> &seeds, unsigned int t=0);

    
	/// Extract the 26-connected component from Seed point. 
	/// NOT IMPELEMENTED YET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int s, unsigned int r, unsigned int c);

	/// Calculate squared radius of largest inscribed sphere in segmented image for voxel positions
	/// of the input_skeleton_image (which acts as a mask).
	int RadiusSquared_Sphere(bdImage &input_skeleton_image, bdImage &input_segmented_image, bdImage &output_image);

	/// Calculate squared radius of largest inscribed circle per 2D slice in segmented image for voxel positions 
	/// of the input_mask_image (which acts as a mask).
	int RadiusSquared_CirclePer2DSlice(bdImage &input_mask_image, bdImage &input_segmented_image, bdImage &output_image);

		
	/// Calculates the 3D radius in pixels of Skeletonized image, where the radius size is searched in Segmented image. This function returns the 
	/// number of pixels for the grown radius (e.g. for squared radius 1 max number of pixels is 6 (6-neighborhood). If any of those 6 pixels is 
	/// a zero value pixel in the segmented image, the radius will not be furhter grown and the output pixel value will be the number of non-zero 
	/// pixels till now, for example in case of only one zero pixel, it will be 5. If all 6 pixels in the segmented image are non-zero, the radius 
	/// is increased to 2 and number of non-zero pixels in this radius is added to existing number of non-zero pixels).
	/// The 'offset_percent' (must be in range [1,100]) relaxes the above mentined condition (e.g. for 50% it means that for squared radius 1 (6-neighborhood)
	/// only 3 pixels need to be non-zero for the radius to continue growing).
	int ProfileVolume_Sphere(bdImage &input_mask_image, bdImage &input_segmented_image, bdImage &output_image, unsigned int offset_percent=100);
	
	//int Calculate3DSquaredRadiusWithNumberOfElementsAndGrayImage_Sphere(Image3D<T> &input_original, Image3D<T> &input_skeleton_image, Image3D<T> &input_segmented_image, Image3D<T> &output_image);
	//void Calculate3DCubeSideHalfSize(Image3D<T> &input_skeleton_image, Image3D<T> &input_segmented_image, Image3D<T> &output_image);
	//void InvertBinary(Image3D<T> &input, Image3D<T> &output, int object_pixel_value = 255);
	
	/// Morphological dilation for a BINARY image.
	int Dilation_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

	/// Morphological dilation for a BINARY image which resizes the image to produce certainly a proper result (not affected by image size).
	int DilationWithImageResize_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

	/// Morphological erosion for a BINARY image.
	int Erosion_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

    /// Morphological erosion for a BINARY image on a masked region using circular SE.
    int Erosion_Circle(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius);
    
    /// Morphological erosion for a BINARY image on a masked region using circular SE. The whole SE is recorded (instead of only the pixel that fullfils requirements for output).
    int ErosionByMapping_Circle(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius);

    /// Morphological erosion for a BINARY image on a masked region using circular SE. The whole SE is recorded (instead of only the pixel that fullfils requirements for output).
    int ErosionByMapping_Circle2(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius, unsigned int percentage);

    
	/// Morphological opening for a BINARY image.
	int Opening_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

	/// Morphological closing for a BINARY image.
	int Closing_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

	/// Morphological closing for a BINARY image which uses dilation with image resizing to always produce accurate result (not affected by image size).
	int ClosingAccurate_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius);

};


//=======================================================================================================






#endif
