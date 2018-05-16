/**************************************************************************
 Template class for processing gray images.

 Author: Danilo Babin
 File name: "bdGIP.h"
***************************************************************************/



#if defined(BD_GIP_USE_SOURCE_CODE)
    #define BD_GIP_API
#else
    #if defined(_MSC_VER) //  Microsoft
        #ifdef BD_GIP_EXPORTS
            #define BD_GIP_API __declspec(dllexport)
        #else
            #define BD_GIP_API __declspec(dllimport)
        #endif
    #else // GCC
        #ifdef BD_GIP_EXPORTS
            #define BD_GIP_API __attribute__((visibility("default")))
        #else
            #define BD_GIP_API
        #endif
    #endif
#endif


#ifndef BD_GIP_DEF
	#define BD_GIP_DEF


#include "bdObject.h"
#include "bdImage.h"
#include "bdGeometry.h"
#include "bdVoxel.h"
#include "bdList.h"
#include "bdRegion2D.h"




//#pragma warning(disable: 4661)



//This is the class of operations on gray 3D images. CAUTION: THIS DOES NOT MEAN THAT THEY PROCESS ONLY SPACE
// DIMENSIONS, THEY CAN ALSO PROCESS 2D TIME SEQUENCES IF THE INPUT IMAGE IS DEFINED IN THAT WAY!!!

class BD_GIP_API bdGIP3D : public bdFunctionObject
{
public:



};


//----------------------------------------------------------------------------------------------------------------------------------------------------------------




class BD_GIP_API bdGIP : public bdGIP3D
{
public:
    
    /// Maximum operator on 2D circle structuring element working on each slice of the image.
    int MaximumCircleFor2DSlice(bdImage &input, bdImage &mask, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask, unsigned int circle_squared_radius, bdImage &output);
    
    /// Mask the input image with the mask image. 'mask' is a single slice image, while input may be a multi-slice image.
    int MaskPer_2D_Slice(bdImage &input, bdImage &mask, bdImage &output);

    /// Mean operator on 2D circle structuring element working on each slice of the image.
    int MeanCircleFor2DSlice(bdImage &input, bdImage &mask, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask, unsigned int circle_squared_radius, bdImage &output);

	/// Invert voxel values in the image (make a 'negative') for the whole image (4D).
	int Negative(bdImage &input, bdImage &output);

	/// Rescales the image voxel value range to the input range. 
	int RescaleWholeRange(bdImage &inut_image, unsigned int lower_value, unsigned int upper_value, bdImage &output_image);

};






#endif
