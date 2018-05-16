/**************************************************************************
 Template bdImageT for images up to 4D data

 Author: Danilo Babin
 File name: "bdImage.h"
***************************************************************************/




//To build as DLL, add:" /D "BD_IMAGE_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_IMAGE_EXPORTS
		#define BD_IMAGE_API __declspec(dllexport) 
	#else
		#define BD_IMAGE_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_IMAGE_EXPORTS
		#define BD_IMAGE_API __attribute__((visibility("default")))
	#else
		#define BD_IMAGE_API
	#endif
#endif




#ifndef BD_IMAGE_DEF
	#define BD_IMAGE_DEF



#include <typeinfo>
#include "bdArray.h"
#include "bdObject.h"
#include "bdRegularGrid.h"
#include "bdImageSliceTag.h"



#pragma warning(disable: 4661)



//template class BD_IMAGE_API bdArray<bdImageSliceTag>;


/// Image template class.

template<class T>
class BD_IMAGE_API bdImageT : public bdRegularGridT<T>
{
protected:

	/// Tags from DICOM files. Number of tags is 'number_of_slices * number_of_time_series'.
	bdArray<bdImageSliceTag> m_tags;
    bdImageSliceTag m_tag_default; // default tag to show if no other exist.

public:
	
	/// Constructor.
	bdImageT();

	/// Destructor.
	~bdImageT();

	/// Copy entire image from input image.
	virtual void CopyFrom(bdImageT<T> &m);

	/// Calls CopyFrom() method.
	bdImageT<T>& operator =(bdImageT<T> &m);

	/// Print type and header info.
	ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	//ostream& PrintInfo(ostream &o);


//	int IsEqualSpacingAs(bdImageT<T> &image, double offset = 0.0001);


	/// Get the full range of voxel values: for 8U it is [0,255], for 16U it is [0, 65535].
	void GetVoxelValueFullRange(unsigned int *full_range_min, unsigned int *full_range_max);

	/// Get the minimum and maximum voxel value in the image. If image is empty, returns fail (0).
	int GetVoxelMinimumAndMaximumValue(unsigned int *min, unsigned int *max);
	

	/// Get voxel value for input world coordinates. If coordinates fall outside the image return 0, else 1.
	int GetVoxelValueForWorldCoordinates(double w_t, double w_s, double w_r, double w_c, T &out_value);
	int GetVoxelValueForWorldCoordinates(double w_s, double w_r, double w_c, T &out_value);

	/// Get tag with index values t,s.
	bdImageSliceTag& Tag(unsigned int t, unsigned int s);

	/// Get tag with index value s (assumes t=0).
	bdImageSliceTag& Tag(unsigned int s);

	/// Set the size of the grid with properties to DEFAULT values: spacing:(1,1,1,1), origin:(0,0,0,0)...
	/// If the grid is already the requested size, does not do anything.
	virtual int SetSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	virtual int SetSizeOf3DGrid(unsigned int s, unsigned int r, unsigned int c);
	virtual int SetSizeOf2DTimeSeries(unsigned int t, unsigned int r, unsigned int c);
	
	///Set only size of the grid, properties remain as before
	virtual void SetSizeOnlyAs(bdImageT<T> &grid);
	
	///Set size and properties of the grid
	virtual void SetSizeAndPropertiesAs(bdImageT<T> &grid);

	/// Releases the data and resets the image properties to the default values.
	/// Releases data (or not) based on '_is_data_created_internally' indicator.
	virtual void Reset();

	/// Count the number of non-zero value voxels in 8-neighborhood of voxel with coordinates (t,s,r,c).
	int NumberOfNonZeroVexelsIn_8_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVexelsIn_8_Neighborhood(unsigned int s, unsigned int r, unsigned int c);

	/// Count the number of non-zero value voxels in 9-neighborhood of voxel with coordinates (t,s,r,c).
	int NumberOfNonZeroVexelsIn_9_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVexelsIn_9_Neighborhood(unsigned int s, unsigned int r, unsigned int c);

	/// Count the number of non-zero value voxels in 26-neighborhood of voxel with coordinates (t,s,r,c).
	int NumberOfNonZeroVexelsIn_26_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVexelsIn_26_Neighborhood(unsigned int s, unsigned int r, unsigned int c);

	/// Count the number of non-zero value voxels in 27-neighborhood of voxel with coordinates (t,s,r,c).
	int NumberOfNonZeroVexelsIn_27_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	int NumberOfNonZeroVexelsIn_27_Neighborhood(unsigned int s, unsigned int r, unsigned int c);


	/// Sort slices by given image slice tag criteria. The order of criteria defines their priority - slices are first sorted using criteria1,
	/// then using criteria2 (but maintining order of criteria1) and finaly using criteria3 (but maintaining order of criteria1 and criteria2).
	//void SortSlices(const char *tag_criteria1, const char *tag_criteria2, const char *tag_criteria3);
    virtual void SortSlices(const char *criteria1_tag_number1, const char *criteria1_tag_number2, const char *criteria2_tag_number1, const char *criteria2_tag_number2, const char *criteria3_tag_number1, const char *criteria3_tag_number2);

	/// Load slice tags from a file. If 'is_loading_image_properties' is true (non-zero) the number of tags (tag array size) must be set prior to calling 
	/// this method (use SetSize or similar method). If it is false (zero), the method sets size and properties of the image.
	/// !!!!! NOT CCOMPLETED !!!!!
	int LoadInfoAndTags(const char *file_name, int is_loading_image_properties = 0);

	/// Save slice tags to a file. 
	void SaveInfoAndTags(const char *file_name);

    /// Save an image to a raw 16 bit unsigned values VTK file.
    int SaveToRawVTK16UFile(const char *file_name);
    
    /// Load an raw 16 bit VTK unsigned values file.
    int LoadRawVTK16UFile(const char *file_name);
};



template<class T> class bdImageT;  // pre-declare the template class itself

//template bdImageT<unsigned char>;
//template bdImageT<unsigned short>;

typedef bdImageT<unsigned char> bdImage8U;
typedef bdImageT<unsigned short> bdImage16U;


typedef bdImage16U bdImage;






#endif
