/**************************************************************************
 Class of 3D region 

 Author: Danilo Babin
 File name: "bdRegion3D.h"
***************************************************************************/



//To build as DLL, add:" /D "BD_REGION_3D_EXPORTS" "
// in command line build options of the project.



#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_REGION_3D_EXPORTS
		#define BD_REGION_3D_API __declspec(dllexport) 
	#else
		#define BD_REGION_3D_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_REGION_3D_EXPORTS
		#define BD_REGION_3D_API __attribute__((visibility("default")))
	#else
		#define BD_REGION_3D_API
	#endif
#endif





#ifndef BD_REGION_3D_DEF
	#define BD_REGION_3D_DEF


#include "bdRegion2D.h"




class BD_REGION_3D_API bdRegionSliceNode
{
public:
	unsigned int m_slice_index;
	bdRegion2D m_slice_region_2D;
};



//-------------------------------------------------------------------------------------------------------------------------------------------------------------



class BD_REGION_3D_API bdRegion3DNode
{
public:
	bdListNode<bdRegionSliceNode> *m_node_of_slice_node; 
	bdRegion2DNode m_region_2D_node; 

	/// Constructor.
	bdRegion3DNode();

	/// Descructor.
	~bdRegion3DNode();

	/// Check if the pointers are valid.
	int IsNonZero();

	/// Get element to which this object points.
	int& GetElement();

	/// Get slice index for the point.
	int GetIndexSlice();

	/// Get row index for the point.
	int GetIndexRow();

	/// Get column index for the point.
	int GetIndexColumn();
};




//-------------------------------------------------------------------------------------------------------------------------------------------------------------



class bdRegion3D;


class BD_REGION_3D_API bdRegion3DIterator : public bdRegion3DNode
{
public:

	/// Constructor.
	bdRegion3DIterator();

	/// Descructor.
	~bdRegion3DIterator();

	/// Set the iterator to the start of the region
	void SetBegin(bdRegion3D *region);
	//void SetBegin(bdRegion3D *region, int start_r_index);
//	void SetBegin(bdRegion3D<T> &region);

	/// Check if the iterator is valid with respect to the set region.
	int IsValid();

	/// Move to the next point in the region.
	void MoveToNext();

};




//-------------------------------------------------------------------------------------------------------------------------------------------------------------




//template class BD_REGION_3D_API bdListNode<bdRegionSliceNode>;
//template class BD_REGION_3D_API bdList<bdRegionSliceNode>;


class BD_REGION_3D_API bdRegion3D
{
private:

	/// Size of the original image, from which the region is created
	unsigned int m_original_image_size_CRS[3];

	/// Number of elements contained in the region.
	unsigned int m_number_of_elements;
	
	/// Region value used like a label.
	int m_region_value;

	/// Data.
	bdList<bdRegionSliceNode> m_slice_list;


public:

	/// Constructor.
	bdRegion3D();

	/// Destructor.
	~bdRegion3D();

	/// Reset the object to its initial state (just after construction).
	void Reset();

	/// Removes all positions from the region, but does not change anything else (original dimensions of image, region value remain).
	void RemoveAllPoints();

	/// Check if the region is empty.
	int IsEmpty();
    
    /// Get the bounds of the region.
    int Bounds(unsigned int &output_s_min, unsigned int &output_s_max, unsigned int &output_r_min, unsigned int &output_r_max, unsigned int &output_c_min, unsigned int &output_c_max);

    /// Get the 6 bounding points of the region (points min and max index values for all dimensions).
    int BoundingPoints(bdDiscreteCoordinates3D &output_s_min, bdDiscreteCoordinates3D &output_s_max, bdDiscreteCoordinates3D &output_r_min, bdDiscreteCoordinates3D &output_r_max, bdDiscreteCoordinates3D &output_c_min, bdDiscreteCoordinates3D &output_c_max);

	/// Return pointer to List_Node of Slice_Node for given input slice index.
	bdListNode<bdRegionSliceNode>* GetNodeOfSliceNode(unsigned int s_index);

	/// Return Region2DNode for given indexes.
	bdRegion3DNode GetNode(unsigned int s_index, unsigned int r_index, unsigned int c_index);

	/// Get the (label) value of the region (single one for the whole region).
	int GetRegionValue(){return m_region_value;};

	/// Get the row_list of the region.
	bdList<bdRegionSliceNode>* GetSliceList(){return &m_slice_list;};

	/// Getnumber of elements contained (entered) in the whole region.
	unsigned int GetNumberOfElements(){return m_number_of_elements;};

	/// Get the size of the original image.
	unsigned int GetNumberOfSlicesOfOriginalImage(){return m_original_image_size_CRS[2];};
	unsigned int GetNumberOfRowsOfOriginalImage(){return m_original_image_size_CRS[1];};
	unsigned int GetNumberOfColumnsOfOriginalImage(){return m_original_image_size_CRS[0];};

	/// Remove a point from the region.
	int RemovePoint(int s_index, int r_index, int c_index);
	int RemovePoint(bdRegion3DNode &node);

	/// Add a point to the region with or without setting element values.
	int AddPoint(unsigned int s, unsigned int r, unsigned int c);
	int AddPoint(unsigned int s, unsigned int r, unsigned int c, int &element);
	int AddPoint(unsigned int s, unsigned int r, unsigned int c, int *p_element);

	/// Add the whole 2d/3d region to this region. In case of 2d region, the slice index must be specified.
	void AddRegion(bdRegion2D &region, unsigned int s);
	void AddRegion(bdRegion3D &region);

	/// Check if the point with given indexes is in the region.
	int IsPointInRegion(unsigned int s_index, unsigned int r_index, unsigned int c_index);

	/// Check if this region overlaps with the input region.
	int IsOverlappingWithRegion(bdRegion3D &region);

	/// Returns pointer to value stored at given indexes of the region. If not in the region, returns NULL.
	int* GetValue(unsigned int s, unsigned int r, unsigned int c);

	/// Make an exact copy.
	bdRegion3D& operator =(bdRegion3D &region);

	/// Set the region label value.
	void SetRegionValue(int value){m_region_value = value;};

	/// Set the size of the original image. Needed later on for reconstruction of region to image of correct size.
	void SetSizeOfOriginalImage(unsigned int slices, unsigned int rows, unsigned int columns);

	//int NumberOfElementsInRow(int r_index);//tested
	//int NumberOfElementsInColumn(int column_index);//tested

	/// Create region from the whole 3D image. Put non-zero pixels in the region and records their values.
	int CreateRegionFromImage(bdImage16U &image, unsigned int t=0);

	/// Create region from the whole 3D image for given voxel value. Put value_for_growing valued pixels in the region.
	int CreateRegionFromImage(bdImage16U &image, unsigned short value_for_growing, unsigned int t=0);

	/// Create region from the seed point by searching 26-neighborhoods. Put non-zero valued pixels in the region and records their values.
	int CreateRegion_26_FromSeed(bdImage16U &image, unsigned int c, unsigned int r, unsigned int s, unsigned int t=0);

	/// Create region from the 27-neighborhood from the central point as seed. Put non-zero valued pixels in the region and records their values.
	int CreateRegionFrom_27_Neighborhood(bdImage16U &image, unsigned int c, unsigned int r, unsigned int s, unsigned int t=0);
	int CreateRegionFrom_27_Neighborhood(bdRegion3D &region, unsigned int c, unsigned int r, unsigned int s);
	
    /// Center of mass of the region.
    int CenterOfMass(unsigned int &output_s, unsigned int &output_r, unsigned int &output_c);
	
//	void PasteRegionToImage(bdImage8U &image);


	//void InnerAndOuterRadiusCalculation(int *inner_radius, int *outer_radius, Point &inner_radius_center_point);
	//void InnerAndOuterRadiusCalculation_fast(int *inner_radius, int *outer_radius, Point &inner_radius_center_point, int area);
	//void CenterOfMassPointCalculation(int *r_center_of_mass, int *c_center_of_mass);
	//void CenterOfMassPointCalculation(Point &center_of_mass);
	//int OverlappingAreaWithImage(Image8U &image);
	//int OverlappingAreaWithRegion(Region<T> &region);

	/// Get List or Region of edge points (internal or external).
//	int Internal_6_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
//	int Internal_6_EdgePointsRegion(bdRegion2D &region_of_egde_points);
	int Internal_26_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
	int Internal_26_EdgePointsRegion(bdRegion3D &region_of_egde_points);
//	int External_6_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
//	int External_6_EdgePointsRegion(bdRegion2D &region_of_egde_points);
//	int External_26_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
//	int External_26_EdgePointsRegion(bdRegion2D &region_of_egde_points);


	/// Make disjoint regions from this one. Returns the number of disjoint regions, (0 if none).
//	int BreakIntoDisjointRegionsOf_6_Neighborhood(bdList<bdRegion3D> &list_of_regions);
	int BreakIntoDisjointRegionsOf_26_Neighborhood(bdList<bdRegion3D> &list_of_regions);


	//void ListOfExternalEdgePoints(DoubleList<Point> &list_of_egde_points);
	//void CreateRegionFromTheListOfDisjointRegions(DoubleList<Region<T>> &list_of_regions);

	//void Erode_4();
	//void Dilate_4();
	//void Opening_4();
	//void Closing_4();

	//friend std::ostream& operator << <>(ostream &o, const Region<T> &r);
		
	
	//int RadiusOfAbsoluteValue(int r, int c, double normalized_difference, unsigned char average_value);
	//int RadiusMedianOfValue(int r, int c, double normalized_difference, unsigned char average_value);
	//int RadiusMedianAbsoluteOfValue(int r, int c, double normalized_difference, unsigned char average_value);
	//void Normalize(int multiplication_value=1);
	//void ConvertToBinary(unsigned char treshold);

	//void GenerateLabelForValue(DoubleList<Point> &label_region, unsigned char label_value);
	//void InnerRadiusValueAndCenterForLabelValue(unsigned char label_value, int *max_radius_value, Point &center);
	
};



#endif