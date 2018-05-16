/**************************************************************************
 Class of 2D region 

 Author: Danilo Babin
 File name: "bdRegion2D.h"
***************************************************************************/



//To build as DLL, add:" /D "BD_REGION_2D_EXPORTS" "
// in command line build options of the project.



#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_REGION_2D_EXPORTS
		#define BD_REGION_2D_API __declspec(dllexport) 
	#else
		#define BD_REGION_2D_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_REGION_2D_EXPORTS
		#define BD_REGION_2D_API __attribute__((visibility("default")))
	#else
		#define BD_REGION_2D_API
	#endif
#endif





#ifndef BD_REGION_2D_DEF
	#define BD_REGION_2D_DEF


#include "bdVoxel.h"
#include "bdGeometry.h"
#include "bdList.h"
#include "bdImage.h"





class BD_REGION_2D_API bdRegionColumnNode
{
public:
	unsigned int m_column_index;
	int m_element;
};


//template class BD_REGION_2D_API bdList< bdRegionColumnNode >;


class BD_REGION_2D_API bdRegionRowNode
{
public:
	bdList< bdRegionColumnNode > m_column_list;
	unsigned int m_row_index;
};



//-------------------------------------------------------------------------------------------------------------------------------------------------------------


//template class BD_REGION_2D_API bdListNode<bdRegionRowNode>;
//template class BD_REGION_2D_API bdListNode<bdRegionColumnNode>;
//template class BD_REGION_2D_API bdList<bdRegionRowNode>;



class BD_REGION_2D_API bdRegion2DNode
{
public:
	bdListNode<bdRegionRowNode> *m_node_of_row_node; 
	bdListNode<bdRegionColumnNode> *m_node_of_column_node; 

	/// Constructor.
	bdRegion2DNode();

	/// Descructor.
	~bdRegion2DNode();

	/// Check if the pointers are valid.
	int IsNonZero();

	/// Get element to which this object points.
	int& GetElement();

	/// Get row index for the point.
	int GetIndexRow();

	/// Get column index for the point.
	int GetIndexColumn();

	/// Move to the next point in the region.
	void MoveToNext();

};



//-------------------------------------------------------------------------------------------------------------------------------------------------------------


class bdRegion2D;


class BD_REGION_2D_API bdRegion2DIterator : public bdRegion2DNode
{
public:

	/// Constructor.
	bdRegion2DIterator();

	/// Descructor.
	~bdRegion2DIterator();

	/// Set the iterator to the start of the region
	void SetBegin(bdRegion2D *region);
	void SetBegin(bdRegion2D *region, unsigned int start_r_index);
//	void SetBegin(bdRegion2D<T> &region);

	/// Check if the iterator is valid with respect to the set region.
	int IsValid();

};



//-------------------------------------------------------------------------------------------------------------------------------------------------------------


//template class BD_REGION_2D_API bdRegularGridT<unsigned short>;//><bdRegionRowNode>;



class BD_REGION_2D_API bdRegion2D
{
protected:

	/// Size of the original image, from which the region is created
	unsigned int m_original_image_size_CR[2];

	/// Number of elements contained in the region.
	unsigned int m_number_of_elements;
	
	/// Region value used like a label.
	int m_region_value;

	/// Data.
	bdList<bdRegionRowNode> m_row_list;


public:
	
	/// Constructor.
	bdRegion2D();

	/// Destructor.
	~bdRegion2D();
	
	/// Reset the object to its initial state (just after construction).
	void Reset();

	/// Removes all positions from the region, but does not change anything else (original dimensions of image, region value remain).
	void RemoveAllPoints();

	/// Check if the region is empty.
	int IsEmpty();

	/// Return pointer to List_Node of Row_Node for given input row index.
	bdListNode<bdRegionRowNode>* GetNodeOfRowNode(unsigned int r_index);

	/// Return Region2DNode for given indexes.
	bdRegion2DNode GetNode(unsigned int r_index, unsigned int c_index);

	/// Return Region2DNode for the fist point in the region.
	bdRegion2DNode GetNodeForFirstPoint();

	/// Get the (label) value of the region (single one for the whole region).
	int GetRegionValue(){return m_region_value;};

	/// Get the row_list of the region.
	bdList<bdRegionRowNode>* GetRowList(){return &m_row_list;};

	/// Getnumber of elements contained (entered) in the whole region.
	unsigned int GetNumberOfElements(){return m_number_of_elements;};

	/// Get the size of the original image.
	unsigned int GetNumberOfRowsOfOriginalImage(){return m_original_image_size_CR[1];};
	unsigned int GetNumberOfColumnsOfOriginalImage(){return m_original_image_size_CR[0];};
	
	/// Remove a point from the region.
	int RemovePoint(unsigned int r_index, unsigned int c_index);
	int RemovePoint(bdRegion2DNode &node);
	
	/// Add a point to the region with or without setting element values.
	int AddPoint(unsigned int r, unsigned int c);
	int AddPoint(unsigned int r, unsigned int c, int element);
	int AddPoint(unsigned int r, unsigned int c, int *p_element);

	/// Add the whole region to this region.
	void AddRegion(bdRegion2D &region);
	
	/// Check if the point with given indexes is in the region.
	int IsPointInRegion(unsigned int r_index, unsigned int c_index);

	/// Check if this region overlaps with the input region. Use to check if regions in neighboring slices of 3D image "touch".
	int IsOverlappingWithRegion(bdRegion2D &region);

	/// Returns pointer to value stored at given indexes of the region. If not in the region, returns NULL.
	int* GetValue(unsigned int r, unsigned int c);

	/// Make an exact copy.
	bdRegion2D& operator =(bdRegion2D &region);

	/// Set the region label value.
	void SetRegionValue(int value){m_region_value = value;};

	/// Set the size of the original image. Needed later on for reconstruction of region to image of correct size.
	void SetSizeOfOriginalImage(unsigned int rows, unsigned int columns);

	/// Get number of points of the region in the given row index.
	unsigned int GetNumberOfPointsInRow(unsigned int r_index);

	/// Get number of elemenpointsts of the region in the given column index. NOTICE: less efficient than GetNumberOfPointsInRow().
	unsigned int GetNumberOfPointsInColumn(unsigned int c_index);

	////Puts the region into the framed image. Framed image is the one with dimensions R_frame+2,C_frame+2,
	//// where R_frame is the biggest region span in terms of rows (C_frame respectively) and is calculated by 
	//// FrameCoordinatesForRegion function. Framed image is used to easily determine the outer edge!
	//// bdRegion2D pixels will have value region_pixel_value, while the background will be 0.
	//void ConvertRegionToFramedImage(Image8U &output_framed_image, unsigned char region_pixel_value);

	/// For existing region insert pixel values from input image.
	int InsertValuesFromImage(bdImage16U &image, unsigned int s, unsigned int t=0);

	/// Find points with the highest value (found in the total region) and record them to output_ridges region.
	int RidgePointsRegion(bdRegion2D &output_ridges);

	/// Create region from the whole image slice. Put non-zero pixels in the region and records their values.
	int CreateRegionFromImage(bdImage16U &image, unsigned int s, unsigned int t=0);

	/// Create region from the whole image slice for given voxel value. Put value_for_growing valued pixels in the region.
	int CreateRegionFromImage(bdImage16U &image, unsigned short value_for_growing, unsigned int s, unsigned int t=0);

	/// Create region from the seed point by searching 8-neighborhoods. Put non-zero valued pixels in the region and records their values.
	int CreateRegion_8_FromSeedPoint(bdImage16U &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	/// Create region from the seed point by searching 8-neighborhoods. Put value_for_growing valued pixels in the region.
	int CreateRegion_8_FromSeedPoint(bdImage16U &image, unsigned short value_for_growing, unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	/// Create region from the seed point by searching 8-neighborhoods. Put pixels with values higher than 'value_to_grow_above' in the region.
	int CreateRegion_8_FromSeedPointAndValueHigherThan(bdImage16U &image, unsigned short value_to_grow_above, unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	/// Calculate the center of mass position (each point in the region is equal) of the region.
	int CenterOfMass(double &center_of_mass_R, double &center_of_mass_C);

	/// Calculate Dice set similarity measure: dice = 2*|region1 & region2| / |region1|+|region2|, where '&' is intersection of regions.
	double Dice(bdRegion2D &region);

	/// Calculate the geometric median position of the region - it is a position IN the region closest to the center of mass.
	int GeometricMedian(double &geometric_median_R, double &geometric_median_C);

	/// Paste region to an already existing image (image with preset size). The value to paste with is given by 'region_value'.
	int PasteRegionToExistingImage(bdImage &img, unsigned short region_value, unsigned int t=0, unsigned int s=0);

	/// Get list of overlapping regions in image 2D slice with this 2D region.
	int ListOfOverlappingRegionsInImageSlice(bdImage &img, bdList<bdRegion2D> &output_list, unsigned int t=0, unsigned int s=0);

	/// Get list of overlapping regions in image 2D slice with this 2D region. A region is a connected component of pixels with the SAME gray value!
	int ListOfOverlappingRegionsInImageSlice_SameIntensityForRegion(bdImage &img, bdList<bdRegion2D> &output_list, unsigned int t=0, unsigned int s=0);

	/// Calculate the overlap area with the input region. Area is calculated in terms of number of positions that overlap.
	unsigned int AreaOfOverlapWith(bdRegion2D &region);
	
	//void PasteRegionToImage(Image16U &image);
	////Pastes the region to the existing image and gives the pixels in the region desired ('region_pixel_value') value
	//void PasteRegionToImage_WithValue(Image8U &image, unsigned char region_pixel_value);

	////Approximates the inner radius size by shrinking (thinning) the whole region. If multiple points are candidates
	//// for the inner radius center point, the one that is the closest to the center of mass point is chosen. The radius
	//// value is equal to the number of erosions needed to erase the whole region.
	//void InnerRadiusCalculation_Approximate(int *inner_radius, Point &inner_radius_center_point);

	////Approximates the outer radius size by checking 8 directions in the coordinates of region frame. This gives 16 intersection
	//// points. The longest distance between these points is calculated and this is the outer radius value. The center point is 
	//// calculated as the point on the longest line (line connecting the points with longest distance) that is the closest 
	//// to the center of mass point.
	//void OuterRadiusCalculation_Approximate(int *outer_radius, Point &outer_radius_center_point);

	////Calculates the outer radius by making outer edge list of points and calculates the longest distance between any two
	//// points in that list. The center point is calculated as the point on the longest line (line connecting the points with
	//// longest distance) that is the closest to the center of mass point.
	//void OuterRadiusCalculation(int *outer_radius, Point &outer_radius_center_point);

	////Checks how circular the region is, 'center_of_mass' is calculated center of mass, 'average' is the average distance of all
	//// internal outer edge pixels from center_of_mass, 'l1' is the L1 norm of the distances compared to average distance, 'l2' is
	//// the L2 norm and min and max are minimum and maximum difference between any distance and average.
	//void CircularityCheck(double *l1, double *l2, double *min, double *max, double *average, Position &center_of_mass);

	////Checks how circular the region is with giving the priority to large regions. This is done by calculating center of mass point
	//// and growing the circle (inscribed circle) until the region edge is touched (if the center of mass point does not fall inside
	//// the region the growing radius value is 0). The returned value is equal to 'area of the inscribed circle' - 'difference between 
	//// the areas of inscribed circle and region'. In other words, it is: 2*'area of the inscribed circle' - 'region area'. Note that the
	//// return value can be negative (for very non-circular objects).
	//int CircularityCheckWithLargeObjectPriority();
	//	
	//void InnerAndOuterRadiusCalculation(int *inner_radius, int *outer_radius, Point &inner_radius_center_point);
	//void InnerAndOuterRadiusCalculation_fast(int *inner_radius, int *outer_radius, Point &inner_radius_center_point, int area);
	//void CenterOfMassPointCalculation(int *r_center_of_mass, int *c_center_of_mass);
	//void CenterOfMassPointCalculation(Point &center_of_mass);
	//void CenterOfMassPositionCalculation(Position &center_of_mass);
	//int OverlappingAreaWithImage(Image8U &image);
	//int OverlappingAreaWithRegion(bdRegion2D<T> &region);
	
	/// Get List or Region of edge points (internal or external).
	int Internal_4_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
	int Internal_4_EdgePointsRegion(bdRegion2D &region_of_egde_points);
	int Internal_8_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
	int Internal_8_EdgePointsRegion(bdRegion2D &region_of_egde_points);
	int External_4_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
	int External_4_EdgePointsRegion(bdRegion2D &region_of_egde_points);
	int External_8_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points);
	int External_8_EdgePointsRegion(bdRegion2D &region_of_egde_points);

	////Make the region convex by connecting all the edge pixels within the 'squared distance' value.
	//void MakeConvex(int distance_squared);

	/// Gets the outer internal_8 edge list. This is important to define the outer shape of the region and is done 
	/// by finding extenal edge and internal edge, and the result is the part of the internal edge that "leans on" 
	/// either the external edge or the image edges.
//	int ListOfOuterInternal_8_EdgePoints(bdList<bdDiscreteCoordinates3D> &list_of_outer_internal_edge_points);

	////Gets the outer internal_4 edge list. This is important to define the outer shape of the region and is done 
	//// by finding extenal edge and internal edge, and the result is the part of the internal edge that "leans on" 
	//// either the external edge or the image edges.
	//void ListOfOuterInternal_4_EdgePoints(bdList<Point> &list_of_outer_internal_edge_points);
	//	
	////For the region fills in all the "holes" between outer internal_8 edge. This means that there will be no inner edges
	//// in the resulting region, just the outer edge remains the same as before.
	//void FillInTheWholeRegionBetweenOuterInternal_8_Edge();

	////For the region fills in all the "holes" between outer internal_4 edge. This means that there will be no inner edges
	//// in the resulting region, just the outer edge remains the same as before.
	//void FillInTheWholeRegionBetweenOuterInternal_4_Edge();


	/// Make disjoint regions from this one. Returns the number of disjoint regions, (0 if none).
	int BreakIntoDisjointRegionsOf_4_Neighborhood(bdList<bdRegion2D> &list_of_regions);
	int BreakIntoDisjointRegionsOf_8_Neighborhood(bdList<bdRegion2D> &list_of_regions);


	//void CreateRegionFromTheListOfDisjointRegions(bdList<bdRegion2D<T>> &list_of_regions);

	//void IntersectionOfTwoRegions(bdRegion2D<T> &region1, bdRegion2D<T> &region2);

	//void Erode_4();
	//void Erode_4_NotAppliedToImageEdges();
	//void Erode_8();
	//void Dilate_4();
	//void Dilate_4_NotAppliedToImageEdges();
	//void Dilate_8();
	//void Opening_4(int radius);
	//void Opening_8(int radius);
	//void Closing_4(int radius);
	//void Closing_8(int radius);

	//void FrameCoordinatesForRegion(int *r_start_included, int *r_end_included, int *c_start_included, int *c_end_included);

	//friend std::ostream& operator << <>(ostream &o, const bdRegion2D<T> &r);
	//	
	//
	////int RadiusOfAbsoluteValue(int r, int c, double normalized_difference, unsigned char average_value);
	////int RadiusMedianOfValue(int r, int c, double normalized_difference, unsigned char average_value);
	////int RadiusMedianAbsoluteOfValue(int r, int c, double normalized_difference, unsigned char average_value);
	////void Normalize(int multiplication_value=1);
	////void ConvertToBinary(unsigned char treshold);

	////void GenerateLabelForValue(bdList<Point> &label_region, unsigned char label_value);
	////void InnerRadiusValueAndCenterForLabelValue(unsigned char label_value, int *max_radius_value, Point &center);
	
};



#endif



