/**************************************************************************
 Geometric shapes and objects for 2D and 3D images

 Author: Danilo Babin
 File name: "bdGeometry.h"
***************************************************************************/


//To build as DLL, add:" /D "BD_GEOMETRY_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_GEOMETRY_EXPORTS
		#define BD_GEOMETRY_API __declspec(dllexport) 
	#else
		#define BD_GEOMETRY_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_GEOMETRY_EXPORTS
		#define BD_GEOMETRY_API __attribute__((visibility("default")))
	#else
		#define BD_GEOMETRY_API
	#endif
#endif



#ifndef BD_GEOMETRY_DEF
	#define BD_GEOMETRY_DEF


#include "bdArray.h"
#include "bdList.h"
#include "bdVoxel.h"
#include "bdMatrix.h"//#include "bdRegularGrid.h"
#include <list>


template class BD_GEOMETRY_API bdArray<unsigned int>;


/// Container for bdGeometry class.

class BD_GEOMETRY_API bdGeometryPointContainer
{
public:

	/// Data.
	bdDiscreteCoordinates3D *m_points;
	
	/// Number of data elements contained in the m_points array.
	unsigned int m_number_of_elements;

	/// We assume that data is read/written in a single continuous memory space. 
	/// 'm_radius' is an array of pointers to different arrays of positions for a single radius value. 
	/// Address of the data element with Size_start(r) is: address = m_points + m_size[r].
	/// Address of the data element with Size_end(r) is: address = m_points + m_size[r+1] - 1.
	bdArray<unsigned int> m_size;
	
	/// Constructor.
	bdGeometryPointContainer();

	/// Destructor.
	~bdGeometryPointContainer();

	/// Check if the container of points is empty.
	int IsEmpty();

	/// Reset to empty default state.
	void Reset();

	/// Set 'm_points' array number of elements (total number of points stored).
	void i_SetPointsArray(unsigned int number_of_points);

	/// Get start index value of the object with size 'r'.
	unsigned int SizeStartIndex(unsigned int r);

	/// Get ending index value of the object with size 'r'.
	unsigned int SizeEndIndex(unsigned int r);

	/// Get point at the given input size index.
	bdDiscreteCoordinates3D& operator[](unsigned int r);
	
};




//-----------------------------------------------------------------------------------------------------------


/// Geometric shapes and objects for 2D and 3D images.

class BD_GEOMETRY_API bdGeometry
{
private:

	/// Size of the image that we calculate geometry for.
	int m_dimensions_CRS[3];

	static bdGeometryPointContainer *m_sphere;
	/// Circle in the 's' direction ('s' values are constant).
	static bdGeometryPointContainer *m_circle;
	/// Circle in the 'r' direction ('r' values are constant).
	static bdGeometryPointContainer *m_circle_R;
	/// Circle in the 'c' direction ('c' values are constant).
	static bdGeometryPointContainer *m_circle_C;
	static bdGeometryPointContainer *m_cube;
	static bdGeometryPointContainer *m_square;

	/// Used in For() and Get() methods to allow easy FOR loop.
	unsigned int m_FOR_loop_current_index;

	/// Used in For() and Get() methods to allow easy FOR loop.
	bdDiscreteCoordinates3D m_FOR_loop_central_coordinates;
    
    /// Line frame pattern for looping over line positions.
    bdArray<bdDiscreteCoordinates3D> m_line_frame_pattern;


public:
	///
	/// Example use:
	/// for(unsigned int i=g.SphereSquaredRadiusStartIndex(squared_radius1); i<=g.SphereSquaredRadiusEndIndex(squared_radius2); i++)
	/// {
	/// 	  if(g.GetCoordinates(p1, g.Sphere(i), p2)) cout<<p2<<" ";
	/// }
	/// ...or with an internal FOR loop index counter:
	/// for(g.ForCoordinates_Sphere(s,r,c,squared_radius1); g.Get_Sphere(squared_radius2,sn,rn,cn); ){... do something...}
	///


	/// Constructor.
	bdGeometry();

	/// Destructor.
	~bdGeometry();

	/// Check is 2D/3D dimensions are set.
	int IsDimensions3DSet();
	int IsDimensions2DSet();

	/// Set dimensions.
	void SetDimensions(unsigned int s_dimension, unsigned int r_dimension, unsigned int c_dimension);
	void SetDimensions(unsigned int r_dimension, unsigned int c_dimension);
	
	/// Checks if the given coordinates are valid in the set dimensions of the object (values set with SetDimensions() method).
	int IsValidPosition(int s, int r, int c);

	/// Use as start of FOR loop. Warning! A single bdGeometry object CAN NOT be used for double loop
	/// with ForCoordinates() and Get_() methods!
	void ForCoordinates(int s, int r, int c, unsigned int start_FOR_loop_index);
	void ForCoordinates(int r, int c, unsigned int start_FOR_loop_index);

	void ForCoordinates_4_Neighbors(int r, int c);
	void ForCoordinates_8_Neighbors(int r, int c);
	void ForCoordinates_9_Neighbors(int r, int c);
	void ForCoordinates_26_Neighbors(int s, int r, int c);
	void ForCoordinates_27_Neighbors(int s, int r, int c);
	void ForCoordinates_Circle(int r, int c, unsigned int start_squared_radius);
	void ForCoordinates_Sphere(int s, int r, int c, unsigned int start_squared_radius);
    void ForCoordinates_Line(int s, int r, int c, int direction_s, int direction_r, int direction_c);

	/// Get the value of the current index for the FOR loop 'm_FOR_loop_current_index'.
	unsigned int GetCurrentFORLoopIndex();

	/// Use as the second argument of FOR loop. Warning! A single bdGeometry object CANNOT be used for double loop
	/// with ForCoordinates() and Get_() methods!
	int Get_4_Neighbors(int &r, int &c);
	int Get_4_Neighbors(bdDiscreteCoordinates3D &output_point);
	int Get_8_Neighbors(int &r, int &c);
	int Get_8_Neighbors(bdDiscreteCoordinates3D &output_point);
	int Get_9_Neighbors(int &r, int &c);
	int Get_9_Neighbors(bdDiscreteCoordinates3D &output_point);
	int Get_26_Neighbors(int &s, int &r, int &c);
	int Get_26_Neighbors(bdDiscreteCoordinates3D &output_point);
	int Get_27_Neighbors(int &s, int &r, int &c);
	int Get_27_Neighbors(bdDiscreteCoordinates3D &output_point);
	int Get_Circle(unsigned int end_square_radius_included, int &r, int &c);
	int Get_Circle(unsigned int end_square_radius_included, bdDiscreteCoordinates3D &output_point);
	int Get_Sphere(unsigned int end_square_radius_included, int &s, int &r, int &c);
	int Get_Sphere(unsigned int end_square_radius_included, bdDiscreteCoordinates3D &output_point);
    int Get_Line(unsigned int end_squared_distance_included, int &s, int &r, int &c);
    int Get_Line(unsigned int end_squared_distance_included, bdDiscreteCoordinates3D &output_point);

    
	/// Get number of elements in Sphere array - equals maximum squared radius that can be used.
	unsigned int GetNumberOfSphereElements();

    /// Get number of elements in Circle array - equals maximum squared radius that can be used.
    unsigned int GetNumberOfCircleElements();

	/// For a fixed 'start_point' and relative position 'relative_point', calculates the output point.
	/// If the output point is within the defined geometry dimensions retuns success 1, else fail 0.
	int GetCoordinates(bdDiscreteCoordinates3D &start_point, bdDiscreteCoordinates3D &relative_point, bdDiscreteCoordinates3D &output_point);
	int GetCoordinates(int s_start, int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, bdDiscreteCoordinates3D &output_point);
	int GetCoordinates(int s_start, int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, int &output_s, int &output_r, int &output_c);
	int GetCoordinates(int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, int &output_r, int &output_c);


	/// For a point with position (r,c) returns how many 9-neighboring pixels fall inside the image. 
	int NumberOfValid_9_Positions(int r, int c);
	
	/// For a point with position (s,r,c) returns how many 27-neighboring voxels fall inside the image. 
	int NumberOfValid_27_Positions(int s, int r, int c);

	/// "Make" functions calculate the positions and store them into an array
	int Make3DCross(bdDiscreteCoordinates3D &starting_point, bdArray<bdDiscreteCoordinates3D> &cross_points_array, int cross_size = 1);
	int Make3DLineSegment(bdDiscreteCoordinates3D &start_point, bdDiscreteCoordinates3D &end_point, bdArray<bdDiscreteCoordinates3D> &line_segment_points_array);

	/// Make a array of line positions to the input_direction position such that the line is 1 pixel wide.
    /// WARNING: The input direction coordinates must all be non-negative (i.e. they must fall in the first quadrant).
	int MakeLineFramePattern(bdDiscreteCoordinates3D &input_direction, bdArray<bdDiscreteCoordinates3D> &output_line_frame_pattern);

	/// For directions in the input array finds the exact same directions or the ones that are the exact opposite 
	/// (e.g. [1,2,-3] and [-1,-2,3]) and deletes these repeated directions from output array.
	int FilterMultipleDirections(bdArray<bdDiscreteCoordinates3D> &input_array_of_directions, bdArray<bdDiscreteCoordinates3D> &output_array_of_directions);

	/// Get points in the SPHERE
	unsigned int SphereSquaredRadiusStartIndex(unsigned int squared_radius);
	unsigned int SphereSquaredRadiusEndIndex(unsigned int squared_radius);
	bdDiscreteCoordinates3D& Sphere(unsigned int index_of_point);

	/// Get points in the CIRCLE ('s' is constant)
	unsigned int CircleSquaredRadiusStartIndex(unsigned int squared_radius);
	unsigned int CircleSquaredRadiusEndIndex(unsigned int squared_radius);
	bdDiscreteCoordinates3D& Circle(unsigned int index_of_point);

	/// Get points in the CIRCLE of R direction ('r' is constant)
	unsigned int Circle_R_SquaredRadiusStartIndex(unsigned int squared_radius);
	unsigned int Circle_R_SquaredRadiusEndIndex(unsigned int squared_radius);
	bdDiscreteCoordinates3D& Circle_R(unsigned int index_of_point);

	/// Get points in the CIRCLE of C direction ('c' is constant)
	unsigned int Circle_C_SquaredRadiusStartIndex(unsigned int squared_radius);
	unsigned int Circle_C_SquaredRadiusEndIndex(unsigned int squared_radius);
	bdDiscreteCoordinates3D& Circle_C(unsigned int index_of_point);

	/// Get points in the CUBE
	unsigned int CubeSizeStartIndex(unsigned int size);
	unsigned int CubeSizeEndIndex(unsigned int size);
	bdDiscreteCoordinates3D& Cube(unsigned int index_of_point);

	/// Get points in the SQUARE
	unsigned int SquareSizeStartIndex(unsigned int size);
	unsigned int SquareSizeEndIndex(unsigned int size);
	bdDiscreteCoordinates3D& Square(unsigned int index_of_point);

	/// Get points in 4-Neighborhood
	unsigned int Neighborhood4StartIndex();
	unsigned int Neighborhood4EndIndex();
	bdDiscreteCoordinates3D& Neighborhood4(unsigned int index_of_point);

	/// Get points in 8-Neighborhood
	unsigned int Neighborhood8StartIndex();
	unsigned int Neighborhood8EndIndex();
	bdDiscreteCoordinates3D& Neighborhood8(unsigned int index_of_point);

	/// Get points in 9-Neighborhood
	unsigned int Neighborhood9StartIndex();
	unsigned int Neighborhood9EndIndex();
	bdDiscreteCoordinates3D& Neighborhood9(unsigned int index_of_point);

	/// Get points in 27-Neighborhood
	unsigned int Neighborhood27StartIndex();
	unsigned int Neighborhood27EndIndex();
	bdDiscreteCoordinates3D& Neighborhood27(unsigned int index_of_point);

	/// Get points in 26-Neighborhood
	unsigned int Neighborhood26StartIndex();
	unsigned int Neighborhood26EndIndex();
	bdDiscreteCoordinates3D& Neighborhood26(unsigned int index_of_point);


	/// Generate arrays of points for various geometric objects/neighboroods.
	void GenerateSphereArray(unsigned int maximum_radius = 30);
	void GenerateCircleArray(unsigned int maximum_radius = 50);
	void GenerateCircle_R_Array(unsigned int maximum_radius = 50);
	void GenerateCircle_C_Array(unsigned int maximum_radius = 50);
	void GenerateCubeArray(unsigned int maximum_size = 30);
	void GenerateSquareArray(unsigned int maximum_radius = 50);

	/// Save arrays of points for various geometric objects/neighboroods.
	int Save(bdGeometryPointContainer &point_container, const char *file_name_root, const char *extension);
	int SaveSphereArray(const char *file_name_root = "bd_geometry_sphere", const char *extension = "bdgm");
	int SaveCircleArray(const char *file_name_root = "bd_geometry_circle", const char *extension = "bdgm");
	int SaveCircle_R_Array(const char *file_name_root = "bd_geometry_circle_R", const char *extension = "bdgm");
	int SaveCircle_C_Array(const char *file_name_root = "bd_geometry_circle_C", const char *extension = "bdgm");
	int SaveCubeArray(const char *file_name_root = "bd_geometry_cube", const char *extension = "bdgm");
	int SaveSquareArray(const char *file_name_root = "bd_geometry_square", const char *extension = "bdgm");
	
	/// Load arrays of points for various geometric objects/neighboroods.
	int Load(bdGeometryPointContainer &point_container, const char *file_name_root, const char *extension);
	int LoadSphereArray(const char *file_name_root = "bd_geometry_sphere", const char *extension = "bdgm");
	int LoadCircleArray(const char *file_name_root = "bd_geometry_circle", const char *extension = "bdgm");
	int LoadCircle_R_Array(const char *file_name_root = "bd_geometry_circle_R", const char *extension = "bdgm");
	int LoadCircle_C_Array(const char *file_name_root = "bd_geometry_circle_C", const char *extension = "bdgm");
	int LoadCubeArray(const char *file_name_root = "bd_geometry_cube", const char *extension = "bdgm");
	int LoadSquareArray(const char *file_name_root = "bd_geometry_square", const char *extension = "bdgm");

};


#endif