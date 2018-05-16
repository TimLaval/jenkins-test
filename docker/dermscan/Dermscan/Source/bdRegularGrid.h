/**************************************************************************
 Template bdRegularGrid for grids up to 4D data

 Author: Danilo Babin
 File name: "bdRegularGrid.h"
***************************************************************************/




//To build as DLL, add:" /D "BD_REGULAR_GRID_EXPORTS" "
// in command line build options of the project.



#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_REGULAR_GRID_EXPORTS
		#define BD_REGULAR_GRID_API __declspec(dllexport) 
	#else
		#define BD_REGULAR_GRID_API __declspec(dllimport)
	#endif
#else // GCC
	#ifdef BD_REGULAR_GRID_EXPORTS
		#define BD_REGULAR_GRID_API __attribute__((visibility("default")))
	#else
		#define BD_REGULAR_GRID_API
	#endif
#endif




#ifndef BD_REGULAR_GRID_DEF
	#define BD_REGULAR_GRID_DEF



#include <typeinfo>
#include "bdArray.h"
#include "bdString.h"
#include "bdObject.h"



#pragma warning(disable: 4661)


template class BD_REGULAR_GRID_API bdArray<int>;
	



template<class T>
class BD_REGULAR_GRID_API bdRegularGridT : public bdDataObject
{
protected:

	/// Number of Columns, Rows, Slices and Time sequences.
	unsigned int m_dimensions_CRST[4];
	/// Based on the dimensions, the data can be of tyes: CRST, CRT, CST, CRS, RST. 
	
	/// Orientation in 3D is saved here. CxCyCz indicate the increasing direction of the columns, 
	/// while RxRyRz indicate the increasing direction of the rows. Sx,Sy,Sz are calculated from input R and C values.
	double m_orientation_CxCyCzRxRyRzSxSySz[9];
	
	/// Spacing: for S,R,C this is the spacing between slices/pixels in mm; for T this is the time interval in ms.
	double m_spacing_CRST[4]; 

	/// Origin: for S,R,C this is the origin of data in 3D space in mm; for T this is the start time in ms.
	double m_origin_CRST[4];

	/// We assume that any 4D data is read/written and a single continuous memory space. 
	/// 'm_pT' is an array of pointers to different 3D grids in the 4D grid. 'm_pS' is an array of RELATIVE pointers to different
	/// 2D grids in the 3D grid. 'm_pR' is an array of RELATIVE pointers to different grid rows in the 2D grid.
	/// Address of the data element with index (t,s,r,c) is: address = m_pT[t] + m_pS[s] + m_pR[r] + c.
	bdArray<int> m_pT, m_pS, m_pR;

	/// Pointer to data
	T* m_data;

	/// Indicates if the _data was internaly created with 'new' function (e.g. in the SetSize method). If so, the data will 
	/// have to be released with 'delete' function.
	int m_is_data_created_internally;

	///// Keeps track of the current index in the grid. This is important for functions which work only on parts of the grid, 
	///// e.g. only on slices. In that case we set m_current_index_ST[0]=s_desired to process the slice with index s_desired.
	//int m_current_index_ST[2];

	//2D Matrix containing tags of every slice. It's size is m_dimensions_CRST[3]*m_dimensions_CRST[2], in other words it is
	// number_of_time_series*number_of_slices.
	//bdMatrix3D<bdImage2DTag> m_tags_ST;

	/// Take the existing data and attach it to this object. Sets default visualization properties: spacing:(1,1,1,1,1), 
	/// origin:(0,0,0,0,0). Be sure that the data stored at data_address was allocated with enough size using 'new' function.
	int GetDataFromAddress(void *data_address, unsigned int n_of_columns, unsigned int n_of_rows, unsigned int n_of_slices, unsigned int n_of_time_series);

	/// Detach the array structures from the data and returns the pointer to data
	T* DetachAndGetPointerToData();


public:
	
	/// Constructor.
	bdRegularGridT();

	/// Destructor.
	~bdRegularGridT();

	/// Copy entire grid from input grid.
	virtual void CopyFrom(bdRegularGridT<T> &m);

	/// Calls the CopyFrom() method.
	bdRegularGridT<T>& operator =(bdRegularGridT<T> &m);

	/// Print type and header info.
	ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	ostream& PrintInfo(ostream &o);

	/// Check if the grid is empty.
	int IsEmpty();

	/// Check if the grid is equal in size as the input grid.
	int IsEqualSizeAs(bdRegularGridT<T> &grid);
	/// Check if the grid is equal in size as the input grid. Works also for 2D+time matrices.
	int IsEqualSizeAs_3D(bdRegularGridT<T> &grid);
	int IsEqualSizeAs_2D(bdRegularGridT<T> &grid);
//	int IsEqualSpacingAs(bdRegularGridT<T> &grid, double offset = 0.0001);

	/// Get pointer to data elements.
	T* GetData();

	/// Get grid dimensions.
	unsigned int* GetDimensions();
	unsigned int GetNumberOfTimeSeries();
	unsigned int GetNumberOfSlices();
	unsigned int GetNumberOfRows();
	unsigned int GetNumberOfColumns();
	unsigned int GetNumberOfDataElements();
	///// If the type of the grid is TRC (2D in time) returns number of TIME SERIES.
	///// For any other type returns number of slices.
	////int GetNumberOfSlicesOrTimeSeries();
	
	/// Get/Set Origin.
	double* GetOrigin();
	double GetOrigin_T();
	double GetOrigin_S();
	double GetOrigin_R();
	double GetOrigin_C();
	virtual void SetOrigin(double t, double s, double r, double c);

	/// Get/Set Spacing.
	double* GetSpacing();
	double GetSpacing_T();
	double GetSpacing_S();
	double GetSpacing_R();
	double GetSpacing_C();
	virtual void SetSpacing(double t, double s, double r, double c);

	/// Get/Set Orientation.
	double* GetOrientation();
	double GetOrientation_Cx();
	double GetOrientation_Cy();
	double GetOrientation_Cz();
	double GetOrientation_Rx();
	double GetOrientation_Ry();
	double GetOrientation_Rz();
    double GetOrientation_Sx();
    double GetOrientation_Sy();
    double GetOrientation_Sz();
	virtual void SetOrientationCxCyCzRxRyRz(double cx, double cy, double cz, double rx, double ry, double rz);

	/// Indexing for a 2D grid (assumes m_current_index values for indexes t,s).
	T& RC(unsigned int r, unsigned int c);

	/// Indexing for a 3D space grid (assumes m_current_index[1] value for index t).
	T& SRC(unsigned int s, unsigned int r, unsigned int c);
	
	/// Indexing for time series of 2D grids (assumes m_current_index[0] value for index s).
	T& TRC(unsigned int t, unsigned int r, unsigned int c);
	
	/// Indexing for time series of 3D grids
	T& TSRC(unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	///// Sets the CURRENT T index (modifies _current_index_ST[1]).
	//bdRegularGridT<T>& SetT(int t);
	//
	///// Sets the CURRENT S index (modifies _current_index_ST[0]).
	//bdRegularGridT<T>& SetS(int s);

	T& operator()(unsigned int r, unsigned int c);
	T& operator()(unsigned int s, unsigned int r, unsigned int c);
	T& operator()(unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	/// Get local indexes for input world coordinates. If coordinates fall outside the grid return 0, else 1.
	int GetIndexesForWorldCoordinates(double w_t, double w_s, double w_r, double w_c, int &out_t, int &out_s, int &out_r, int &out_c);
	int GetIndexesForWorldCoordinates(double w_s, double w_r, double w_c, int &out_s, int &out_r, int &out_c);

    /// Get world coordinates for given indexes.
    int GetWorldCoordinatesForIndexes(unsigned int t, unsigned int s, unsigned int r, unsigned int c, double &out_w_t, double &out_w_z, double &out_w_y, double &out_w_x);
    int GetWorldCoordinatesForIndexes(unsigned int s, unsigned int r, unsigned int c, double &out_w_z, double &out_w_y, double &out_w_x);

    /// Calculate "center of mass" for the 3D part of the image.
    int CenterOfMass(double &out_world_x, double &out_world_y, double &out_world_z);
    
    /// Bounds of the regular grid in world coordinates.
    int Bounds_3D_WorldCoordinates(double &x_min, double &x_max, double &y_min, double &y_max, double &z_min, double &z_max);
    
	/// Set the size of the grid with properties to DEFAULT values: spacing:(1,1,1,1), origin:(0,0,0,0)...
	/// If the grid is already the requested size, does not do anything.
	virtual int SetSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c);
	virtual int SetSizeOf3DGrid(unsigned int s, unsigned int r, unsigned int c);
	virtual int SetSizeOf2DTimeSeries(unsigned int t, unsigned int r, unsigned int c);
	virtual int SetSizeOf2DGrid(unsigned int r, unsigned int c);
	
	///Set only size of the grid, properties remain as before
	void SetSizeOnlyAs(bdRegularGridT<T> &grid);
	
	///Set size and properties of the grid
	void SetSizeAndPropertiesAs(bdRegularGridT<T> &grid);
	void SetVisualizationPropertiesToMatchInput(bdRegularGridT<T> &grid);

	/// Releases the data and resets the grid properties to the default values.
	/// Releases data (or not) based on '_is_data_created_internally' indicator.
	virtual void Reset();

	/// For an existing grid, this method will re-arrange the pointers to get a desired number of time_series 't', slices 's',
	/// rows 'r' and columns 'c'. ATTENTION! The number of total data elements MUST REMAIN THE SAME: t*s*r*c must be
	/// equal to existing total number of elements, if not the conversion will abort and the method will return FAIL (0).
	/// Use this method to convert the externally loaded 3D data to desired type of grid.
	/// It is assumed that original grid is organized as: Time0{slice0,...,sliceN},...,TimeM{slice0,...,sliceN}.
	virtual int ConvertToSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c);

	/// Fill in the whole grid with input value.
	void FillInWith(T value);

	/// Fill in the slice with indexes (t,s) with input value.
	void FillInSliceWith(unsigned int t, unsigned int s, T value);


	/// Print properties from grid.
	bdString PrintSpacing();
	bdString PrintOrigin();
	bdString PrintDimensions();
	bdString PrintOrientation();

};




//template class BD_REGULAR_GRID_API bdRegularGridT<unsigned char>;
//template class BD_REGULAR_GRID_API bdRegularGridT<unsigned short>;





#endif
