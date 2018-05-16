/**************************************************************************
 Voxel class 

 Author: Danilo Babin
 File name: "bdVoxel.h"
***************************************************************************/


//To build as DLL, add:" /D "BD_VOXEL_EXPORTS" "
// in command line build options of the project.



#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_VOXEL_EXPORTS
		#define BD_VOXEL_API __declspec(dllexport) 
	#else
		#define BD_VOXEL_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_VOXEL_EXPORTS
		#define BD_VOXEL_API __attribute__((visibility("default")))
	#else
		#define BD_VOXEL_API
	#endif
#endif



#ifndef BD_VOXEL_DEF
	#define BD_VOXEL_DEF


#include <iostream>
#include <math.h>
#include "bdString.h"





class BD_VOXEL_API bdDiscreteCoordinates3D
{
protected:

	/// Coordinates: C, R, S.
	int m_coordinates[3];

public:

	/// Constructor / Destructor
	bdDiscreteCoordinates3D();
	bdDiscreteCoordinates3D(int s_coordinate, int r_coordinate, int c_coordinate);
	~bdDiscreteCoordinates3D();

	/// Set Coordinates of the voxel
	bdDiscreteCoordinates3D& operator ()(int r_coordinate, int c_coordinate);
	bdDiscreteCoordinates3D& operator ()(int s_coordinate, int r_coordinate, int c_coordinate);
	void SetCoordinates3D(int *xyz);

	/// Get coordinates
	int* GetCoordinates();
	int& S();
	int& R();
	int& C();
	
	/// Make an exact copy of the input voxel. operator= and CopyFrom do exactly the same.
	bdDiscreteCoordinates3D& operator =(const bdDiscreteCoordinates3D &p);//LEGACY
	void CopyFrom(const bdDiscreteCoordinates3D &p);
	
	/// Check if the voxels are equal based on different criteria.
	int HasEqual_2D_CoordinatesAs(const bdDiscreteCoordinates3D &p);// Checks only 2D coordinates of voxels.
	int HasEqualCoordinatesAs(const bdDiscreteCoordinates3D &p);// Checks only 3D coordinates of voxels.
	int IsEqualAs(const bdDiscreteCoordinates3D &p);// Checks coordinates and values stored in voxels.
	
	/// Calculate distances between voxels.
	double DistanceEuclidean(bdDiscreteCoordinates3D &p);
	unsigned int DistanceEuclideanSquared(bdDiscreteCoordinates3D &p);
	unsigned int DistanceCityBlock(bdDiscreteCoordinates3D &p);// L1 norm distance: d = ds+dr+dc.
	unsigned int DistanceMaximumNorm(bdDiscreteCoordinates3D &p);// d = max{|ds|,|dr|,|dc|}.

	/// Print out the coordinates and value of the voxel.
	BD_VOXEL_API friend std::ostream& operator << (ostream &o, bdDiscreteCoordinates3D &p);
	
	/// Read in the voxel by reading the input string.
	int ReadFromCharString(char *string_with_position);

};




//-----------------------------------------------------------------------------------------------------------




class BD_VOXEL_API bdVoxel3D : public bdDiscreteCoordinates3D
{
protected:
	unsigned int m_value;

public:

	/// Constructor / Destructor
	bdVoxel3D();
	bdVoxel3D(int s_coordinate, int r_coordinate, int c_coordinate);
	~bdVoxel3D();

	/// Set Coordinates of the voxel
	bdVoxel3D& operator ()(int r_coordinate, int c_coordinate);
	bdVoxel3D& operator ()(int s_coordinate, int r_coordinate, int c_coordinate);

	/// Set/Get the value of the voxel
	unsigned int& V();
	void SetValue(unsigned int value);
	
	/// Make an exact copy of the input voxel. operator= and CopyFrom do exactly the same.
	bdVoxel3D& operator =(const bdVoxel3D &voxel);///LEGACY
	bdVoxel3D& operator =(bdDiscreteCoordinates3D &p);//LEGACY
	void CopyFrom(const bdVoxel3D &voxel);
	
	/// Check if the voxels are equal based on different criteria.
	int IsEqualAs(bdVoxel3D &voxel);// Checks coordinates and values stored in voxels.
	
	/// Print out the coordinates and value of the voxel.
	BD_VOXEL_API friend std::ostream& operator << (ostream &o, bdVoxel3D &voxel);
	
	/// Read in the voxel by reading the input string.
	int ReadFromCharString(char *string_with_position);

};




typedef bdVoxel3D bdVoxel;


#endif