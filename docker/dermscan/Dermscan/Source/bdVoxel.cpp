/**************************************************************************
 Voxel class 

 Author: Danilo Babin
 File name: "bdVoxel.cpp"
***************************************************************************/



#include "bdVoxel.h"






//class bdDiscreteCoordinates3D
//{
//protected:
//	int *m_coordinates;// Coordinates: C, R, S.
//
//public:
//
//	// Constructor / Destructor
//	bdDiscreteCoordinates3D();
//	bdDiscreteCoordinates3D(int s_coordinate, int r_coordinate, int c_coordinate);
//	~bdDiscreteCoordinates3D();
//
//	// Set Coordinates of the voxel
//	bdDiscreteCoordinates3D& operator ()(int r_coordinate, int c_coordinate);
//	bdDiscreteCoordinates3D& operator ()(int s_coordinate, int r_coordinate, int c_coordinate);
//	void SetCoordinates3D(int *xyz);
//
//	// Get coordinates
//	int* GetCoordinates();
//	int& S();
//	int& R();
//	int& C();
//	
//	// Make an exact copy of the input voxel. operator= and CopyFrom do exactly the same.
//	bdDiscreteCoordinates3D& operator =(const bdDiscreteCoordinates3D &p);//LEGACY
//	void CopyFrom(const bdDiscreteCoordinates3D &p);
//	
//	// Check if the voxels are equal based on different criteria.
//	int HasEqual_2D_CoordinatesAs(const bdDiscreteCoordinates3D &p);// Checks only 2D coordinates of voxels.
//	int HasEqualCoordinatesAs(const bdDiscreteCoordinates3D &p);// Checks only 3D coordinates of voxels.
//	int IsEqualAs(const bdDiscreteCoordinates3D &p);// Checks coordinates and values stored in voxels.
//	
//	// Calculate distances between voxels.
//	double DistanceEuclidean(bdDiscreteCoordinates3D &p);
//	unsigned int DistanceEuclideanSquared(bdDiscreteCoordinates3D &p);
//	unsigned int DistanceCityBlock(bdDiscreteCoordinates3D &p);// L1 norm distance: d = ds+dr+dc.
//	unsigned int DistanceMaximumNorm(bdDiscreteCoordinates3D &p);// d = max{|ds|,|dr|,|dc|}.
//
//	
//	// Print out the coordinates and value of the voxel.
//	friend std::ostream& operator << (ostream &o, bdDiscreteCoordinates3D &p);
//	
//	// Read in the voxel by reading the input string.
//	int ReadFromCharString(char *string_with_position);
//
//};





//-----------------------------------------------------------------------------------------------------------



//class bdVoxel4D : public bdVoxel3D
//{
//public:
//
//	// Constructor / Destructor
//	bdVoxel4D();
//	bdVoxel4D(int t_coordinate, int s_coordinate, int r_coordinate, int c_coordinate);
//	~bdVoxel4D();
//
//	// Set Coordinate values of the voxel
//	bdVoxel4D& operator ()(int t_coordinate, int s_coordinate, int r_coordinate, int c_coordinate);
//	void SetCoordinates4D(int *xyzt);
//
//
//	// Get coordinates
//	int& T();
//	
//	// Make an exact copy of the input voxel. operator= and CopyFrom do exactly the same.
//	bdVoxel4D& operator =(bdVoxel4D &voxel);//LEGACY
//	void CopyFrom(const bdVoxel4D &voxel);
//	
//	// Check if the voxels are equal based on different criteria.
//	int HasEqual4DCoordinatesAs(const bdVoxel4D &voxel);
//	int IsEqualAs(bdVoxel4D &voxel);
//		
//	// Print out the coordinates and value of the voxel.
//	friend std::ostream& operator << (ostream &o, bdVoxel4D &voxel);
//	
////	// Read in the voxel by reading the input string.
////	int ReadFromCharString(char *string_with_postion);
//
//};
//


//class bdVoxel3D
//{
//private:
//	//int m_CRS[3];//
//	int *m_coordinates;
//	unsigned int m_value;
//
//public:
//
//	bdVoxel3D(){ m_CRS[0] = 0; m_CRS[1] = 0; m_CRS[2] = 0; m_value = 0; };
//	bdVoxel3D(int s_coordinate, int r_coordinate, int c_coordinate);
//	bdVoxel3D& operator ()(int r_coordinate, int c_coordinate);
//	bdVoxel3D& operator ()(int s_coordinate, int r_coordinate, int c_coordinate);
//	bdVoxel3D& operator ()(int s_coordinate, int r_coordinate, int c_coordinate, int _value);
//	
//	//bdVoxel3D& operator =(bdVoxel3D &p);
//	//bdVoxel3D& operator =(const bdVoxel3D &p);
//	void CopyFrom(const bdVoxel3D &p);
//	//bdVoxel3D& operator =(int value_);
//	int operator ==(bdVoxel3D &p);
//	int HasSame_3D_CoordinatesAs(const bdVoxel3D &p);
//	double Distance_Euclidean(const bdVoxel3D &p2);
//	int Distance_EuclideanSquared(const bdVoxel3D &p2);
//	int Distance_CityBlock(const bdVoxel3D &p2);
//	friend std::ostream& operator << (ostream &o, const bdVoxel3D &p);
//	int ReadFromCharString(char *string_with_postion);
//};







bdDiscreteCoordinates3D::bdDiscreteCoordinates3D()
{
//	m_coordinates = new int [3];
	m_coordinates[0] = 0;
	m_coordinates[1] = 0;
	m_coordinates[2] = 0;
}


bdDiscreteCoordinates3D::bdDiscreteCoordinates3D(int s_coordinate, int r_coordinate, int c_coordinate)
{
//	m_coordinates = new int [3];
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	m_coordinates[2] = s_coordinate;
}


bdDiscreteCoordinates3D::~bdDiscreteCoordinates3D()
{
//	delete [] m_coordinates;
}


bdDiscreteCoordinates3D& bdDiscreteCoordinates3D::operator()(int r_coordinate, int c_coordinate)
{
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	return *this;
}


bdDiscreteCoordinates3D& bdDiscreteCoordinates3D::operator()(int s_coordinate, int r_coordinate, int c_coordinate)
{
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	m_coordinates[2] = s_coordinate;
	return *this;
}


void bdDiscreteCoordinates3D::SetCoordinates3D(int *xyz)
{
	m_coordinates[0] = xyz[0];
	m_coordinates[1] = xyz[1];
	m_coordinates[2] = xyz[2];
}


int* bdDiscreteCoordinates3D::GetCoordinates()
{
	return m_coordinates;
}


int& bdDiscreteCoordinates3D::S()
{
	return m_coordinates[2];
}


int& bdDiscreteCoordinates3D::R()
{
	return m_coordinates[1];
}


int& bdDiscreteCoordinates3D::C()
{
	return m_coordinates[0];
}


bdDiscreteCoordinates3D& bdDiscreteCoordinates3D::operator =(const bdDiscreteCoordinates3D &p)
{
    if(&p==this) return *this;
	m_coordinates[0] = p.m_coordinates[0];
	m_coordinates[1] = p.m_coordinates[1];
	m_coordinates[2] = p.m_coordinates[2];
    return *this;
}


void bdDiscreteCoordinates3D::CopyFrom(const bdDiscreteCoordinates3D &p)
{
    if(&p==this) return;
	m_coordinates[0] = p.m_coordinates[0];
	m_coordinates[1] = p.m_coordinates[1];
	m_coordinates[2] = p.m_coordinates[2];
}


int bdDiscreteCoordinates3D::IsEqualAs(const bdDiscreteCoordinates3D &p)
{
	return this->HasEqualCoordinatesAs(p);
}


int bdDiscreteCoordinates3D::HasEqualCoordinatesAs(const bdDiscreteCoordinates3D &p)
{
	if(m_coordinates[0] != p.m_coordinates[0]) return 0;
	if(m_coordinates[1] != p.m_coordinates[1]) return 0;
	if(m_coordinates[2] != p.m_coordinates[2]) return 0;
	return 1;
}


int bdDiscreteCoordinates3D::HasEqual_2D_CoordinatesAs(const bdDiscreteCoordinates3D &p)
{
	if(m_coordinates[0] != p.m_coordinates[0]) return 0;
	if(m_coordinates[1] != p.m_coordinates[1]) return 0;
	return 1;
}


double bdDiscreteCoordinates3D::DistanceEuclidean(bdDiscreteCoordinates3D &p)
{
	return(sqrt((double)( (S()-p.S())*(S()-p.S()) + (R()-p.R())*(R()-p.R()) + (C()-p.C())*(C()-p.C()) )));
}


unsigned int bdDiscreteCoordinates3D::DistanceEuclideanSquared(bdDiscreteCoordinates3D &p)
{
	return( (S()-p.S())*(S()-p.S()) + (R()-p.R())*(R()-p.R()) + (C()-p.C())*(C()-p.C()) );
}


unsigned int bdDiscreteCoordinates3D::DistanceCityBlock(bdDiscreteCoordinates3D &p)
{
	int ds = S()-p.S(); if(ds<0) ds = -ds;
	int dr = R()-p.R(); if(dr<0) dr = -dr;
	int dc = C()-p.C(); if(dc<0) dc = -dc;
	return(ds+dr+dc);
}


unsigned int bdDiscreteCoordinates3D::DistanceMaximumNorm(bdDiscreteCoordinates3D &p)
{
	int ds = S()-p.S(); if(ds<0) ds = -ds;
	int dr = R()-p.R(); if(dr<0) dr = -dr;
	int dc = C()-p.C(); if(dc<0) dc = -dc;
	if(ds>=dr && ds>=dc) return ds;
	else
	{
		if(dr>=ds && dr>=dc) return dr;
		else return dc;
	}
}


ostream& operator<< (ostream &o, bdDiscreteCoordinates3D &p)
{ 
	o<<"["<<p.S()<<","<<p.R()<<","<<p.C()<<"]";
    return o;
}


int bdDiscreteCoordinates3D::ReadFromCharString(char *string_with_position)
{
	bdString s, s_coord;
	s.Append(string_with_position);

	//Read the coordinates to list
	s.ExtractStringBetweenCharacters('[', ']', s_coord);
	std::list<int> list_of_coordinates;
	s_coord.ExtractAllIntNumbersToList(list_of_coordinates);
	if(list_of_coordinates.size() != 3)
	{
		std::cout<<"bdPosition3D::ReadFromCharString(char *): Error reading coordinates!"<<endl;
		return 0;
	}

	//If no errors have been encountered, put the coordinates and values into the bdPosition
	int i = (int) (list_of_coordinates.size()-1);
	for(std::list<int>::iterator it = list_of_coordinates.begin(); it != list_of_coordinates.end(); it++)
	{
		m_coordinates[i] = (*it);
		i--;
	}

	return 1;
}



//-----------------------------------------------------------------------------------------------------------




bdVoxel3D::bdVoxel3D()
{
//	m_coordinates = new int [3];
	m_coordinates[0] = 0;
	m_coordinates[1] = 0;
	m_coordinates[2] = 0;
	m_value = 0;
}


bdVoxel3D::bdVoxel3D(int s_coordinate, int r_coordinate, int c_coordinate)
{
//	m_coordinates = new int [3];
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	m_coordinates[2] = s_coordinate;
	m_value = 0;
}


bdVoxel3D::~bdVoxel3D()
{
//	delete [] m_coordinates;
	m_value = 0;
}


bdVoxel3D& bdVoxel3D::operator()(int r_coordinate, int c_coordinate)
{
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	return *this;
}


bdVoxel3D& bdVoxel3D::operator()(int s_coordinate, int r_coordinate, int c_coordinate)
{
	m_coordinates[0] = c_coordinate;
	m_coordinates[1] = r_coordinate;
	m_coordinates[2] = s_coordinate;
	return *this;
}

//
//void bdVoxel3D::SetCoordinates3D(double *xyz);
//{
//	m_coordinates[0] = xyz[0];
//	m_coordinates[1] = xyz[1];
//	m_coordinates[2] = xyz[2];
//	return *this;
//}


unsigned int& bdVoxel3D::V()
{
	return m_value;
}


void bdVoxel3D::SetValue(unsigned int value)
{
	m_value = value;
}

//
//int* bdVoxel3D::GetCoordinates()
//{
//	return m_coordinates;
//}
//
//
//int bdVoxel3D::S()
//{
//	return m_coordinates[2];
//}
//
//
//int bdVoxel3D::R()
//{
//	return m_coordinates[1];
//}
//
//
//int bdVoxel3D::C()
//{
//	return m_coordinates[0];
//}


bdVoxel3D& bdVoxel3D::operator =(const bdVoxel3D &voxel)
{
    if(&voxel==this) return *this;
	m_coordinates[0] = voxel.m_coordinates[0];
	m_coordinates[1] = voxel.m_coordinates[1];
	m_coordinates[2] = voxel.m_coordinates[2];
	m_value = voxel.m_value;
    return *this;
}


bdVoxel3D& bdVoxel3D::operator =(bdDiscreteCoordinates3D &p)//bdVoxel3D& bdVoxel3D::operator =(const bdDiscreteCoordinates3D &p)
{
    //if(&p==this) return *this;
	m_coordinates[0] = p.C();
	m_coordinates[1] = p.R();
	m_coordinates[2] = p.S();
    return *this;
}


void bdVoxel3D::CopyFrom(const bdVoxel3D &voxel)
{
    if(&voxel==this) return;
	m_coordinates[0] = voxel.m_coordinates[0];
	m_coordinates[1] = voxel.m_coordinates[1];
	m_coordinates[2] = voxel.m_coordinates[2];
	m_value = voxel.m_value;
}


int bdVoxel3D::IsEqualAs(bdVoxel3D &voxel)
{
	if(m_coordinates[0] != voxel.m_coordinates[0]) return 0;
	if(m_coordinates[1] != voxel.m_coordinates[1]) return 0;
	if(m_coordinates[2] != voxel.m_coordinates[2]) return 0;
	if(m_value != voxel.m_value) return 0;
	return 1;
}
//
//
//int bdVoxel3D::HasEqual3DCoordinatesAs(const bdVoxel3D &voxel)
//{
//	if(m_coordinates[0] != voxel.m_coordinates[0]) return 0;
//	if(m_coordinates[1] != voxel.m_coordinates[1]) return 0;
//	if(m_coordinates[2] != voxel.m_coordinates[2]) return 0;
//	return 1;
//}
//
//
//int bdVoxel3D::HasEqual2DCoordinatesAs(const bdVoxel3D &voxel)
//{
//	if(m_coordinates[0] != voxel.m_coordinates[0]) return 0;
//	if(m_coordinates[1] != voxel.m_coordinates[1]) return 0;
//	return 1;
//}
//
//
//double bdVoxel3D::DistanceEuclidean(const bdVoxel3D &v)
//{
//	return(sqrt((double)((S()-v.S())*(S()-v.S())+(R()-v.R())*(R()-v.R())+(C()-v.C())*(C()-v.C()))));
//}
//
//
//unsigned int bdVoxel3D::DistanceEuclideanSquared(const bdVoxel3D &v)
//{
//	return( (S()-v.S())*(S()-v.S()) + (R()-v.R())*(R()-v.R()) + (C()-v.C())*(C()-v.C()) );
//}
//
//
//unsigned int bdVoxel3D::DistanceCityBlock(const bdVoxel3D &v)
//{
//	int ds = S()-v.S(); if(ds<0) ds = -ds;
//	int dr = R()-v.R(); if(dr<0) dr = -dr;
//	int dc = C()-v.C(); if(dc<0) dc = -dc;
//	return(ds+dr+dc);
//}


ostream& operator<< (ostream &o, bdVoxel3D &voxel)
{ 
	o<<"["<<voxel.S()<<","<<voxel.R()<<","<<voxel.C()<<"]{"<<voxel.m_value<<"}";
    return o;
}


int bdVoxel3D::ReadFromCharString(char *string_with_position)
{

    bdString bds, bds_coordinates, bds_value;
    bds.Append(string_with_position);
    bds.ExtractStringBetweenCharacters('[', ']', bds_coordinates);
    bdList<int> l;
    bds_coordinates.ExtractAllIntNumbersToList(l);
    if(l.GetNumberOfElements()!=3) return 0;
    this->S() = l.GetLeftEnd();
    this->R() = l[1];
    this->C() = l.GetRightEnd();
    bds.ExtractStringBetweenCharacters('{', '}', bds_value);
    bds_value.ExtractAllIntNumbersToList(l);
    if(l.GetNumberOfElements()!=1) return 0;
    this->V() = l.GetLeftEnd();
    
    
	return 1;
}



//-----------------------------------------------------------------------------------------------------------




//bdVoxel4D::bdVoxel4D()
//{
//	m_coordinates = new int [4];
//	m_coordinates[0] = 0;
//	m_coordinates[1] = 0;
//	m_coordinates[2] = 0;
//	m_coordinates[3] = 0;
//	m_value = 0;
//}
//
//
//bdVoxel4D::bdVoxel4D(int t_coordinate, int s_coordinate, int r_coordinate, int c_coordinate)
//{
//	m_coordinates = new int [4];
//	m_coordinates[0] = c_coordinate;
//	m_coordinates[1] = r_coordinate;
//	m_coordinates[2] = s_coordinate;
//	m_coordinates[3] = t_coordinate;
//	m_value = 0;
//}
//
//
//bdVoxel4D::~bdVoxel4D()
//{
//	delete [] m_coordinates;
//	m_value = 0;
//}
//
//
//bdVoxel4D& bdVoxel4D::operator()(int t_coordinate, int s_coordinate, int r_coordinate, int c_coordinate)
//{
//	m_coordinates[0] = c_coordinate;
//	m_coordinates[1] = r_coordinate;
//	m_coordinates[2] = s_coordinate;
//	m_coordinates[3] = t_coordinate;
//	return *this;
//}
//
//
//void bdVoxel4D::SetCoordinates4D(int *xyzt)
//{
//	m_coordinates[0] = xyzt[0];
//	m_coordinates[1] = xyzt[1];
//	m_coordinates[2] = xyzt[2];
//	m_coordinates[3] = xyzt[3];
//}
//
//
//int& bdVoxel4D::T()
//{
//	return m_coordinates[3];
//}
//
//
//bdVoxel4D& bdVoxel4D::operator =(bdVoxel4D &voxel)
//{
//    if(&voxel==this) return *this;
//	m_coordinates[0] = voxel.m_coordinates[0];
//	m_coordinates[1] = voxel.m_coordinates[1];
//	m_coordinates[2] = voxel.m_coordinates[2];
//	m_coordinates[3] = voxel.m_coordinates[3];
//	m_value = voxel.m_value;
//    return *this;
//}
//
//
//void bdVoxel4D::CopyFrom(const bdVoxel4D &voxel)
//{
//    if(&voxel==this) return;
//	m_coordinates[0] = voxel.m_coordinates[0];
//	m_coordinates[1] = voxel.m_coordinates[1];
//	m_coordinates[2] = voxel.m_coordinates[2];
//	m_coordinates[3] = voxel.m_coordinates[3];
//	m_value = voxel.m_value;
//}
//
//
//int bdVoxel4D::IsEqualAs(bdVoxel4D &voxel)
//{
//	if(m_coordinates[0] != voxel.m_coordinates[0]) return 0;
//	if(m_coordinates[1] != voxel.m_coordinates[1]) return 0;
//	if(m_coordinates[2] != voxel.m_coordinates[2]) return 0;
//	if(m_coordinates[3] != voxel.m_coordinates[3]) return 0;
//	if(m_value != voxel.m_value) return 0;
//	return 1;
//}
//
//
//int bdVoxel4D::HasEqual4DCoordinatesAs(const bdVoxel4D &voxel)
//{
//	if(m_coordinates[0] != voxel.m_coordinates[0]) return 0;
//	if(m_coordinates[1] != voxel.m_coordinates[1]) return 0;
//	if(m_coordinates[2] != voxel.m_coordinates[2]) return 0;
//	if(m_coordinates[3] != voxel.m_coordinates[3]) return 0;
//	return 1;
//}
//
//
//ostream& operator<< (ostream &o, bdVoxel4D &voxel)
//{ 
//	o<<"["<<voxel.T()<<","<<voxel.S()<<","<<voxel.R()<<","<<voxel.C()<<"]{"<<voxel.m_value<<"}";
//    return o;
//}
//

