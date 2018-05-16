/**************************************************************************
 Geometric shapes and objects for 2D and 3D images

 Author: Danilo Babin
 File name: "bdGeometry.cpp"
***************************************************************************/



#include "bdGeometry.h"


//Initialize static variables
bdGeometryPointContainer* bdGeometry::m_sphere = 0;
bdGeometryPointContainer* bdGeometry::m_circle = 0;
bdGeometryPointContainer* bdGeometry::m_circle_R = 0;
bdGeometryPointContainer* bdGeometry::m_circle_C = 0;
bdGeometryPointContainer* bdGeometry::m_cube = 0;
bdGeometryPointContainer* bdGeometry::m_square = 0;




//-----------------------------------------------------------------------------------------------------------




bdGeometryPointContainer::bdGeometryPointContainer()
{
	m_number_of_elements = 0;
	m_points = NULL;
}


bdGeometryPointContainer::~bdGeometryPointContainer()
{
	delete [] m_points;
	m_points = NULL;
	m_number_of_elements = 0;
	m_size.Reset();
}


int bdGeometryPointContainer::IsEmpty()
{
	if(m_number_of_elements<=0) return 1;
	else return 0;
}


void bdGeometryPointContainer::Reset()
{
	delete [] m_points;
	m_points = NULL;
	m_number_of_elements = 0;
	m_size.Reset();
}


void bdGeometryPointContainer::i_SetPointsArray(unsigned int number_of_points)
{
	if(!number_of_points) return;
	delete [] m_points;
	m_points = NULL;
	m_number_of_elements = number_of_points;
	m_points = new bdDiscreteCoordinates3D[m_number_of_elements];
}


unsigned int bdGeometryPointContainer::SizeStartIndex(unsigned int r)
{
	if(r >= m_size.GetNumberOfElements()) { return m_number_of_elements-1; }
	else { return m_size[r]; }
}


unsigned int bdGeometryPointContainer::SizeEndIndex(unsigned int r)
{
	// We check here for 'r+1' because we need to get the END of the size 'r'.
	if(r+1 >= m_size.GetNumberOfElements()) { return m_number_of_elements-1; }
	else { return ( m_size[r+1]-1 ); }
}


bdDiscreteCoordinates3D& bdGeometryPointContainer::operator[](unsigned int r)
{
	return m_points[r];
}



//bdDiscreteCoordinates3D* bdGeometryPointContainer::Size_start(unsigned int r)
//{
//	if(r >= m_size.GetNumberOfElements()) { return &(m_points[m_number_of_elements-1]); }
//	else { return &(m_points[ (m_size[r]) ]); }
//}
//
//
//bdDiscreteCoordinates3D* bdGeometryPointContainer::Size_end(unsigned int r)
//{
//	// We check here for 'r+1' because we need to get the END of the size 'r'.
//	if(r+1 >= m_size.GetNumberOfElements()) { return &(m_points[m_number_of_elements-1]); }
//	else { return &(m_points[ ( m_size[r+1]-1 ) ]); }
//}


//bdDiscreteCoordinates3D* bdGeometryPointContainer::GetPoint(unsigned int r)
//{
//	//if(r >= m_number_of_elements) return NULL;
//	return &(m_points[r]);
//}




//-----------------------------------------------------------------------------------------------------------




bdGeometry::bdGeometry()
{
	m_dimensions_CRS[0] = m_dimensions_CRS[1] = m_dimensions_CRS[2] = 0;
	m_FOR_loop_current_index = 0;
	m_FOR_loop_central_coordinates(0,0,0);
    m_line_frame_pattern.Set(1);
	if(!m_sphere)
	{
		m_sphere = new bdGeometryPointContainer();
		if(!this->LoadSphereArray())// Try loading the array.
		{
			this->GenerateSphereArray();// If loading fails, generate the array...
			this->SaveSphereArray();// ...and save it!
		}
	}
	if(!m_cube)
	{
		m_cube = new bdGeometryPointContainer();
		if(!this->LoadCubeArray())// Try loading the array.
		{
			this->GenerateCubeArray();// If loading fails, generate the array...
			this->SaveCubeArray();// ...and save it!
		}
	}
	if(!m_circle)
	{
		m_circle = new bdGeometryPointContainer();
		if(!this->LoadCircleArray())// Try loading the array.
		{
			this->GenerateCircleArray();// If loading fails, generate the array...
			this->SaveCircleArray();// ...and save it!
		}
	}
	if(!m_circle_R)
	{
		m_circle_R = new bdGeometryPointContainer();
		if(!this->LoadCircle_R_Array())// Try loading the array.
		{
			this->GenerateCircle_R_Array();// If loading fails, generate the array...
			this->SaveCircle_R_Array();// ...and save it!
		}
	}
	if(!m_circle_C)
	{
		m_circle_C = new bdGeometryPointContainer();
		if(!this->LoadCircle_C_Array())// Try loading the array.
		{
			this->GenerateCircle_C_Array();// If loading fails, generate the array...
			this->SaveCircle_C_Array();// ...and save it!
		}
	}
	if(!m_square)
	{
		m_square = new bdGeometryPointContainer();
		if(!this->LoadSquareArray())// Try loading the array.
		{
			this->GenerateSquareArray();// If loading fails, generate the array...
			this->SaveSquareArray();// ...and save it!
		}
	}
}


bdGeometry::~bdGeometry()
{
}


int bdGeometry::IsDimensions2DSet()
{
	if(m_dimensions_CRS[0]==0 || m_dimensions_CRS[1]==0) return 0;
	else return 1;
}


int bdGeometry::IsDimensions3DSet()
{
	if(m_dimensions_CRS[0]==0 || m_dimensions_CRS[1]==0 || m_dimensions_CRS[2]==0) return 0;
	else return 1;
}


void bdGeometry::SetDimensions(unsigned int s_dimension, unsigned int r_dimension, unsigned int c_dimension)
{
	m_dimensions_CRS[0] = c_dimension;
	m_dimensions_CRS[1] = r_dimension;
	m_dimensions_CRS[2] = s_dimension;
}


void bdGeometry::SetDimensions(unsigned int r_dimension, unsigned int c_dimension)
{
	m_dimensions_CRS[0] = c_dimension;
	m_dimensions_CRS[1] = r_dimension;
//	m_dimensions_CRS[2] = 1;
}



//void bdGeometry::ListOfClockwisePositionsOfSquareRadius(DoubleList<bdDiscreteCoordinates3D> &list, int r, int c, int square_r0, int square_r1)
//{
//	assert(Dimensions2DAreSet());
//	assert(square_r1>square_r0 && square_r0>=0);
//	assert(r>=0 && c>=0);
//
//	if(!list.IsEmpty()) list.DeleteAll();
//
//	DoubleList<bdDiscreteCoordinates3D> quadrant1, quadrant2, quadrant3, quadrant4;
//	bdDiscreteCoordinates3D point;
//
//	for(ForPoint(r,c); GetPositionsOfRadiusSquare(square_r0,square_r1,&point.r,&point.c); )
//	{
//		//If the point is in the first quadrant
//		if(point.r<r && point.c>=c)
//		{
//			quadrant1.AddToRightEnd(point);
//		}
//		//If the point is in the forth quadrant
//		if(point.r>=r && point.c>c)
//		{
//			quadrant4.AddToRightEnd(point);
//		}
//		//If the point is in the third quadrant
//		if(point.r>r && point.c<=c)
//		{
//			quadrant3.AddToLeftEnd(point);
//		}
//		//If the point is in the second quadrant
//		if(point.r<=r && point.c<c)
//		{
//			quadrant2.AddToLeftEnd(point);
//		}
//	}
//
//	//Prvi i treci kvadrant ne treba sortirati (s tim sto je treci ispunjen suprotnim smerom)!!!
//
//	//Sorting the points in the forth quadrant
//	if(!(quadrant4.IsEmpty() || quadrant4.NumberOfElements()==1))
//	{
//		int change_made = 1;
//		while(change_made)
//		{
//			change_made = 0;
//			for(int i=0; i<quadrant4.NumberOfElements()-1; i++)
//			{
//				if(quadrant4[i].r==quadrant4[i+1].r && quadrant4[i].c<quadrant4[i+1].c)
//				{
//					bdDiscreteCoordinates3D temp;
//					temp = quadrant4[i];
//					quadrant4[i] = quadrant4[i+1];
//					quadrant4[i+1] = temp;
//					change_made = 1;
//				}
//			}	
//		}
//	}
//
//	//Sorting the points in the second quadrant
//	if(!(quadrant2.IsEmpty() || quadrant2.NumberOfElements()==1))
//	{
//		int change_made = 1;
//		while(change_made)
//		{
//			change_made = 0;
//			for(int i=0; i<quadrant2.NumberOfElements()-1; i++)
//			{
//				if(quadrant2[i].r==quadrant2[i+1].r && quadrant2[i].c>quadrant2[i+1].c)
//				{
//					bdDiscreteCoordinates3D temp;
//					temp = quadrant2[i];
//					quadrant2[i] = quadrant2[i+1];
//					quadrant2[i+1] = temp;
//					change_made = 1;
//				}
//			}	
//		}
//	}
//
//	while(!quadrant1.IsEmpty())
//	{
//		list.AddToRightEnd(quadrant1.GetLeftEnd());
//		quadrant1.DeleteLeftEnd();
//	}
//	while(!quadrant4.IsEmpty())
//	{
//		list.AddToRightEnd(quadrant4.GetLeftEnd());
//		quadrant4.DeleteLeftEnd();
//	}
//	while(!quadrant3.IsEmpty())
//	{
//		list.AddToRightEnd(quadrant3.GetLeftEnd());
//		quadrant3.DeleteLeftEnd();
//	}
//	while(!quadrant2.IsEmpty())
//	{
//		list.AddToRightEnd(quadrant2.GetLeftEnd());
//		quadrant2.DeleteLeftEnd();
//	}
//}


int bdGeometry::IsValidPosition(int s, int r, int c)
{
	if(s<0 || s>=this->m_dimensions_CRS[2]) return 0;
	if(r<0 || r>=this->m_dimensions_CRS[1]) return 0;
	if(c<0 || c>=this->m_dimensions_CRS[0]) return 0;
	return 1;
}


void bdGeometry::ForCoordinates(int s, int r, int c, unsigned int start_FOR_loop_index)
{
	m_FOR_loop_central_coordinates(s,r,c);
	m_FOR_loop_current_index = start_FOR_loop_index;
}


void bdGeometry::ForCoordinates(int r, int c, unsigned int start_FOR_loop_index)
{
	this->ForCoordinates(0,r,c,start_FOR_loop_index);
}


void bdGeometry::ForCoordinates_4_Neighbors(int r, int c)
{
	this->ForCoordinates(0,r,c,1);
}


void bdGeometry::ForCoordinates_8_Neighbors(int r, int c)
{
	this->ForCoordinates(0,r,c,1);
}


void bdGeometry::ForCoordinates_9_Neighbors(int r, int c)
{
	this->ForCoordinates(0,r,c,0);
}


void bdGeometry::ForCoordinates_26_Neighbors(int s, int r, int c)
{
	this->ForCoordinates(s,r,c,1);
}


void bdGeometry::ForCoordinates_27_Neighbors(int s, int r, int c)
{
	this->ForCoordinates(s,r,c,0);
}


void bdGeometry::ForCoordinates_Circle(int r, int c, unsigned int start_squared_radius)
{
	this->ForCoordinates(0,r,c,this->CircleSquaredRadiusStartIndex(start_squared_radius));
}


void bdGeometry::ForCoordinates_Sphere(int s, int r, int c, unsigned int start_squared_radius)
{
	this->ForCoordinates(s,r,c,this->SphereSquaredRadiusStartIndex(start_squared_radius));
}


void bdGeometry::ForCoordinates_Line(int s, int r, int c, int direction_s, int direction_r, int direction_c)
{
    this->ForCoordinates(s,r,c,0);
    bdDiscreteCoordinates3D direction;
    direction.S() = direction_s - s;
    if(direction.S()<0) direction.S() = -direction.S();
    direction.R() = direction_r - r;
    if(direction.R()<0) direction.R() = -direction.R();
    direction.C() = direction_c - c;
    if(direction.C()<0) direction.C() = -direction.C();
    
    this->MakeLineFramePattern(direction, m_line_frame_pattern);
    
    for(unsigned int i=0; i<m_line_frame_pattern.GetNumberOfElements(); i++)
    {
        if(direction_s-s<0) m_line_frame_pattern[i].S() = -m_line_frame_pattern[i].S();
        if(direction_r-r<0) m_line_frame_pattern[i].R() = -m_line_frame_pattern[i].R();
        if(direction_c-c<0) m_line_frame_pattern[i].C() = -m_line_frame_pattern[i].C();
    }
}


unsigned int bdGeometry::GetCurrentFORLoopIndex()
{
	return m_FOR_loop_current_index;
}


int bdGeometry::Get_4_Neighbors(int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->Neighborhood4EndIndex())
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.R(),m_FOR_loop_central_coordinates.C(),
			this->Neighborhood4(m_FOR_loop_current_index), r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_4_Neighbors(bdDiscreteCoordinates3D &output_point)
{
	return this->Get_4_Neighbors(output_point.R(),output_point.C());
}


int bdGeometry::Get_8_Neighbors(int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->Neighborhood8EndIndex())
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.R(),m_FOR_loop_central_coordinates.C(),
			this->Neighborhood8(m_FOR_loop_current_index), r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_8_Neighbors(bdDiscreteCoordinates3D &output_point)
{
	return this->Get_8_Neighbors(output_point.R(),output_point.C());
}


int bdGeometry::Get_9_Neighbors(int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->Neighborhood9EndIndex())
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.R(),m_FOR_loop_central_coordinates.C(),
			this->Neighborhood9(m_FOR_loop_current_index), r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_9_Neighbors(bdDiscreteCoordinates3D &output_point)
{
	return this->Get_9_Neighbors(output_point.R(),output_point.C());
}


int bdGeometry::Get_26_Neighbors(int &s, int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->Neighborhood26EndIndex())
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.S(),m_FOR_loop_central_coordinates.R(),
			m_FOR_loop_central_coordinates.C(), this->Neighborhood26(m_FOR_loop_current_index), s,r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_26_Neighbors(bdDiscreteCoordinates3D &output_point)
{
	return this->Get_26_Neighbors(output_point.S(), output_point.R(),output_point.C());
}


int bdGeometry::Get_27_Neighbors(int &s, int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->Neighborhood27EndIndex())
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.S(),m_FOR_loop_central_coordinates.R(),
			m_FOR_loop_central_coordinates.C(), this->Neighborhood27(m_FOR_loop_current_index), s,r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_27_Neighbors(bdDiscreteCoordinates3D &output_point)
{
	return this->Get_27_Neighbors(output_point.S(), output_point.R(),output_point.C());
}


int bdGeometry::Get_Circle(unsigned int end_square_radius_included, int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->CircleSquaredRadiusEndIndex(end_square_radius_included))
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.R(),m_FOR_loop_central_coordinates.C(), this->Circle(m_FOR_loop_current_index), r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_Circle(unsigned int end_square_radius_included, bdDiscreteCoordinates3D &output_point)
{
	return this->Get_Circle(end_square_radius_included, output_point.R(),output_point.C());
}


int bdGeometry::Get_Sphere(unsigned int end_square_radius_included, int &s, int &r, int &c)
{
	while(m_FOR_loop_current_index <= this->SphereSquaredRadiusEndIndex(end_square_radius_included))
	{
		if(this->GetCoordinates(m_FOR_loop_central_coordinates.S(),m_FOR_loop_central_coordinates.R(),
			m_FOR_loop_central_coordinates.C(), this->Sphere(m_FOR_loop_current_index), s,r,c) )
		{
			m_FOR_loop_current_index++;
			return 1;
		}
		else m_FOR_loop_current_index++;
	}
	return 0;
}


int bdGeometry::Get_Sphere(unsigned int end_square_radius_included, bdDiscreteCoordinates3D &output_point)
{
	return this->Get_Sphere(end_square_radius_included, output_point.S(),output_point.R(),output_point.C());
}


int bdGeometry::Get_Line(unsigned int end_squared_distance_included, int &s, int &r, int &c)
{
    unsigned int line_frame_index = m_FOR_loop_current_index % m_line_frame_pattern.GetNumberOfElements();
    unsigned int line_segment_index = m_FOR_loop_current_index / m_line_frame_pattern.GetNumberOfElements();
    
    s = m_FOR_loop_central_coordinates.S() + line_segment_index*(m_line_frame_pattern[(m_line_frame_pattern.GetNumberOfElements()-1)].S()) + m_line_frame_pattern[line_frame_index].S();
    r = m_FOR_loop_central_coordinates.R() + line_segment_index*(m_line_frame_pattern[(m_line_frame_pattern.GetNumberOfElements()-1)].R()) + m_line_frame_pattern[line_frame_index].R();
    c = m_FOR_loop_central_coordinates.C() + line_segment_index*(m_line_frame_pattern[(m_line_frame_pattern.GetNumberOfElements()-1)].C()) + m_line_frame_pattern[line_frame_index].C();
    
    unsigned int d = (s - m_FOR_loop_central_coordinates.S())*(s - m_FOR_loop_central_coordinates.S()) +
        (r - m_FOR_loop_central_coordinates.R())*(r - m_FOR_loop_central_coordinates.R()) +
        (c - m_FOR_loop_central_coordinates.C())*(c - m_FOR_loop_central_coordinates.C());
    
    //cout<<" src ="<<s<<","<<r<<","<<c<<", d="<<d<<endl;
    if(d>end_squared_distance_included) return 0;
    
    m_FOR_loop_current_index++;
    
    return this->IsValidPosition(s,r,c);
}


int bdGeometry::Get_Line(unsigned int end_squared_distance_included, bdDiscreteCoordinates3D &output_point)
{
    return this->Get_Line(end_squared_distance_included, output_point.S(),output_point.R(),output_point.C());
}



unsigned int bdGeometry::GetNumberOfCircleElements()
{
    return m_circle->m_number_of_elements;
}


unsigned int bdGeometry::GetNumberOfSphereElements()
{
	return m_sphere->m_number_of_elements;
}


int bdGeometry::GetCoordinates(bdDiscreteCoordinates3D &start_point, bdDiscreteCoordinates3D &relative_point, bdDiscreteCoordinates3D &output_point)
{
	output_point.S() = start_point.S() + relative_point.S();
	output_point.R() = start_point.R() + relative_point.R();
	output_point.C() = start_point.C() + relative_point.C();
	if(output_point.S()>=0 && output_point.S()<m_dimensions_CRS[2] && 
		output_point.R()>=0 && output_point.R()<m_dimensions_CRS[1] && 
			output_point.C()>=0 && output_point.C()<m_dimensions_CRS[0]) return 1;
	else return 0;
}


int bdGeometry::GetCoordinates(int s_start, int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, bdDiscreteCoordinates3D &output_point)
{
	output_point.S() = s_start + relative_point.S();
	output_point.R() = r_start + relative_point.R();
	output_point.C() = c_start + relative_point.C();
	if(output_point.S()>=0 && output_point.S()<m_dimensions_CRS[2] && 
		output_point.R()>=0 && output_point.R()<m_dimensions_CRS[1] && 
			output_point.C()>=0 && output_point.C()<m_dimensions_CRS[0]) return 1;
	else return 0;
}


int bdGeometry::GetCoordinates(int s_start, int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, int &output_s, int &output_r, int &output_c)
{
	output_s = s_start + relative_point.S();
	output_r = r_start + relative_point.R();
	output_c = c_start + relative_point.C();
	if(output_s>=0 && output_s<m_dimensions_CRS[2] && 
		output_r>=0 && output_r<m_dimensions_CRS[1] && 
			output_c>=0 && output_c<m_dimensions_CRS[0]) return 1;
	else return 0;
}


int bdGeometry::GetCoordinates(int r_start, int c_start, bdDiscreteCoordinates3D &relative_point, int &output_r, int &output_c)
{
	output_r = r_start + relative_point.R();
	output_c = c_start + relative_point.C();
	if(output_r>=0 && output_r<m_dimensions_CRS[1] && 
		output_c>=0 && output_c<m_dimensions_CRS[0]) return 1;
	else return 0;
}


int bdGeometry::NumberOfValid_9_Positions(int r, int c)
{
	if(!this->IsDimensions2DSet()) return 0;

	//DON'T CHANGE THE ORDER OF THE FOLLOWING CONDITIONS!!!
	int number_of_pixels = 9;
	int number_of_reductions = 0;
	if(r==0 || r==m_dimensions_CRS[1]-1) 
	{
		number_of_pixels-=3;
		number_of_reductions++;
	}
	if(c==0 || c==m_dimensions_CRS[0]-1) 
	{
		if(number_of_reductions==0)
		{
			number_of_pixels-=3;
		}
		if(number_of_reductions==1)
		{
			number_of_pixels-=2;
		}
	}
	return number_of_pixels;
}


int bdGeometry::NumberOfValid_27_Positions(int s, int r, int c)
{
	if(!this->IsDimensions3DSet()) return 0;

	//DON'T CHANGE THE ORDER OF THE FOLLOWING CONDITIONS!!!
	int number_of_voxels = 27;
	int number_of_reductions = 0;
	if(s==0 || s==m_dimensions_CRS[2]-1) 
	{
		number_of_voxels-=9;
		number_of_reductions++;
	}
	if(r==0 || r==m_dimensions_CRS[1]-1) 
	{
		if(number_of_reductions==0)
		{
			number_of_voxels-=9;
			number_of_reductions++;
		}
		else
		{
			if(number_of_reductions==1)
			{
				number_of_voxels-=6;
				number_of_reductions++;						
			}
		}
	}
	if(c==0 || c==m_dimensions_CRS[0]-1) 
	{
		if(number_of_reductions==0) { number_of_voxels-=9; }
		if(number_of_reductions==1) { number_of_voxels-=6; }
		if(number_of_reductions==2) { number_of_voxels-=4; }
	}

	return number_of_voxels;
}


int bdGeometry::Make3DCross(bdDiscreteCoordinates3D &starting_point, bdArray<bdDiscreteCoordinates3D> &cross_points_array, int cross_size)
{
	if(!this->IsDimensions3DSet()) return 0;
	if(cross_size>=6 || cross_size<=0) return 0;

	std::list<bdDiscreteCoordinates3D> l;
	bdDiscreteCoordinates3D s,r,o,p;
	//s(starting_point.s,starting_point.r,starting_point.c);
	s = starting_point;

	r(0,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	
	if(cross_size>=1)
	{
		r(1,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(-1,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,1,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,-1,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,1); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,-1); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	}
	if(cross_size>=2)
	{
		r(2,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(-2,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,2,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,-2,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,2); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,-2); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	}
	if(cross_size>=3)
	{
		r(3,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(-3,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,3,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,-3,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,3); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,-3); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	}
	if(cross_size>=4)
	{
		r(4,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(-4,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,4,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,-4,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,4); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,-4); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	}
	if(cross_size>=5)
	{
		r(5,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(-5,0,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,5,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,-5,0); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,5); if(this->GetCoordinates(s,r,o)) l.push_front(o);
		r(0,0,-5); if(this->GetCoordinates(s,r,o)) l.push_front(o);
	}
	
	cross_points_array.Set((unsigned int)(l.size()));
	int i=0;
	for(std::list<bdDiscreteCoordinates3D>::iterator it = l.begin(); it != l.end(); it++)
	{
		cross_points_array[i] = (*it);
		i++;
	}
	return 1;
}




int bdGeometry::Make3DLineSegment(bdDiscreteCoordinates3D &start_point, bdDiscreteCoordinates3D &end_point, bdArray<bdDiscreteCoordinates3D> &line_segment_points_array)
{
	//cout<<"start_point="<<start_point<<"  end_point="<<end_point<<endl;
	if(!this->IsDimensions3DSet()) return 0;
	if(start_point.S()<0 || start_point.S()>=this->m_dimensions_CRS[2]) return 0;
	if(start_point.R()<0 || start_point.R()>=this->m_dimensions_CRS[1]) return 0;
	if(start_point.C()<0 || start_point.C()>=this->m_dimensions_CRS[0]) return 0;

	line_segment_points_array.Reset();

	//c = c_0 + a_c*t
	//r = r_0 + a_r*t
	//s = s_0 + a_s*t
	// for t=[0,1]: a = 'end_coordinate' - 'start_coordinate'; example: a_s = end_point.s - start_point.s

	//Number of steps for 't' will be the distance in 3D between the start and end point
	double distance = start_point.DistanceEuclidean(end_point);
	//cout<<"distance="<<distance<<" ";
	int number_of_steps = (int)(distance+1);//Make it distance+1, so that it's always rounded to the first greater value
	//cout<<"number_of_steps="<<number_of_steps<<" ";


	double a_s = end_point.S() - start_point.S();
	double a_r = end_point.R() - start_point.R();
	double a_c = end_point.C() - start_point.C();

	double d_s = a_s/((double)(number_of_steps));
	double d_r = a_r/((double)(number_of_steps));
	double d_c = a_c/((double)(number_of_steps));

	//cout<<"a_s="<<a_s<<" a_r="<<a_r<<" a_c="<<a_c<<endl;
	//cout<<"d_s="<<d_s<<" d_r="<<d_r<<" d_c="<<d_c<<endl;

	//Initialize a temporary array
	bdArray<bdDiscreteCoordinates3D> temp;
	temp.Set(number_of_steps+1);//!!!!!!!!!!!!!!!!!!!!OVDE DODAO 1  - izgleda da radi!!!
	temp[0] = start_point;

	bdDiscreteCoordinates3D p = start_point;

	//Because of a small step of sampling the line, there will be positions that are the same for 2 or more calculations - we want to
	// know how many of different positions exist for setting the length of the output array
	int number_of_points_that_are_NOT_the_same = 1;
	bdDiscreteCoordinates3D previous_point = start_point;
	int i=1;
	while((p.S()!=end_point.S() || p.R()!=end_point.R() || p.C()!=end_point.C()) && i<number_of_steps+1)
	{
		p.S() = (int)(start_point.S() + ((double)(i))*d_s);
		p.R() = (int)(start_point.R() + ((double)(i))*d_r);
		p.C() = (int)(start_point.C() + ((double)(i))*d_c);

		if(p.S()!=previous_point.S() || p.R()!=previous_point.R() || p.C()!=previous_point.C())
		{
			temp[number_of_points_that_are_NOT_the_same] = p;
			number_of_points_that_are_NOT_the_same++;
			previous_point = p;
		}

		i++;
	}

	line_segment_points_array.Set(number_of_points_that_are_NOT_the_same);
	
	for(i=0; i<number_of_points_that_are_NOT_the_same; i++)
	{
		line_segment_points_array[i] = temp[i];
	}	

	return 1;
}


int bdGeometry::MakeLineFramePattern(bdDiscreteCoordinates3D &input_direction, bdArray<bdDiscreteCoordinates3D> &output_line_frame_pattern)
{
	if(input_direction.S()==0 && input_direction.R()==0 && input_direction.C()==0) return 0;

    //cout<<" 1 ";
    
	bdGeometry g;
	g.SetDimensions(200,200,200);//these dimensions should be enough...

	bdDiscreteCoordinates3D p1, p2;
	p1(100,100,100); p2(100+input_direction.S(), 100+input_direction.R(), 100+input_direction.C());//start and end of Vector for direction
	bdArray<bdDiscreteCoordinates3D> line_segment;
	g.Make3DLineSegment(p1,p2,line_segment);//Make a simple line segment

	//The line segment has to be refined to be always one pixel wide. We do this 
	// with tracking the highest distance voxels from the start voxel

    //cout<<" 2 ";

	bdMatrix<int> line_segment_image;//bdRegularGridT<int> line_segment_image;
	//Find the size of the line segment image
	int ds = input_direction.S(); if(ds<0) ds = -ds;//take the absolute value
	int dr = input_direction.R(); if(dr<0) dr = -dr;
	int dc = input_direction.C(); if(dc<0) dc = -dc;
	line_segment_image.SetSize(ds+1,dr+1,dc+1);//line_segment_image.SetSizeOf3DGrid(ds+1,dr+1,dc+1);
	line_segment_image.FillInWith(0);
	for(unsigned int i=0; i<line_segment.GetNumberOfElements(); i++)
	{ line_segment_image((line_segment[i].S()-100),(line_segment[i].R()-100),(line_segment[i].C()-100)) = 255; }
	
	//To find the shortest path it is enough to find (squared) Euclidean distances of each
	// point from the end point and do backtracking on these values. This will be enough 
	// because all voxels in line segment follow the same direction (greedy algorithm can work here).
	//Do the region growing with backtracking to find the shortest path.
	
	bdMatrix<int> line_segment_distances_image;
	line_segment_distances_image.SetSize(ds+1,dr+1,dc+1);//line_segment_distances_image.SetSizeOf3DGrid(ds+1,dr+1,dc+1);
	line_segment_distances_image.FillInWith(0);

    //cout<<" 3 ";
    
	//Calculate the squared Euclidean distances
	for(unsigned int i=0; i<line_segment.GetNumberOfElements(); i++)
	{		
		int d = p1.DistanceEuclideanSquared(line_segment[i]);
		line_segment_distances_image((line_segment[i].S()-100),(line_segment[i].R()-100),(line_segment[i].C()-100)) = d;
	}
    
    
    //cout<<" 4 ";

	//Start from start voxel and take neighboring voxel with higher squared Euclidean distance as the
	// part of the shortest path
	bdList<bdDiscreteCoordinates3D> l;//This is the output in the form of list - we will later copy it to the array
	bdDiscreteCoordinates3D p, end_point;
	end_point(p2.S()-100,p2.R()-100,p2.C()-100);
	l.AddToLeftEnd(p(0,0,0));//This point will be in the output for sure.
	g.SetDimensions(line_segment_image.GetNumberOfSlices(), line_segment_image.GetNumberOfRows(), line_segment_image.GetNumberOfColumns());
	//while(!l.IsEmpty())
	while(!p.HasEqualCoordinatesAs(end_point))
	{
		int max_distance = 0;
		bdDiscreteCoordinates3D max_distance_voxel;
		max_distance_voxel(0,0,0);
		int sn, rn, cn;
		for(g.ForCoordinates_26_Neighbors(p.S(),p.R(),p.C()); g.Get_26_Neighbors(sn,rn,cn); )
		{			
			if(line_segment_image(sn,rn,cn)!=0)//voxel is a part of the line segment
			{
				//...check if it has the highest distance value...
				if(line_segment_distances_image(sn,rn,cn)>max_distance) 
				{
					//... if so, record the max distance and its voxel.
					max_distance = line_segment_distances_image(sn,rn,cn);
					max_distance_voxel(sn,rn,cn);
				}
			}
		}
		l.AddToRightEnd(max_distance_voxel);
		p = max_distance_voxel;
	}
    
    
    //cout<<" 5 ";

	//At this moment the refined line segment indexes are in the list 'l', so copy the list to the output
	output_line_frame_pattern.Set(l.GetNumberOfElements());
	//bdDiscreteCoordinates3D *pp;
	int i=0;
	bdListIterator<bdDiscreteCoordinates3D> it;
	for(it.SetLeftEnd(l); it.IsValid(); it.MoveRight())
	{ 
		output_line_frame_pattern[i] = it.GetElement();
		i++;
	}
	//for(l.For_StartFromLeftEnd(); l.GetAllPointersToElemetsByMovingToRight(&pp); i++)
	//	{ output_line_frame_pattern[i] = (*pp); }
    
    //cout<<" 6 ";

	return 1;
}


int bdGeometry::FilterMultipleDirections(bdArray<bdDiscreteCoordinates3D> &input_array_of_directions, bdArray<bdDiscreteCoordinates3D> &output_array_of_directions)
{
	if(input_array_of_directions.IsEmpty()) return 0;

	bdList<bdDiscreteCoordinates3D> l;
	//for(unsigned int i=0; i<input_array_of_directions.GetNumberOfElements()-1; i++)
	for(unsigned int i=0; i<input_array_of_directions.GetNumberOfElements(); i++)
	{
		int is_duplicate_found = 0;
		//for(unsigned int j=i+1; j<input_array_of_directions.GetNumberOfElements() && !is_duplicate_found; j++)
		//{
		//	//check if the directions are the same
		//	if( (input_array_of_directions[i].S()==input_array_of_directions[j].S()
		//		&& input_array_of_directions[i].R()==input_array_of_directions[j].R()
		//		&& input_array_of_directions[i].C()==input_array_of_directions[j].C()) 
		//		|| (input_array_of_directions[i].S()== - input_array_of_directions[j].S()
		//		&& input_array_of_directions[i].R()== - input_array_of_directions[j].R()
		//		&& input_array_of_directions[i].C()== - input_array_of_directions[j].C()) ) 
		//			is_duplicate_found = 1;
		//}

		bdListIterator<bdDiscreteCoordinates3D> it;
		for(it.SetLeftEnd(l); it.IsValid(); it.MoveRight())
		{
			bdDiscreteCoordinates3D *pc = it.GetElementPointer();
			//check if the directions are the same
			if( (input_array_of_directions[i].S()==pc->S() && input_array_of_directions[i].R()==pc->R() && input_array_of_directions[i].C()==pc->C()) 
				|| (input_array_of_directions[i].S()== - pc->S() && input_array_of_directions[i].R()== - pc->R() && input_array_of_directions[i].C()== - pc->C()) ) 
					is_duplicate_found = 1;
		}

		//If a duplicate for a direction is not found, copy the direction with index 'i' into the list.
		if(!is_duplicate_found) l.AddToRightEnd(input_array_of_directions[i]);
	}
	////Add the last direction in the array to the list
	//l.AddToRightEnd(input_array_of_directions[(input_array_of_directions.GetNumberOfElements()-1)]);

	//Copy the elements from the list to output array
	output_array_of_directions.Set(l.GetNumberOfElements());
	for(unsigned int i=0; i<output_array_of_directions.GetNumberOfElements(); i++)
	{
		output_array_of_directions[i] = l.GetLeftEnd();
		l.DeleteLeftEnd();
	}

	return 1;
}




//void bdGeometry::GenerateArraysOfPositionsForDoubleRadiusValues(int maximum_radius, bdArray<bdArray<bdDiscreteCoordinates3D>> &array_of_s_r_c_coordinates)
//{
//	int radius_max = maximum_radius*maximum_radius;//Square of the input radius 
//
//	
//	bdArray<DoubleList<bdDiscreteCoordinates3D>> coordinates_array;
//	coordinates_array.Set(radius_max+2);
//
//	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
//	p_zero(0,0,0);
//	int distance;
//	for(int i=0; i<coordinates_array.NumberOfElements(); i++)
//	{
//		for(int s=-i; s<i+1; s++)
//		{
//			for(int r=-i; r<i+1; r++)
//			{
//				for(int c=-i; c<i+1; c++)
//				{
//					distance = p_zero.Distance_EuclideanSquared(p_relative(s,r,c));
//					if(distance==i)
//					{
//						pp = coordinates_array[i].AddNewToRightEnd();
//						(*pp)(s,r,c);//Set the coordinates
//						pp->value = distance;//Set the squared Euclidean distance from center point
//					}
//				}
//			}
//		}
//		cout<<" "<<coordinates_array[i].NumberOfElements();
//	}
//
//	//Now analyse the coordinates_array, 'n' will be the number of non-empty lists in the array
//	int n=0;
//	for(int i=0; i<coordinates_array.NumberOfElements(); i++)
//	{
//		if(coordinates_array[i].NumberOfElements()!=0) n++;
//	}
//
//	cout<<endl<<"n="<<n<<endl;
//	//Change the input array based on the 'n' and coordinates_array
//	array_of_s_r_c_coordinates.Set(n);
//	int j=0;
//	for(int i=0; i<coordinates_array.NumberOfElements() && j<n; i++)
//	{
//		if(coordinates_array[i].NumberOfElements()!=0)
//		{
//			array_of_s_r_c_coordinates[j].Set(coordinates_array[i].NumberOfElements());
//			for(int k=0; k<array_of_s_r_c_coordinates[j].NumberOfElements(); k++)
//			{
//				array_of_s_r_c_coordinates[j][k] = coordinates_array[i].GetLeftEnd();
//				coordinates_array[i].DeleteLeftEnd();//IsEmpty the list
//			}
//			j++;
//		}
//	}
//}


unsigned int bdGeometry::SphereSquaredRadiusStartIndex(unsigned int squared_radius)
{
	return (m_sphere->SizeStartIndex(squared_radius));
}


unsigned int bdGeometry::SphereSquaredRadiusEndIndex(unsigned int squared_radius)
{
	return (m_sphere->SizeEndIndex(squared_radius));
}


bdDiscreteCoordinates3D& bdGeometry::Sphere(unsigned int index_of_point)
{
	return (*m_sphere)[index_of_point];
}


unsigned int bdGeometry::CircleSquaredRadiusStartIndex(unsigned int squared_radius)
{
	return (m_circle->SizeStartIndex(squared_radius));
}


unsigned int bdGeometry::CircleSquaredRadiusEndIndex(unsigned int squared_radius)
{
	return (m_circle->SizeEndIndex(squared_radius));
}


bdDiscreteCoordinates3D& bdGeometry::Circle(unsigned int index_of_point)
{
	return (*m_circle)[index_of_point];
}


unsigned int bdGeometry::Circle_R_SquaredRadiusStartIndex(unsigned int squared_radius)
{
	return (m_circle_R->SizeStartIndex(squared_radius));
}


unsigned int bdGeometry::Circle_R_SquaredRadiusEndIndex(unsigned int squared_radius)
{
	return (m_circle_R->SizeEndIndex(squared_radius));
}


bdDiscreteCoordinates3D& bdGeometry::Circle_R(unsigned int index_of_point)
{
	return (*m_circle_R)[index_of_point];
}


unsigned int bdGeometry::Circle_C_SquaredRadiusStartIndex(unsigned int squared_radius)
{
	return (m_circle_C->SizeStartIndex(squared_radius));
}


unsigned int bdGeometry::Circle_C_SquaredRadiusEndIndex(unsigned int squared_radius)
{
	return (m_circle_C->SizeEndIndex(squared_radius));
}


bdDiscreteCoordinates3D& bdGeometry::Circle_C(unsigned int index_of_point)
{
	return (*m_circle_C)[index_of_point];
}


unsigned int bdGeometry::CubeSizeStartIndex(unsigned int size)
{
	return (m_cube->SizeStartIndex(size));
}


unsigned int bdGeometry::CubeSizeEndIndex(unsigned int size)
{
	return (m_cube->SizeEndIndex(size));
}


bdDiscreteCoordinates3D& bdGeometry::Cube(unsigned int index_of_point)
{
	return (*m_cube)[index_of_point];
}


unsigned int bdGeometry::SquareSizeStartIndex(unsigned int size)
{
	return (m_square->SizeStartIndex(size));
}


unsigned int bdGeometry::SquareSizeEndIndex(unsigned int size)
{
	return (m_square->SizeEndIndex(size));
}


bdDiscreteCoordinates3D& bdGeometry::Square(unsigned int index_of_point)
{
	return (*m_square)[index_of_point];
}


unsigned int bdGeometry::Neighborhood4StartIndex()
{
	return (m_circle->SizeStartIndex(1));
}


unsigned int bdGeometry::Neighborhood4EndIndex()
{
	return (m_circle->SizeEndIndex(1));
}


bdDiscreteCoordinates3D& bdGeometry::Neighborhood4(unsigned int index_of_point)
{
	return (*m_circle)[index_of_point];
}


unsigned int bdGeometry::Neighborhood8StartIndex()
{
	return (m_square->SizeStartIndex(1));
}


unsigned int bdGeometry::Neighborhood8EndIndex()
{
	return (m_square->SizeEndIndex(1));
}


bdDiscreteCoordinates3D& bdGeometry::Neighborhood8(unsigned int index_of_point)
{
	return (*m_square)[index_of_point];
}


unsigned int bdGeometry::Neighborhood9StartIndex()
{
	return (m_square->SizeStartIndex(0));
}


unsigned int bdGeometry::Neighborhood9EndIndex()
{
	return (m_square->SizeEndIndex(1));
}


bdDiscreteCoordinates3D& bdGeometry::Neighborhood9(unsigned int index_of_point)
{
	return (*m_square)[index_of_point];
}


unsigned int bdGeometry::Neighborhood27StartIndex()
{
	return (m_cube->SizeStartIndex(0));
}


unsigned int bdGeometry::Neighborhood27EndIndex()
{
	return (m_cube->SizeEndIndex(1));
}


bdDiscreteCoordinates3D& bdGeometry::Neighborhood27(unsigned int index_of_point)
{
	return (*m_cube)[index_of_point];
}


unsigned int bdGeometry::Neighborhood26StartIndex()
{
	return (m_cube->SizeStartIndex(1));
}


unsigned int bdGeometry::Neighborhood26EndIndex()
{
	return (m_cube->SizeEndIndex(1));
}


bdDiscreteCoordinates3D& bdGeometry::Neighborhood26(unsigned int index_of_point)
{
	return (*m_cube)[index_of_point];
}


void bdGeometry::GenerateSphereArray(unsigned int maximum_radius)
{
	int radius_max = maximum_radius*maximum_radius;// Square of the input radius 

	m_sphere->Reset();
	m_sphere->m_size.Set(radius_max+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(radius_max+1);

	unsigned int limit = 2*maximum_radius+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(maximum_radius+1, maximum_radius+1, maximum_radius+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_sphere->m_number_of_elements = 0;
	for(unsigned int s=0; s<limit; s++)
	{
		for(unsigned int r=0; r<limit; r++)
		{
			for(unsigned int c=0; c<limit; c++)
			{
				unsigned int d = p_zero.DistanceEuclideanSquared(p_relative(s,r,c));
				if(d<coordinates_array.GetNumberOfElements())
				{
					pp = coordinates_array[d].AddNewToRightEnd();
					pp->S() = p_zero.S()-s; 
					pp->R() = p_zero.R()-r; 
					pp->C() = p_zero.C()-c;
					m_sphere->m_number_of_elements++;// Increase the total number of points.
				}
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_sphere->m_points = new bdDiscreteCoordinates3D[m_sphere->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_sphere->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_sphere->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}



void bdGeometry::GenerateCircleArray(unsigned int maximum_radius)
{
	int radius_max = maximum_radius*maximum_radius;// Square of the input radius 

	m_circle->Reset();
	m_circle->m_size.Set(radius_max+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(radius_max+1);

	unsigned int limit = 2*maximum_radius+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(0, maximum_radius+1, maximum_radius+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_circle->m_number_of_elements = 0;

	for(unsigned int r=0; r<limit; r++)
	{
		for(unsigned int c=0; c<limit; c++)
		{
			unsigned int d = p_zero.DistanceEuclideanSquared(p_relative(0,r,c));
			if(d<coordinates_array.GetNumberOfElements())
			{
				pp = coordinates_array[d].AddNewToRightEnd();
				pp->S() = 0; 
				pp->R() = p_zero.R()-r; 
				pp->C() = p_zero.C()-c;
				m_circle->m_number_of_elements++;// Increase the total number of points.
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_circle->m_points = new bdDiscreteCoordinates3D[m_circle->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_circle->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_circle->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}


void bdGeometry::GenerateCircle_R_Array(unsigned int maximum_radius)
{
	int radius_max = maximum_radius*maximum_radius;// Square of the input radius 

	m_circle_R->Reset();
	m_circle_R->m_size.Set(radius_max+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(radius_max+1);

	unsigned int limit = 2*maximum_radius+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(0, maximum_radius+1, maximum_radius+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_circle_R->m_number_of_elements = 0;

	for(unsigned int r=0; r<limit; r++)
	{
		for(unsigned int c=0; c<limit; c++)
		{
			unsigned int d = p_zero.DistanceEuclideanSquared(p_relative(0,r,c));
			if(d<coordinates_array.GetNumberOfElements())
			{
				pp = coordinates_array[d].AddNewToRightEnd();
				pp->S() = p_zero.R()-r;  
				pp->R() = 0;
				pp->C() = p_zero.C()-c;
				m_circle_R->m_number_of_elements++;// Increase the total number of points.
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_circle_R->m_points = new bdDiscreteCoordinates3D[m_circle_R->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_circle_R->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_circle_R->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}


void bdGeometry::GenerateCircle_C_Array(unsigned int maximum_radius)
{
	int radius_max = maximum_radius*maximum_radius;// Square of the input radius 

	m_circle_C->Reset();
	m_circle_C->m_size.Set(radius_max+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(radius_max+1);

	unsigned int limit = 2*maximum_radius+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(0, maximum_radius+1, maximum_radius+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_circle_C->m_number_of_elements = 0;

	for(unsigned int r=0; r<limit; r++)
	{
		for(unsigned int c=0; c<limit; c++)
		{
			unsigned int d = p_zero.DistanceEuclideanSquared(p_relative(0,r,c));
			if(d<coordinates_array.GetNumberOfElements())
			{
				pp = coordinates_array[d].AddNewToRightEnd();
				pp->S() = p_zero.R()-r;  
				pp->R() = p_zero.C()-c;
				pp->C() = 0;
				m_circle_C->m_number_of_elements++;// Increase the total number of points.
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_circle_C->m_points = new bdDiscreteCoordinates3D[m_circle_C->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_circle_C->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_circle_C->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}


void bdGeometry::GenerateCubeArray(unsigned int maximum_size)
{
	//int radius_max = maximum_radius*maximum_radius;// Square of the input radius 

	m_cube->Reset();
	m_cube->m_size.Set(maximum_size+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(maximum_size+1);

	unsigned int limit = 2*maximum_size+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(maximum_size+1, maximum_size+1, maximum_size+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_cube->m_number_of_elements = 0;
	for(unsigned int s=0; s<limit; s++)
	{
		for(unsigned int r=0; r<limit; r++)
		{
			for(unsigned int c=0; c<limit; c++)
			{
				unsigned int d = p_zero.DistanceMaximumNorm(p_relative(s,r,c));
				if(d<coordinates_array.GetNumberOfElements())
				{
					pp = coordinates_array[d].AddNewToRightEnd();
					pp->S() = p_zero.S()-s; 
					pp->R() = p_zero.R()-r; 
					pp->C() = p_zero.C()-c;
					m_cube->m_number_of_elements++;// Increase the total number of points.
				}
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_cube->m_points = new bdDiscreteCoordinates3D[m_cube->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_cube->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_cube->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}


void bdGeometry::GenerateSquareArray(unsigned int maximum_size)
{
	m_square->Reset();
	m_square->m_size.Set(maximum_size+1);
	
	bdArray< bdList<bdDiscreteCoordinates3D> > coordinates_array;
	coordinates_array.Set(maximum_size+1);

	unsigned int limit = 2*maximum_size+3;
	bdDiscreteCoordinates3D p_zero, p_relative, *pp;
	p_zero(0, maximum_size+1, maximum_size+1);

	// Generate all distances and enter them to the array of lists 'coordinates_array'.
	m_square->m_number_of_elements = 0;
	
	for(unsigned int r=0; r<limit; r++)
	{
		for(unsigned int c=0; c<limit; c++)
		{
			unsigned int d = p_zero.DistanceMaximumNorm(p_relative(0,r,c));
			if(d<coordinates_array.GetNumberOfElements())
			{
				pp = coordinates_array[d].AddNewToRightEnd();
				pp->S() = 0; 
				pp->R() = p_zero.R()-r; 
				pp->C() = p_zero.C()-c;
				m_square->m_number_of_elements++;// Increase the total number of points.
			}
		}
	}

	// Allocate the array to store the found number of elements
	m_square->m_points = new bdDiscreteCoordinates3D[m_square->m_number_of_elements];

	//Copy the points into the container array and set its 'size' index array.
	int current_number_of_points = 0;
	for(unsigned int c_size=0; c_size<coordinates_array.GetNumberOfElements(); c_size++)
	{
		m_square->m_size[c_size] = current_number_of_points;
		while(!coordinates_array[c_size].IsEmpty())
		{
			m_square->m_points[current_number_of_points] = coordinates_array[c_size].GetLeftEnd();
			coordinates_array[c_size].DeleteLeftEnd();// Empty the list
			current_number_of_points++;
		}	
	}
}


int bdGeometry::Save(bdGeometryPointContainer &point_container, const char *file_name_root, const char *extension)
{
	bdString bds_name;
	bds_name.Append(file_name_root); bds_name.Append("."); bds_name.Append(extension);
	std::ofstream file(bds_name.C_String(),ios::binary);
	if(!file)
	{
		std::cout<<"Error: bdGeometry::Save(bdGeometryPointContainer&, char*, char*): Unable to open "<<bds_name.C_String()<<" !"<<endl;
		return 0;
	}

	// Save the size of 'size' array
	unsigned int n = point_container.m_size.GetNumberOfElements();
	file.write((char*) &n, sizeof(int));

	// Save the size of 'points' array
	n = point_container.m_number_of_elements;
	file.write((char*) &n, sizeof(int));
	
	//Save the indexes of 'size' array
	for(unsigned int i=0; i<point_container.m_size.GetNumberOfElements(); i++)
	{
		n = point_container.m_size[i];
		file.write((char*) &n, sizeof(int));
	}

	//Save all the points consisting of s,r,c coordinates
	for(unsigned int i=0; i<point_container.m_number_of_elements; i++)
	{
		file.write((char*) &(point_container.m_points[i].S()), sizeof(int));
		file.write((char*) &(point_container.m_points[i].R()), sizeof(int));
		file.write((char*) &(point_container.m_points[i].C()), sizeof(int));
		//file.write((char*) &(sphere[i][j].value), sizeof(int));
	}
	
	////Print the size arrays
	//file<<"sphere.Set("<<array_of_s_r_c_coordinates.NumberOfElements()<<");"<<endl;
	//for(int i=0; i<array_of_s_r_c_coordinates.NumberOfElements(); i++)
	//{
	//	if(array_of_s_r_c_coordinates[i].NumberOfElements()!=0)
	//	{
	//		file<<"sphere["<<i<<"].Set("<<array_of_s_r_c_coordinates[i].NumberOfElements()<<");"<<endl;
	//	}
	//}
	//
	//for(int i=0; i<array_of_s_r_c_coordinates.NumberOfElements(); i++)
	//{
	//	for(int j=0; j<array_of_s_r_c_coordinates[i].NumberOfElements(); j++)
	//	{
	//		file<<"sphere["<<i<<"]["<<j<<"]="<<array_of_s_r_c_coordinates[i][j]<<";";
	//	}
	//	file<<endl;
	//}

	file.close();
	return 1;
}


int bdGeometry::SaveSphereArray(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_sphere,file_name_root,extension));
}


int bdGeometry::SaveCircleArray(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_circle,file_name_root,extension));
}


int bdGeometry::SaveCircle_R_Array(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_circle_R,file_name_root,extension));
}


int bdGeometry::SaveCircle_C_Array(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_circle_C,file_name_root,extension));
}


int bdGeometry::SaveCubeArray(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_cube,file_name_root,extension));
}


int bdGeometry::SaveSquareArray(const char *file_name_root, const char *extension)
{
	return (this->Save(*m_square,file_name_root,extension));
}


int bdGeometry::Load(bdGeometryPointContainer &point_container, const char *file_name_root, const char *extension)
{
	point_container.Reset();
	
	bdString bds_name;
	bds_name.Append(file_name_root); bds_name.Append("."); bds_name.Append(extension);
	std::ifstream file(bds_name.C_String(),ios::binary);
	if(!file)
	{
		std::cout<<"bdGeometry::Load(bdGeometryPointContainer&, char*, char*): Unable to open "<<bds_name.C_String()<<" !"<<endl;
		return 0;
	}
	
	// Load the size of the 'size' array
	int n;
	file.read((char*) &n, sizeof(int));
	point_container.m_size.Set(n);

	// Load the size of 'points' array
	file.read((char*) &n, sizeof(int));
	point_container.i_SetPointsArray(n);
	
	// Load the indexes of 'size' array
	for(unsigned int i=0; i<point_container.m_size.GetNumberOfElements(); i++)
	{
		file.read((char*) &n, sizeof(int));
		point_container.m_size[i] = n;
	}

	// Load all the points consisting of s,r,c coordinates
	for(unsigned int i=0; i<point_container.m_number_of_elements; i++)
	{
		file.read((char*) &n, sizeof(int));
		point_container.m_points[i].S() = n;
		file.read((char*) &n, sizeof(int));
		point_container.m_points[i].R() = n;
		file.read((char*) &n, sizeof(int));
		point_container.m_points[i].C() = n;
	}
	file.close();
	return 1;
}


int bdGeometry::LoadSphereArray(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_sphere,file_name_root,extension));
}


int bdGeometry::LoadCircleArray(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_circle,file_name_root,extension));
}


int bdGeometry::LoadCircle_R_Array(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_circle_R,file_name_root,extension));
}


int bdGeometry::LoadCircle_C_Array(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_circle_C,file_name_root,extension));
}


int bdGeometry::LoadCubeArray(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_cube,file_name_root,extension));
}


int bdGeometry::LoadSquareArray(const char *file_name_root, const char *extension)
{
	return (this->Load(*m_square,file_name_root,extension));
}





//void GenerateArraysOfPositionsForDoubleRadiusValues(int maximum_radius, bdArray<bdArray<bdDiscreteCoordinates3D>> &array_of_s_r_c_coordinates);


//	void ListOfClockwisePositionsOfSquareRadius(DoubleList<bdDiscreteCoordinates3D> &list, int r, int c, int square_r0, int square_r1);
