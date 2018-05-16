/**************************************************************************
 Class of 2D region 

 Author: Danilo Babin
 File name: "bdRegion2D.cpp"
***************************************************************************/


#include "bdRegion2D.h"





bdRegion2DNode::bdRegion2DNode()
{
	m_node_of_column_node = NULL;
	m_node_of_row_node = NULL;
}


bdRegion2DNode::~bdRegion2DNode()
{
	m_node_of_column_node = NULL;
	m_node_of_row_node = NULL;
}


int bdRegion2DNode::IsNonZero()
{
	if((!m_node_of_column_node) || (!m_node_of_row_node)) return 0;
	return 1;
}


int& bdRegion2DNode::GetElement()
{
	return m_node_of_column_node->GetElementPointer()->m_element;
}


int bdRegion2DNode::GetIndexRow()
{
	return m_node_of_row_node->GetElementPointer()->m_row_index;
}
	

int bdRegion2DNode::GetIndexColumn()
{
	return m_node_of_column_node->GetElementPointer()->m_column_index;
}


void bdRegion2DNode::MoveToNext()
{
	m_node_of_column_node = m_node_of_column_node->GetRight();
	if(!m_node_of_column_node)
	{
		m_node_of_row_node = m_node_of_row_node->GetRight();
		if(m_node_of_row_node) m_node_of_column_node = m_node_of_row_node->GetElementPointer()->m_column_list.GetLeftEndNodePointer();
	}
}




//-------------------------------------------------------------------------------------------------------------------------------------------------------------





bdRegion2DIterator::bdRegion2DIterator()
{
	m_node_of_column_node = NULL;
	m_node_of_row_node = NULL;
}


bdRegion2DIterator::~bdRegion2DIterator()
{
	m_node_of_column_node = NULL;
	m_node_of_row_node = NULL;
}


void bdRegion2DIterator::SetBegin(bdRegion2D *region)
{
	if(!region)
	{
		m_node_of_column_node = NULL;
		m_node_of_row_node = NULL;
		return;
	}
	if(region->IsEmpty())
	{
		m_node_of_column_node = NULL;
		m_node_of_row_node = NULL;
		return;
	}
	m_node_of_row_node = region->GetRowList()->GetLeftEndNodePointer();
	m_node_of_column_node = m_node_of_row_node->GetElementPointer()->m_column_list.GetLeftEndNodePointer();
}


void bdRegion2DIterator::SetBegin(bdRegion2D *region, unsigned int start_r_index)
{
	if(!region)
	{
		m_node_of_column_node = NULL;
		m_node_of_row_node = NULL;
		return;
	}
	bdListIterator<bdRegionRowNode> it;
	for(it.SetLeftEnd(region->GetRowList()); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_row_index==start_r_index)
		{
			m_node_of_row_node = it.GetNodeAddress();
			m_node_of_column_node = m_node_of_row_node->GetElementPointer()->m_column_list.GetLeftEndNodePointer();
			return;
		}
	}
	//If gets here, no start_r_index was found, set pointers to NULLs.
	m_node_of_row_node = NULL;
	m_node_of_column_node = NULL;
}


//template<class T>
//void bdRegion2DIterator<T>::SetBegin(bdRegion2D<T> &region)
//{
//	if(region.IsEmpty())
//	{
//		m_node_of_column_node = NULL;
//		m_node_of_row_node = NULL;
//		return;
//	}
//	m_node_of_row_node = region.m_row_list.GetLeftEndPointerToElement();
//	m_node_of_column_node = m_node_of_row_node->m_column_list.GetLeftEndPointerToElement();
//}


int bdRegion2DIterator::IsValid()
{
	if(m_node_of_column_node && m_node_of_row_node) return 1;
	return 0;
}








//-------------------------------------------------------------------------------------------------------------------------------------------------------------





bdRegion2D::bdRegion2D()
{
	m_original_image_size_CR[1] = m_original_image_size_CR[0] = 0;
	m_region_value = 0;
	m_number_of_elements = 0;
}


bdRegion2D::~bdRegion2D()
{
	this->Reset();
}


void bdRegion2D::Reset()
{
	m_row_list.Reset();
	m_number_of_elements = 0;
	m_region_value = 0;
	m_original_image_size_CR[1] = m_original_image_size_CR[0] = 0;
}


void bdRegion2D::RemoveAllPoints()
{
	m_row_list.Reset();
	m_number_of_elements = 0;
}


int bdRegion2D::IsEmpty()
{
	if(m_number_of_elements!=0) return 0;
	return 1;
}


bdListNode<bdRegionRowNode>* bdRegion2D::GetNodeOfRowNode(unsigned int r_index)
{
	bdListIterator<bdRegionRowNode> it;
	for(it.SetLeftEnd(m_row_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_row_index == r_index)
			return it.GetNodeAddress();
		if(it.GetElementPointer()->m_row_index > r_index) return NULL;
	}
	return NULL;
}


bdRegion2DNode bdRegion2D::GetNode(unsigned int r_index, unsigned int c_index)
{
	bdRegion2DNode node;
	node.m_node_of_row_node = NULL;
	node.m_node_of_column_node = NULL;

	//If empty, local index not found
	if(this->IsEmpty()) return node;

	bdListNode<bdRegionRowNode> *row_node = this->GetNodeOfRowNode(r_index);
	if(!row_node) return node;

	bdListIterator<bdRegionColumnNode> it;
	for(it.SetLeftEnd(row_node->GetElementPointer()->m_column_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_column_index == c_index)
		{
			node.m_node_of_row_node = row_node;
			node.m_node_of_column_node = it.GetNodeAddress();
			return node;
		}
	}
	return node;
}


bdRegion2DNode bdRegion2D::GetNodeForFirstPoint()
{
	bdRegion2DNode node;

	if(this->IsEmpty())
	{
		node.m_node_of_row_node = NULL;
		node.m_node_of_column_node = NULL;
	}
	else
	{
		node.m_node_of_row_node = m_row_list.GetLeftEndNodePointer();
		node.m_node_of_column_node = m_row_list.GetLeftEndPointerToElement()->m_column_list.GetLeftEndNodePointer();
	}
	return node;
}


int bdRegion2D::RemovePoint(unsigned int r_index, unsigned int c_index)
{
	bdRegion2DNode node = this->GetNode(r_index, c_index);
	return this->RemovePoint(node);
}


int bdRegion2D::RemovePoint(bdRegion2DNode &node)
{
	if(this->IsEmpty()) return 1;

	if(!node.IsNonZero()) return 0;

	node.m_node_of_row_node->GetElementPointer()->m_column_list.DeleteNode(node.m_node_of_column_node);
	m_number_of_elements--;

	//If a column list had only one element and is now empty, remove it from the row list
	if(node.m_node_of_row_node->GetElementPointer()->m_column_list.IsEmpty())
		{ this->m_row_list.DeleteNode(node.m_node_of_row_node); }

	return 1;
}


int bdRegion2D::AddPoint(unsigned int r, unsigned int c, int element)
{
	if(this->IsEmpty()) 
	{
		//cout<<"--> if(Empty())"<<endl;
		bdRegionColumnNode cn;
		cn.m_column_index = c;
		cn.m_element = element;
		bdRegionRowNode rn;
		rn.m_row_index = r;
		rn.m_column_list.AddToRightEnd(cn);
		m_row_list.AddToRightEnd(rn);
		m_number_of_elements++;
		return 1;
	}

	//If the region is not empty
	//cout<<"--> find pointer location for new coordinates"<<endl;

	// Find the node of row_node for which the index is the last one smaller or equal to 'r'.
	bdListIterator<bdRegionRowNode> it;
	bdListNode<bdRegionRowNode> *recorded_node_of_row_node = m_row_list.GetLeftEndNodePointer();
	for(it.SetLeftEnd(m_row_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_row_index<=r) recorded_node_of_row_node = it.GetNodeAddress();
		else { break; }
	}

	// There are 3 possible cases:

	//--- 1. The exact row index was found ---
	if(recorded_node_of_row_node->GetElementPointer()->m_row_index == r)
	{
		//... add a new element to the column list.
		// Find the node of column_node for which the index is the last one smaller or equal to 'c'.
		bdListIterator<bdRegionColumnNode> itc;
		bdListNode<bdRegionColumnNode> *recorded_node_of_column_node = recorded_node_of_row_node->GetElementPointer()->m_column_list.GetLeftEndNodePointer();
		for(itc.SetLeftEnd(recorded_node_of_row_node->GetElementPointer()->m_column_list); itc.IsValid(); itc.MoveRight())
		{
			if(itc.GetElementPointer()->m_column_index<=c) recorded_node_of_column_node = itc.GetNodeAddress();
			else { break; }
		}
		
		// There are 3 possible cases:

		//--- 1a. The exact column index was found ---
		if(recorded_node_of_column_node->GetElementPointer()->m_column_index == c)
		{
			//The index was already entered, just exit.
			return 1;
		}
		//--- 1b. The smaller column index was found ---
		if(recorded_node_of_column_node->GetElementPointer()->m_column_index < c)
		{
			bdRegionColumnNode cn;
			cn.m_column_index = c;
			cn.m_element = element;
			recorded_node_of_row_node->GetElementPointer()->m_column_list.AddToRightOfNode(recorded_node_of_column_node,cn);
			m_number_of_elements++;
			return 1;
		}
		//--- 1c. The higher column index was found ---
		if(recorded_node_of_column_node->GetElementPointer()->m_column_index > c)
		{
			bdRegionColumnNode cn;
			cn.m_column_index = c;
			cn.m_element = element;
			recorded_node_of_row_node->GetElementPointer()->m_column_list.AddToLeftOfNode(recorded_node_of_column_node,cn);
			m_number_of_elements++;
			return 1;
		}

	}
	//--- 2. Lower row index was found ---
	if(recorded_node_of_row_node->GetElementPointer()->m_row_index < r)
	{
		//... make a new column list with given element and add to right of the recorded_node_of_row_node.
		bdRegionColumnNode cn;
		cn.m_column_index = c;
		cn.m_element = element;
		bdRegionRowNode rn;
		rn.m_row_index = r;
		rn.m_column_list.AddToRightEnd(cn);
		m_row_list.AddToRightOfNode(recorded_node_of_row_node,rn); //m_row_list.AddToRightEnd(rn);
		m_number_of_elements++;
		return 1;

	}
	//--- 3. Higher row index was found ---
	if(recorded_node_of_row_node->GetElementPointer()->m_row_index > r)
	{
		//... make a new column list with given element and add to right of the recorded_node_of_row_node.
		bdRegionColumnNode cn;
		cn.m_column_index = c;
		cn.m_element = element;
		bdRegionRowNode rn;
		rn.m_row_index = r;
		rn.m_column_list.AddToRightEnd(cn);
		m_row_list.AddToLeftOfNode(recorded_node_of_row_node,rn); //m_row_list.AddToLeftEnd(rn);
		m_number_of_elements++;
		return 1;
	}

	return 1;//This means that adding succeeded
}


int bdRegion2D::AddPoint(unsigned int r, unsigned int c, int *p_element)
{
	if(!p_element) return 0;
	return (this->AddPoint(r,c,*p_element));
}


int bdRegion2D::AddPoint(unsigned int r, unsigned int c)
{
	int el = 0;
	return (this->AddPoint(r,c,el));
}


void bdRegion2D::AddRegion(bdRegion2D &region)
{
	if(this->IsEmpty()) 
	{
		(*this) = region;
		return;
	}

	if(region.IsEmpty()) return;

	bdRegion2DIterator it;
	for(it.SetBegin(&region); it.IsValid(); it.MoveToNext())
	{
		this->AddPoint(it.GetIndexRow(),it.GetIndexColumn(),it.GetElement());
	}
}


int bdRegion2D::IsPointInRegion(unsigned int r_index, unsigned int c_index)
{
	//If empty, local index not found
	if(this->IsEmpty()) return 0;

	bdListNode<bdRegionRowNode> *row_node = this->GetNodeOfRowNode(r_index);
	if(!row_node) return 0;

	bdListIterator<bdRegionColumnNode> it;
	for(it.SetLeftEnd(row_node->GetElementPointer()->m_column_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_column_index == c_index) { return 1; }
	}
	return 0;
}


int bdRegion2D::IsOverlappingWithRegion(bdRegion2D &region)
{
	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		if(region.IsPointInRegion(it.GetIndexRow(),it.GetIndexColumn())) { return 1; }
	}
	return 0;
}


int* bdRegion2D::GetValue(unsigned int r, unsigned int c)
{
	bdRegion2DNode node = this->GetNode(r,c);
	if(node.IsNonZero()) return (&(node.GetElement()));
	return NULL;
}


bdRegion2D& bdRegion2D::operator =(bdRegion2D &region)
{
	if (&region==this) return *this;

	this->Reset();

	m_row_list = region.m_row_list;
	m_region_value = region.m_region_value;
	m_number_of_elements = region.m_number_of_elements;
	m_original_image_size_CR[1] = region.m_original_image_size_CR[1];
	m_original_image_size_CR[0] = region.m_original_image_size_CR[0];

	return *this;
}


void bdRegion2D::SetSizeOfOriginalImage(unsigned int rows, unsigned int columns)
{
	m_original_image_size_CR[1] = rows;
	m_original_image_size_CR[0] = columns;
}


unsigned int bdRegion2D::GetNumberOfPointsInRow(unsigned int r_index)
{
	bdListNode<bdRegionRowNode> *node_of_row_node = this->GetNodeOfRowNode(r_index);
	if(!node_of_row_node) return 0;
	return (node_of_row_node->GetElementPointer()->m_column_list.GetNumberOfElements());
}


unsigned int bdRegion2D::GetNumberOfPointsInColumn(unsigned int c_index)
{
	unsigned int number_of_elements_in_column = 0;
	
	bdListIterator<bdRegionRowNode> itr;
	for(itr.SetLeftEnd(m_row_list); itr.IsValid(); itr.MoveRight())
	{
		bdListIterator<bdRegionColumnNode> itc;
		for(itc.SetLeftEnd(itr.GetElementPointer()->m_column_list); itc.IsValid(); itc.MoveRight())
		{
			if(itc.GetElementPointer()->m_column_index == c_index)
			{
				number_of_elements_in_column++;
				break;
			}
			if(itc.GetElementPointer()->m_column_index > c_index) break;
		}
	}

	return number_of_elements_in_column;
}




//template<class T>
//void bdRegion2D<T>::ConvertRegionToFramedImage(Image8U &output_framed_image, unsigned char region_pixel_value)
//{
//	int r_start_included, r_end_included, c_start_included, c_end_included;
//
//	this->FrameCoordinatesForRegion(&r_start_included, &r_end_included, &c_start_included, &c_end_included);
//
//	int R_frame = (r_end_included - r_start_included + 1) + 2;
//	int C_frame = (c_end_included - c_start_included + 1) + 2;
//
//	output_framed_image.Set(R_frame,C_frame);
//	output_framed_image.FillInWith(0);
//
//	output_framed_image.SetPositionIndexes(r_start_included-1, c_start_included-1);
//
//	int rn, cn;
//	for(this->For_TheWholeRegion(); this->GetAllPoints(&rn, &cn); )
//	{
//		output_framed_image(rn,cn) = region_pixel_value;
//	}
//}




int bdRegion2D::InsertValuesFromImage(bdImage16U &image, unsigned int s, unsigned int t)
{
	if(this->IsEmpty()) return 0;
	if(image.IsEmpty()) return 0;

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		unsigned int r = it.GetIndexRow();
		unsigned int c = it.GetIndexColumn();

		if(r<image.GetNumberOfRows() && c<image.GetNumberOfColumns()) { it.GetElement() = image(t,s,r,c); }
	}
	return 1;
}


int bdRegion2D::RidgePointsRegion(bdRegion2D &output_ridges)
{
	if(this->IsEmpty()) return 0;

	output_ridges.Reset();

	int max = -1;
	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		if(it.GetElement()>=max)
		{
			if(it.GetElement()>max) 
			{
				output_ridges.RemoveAllPoints();
				max = it.GetElement();
			}
			output_ridges.AddPoint(it.GetIndexRow(),it.GetIndexColumn(),it.GetElement());
		}
	}
	return 1;
}


int bdRegion2D::CreateRegionFromImage(bdImage16U &image, unsigned int s, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CR[1] = image.GetNumberOfRows();
	m_original_image_size_CR[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(255);

	for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
	{
		for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
		{
			if(image(t,s,r,c)!=0)
			{
				this->AddPoint(r,c,(int&)image(t,s,r,c));
			}
		}
	}
	return 1;
}


int bdRegion2D::CreateRegionFromImage(bdImage16U &image, unsigned short value_for_growing, unsigned int s, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CR[1] = image.GetNumberOfRows();
	m_original_image_size_CR[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(value_for_growing);

	for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
	{
		for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
		{
			if(image(t,s,r,c)==value_for_growing)
			{
				this->AddPoint(r,c,(int&)value_for_growing);
			}
		}
	}
	return 1;
}


int bdRegion2D::CreateRegion_8_FromSeedPoint(bdImage16U &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(image.IsEmpty()) return 0;
	if(c>=image.GetNumberOfColumns() || r>=image.GetNumberOfRows() || s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CR[1] = image.GetNumberOfRows();
	m_original_image_size_CR[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(image(t,s,r,c));

	//'temp_image' contains pixels that are in the region
	bdImage temp_image;
	temp_image.SetSize(1,1,image.GetNumberOfRows(),image.GetNumberOfColumns());
	temp_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(1, image.GetNumberOfRows(), image.GetNumberOfColumns());

	bdDiscreteCoordinates3D p(s,r,c);
	bdList<bdDiscreteCoordinates3D> list_of_points;
	list_of_points.AddToRightEnd(p);

	this->AddPoint(p.R(),p.C(),(int&)image(t,s,r,c));
	temp_image(0,0,p.R(),p.C()) = 1;

	int rn,cn;
	while(!list_of_points.IsEmpty())
	{
		p = list_of_points.GetLeftEnd();
		for(g.ForCoordinates_8_Neighbors(p.R(),p.C()); g.Get_8_Neighbors(rn,cn); )
		{
			if(image(t,s,rn,cn)!=0)
			{
				//If the point is not in the region
				if(temp_image(0,0,rn,cn)==0) 
				{
					list_of_points.AddToRightEnd(p(s,rn,cn));
					this->AddPoint(rn,cn,(int&)image(t,s,rn,cn));
					temp_image(0,0,rn,cn) = 1;
				}
			}
		}
		list_of_points.DeleteLeftEnd();
	}

	return 1;
}


int bdRegion2D::CreateRegion_8_FromSeedPoint(bdImage16U &image, unsigned short value_for_growing, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(image.IsEmpty()) return 0;
	if(c>=image.GetNumberOfColumns() || r>=image.GetNumberOfRows() || s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CR[1] = image.GetNumberOfRows();
	m_original_image_size_CR[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(image(t,s,r,c));

	//'temp_image' contains pixels that are in the region
	bdImage temp_image;
	temp_image.SetSize(1,1,image.GetNumberOfRows(),image.GetNumberOfColumns());
	temp_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(1, image.GetNumberOfRows(), image.GetNumberOfColumns());

	bdDiscreteCoordinates3D p(s,r,c);
	bdList<bdDiscreteCoordinates3D> list_of_points;
	list_of_points.AddToRightEnd(p);

	this->AddPoint(p.R(),p.C(),(int&)image(t,s,r,c));
	temp_image(0,0,p.R(),p.C()) = 1;

	int rn,cn;
	while(!list_of_points.IsEmpty())
	{
		p = list_of_points.GetLeftEnd();
		for(g.ForCoordinates_8_Neighbors(p.R(),p.C()); g.Get_8_Neighbors(rn,cn); )
		{
			if(image(t,s,rn,cn)==value_for_growing)
			{
				//If the point is not in the region
				if(temp_image(0,0,rn,cn)==0) 
				{
					list_of_points.AddToRightEnd(p(s,rn,cn));
					this->AddPoint(rn,cn,image(t,s,rn,cn)); //this->AddPoint(rn,cn,(int&)image(t,s,rn,cn));
					temp_image(0,0,rn,cn) = 1;
				}
			}
		}
		list_of_points.DeleteLeftEnd();
	}

	return 1;
}


int bdRegion2D::CreateRegion_8_FromSeedPointAndValueHigherThan(bdImage16U &image, unsigned short value_to_grow_above, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(image.IsEmpty()) return 0;
	if(c>=image.GetNumberOfColumns() || r>=image.GetNumberOfRows() || s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CR[1] = image.GetNumberOfRows();
	m_original_image_size_CR[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(image(t,s,r,c));

	//'temp_image' contains pixels that are in the region
	bdImage temp_image;
	temp_image.SetSize(1,1,image.GetNumberOfRows(),image.GetNumberOfColumns());
	temp_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(1, image.GetNumberOfRows(), image.GetNumberOfColumns());

	bdDiscreteCoordinates3D p(s,r,c);
	bdList<bdDiscreteCoordinates3D> list_of_points;
	list_of_points.AddToRightEnd(p);

	this->AddPoint(p.R(),p.C(),(int&)image(t,s,r,c));
	temp_image(0,0,p.R(),p.C()) = 1;

	int rn,cn;
	while(!list_of_points.IsEmpty())
	{
		p = list_of_points.GetLeftEnd();
		for(g.ForCoordinates_8_Neighbors(p.R(),p.C()); g.Get_8_Neighbors(rn,cn); )
		{
			if(image(t,s,rn,cn)>value_to_grow_above)
			{
				//If the point is not in the region
				if(temp_image(0,0,rn,cn)==0) 
				{
					list_of_points.AddToRightEnd(p(s,rn,cn));
					this->AddPoint(rn,cn,(int&)image(t,s,rn,cn));
					temp_image(0,0,rn,cn) = 1;
				}
			}
		}
		list_of_points.DeleteLeftEnd();
	}

	return 1;
}


int bdRegion2D::CenterOfMass(double &center_of_mass_R, double &center_of_mass_C)
{
	if(this->IsEmpty()) return 0;

	center_of_mass_R = center_of_mass_C = 0;

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		center_of_mass_R += it.GetIndexRow();
		center_of_mass_C += it.GetIndexColumn();
	}

	center_of_mass_R = center_of_mass_R / ((double) this->GetNumberOfElements());
	center_of_mass_C = center_of_mass_C / ((double) this->GetNumberOfElements());

	return 1;
}


double bdRegion2D::Dice(bdRegion2D &region)
{
	unsigned int overlap = this->AreaOfOverlapWith(region);
	return ( (((double)overlap)*2.0)  /  ((double)(this->GetNumberOfElements()+region.GetNumberOfElements())) );
}


int bdRegion2D::GeometricMedian(double &geometric_median_R, double &geometric_median_C)
{
	double 	center_of_mass_R, center_of_mass_C;
	if(!this->CenterOfMass(center_of_mass_R, center_of_mass_C)) return 0;

	bdRegion2DIterator it;
	double min_distance = 65535;//Set this to a unlikely high value.
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int r =  it.GetIndexRow();
		int c =  it.GetIndexColumn();
		double d = (r-center_of_mass_R)*(r-center_of_mass_R) + (c-center_of_mass_C)*(c-center_of_mass_C);
		if(d<min_distance)
		{
			geometric_median_R = r; geometric_median_C = c;
			min_distance = d;
		}
	}

	return 1;
}


int bdRegion2D::PasteRegionToExistingImage(bdImage &img, unsigned short region_value, unsigned int t, unsigned int s)
{
	if(img.IsEmpty()) return 0;

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		img(t,s,it.GetIndexRow(),it.GetIndexColumn()) = region_value;
	}

	return 1;
}


int bdRegion2D::ListOfOverlappingRegionsInImageSlice(bdImage &img, bdList<bdRegion2D> &output_list, unsigned int t, unsigned int s)
{
	if(img.IsEmpty()) return 0;
	if(this->IsEmpty()) return 0;

	output_list.Reset();

	bdImage temp;
	temp.SetSize(1,1,img.GetNumberOfRows(),img.GetNumberOfColumns());
	temp.FillInWith(0);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int r = it.GetIndexRow();
		int c = it.GetIndexColumn();

		if(temp(r,c)==0 && img(t,s,r,c)!=0)
		{
			bdRegion2D *reg = output_list.AddNewToRightEnd();
			reg->CreateRegion_8_FromSeedPoint(img,t,s,r,c);
			reg->PasteRegionToExistingImage(temp,255);
		}
	}

	return 1;
}


int bdRegion2D::ListOfOverlappingRegionsInImageSlice_SameIntensityForRegion(bdImage &img, bdList<bdRegion2D> &output_list, unsigned int t, unsigned int s)
{
	if(img.IsEmpty()) return 0;
	if(this->IsEmpty()) return 0;

	output_list.Reset();

	bdImage temp;
	temp.SetSize(1,1,img.GetNumberOfRows(),img.GetNumberOfColumns());
	temp.FillInWith(0);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int r = it.GetIndexRow();
		int c = it.GetIndexColumn();

		if(temp(r,c)==0 && img(t,s,r,c)!=0)
		{
			bdRegion2D *reg = output_list.AddNewToRightEnd();
			reg->CreateRegion_8_FromSeedPoint(img,img(t,s,r,c),t,s,r,c);
			reg->PasteRegionToExistingImage(temp,255);
		}
	}

	return 1;
}


unsigned int bdRegion2D::AreaOfOverlapWith(bdRegion2D &region)
{
	if(this->IsEmpty()) return 0;
	if(region.IsEmpty()) return 0;

	bdRegion2D *reg1, *reg2;
	if(this->GetNumberOfElements()<region.GetNumberOfElements()) { reg1 = this; reg2 = &region; }
	else { reg1 = &region; reg2 = this; }

	unsigned int overlap_area = 0;
	bdRegion2DIterator it;
	for(it.SetBegin(reg1); it.IsValid(); it.MoveToNext())
	{
		if(reg2->IsPointInRegion(it.GetIndexRow(),it.GetIndexColumn())) overlap_area++;
	}

	return overlap_area;
}




//
//template <class T>
//void bdRegion2D<T>::InnerRadiusCalculation_Approximate(int *inner_radius, Point &inner_radius_center_point)
//{
//	bdList<Point> l1, l2, l3, *pl1, *pl2, *pl3;
//	Image8U framed_image;
//
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	Geometry g;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//
//	Point p, center_of_mass_point;
//	center_of_mass_point(0,0);
//
//	//Find the edge points
//	int rn,cn;
//	for(int r=0; r<framed_image.NumberOfRows(); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns(); c++)
//		{
//			if(framed_image[r][c]==255)
//			{
//				center_of_mass_point.r = center_of_mass_point.r + r;
//				center_of_mass_point.c = center_of_mass_point.c + c;
//				for(g.ForPoint(r,c); g.Get_4_ClosestPositions(&rn,&cn); )
//				{
//					if(framed_image[rn][cn]==0)
//					{
//						framed_image[r][c] = 1;
//						l1.AddToRightEnd(p(r,c));
//						g.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//	}
//	//At this point all the edges are in list l1 and are marked with value 1 in framed_image
//
//	//Center of mass point calculation
//	center_of_mass_point.r = center_of_mass_point.r / (this->NumberOfElements());
//	center_of_mass_point.c = center_of_mass_point.c / (this->NumberOfElements());
//
//	pl1 = &l1;
//	pl2 = &l2;
//	pl3 = &l3;
//	(*pl2) = (*pl1);
//	int counter_for_number_of_erosions = 1;
//
//	//Find the list of inner radius center point candidates - store it to pl2 list
//	while(!pl1->Empty())
//	{
//		counter_for_number_of_erosions++;
//
//		while(!pl1->Empty())
//		{
//			p = pl1->GetLeftEnd();
//			for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//			{
//				if(framed_image[rn][cn]==255)
//				{
//					framed_image[rn][cn] = 1;
//					//framed_image[rn][cn] = counter_for_number_of_erosions;
//					pl3->AddToRightEnd(p(rn,cn));
//				}
//			}
//			pl1->DeleteLeftEnd();
//		}
//		//At this point list *pl1 is empty and all the next edge points are in *pl3
//
//		//If *pl3 is not empty, store it to *pl1. List *pl2 should be emptied and made the same as *pl3 
//		if(!pl3->Empty())
//		{
//			//Make pl2 list the same as pl3...
//			//pl2->Delete();
//			(*pl2) = (*pl3);
//			//... and switch pl1 and pl3
//			bdList<Point> *pltemp;
//			pltemp = pl1;
//			pl1 = pl3;
//			pl3 = pltemp;
//			//pltemp = NULL;
//		}
//		//If *pl3 is empty it means that the previous edge was the last one. It is still stored in *pl2. Loop exits.
//	}
//	//At this point *pl2 has all the candidate points for the inner circle center point. 
//	//The inner circle radius is determined by growing the edge stored in *pl2 till it gets to a zero pixel in the 
//	// original region. The inner radius value is the distance from the center to that zero pixel. Center point is 
//	// determined as a center of mass for the points in list *pl2.
//
//	//Calculate the center point of inner radius (and put value 2 in the framed_image for these points)
//	Point center_point;
//	center_point(0,0);
//	for(pl2->For_StartFromLeftEnd(); pl2->GetAllElemetsByMovingToRight(&p); )
//	{
//		center_point(center_point.r+p.r,center_point.c+p.c);
//		framed_image[p.r][p.c] = 2;
//	}
//	center_point((center_point.r/(pl2->NumberOfElements())), (center_point.c/(pl2->NumberOfElements())));
//
//	//This is the center point (actually a mass center point) calculated for the remaining edge of the region. The 
//	// inner_radius_center_point is the point in the pl2 list closest to center_point.
//	
//	//Determine the inner_radius_center_point
//	double distance = -1;
//	double d;
//	for(pl2->For_StartFromLeftEnd(); pl2->GetAllElemetsByMovingToRight(&p); )
//	{
//		d = center_point.Distance_Euclidean(p);
//		if(distance<0 || d<distance)
//		{
//			distance = d;
//			inner_radius_center_point = p;
//		}
//	}
//	
//	//Now that the inner_radius_center_point is determined, grow region *pl2 and stop at first non-zero pixels
//	int inner_radius_determined = 0;
//	pl1->Delete();
//	pl3->Delete();
//	while(!inner_radius_determined)
//	{
//		while(!pl2->Empty())
//		{
//			p = pl2->GetLeftEnd();
//
//			for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//			{
//				//This is the new edge point, store to pl1
//				if(framed_image[rn][cn]==1)
//				{
//					framed_image[rn][cn] = 2;
//					pl1->AddToRightEnd(p(rn,cn));
//				}
//				else
//				{
//					//If the zero pixel was found
//					if(framed_image[rn][cn]==0)
//					{
//						//Store into the list of zero pixels
//						pl3->AddToRightEnd(p(rn,cn));
//						inner_radius_determined = 1;
//					}
//				}
//			}
//			pl2->Delete();
//		}
//		//If the radius was not determined give new points to pl2
//		if(!inner_radius_determined)
//		{
//			bdList<Point> *pltemp;
//			pltemp = pl1;
//			pl1 = pl2;
//			pl2 = pltemp;
//		}
//	}
//
//	//At this point all the zero-pixels that were first encountered while growing the region are in pl3
//	// We have to determine the shortest distance of all these points from inner_radius_determined
//
//	//Determine the inner_radius_center_point
//	distance = -1;
//	for(pl3->For_StartFromLeftEnd(); pl3->GetAllElemetsByMovingToRight(&p); )
//	{
//		d = inner_radius_center_point.Distance_Euclidean(p);
//		if(distance<0 || d<distance)
//		{
//			cout<<p<<endl;
//			distance = d;
//		}
//	}
//	*inner_radius = (int) distance;
//
//	//Convert the inner radius center position to original positions
//	framed_image.ConvertIndexToPosition(inner_radius_center_point.r,inner_radius_center_point.c,&inner_radius_center_point.r,&inner_radius_center_point.c);
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::OuterRadiusCalculation_Approximate(int *outer_radius, Point &outer_radius_center_point)
//{
//	Image8U framed_image;
//
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	Geometry g;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//	int rn,cn;
//
//	Point p, center_of_mass_point;
//	center_of_mass_point(0,0);
//
//	//Calculate the center of mass point 
//	for(int r=0; r<framed_image.NumberOfRows(); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns(); c++)
//		{
//			if(framed_image[r][c]==255)
//			{
//				center_of_mass_point.r = center_of_mass_point.r + r;
//				center_of_mass_point.c = center_of_mass_point.c + c;
//			}
//		}
//	}
//	center_of_mass_point.r = center_of_mass_point.r / (this->NumberOfElements());
//	center_of_mass_point.c = center_of_mass_point.c / (this->NumberOfElements());
//
//	Array<Point> array_of_intersection_points;
//	array_of_intersection_points.Set(16);
//
//	//For 8 directions we have to check if the row index of points is smaller or greater than the initial value
//
//	Point first, second;//these are the two intersection points
//
//	//First 4 directions
//	Array<Point> start_points, end_points;
//	start_points.Set(4);
//	end_points.Set(4);
//	start_points[0].r = 0; start_points[0].c = 0; 
//	end_points[0].r = framed_image.NumberOfRows()-1; end_points[0].c = framed_image.NumberOfColumns()-1; 
//	
//	start_points[1].r = 0; start_points[1].c = framed_image.NumberOfColumns()/4; 
//	end_points[1].r = framed_image.NumberOfRows()-1; end_points[1].c = (framed_image.NumberOfColumns()*3)/4; 
//	
//	start_points[2].r = 0; start_points[2].c = framed_image.NumberOfColumns()/2; 
//	end_points[2].r = framed_image.NumberOfRows()-1; end_points[2].c = framed_image.NumberOfColumns()/2; 
//	
//	start_points[3].r = 0; start_points[3].c = (framed_image.NumberOfColumns()*3)/4; 
//	end_points[3].r = framed_image.NumberOfRows()-1; end_points[3].c = framed_image.NumberOfColumns()/4; 
//	for(int i=0; i<4; i++)
//	{
//		first((framed_image.NumberOfRows()/2),(framed_image.NumberOfColumns()/2));
//		second((framed_image.NumberOfRows()/2),(framed_image.NumberOfColumns()/2));
//		for(g.ForPoint(start_points[i].r,start_points[i].c); g.GetPositionsInLineSegment(end_points[i].r,end_points[i].c, &rn, &cn); )
//		{
//			//If the point is part of the region
//			if(framed_image[rn][cn]==255)
//			{
//				//If the row value is smaller than the row value of the first point or greater than the row index of the second point
//				if(rn<first.r) first(rn,cn);
//				if(rn>second.r) second(rn,cn);
//			}			
//		}
//		array_of_intersection_points[2*i] = first;
//		array_of_intersection_points[2*i+1] = second;
//	}
//
//	//Second 4 directions - compare columns
//	start_points[0].r = framed_image.NumberOfRows()/4; start_points[0].c = 0; 
//	end_points[0].r = (framed_image.NumberOfRows()*3)/4; end_points[0].c = framed_image.NumberOfColumns()-1; 
//	
//	start_points[1].r = framed_image.NumberOfRows()/2; start_points[1].c = 0;
//	end_points[1].r = framed_image.NumberOfRows()/2; end_points[1].c = framed_image.NumberOfColumns()-1; 
//	
//	start_points[2].r = (framed_image.NumberOfRows()*3)/4; start_points[2].c = 0;
//	end_points[2].r = framed_image.NumberOfRows()/4; end_points[2].c = framed_image.NumberOfColumns()-1; 
//	
//	start_points[3].r = framed_image.NumberOfRows()-1; start_points[3].c = 0;
//	end_points[3].r = 0; end_points[3].c = framed_image.NumberOfColumns()-1; 
//	for(int i=0; i<4; i++)
//	{
//		first((framed_image.NumberOfRows()/2),(framed_image.NumberOfColumns()/2));
//		second((framed_image.NumberOfRows()/2),(framed_image.NumberOfColumns()/2));
//		for(g.ForPoint(start_points[i].r,start_points[i].c); g.GetPositionsInLineSegment(end_points[i].r,end_points[i].c, &rn, &cn); )
//		{
//			//If the point is part of the region
//			if(framed_image[rn][cn]==255)
//			{
//				//If the row value is smaller than the row value of the first point or greater than the row index of the second point
//				if(cn<first.c) first(rn,cn);
//				if(cn>second.c) second(rn,cn);
//			}
//		}
//		array_of_intersection_points[8+2*i] = first;
//		array_of_intersection_points[8+2*i+1] = second;
//	}
//
//	//At this point we have 16 points. Find the biggest distance between these.
//	int distance = 0;
//	int d;	
//	Point distance_point1, distance_point2;
//	for(int i=0; i<array_of_intersection_points.NumberOfElements()-1; i++)
//	{
//		for(int j=i+1; j<array_of_intersection_points.NumberOfElements(); j++)
//		{			
//			d = (int)(array_of_intersection_points[i].Distance_Euclidean(array_of_intersection_points[j]));			
//			if(d>distance) 
//			{
//				distance = d;
//				distance_point1 = array_of_intersection_points[i];
//				distance_point2 = array_of_intersection_points[j];				
//			}
//		}
//	}
//	*outer_radius = distance/2; 
//
//	//Find the point on the line connecting largest distance points that is the closest to the center of mass point (it does no need
//	// to be in the region)
//	d = 0;
//	distance = -1;
//	int r_pos,c_pos;
//	for(g.ForPoint(distance_point1.r,distance_point1.c); g.GetPositionsInLineSegment(distance_point2.r,distance_point2.c, &rn, &cn); )
//	{
//		d = center_of_mass_point.Distance_Euclidean(p(rn,cn));
//		if(distance==-1 || d<distance)
//		{
//			distance = d;
//			framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//			outer_radius_center_point(r_pos,c_pos);
//		}
//	}
//	//framed_image.SaveToFile("framed_intersection_points.png");
//	//cout<<distance<<" "<<distance_point1<<" "<<distance_point2<<endl;
//}
//
//
//
//
//template<class T>
//void bdRegion2D<T>::OuterRadiusCalculation(int *outer_radius, Point &outer_radius_center_point)
//{
//	bdList<Point> l;
//	this->ListOfOuterInternal_4_EdgePoints(l);
//	//Copy the list into the array
//	Array<Point> edge_array;
//	edge_array.Set(l.NumberOfElements());
//	Point p, center_of_mass_point;
//	int i = 0;
//	for(l.For_StartFromLeftEnd(); l.GetAllElemetsByMovingToRight(&p); )
//	{
//		edge_array[i] = p;
//		i++;
//	}
//	this->CenterOfMassPointCalculation(center_of_mass_point);
//
//	//Find the biggest distance between the edge points.
//	int distance = 0;
//	int d;	
//	Point distance_point1, distance_point2;
//	for(int i=0; i<edge_array.NumberOfElements()-1; i++)
//	{
//		for(int j=i+1; j<edge_array.NumberOfElements(); j++)
//		{			
//			d = (int)(edge_array[i].Distance_Euclidean(edge_array[j]));			
//			if(d>distance) 
//			{
//				distance = d;
//				distance_point1 = edge_array[i];
//				distance_point2 = edge_array[j];				
//			}
//		}
//	}
//	*outer_radius = distance/2; 
//
//	//Find the point on the line connecting largest distance points that is the closest to the center of mass point (it does no need
//	// to be in the region)
//	d = 0;
//	distance = -1;
//	int rn,cn;
//	Geometry g;
//	g.SetDimensions(0,this->NumberOfRowsOfOriginalImage(),this->NumberOfColumnsOfOriginalImage());
//	//int r_pos,c_pos;
//	for(g.ForPoint(distance_point1.r,distance_point1.c); g.GetPositionsInLineSegment(distance_point2.r,distance_point2.c, &rn, &cn); )
//	{
//		d = center_of_mass_point.Distance_Euclidean(p(rn,cn));
//		if(distance==-1 || d<distance)
//		{
//			distance = d;
//			//framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//			outer_radius_center_point(rn,cn);
//		}
//	}
//
//}
//
//
//
//
//
//template<class T>
//void bdRegion2D<T>::CircularityCheck(double *l1, double *l2, double *min, double *max, double *average, Position &center_of_mass)
//{
//	assert(!Empty());
//
//	//Get the list of outer internal edge points
//	bdList<Point> l;
//	this->ListOfOuterInternal_4_EdgePoints(l);
//
//	//Get the center of mass position
//	this->CenterOfMassPositionCalculation(center_of_mass);
//	
//	//Calculate the average distance value
//	Point p;
//	*average = 0;
//	Array<double> distance_array;
//	distance_array.Set(l.NumberOfElements());
//	int i=0;
//	for(l.For_StartFromLeftEnd(); l.GetAllElemetsByMovingToRight(&p); )
//	{
//		//Calculate the distance
//		distance_array[i] =(((double)(p.r))-center_of_mass.r)*(((double)(p.r))-center_of_mass.r)+(((double)(p.c))-center_of_mass.c)*(((double)(p.c))-center_of_mass.c);
//		*average += distance_array[i];
//		i++;
//	}
//	*average = (*average)/ ((double)(l.NumberOfElements()));
//	
//	//Calculate min, max, l1 and l2 norms.
//	*max = 0;
//	*min = -1;//set to irregular value
//	*l1 = 0;
//	*l2 = 0;
//	//for(l.For_StartFromLeftEnd(); l.GetAllElemetsByMovingToRight(&p); )
//	for(i=0; i<distance_array.NumberOfElements(); i++)
//	{
//		//Calculate the distance
//		//double d =(((double)(p.r))-center_of_mass.r)*(((double)(p.r))-center_of_mass.r)+(((double)(p.c))-center_of_mass.c)*(((double)(p.c))-center_of_mass.c);
//		
//		double diff = distance_array[i] - (*average);
//		if(diff<0) diff = -diff;
//
//		if((*max)<diff) *max = diff;
//		if((*min)>diff || (*min)==-1) (*min) = diff;
//
//		*l1 += diff;
//		*l2 += diff*diff;
//	}
//	*l1 = (*l1)/ ((double)(l.NumberOfElements()));
//	*l2 = (*l2)/ ((double)(l.NumberOfElements()));
//}
//
//
//
//
//template<class T>
//int bdRegion2D<T>::CircularityCheckWithLargeObjectPriority()
//{
//	assert(!Empty());
//
//	Point p;
//	this->CenterOfMassPointCalculation(p);
//
//	//If the center of mass point does not belong to the region, the area of the inscibed circle is 0:
//	if(!this->PointIsInTheRegion(p.r,p.c))
//	{
//		return (-(NumberOfElements()));
//	}
//
//	Geometry g;
//	g.SetDimensions(0, this->m_original_image_size_CR[1], this->m_original_image_size_CR[0]);
//
//	//Grow the circle and check if all the pixels belong to the region (growing inscribed circle)
//	int is_edge_found = 0;
//	int rn, cn;
//	int inscribed_circle_area = 1;//starts with 1 because of the center point (center of mass point belongs to the region)
//	for(int radius = 0; !is_edge_found; radius++)
//	{
//		int area_of_current_radius_ring = 0;
//		for(g.ForPoint(p.r,p.c); g.GetPositionsOfRadius(radius, radius+1, &rn, &cn); )
//		{
//			//If the point is not in the region, stop searching
//			if(!PointIsInTheRegion(rn,cn))
//			{
//				g.ExitLoopAtNextFunctionCall();
//				is_edge_found = 1;
//			}
//			else
//			{
//				area_of_current_radius_ring++;
//			}
//		}
//
//		if(!is_edge_found)
//		{
//			inscribed_circle_area += area_of_current_radius_ring;
//		}
//	}
//
//	//Calculate circularity coefficient
//	//return((2*inscribed_circle_area) - NumberOfElements());
//	return((3*inscribed_circle_area) - 2*NumberOfElements());
//}
//
//
//
//template<class T>
//void bdRegion2D<T>::InnerAndOuterRadiusCalculation(int *inner_radius, int *outer_radius, Point &inner_radius_center_point)
//{
//	assert(!Empty());
//
//	Geometry g;
//	g.SetDimensions(0,m_original_image_size_CR[1],m_original_image_size_CR[0]);
//
//	int r,c,rn,cn;
//
//	int max_inner_radius = 0;
//	int min_outer_radius = NumberOfElements();//high value for initialization
//
//
//	Image8U temp_inner_radius_image;
//	temp_inner_radius_image.Set(m_original_image_size_CR[1],m_original_image_size_CR[0]);
//	temp_inner_radius_image.FillInWith(0);
//
//	
//
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		int inner_radius_found = 0;
//		int outer_radius_found = 0;
//		int radius = 0;		
//
//		while(!(inner_radius_found && outer_radius_found))
//		{
//			radius++;
//			outer_radius_found = 1;
//
//			for(g.ForPoint(r,c); g.GetPositionsOfRadius(radius-1,radius,&rn,&cn);)
//			{
//				if(!inner_radius_found)
//				{
//					if(!PointIsInTheRegion(rn,cn))
//					{
//						inner_radius_found = 1;
//						if(radius>max_inner_radius) max_inner_radius = radius;
//
//						temp_inner_radius_image[r][c] = radius;
//
//						//*inner_radius = radius;
//					}
//				}
//
//				if(outer_radius_found)
//				{
//					if(PointIsInTheRegion(rn,cn))
//					{
//						outer_radius_found = 0;
//					}
//				}
//			}
//		}
//
//		if(radius<min_outer_radius) min_outer_radius = radius;
//	}
//
//	*inner_radius = max_inner_radius;
//	*outer_radius = min_outer_radius;
//
//	//Find the center - it is the point with value of the maximum radius which is surrounded with positions
//	// that have the highest radius values (we sum these values and obtain the center position)
//
//	//List of points that have maximum radius value
//	bdList<Point> list_of_maximum_radius_points;
//
//	Point point;
//
//	//Fill in the list of points that have maximum radius value
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		if(temp_inner_radius_image[r][c]==max_inner_radius) list_of_maximum_radius_points.AddToRightEnd(point(0,r,c));
//	}
//
//	//If there is only one element in the list, no need for further checks
//	if(list_of_maximum_radius_points.NumberOfElements()==1)
//	{
//		inner_radius_center_point = list_of_maximum_radius_points.GetLeftEnd();
//		return;
//	}
//
//	//Center point is the one that has max value and is surrounded by highest values of radiuses in 8-neighborhood
//	unsigned int radius_sum, max_radius_sum = 0;
//	for(list_of_maximum_radius_points.For_StartFromLeftEnd(); list_of_maximum_radius_points.GetAllElemetsByMovingToRight(&point); )
//	{
//		radius_sum = 0;
//		for(g.ForPoint(point.r,point.c); g.GetPositionsInCircleOfSquareRadius(max_inner_radius*max_inner_radius,&rn,&cn);)
//		{
//			radius_sum += (unsigned int)(temp_inner_radius_image[rn][cn]);
//		}
//		if(max_radius_sum<radius_sum)
//		{
//			max_radius_sum = radius_sum;
//			inner_radius_center_point = point;
//		}
//	}
//
//}
//
//
//
//
//template<class T>
//void bdRegion2D<T>::InnerAndOuterRadiusCalculation_fast(int *inner_radius, int *outer_radius, Point &inner_radius_center_point, int area)
//{
//	assert(!Empty());
//
//	bdRegion2D<T> region;
//	region = (*this);
//
//
//	Geometry g;
//	g.SetDimensions(0,m_original_image_size_CR[1],m_original_image_size_CR[0]);
//
//
//	//If number of elements is larger than desired, shrink the region and then find the inner radius
//	if(NumberOfElements()>area)
//	{
//		bdList<Point> list1, list2, *p_list1, *p_list2, *p_temp;
//		//region.ListOfEdgePoints(list1);
//		region.ListOfInternal_8_EdgePoints(list1);//PROVERI DA LI JE ISPRAVNA FUNCKIJA POZVANA!!!!!!!!!!!!!!!!!!!!!!!!!!!
//		p_list1 = &list1;
//		p_list2 = &list2;
//		
//		Point point, p;
//		int r, c;
//
//		//Repeat untill the region is not of desired area value
//		while(region.NumberOfElements()>area)
//		{
//			//cout<<"region area = "<<region.NumberOfElements();
//
//			p_list2->DeleteAll();
//
//			//For each edge point 
//			for(p_list1->For_StartFromLeftEnd(); p_list1->GetAllElemetsByMovingToRight(&point); )
//			{
//				if(region.PointIsInTheRegion(point.r,point.c))
//				{
//					region.DeletePoint(point.r,point.c);
//				}
//			}
//
//			//For each edge point 
//			for(p_list1->For_StartFromLeftEnd(); p_list1->GetAllElemetsByMovingToRight(&point); )
//			{
//				//Check all its 4-neighbors
//				for(g.ForPoint(point.r,point.c); g.Get_4_ClosestPositions(&r,&c); )
//				{
//					if(region.PointIsInTheRegion(r,c)) p_list2->AddToRightEnd(p(0,r,c));
//				}
//			}
//
//			p_list1->DeleteAll();
//			p_temp = p_list1;
//			p_list1 = p_list2;
//			p_list2 = p_temp;
//		}
//
//		//If we have the region od adecquate area, empty the last edge list (pointer 1)
//		p_list1->DeleteAll();
//	}
//
//	//After this part we have a region of desired area for which we calculate radiuses based on the original region
//
//	//cout<<"stigao do racunanja"<<endl;
//	int r,c,rn,cn;
//
//	int max_inner_radius = 0;
//	int min_outer_radius = NumberOfElements();//high value for initialization
//
//
//	Image8U temp_inner_radius_image;
//	temp_inner_radius_image.Set(m_original_image_size_CR[1],m_original_image_size_CR[0]);
//	temp_inner_radius_image.FillInWith(0);
//
//	
//
//	for(region.For_TheWholeRegion(); region.GetAllPoints(&r,&c); )
//	{
//		int inner_radius_found = 0;
//		int outer_radius_found = 0;
//		int radius = 0;		
//
//		while(!(inner_radius_found && outer_radius_found))
//		{
//			radius++;
//			outer_radius_found = 1;
//
//			for(g.ForPoint(r,c); g.GetPositionsOfRadius(radius-1,radius,&rn,&cn);)
//			{
//				if(!inner_radius_found)
//				{
//					if(!PointIsInTheRegion(rn,cn))
//					{
//						inner_radius_found = 1;
//						if(radius>max_inner_radius) max_inner_radius = radius;
//
//						temp_inner_radius_image[r][c] = radius;
//
//						//*inner_radius = radius;
//					}
//				}
//
//				if(outer_radius_found)
//				{
//					if(PointIsInTheRegion(rn,cn))
//					{
//						outer_radius_found = 0;
//					}
//				}
//			}
//		}
//
//		if(radius<min_outer_radius) min_outer_radius = radius;
//	}
//
//	*inner_radius = max_inner_radius;
//	*outer_radius = min_outer_radius;
//
//	//Find the center - it is the point with value of the maximum radius which is surrounded with positions
//	// that have the highest radius values (we sum these values and obtain the center position)
//
//	//List of points that have maximum radius value
//	bdList<Point> list_of_maximum_radius_points;
//
//	Point point;
//
//	//Fill in the list of points that have maximum radius value
//	for(region.For_TheWholeRegion(); region.GetAllPoints(&r,&c); )
//	{
//		if(temp_inner_radius_image[r][c]==max_inner_radius) list_of_maximum_radius_points.AddToRightEnd(point(0,r,c));
//	}
//
//	//If there is only one element in the list, no need for further checks
//	if(list_of_maximum_radius_points.NumberOfElements()==1)
//	{
//		inner_radius_center_point = list_of_maximum_radius_points.GetLeftEnd();
//		return;
//	}
//
//	//Center point is the one that has max value and is surrounded by highest values of radiuses in 8-neighborhood
//	unsigned int radius_sum, max_radius_sum = 0;
//	for(list_of_maximum_radius_points.For_StartFromLeftEnd(); list_of_maximum_radius_points.GetAllElemetsByMovingToRight(&point); )
//	{
//		radius_sum = 0;
//		for(g.ForPoint(point.r,point.c); g.GetPositionsInCircleOfSquareRadius(max_inner_radius*max_inner_radius,&rn,&cn);)
//		{
//			radius_sum += (unsigned int)(temp_inner_radius_image[rn][cn]);
//		}
//		if(max_radius_sum<radius_sum)
//		{
//			max_radius_sum = radius_sum;
//			inner_radius_center_point = point;
//		}
//	}
//
//}
//
//
//
//
//
//
//template<class T>
//void bdRegion2D<T>::CenterOfMassPointCalculation(int *r_center_of_mass, int *c_center_of_mass)
//{
//	assert(!Empty());
//
//	long int r_sum = 0;
//	long int c_sum = 0;
//
//	int r,c;
//	
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		r_sum += r;
//		c_sum += c;
//	}
//
//	*r_center_of_mass = (int) (r_sum / m_number_of_elements);
//	*c_center_of_mass = (int) (c_sum / m_number_of_elements);
//}
//
//
//
//
//template<class T>
//void bdRegion2D<T>::CenterOfMassPointCalculation(Point &center_of_mass)
//{
//	assert(!Empty());
//
//	long int r_sum = 0;
//	long int c_sum = 0;
//
//	int r,c;
//	
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		r_sum += r;
//		c_sum += c;
//	}
//
//	center_of_mass.r = (int) (r_sum / m_number_of_elements);
//	center_of_mass.c = (int) (c_sum / m_number_of_elements);
//}
//
//
//
//
//template<class T>
//void bdRegion2D<T>::CenterOfMassPositionCalculation(Position &center_of_mass)
//{
//	assert(!Empty());
//
//	long int r_sum = 0;
//	long int c_sum = 0;
//
//	int r,c;
//	
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		r_sum += r;
//		c_sum += c;
//	}
//
//	center_of_mass.r = ((double)(r_sum)) / ((double)(m_number_of_elements));
//	center_of_mass.c = ((double)(c_sum)) / ((double)(m_number_of_elements));
//}
//
//
//
//
//
//
//template<class T>
//int bdRegion2D<T>::OverlappingAreaWithImage(Image8U &image)
//{
//	if(Empty() || image.Empty()) return 0;
//
//	int overlapping_area = 0;
//
//	int r,c;
//	
//	for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//	{
//		if(image[r][c]==m_region_value) overlapping_area++;
//	}
//	
//
//	return overlapping_area;
//}
//
//
//
//template<class T>
//int bdRegion2D<T>::OverlappingAreaWithRegion(bdRegion2D<T> &region)
//{
//	if(Empty() || region.Empty()) return 0;
//
//	int overlapping_area = 0;
//
//	int r,c;
//
//	//If this region is smaller or equal to region2
//	if(NumberOfElements()<=region.NumberOfElements())
//	{
//		for(For_TheWholeRegion(); GetAllPoints(&r,&c); )
//		{
//			if(region.PointIsInTheRegion(r,c)) overlapping_area++;
//		}
//	}
//	//If region2 is smaller than this region
//	else
//	{
//		for(region.For_TheWholeRegion(); region.GetAllPoints(&r,&c); )
//		{
//			if(PointIsInTheRegion(r,c)) overlapping_area++;
//		}
//	}
//
//	return overlapping_area;
//}
//
//


int bdRegion2D::Internal_4_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	list_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		int n_of_neighbors = 0;//we need this to detect if the edge of a region is on the edge of an image!
		bdDiscreteCoordinates3D p;
		for(g.ForCoordinates_4_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_4_Neighbors(rn,cn); )
		{
			n_of_neighbors++;

			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point. 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				list_of_egde_points.AddToLeftEnd(p(0,it.GetIndexRow(),it.GetIndexColumn()));
				n_of_neighbors = 4; //set the indicator, so that the next IF outside the FOR loop is not entered
				break;
			}
		}

		//If the number of neighbors for the point was less than 4, this point is on the edge of the image, so it IS an internal edge point!
		if(n_of_neighbors<4) list_of_egde_points.AddToLeftEnd(p(0,it.GetIndexRow(),it.GetIndexColumn()));
	}
	return 1;
}


int bdRegion2D::Internal_4_EdgePointsRegion(bdRegion2D &region_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	region_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		int n_of_neighbors = 0;//we need this to detect if the edge of a region is on the edge of an image!
		bdDiscreteCoordinates3D p;
		for(g.ForCoordinates_4_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_4_Neighbors(rn,cn); )
		{
			n_of_neighbors++;

			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point. 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				region_of_egde_points.AddPoint(it.GetIndexRow(),it.GetIndexColumn(),(&(it.GetElement())));
				n_of_neighbors = 4; //set the indicator, so that the next IF outside the FOR loop is not entered
				break;
			}
		}

		//If the number of neighbors for the point was less than 4, this point is on the edge of the image, so it IS an internal edge point!
		if(n_of_neighbors<4) region_of_egde_points.AddPoint(it.GetIndexRow(),it.GetIndexColumn(),(&(it.GetElement())));
	}
	return 1;
}


int bdRegion2D::Internal_8_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	list_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
        //cout<<" R ";
		int rn, cn;
		int n_of_neighbors = 0;//we need this to detect if the edge of a region is on the edge of an image!
		bdDiscreteCoordinates3D p;
		for(g.ForCoordinates_8_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_8_Neighbors(rn,cn); )
		{
			n_of_neighbors++;

			//If the 8-neighborhood point doesn't belong to the region, this point is an edge point. 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				list_of_egde_points.AddToLeftEnd(p(0,it.GetIndexRow(),it.GetIndexColumn()));
				n_of_neighbors = 8; //set the indicator, so that the next IF outside the FOR loop is not entered
				break;
			}
		}

		//If the number of neighbors for the point was less than 8, this point is on the edge of the image, so it IS an internal edge point!
		if(n_of_neighbors<8) list_of_egde_points.AddToLeftEnd(p(0,it.GetIndexRow(),it.GetIndexColumn()));
	}
	return 1;
}


int bdRegion2D::Internal_8_EdgePointsRegion(bdRegion2D &region_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	region_of_egde_points.Reset();
	region_of_egde_points.SetSizeOfOriginalImage(m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		int n_of_neighbors = 0;//we need this to detect if the edge of a region is on the edge of an image!
		for(g.ForCoordinates_8_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_8_Neighbors(rn,cn); )
		{
			n_of_neighbors++;

			//If the 8-neighborhood point doesn't belong to the region, this point is an edge point. 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				region_of_egde_points.AddPoint(it.GetIndexRow(),it.GetIndexColumn(),(&(it.GetElement())));
				n_of_neighbors = 8; //set the indicator, so that the next IF outside the FOR loop is not entered
				break;
			}
		}

		//If the number of neighbors for the point was less than 4, this point is on the edge of the image, so it IS an internal edge point!
		if(n_of_neighbors<8) region_of_egde_points.AddPoint(it.GetIndexRow(),it.GetIndexColumn(),(&(it.GetElement())));
	}
	return 1;
}


//int bdRegion2D::External_4_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
//{
//	if(!this->IsEmpty()) return 0;
//
//	list_of_egde_points.Reset();
//
//	Point point;
//
//	int r,c;
//
//	bdRegion2D temp;
//
//	if(!this-RegionOfExternal_4_EdgePoints(temp)) return 0;
//
//	for(temp.For_TheWholeRegion(); temp.GetAllPoints(&r,&c);)
//	{
//		list_of_egde_points.AddToLeftEnd(point(0,r,c));
//	}
//
//	temp.Delete();
//}


int bdRegion2D::External_4_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	list_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		for(g.ForCoordinates_4_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_4_Neighbors(rn,cn); )
		{
			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				bdDiscreteCoordinates3D p;
				list_of_egde_points.AddToLeftEnd(p(rn,cn));//We first add to the region because multiple points are eliminated, and then copy to the list
			}
		}
	}
	return 1;
}


int bdRegion2D::External_4_EdgePointsRegion(bdRegion2D &region_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	region_of_egde_points.Reset();
	region_of_egde_points.SetSizeOfOriginalImage(m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		for(g.ForCoordinates_4_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_4_Neighbors(rn,cn); )
		{
			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				region_of_egde_points.AddPoint(rn,cn);//We first add to the region because multiple points are eliminated, and then copy to the list
			}
		}
	}
	return 1;
}


int bdRegion2D::External_8_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	list_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		for(g.ForCoordinates_8_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_8_Neighbors(rn,cn); )
		{
			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				bdDiscreteCoordinates3D p;
				list_of_egde_points.AddToLeftEnd(p(rn,cn));//We first add to the region because multiple points are eliminated, and then copy to the list
			}
		}
	}
	return 1;
}


int bdRegion2D::External_8_EdgePointsRegion(bdRegion2D &region_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	region_of_egde_points.Reset();
	region_of_egde_points.SetSizeOfOriginalImage(m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdGeometry g;
	g.SetDimensions(1, m_original_image_size_CR[1], m_original_image_size_CR[0]);

	bdRegion2DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int rn, cn;
		for(g.ForCoordinates_8_Neighbors(it.GetIndexRow(),it.GetIndexColumn()); g.Get_8_Neighbors(rn,cn); )
		{
			//If the 4-neighborhood point doesn't belong to the region, this point is an edge point 
			if(!(this->IsPointInRegion(rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				region_of_egde_points.AddPoint(rn,cn);//We first add to the region because multiple points are eliminated, and then copy to the list
			}
		}
	}
	return 1;
}







//template<class T>
//void bdRegion2D<T>::MakeConvex(int distance_squared)
//{
//	bdList<Point> list_of_edge_points;
//	this->ListOfInternal_4_EdgePoints(list_of_edge_points);
//
//	Geometry g, g2;
//	g.SetDimensions(0,this->NumberOfRowsOfOriginalImage(),this->NumberOfColumnsOfOriginalImage());
//	g2.SetDimensions(0,this->NumberOfRowsOfOriginalImage(),this->NumberOfColumnsOfOriginalImage());
//
//	Image8U temp_image;
//	temp_image.Set(this->NumberOfRowsOfOriginalImage(),this->NumberOfColumnsOfOriginalImage());
//	temp_image.FillInWith(0);
//	//Paste the list of edge points to the temp image
//	Point p,p2;
//	for(list_of_edge_points.For_StartFromLeftEnd(); list_of_edge_points.GetAllElemetsByMovingToRight(&p); )
//	{
//		temp_image[p.r][p.c] = 255;
//	}
//
//	//For all the points in the list, find the edge points with the given distance, connect them with line and add to the region
//	for(list_of_edge_points.For_StartFromLeftEnd(); list_of_edge_points.GetAllElemetsByMovingToRight(&p); )
//	{
//		int rn,cn;
//		for(g.ForPoint(p.r,p.c); g.GetPositionsOfRadiusSquare(0,distance_squared,&rn,&cn); )
//		{
//			//If the point is in the temp_image, draw a line to the point, add all the line points to the region
//			if(temp_image[rn][cn]!=0)
//			{
//				Array<Point> line_points;
//				g2.Make3DLineSegment(p,p2(0,rn,cn),line_points);
//				for(int i=0; i<line_points.NumberOfElements(); i++)
//				{
//					this->AddPointToRegion(line_points[i].r,line_points[i].c);
//				}
//			}
//		}
//	}
//	
//	//Fill in the region
//	this->FillInTheWholeRegionBetweenOuterInternal_4_Edge();
//}



//int bdRegion2D::ListOfOuterInternal_8_EdgePoints(bdList<bdDiscreteCoordinates3D> &list_of_outer_internal_edge_points)
//{
//	assert(!Empty());
//	
//	list_of_outer_internal_edge_points.Delete();
//
//	Image8U framed_image;
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	// At this point all the region edges are NOT leaning on image edges (that is how the frame image is made)
//	// The idea is to start from the beginning of the image, detect the first edge and follow it internaly and
//	// externaly. The intersection of both edges give the outer edge.
//
//	int edge_found = 0;
//	int r_start,c_start;
//	for(int r=0; r<framed_image.NumberOfRows() && (!edge_found); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns()-1 && (!edge_found); c++)
//		{
//			//Work with image indexes, but enter positions to regions. This is because the position might be
//			// negative value, while index can not!
//
//			//Find a first edge and exit the loop
//			if(framed_image[r][c]==0 && framed_image[r][c+1]!=0) 
//			{
//				r_start = r;
//				c_start = c;
//				edge_found = 1;
//			}
//		}
//	}
//	//The first point in the internal edge has coordinates [r][c+1] and the external edge [r][c]
//
//	Geometry g, g2;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//	g2.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//
//	Point p;
//	bdList<Point> l;//, external_edge, internal_edge;
//	l.AddToLeftEnd(p(r_start,c_start));
//
//
//	//Get the external edge
//	//external_edge.AddToLeftEnd(p(r_start,c_start));
//	framed_image[r_start][c_start] = 1;
//	int rn,cn,rnn,cnn,r_pos,c_pos;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			if(framed_image[rn][cn]==0)
//			{				
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==255) 
//					{
//						framed_image[rn][cn] = 1;
//						//external_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Get the internal edge
//	l.AddToLeftEnd(p(r_start,c_start+1));
//	framed_image.ConvertIndexToPosition(r_start,c_start+1,&r_pos,&c_pos);
//	list_of_outer_internal_edge_points.AddToLeftEnd(p(r_pos,c_pos));
//	framed_image[r_start][c_start+1] = 2;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			//If the pixel might be the internal edge
//			if(framed_image[rn][cn]==255)
//			{
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==1) 
//					{
//						framed_image[rn][cn] = 2;
//						framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//						list_of_outer_internal_edge_points.AddToRightEnd(p(r_pos,c_pos));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//}




//template<class T>
//void bdRegion2D<T>::ListOfOuterInternal_4_EdgePoints(bdList<Point> &list_of_outer_internal_edge_points)
//{
//	assert(!Empty());
//	
//	list_of_outer_internal_edge_points.Delete();
//
//	Image8U framed_image;
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	//At this point all the region edges are NOT leaning on image edges (that is hoe the frame image is made)
//	// The idea is to start from the beginning of the image, detect the first edge and follow it internaly and
//	// externaly. The intersection of both edges giver the outer edge.
//
//	int edge_found = 0;
//	int r_start,c_start;
//	for(int r=0; r<framed_image.NumberOfRows() && (!edge_found); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns()-1 && (!edge_found); c++)
//		{
//			//Work with image indexes, but enter positions to regions. This is because the position might be
//			// negative value, while index can not!
//
//			//Find a first edge and exit the loop
//			if(framed_image[r][c]==0 && framed_image[r][c+1]!=0) 
//			{
//				r_start = r;
//				c_start = c;
//				edge_found = 1;
//			}
//		}
//	}
//	//The first point in the internal edge has coordinates [r][c+1] and the external edge [r][c]
//
//	Geometry g, g2;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//	g2.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//
//	Point p;
//	bdList<Point> l;//, external_edge, internal_edge;
//	l.AddToLeftEnd(p(r_start,c_start));
//
//
//	//Get the external edge
//	//external_edge.AddToLeftEnd(p(r_start,c_start));
//	framed_image[r_start][c_start] = 1;
//	int rn,cn,rnn,cnn,r_pos,c_pos;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//		{
//			if(framed_image[rn][cn]==0)
//			{				
//				for(g2.ForPoint(rn,cn); g2.Get_8_Positions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==255) 
//					{
//						framed_image[rn][cn] = 1;
//						//external_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Get the internal edge
//	l.AddToLeftEnd(p(r_start,c_start+1));
//	framed_image.ConvertIndexToPosition(r_start,c_start+1,&r_pos,&c_pos);
//	list_of_outer_internal_edge_points.AddToLeftEnd(p(r_pos,c_pos));
//	framed_image[r_start][c_start+1] = 2;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			//If the pixel might be the internal edge
//			if(framed_image[rn][cn]==255)
//			{
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==1) 
//					{
//						framed_image[rn][cn] = 2;
//						framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//						list_of_outer_internal_edge_points.AddToRightEnd(p(r_pos,c_pos));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//}
//
//
//
//
//
//
//template<class T> 
//void bdRegion2D<T>::FillInTheWholeRegionBetweenOuterInternal_8_Edge()
//{
//	assert(!Empty());
//
//	Image8U framed_image;
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	//At this point all the region edges are NOT leaning on image edges (that is hoe the frame image is made)
//	// The idea is to start from the beginning of the image, detect the first edge and follow it internaly and
//	// externaly. The intersection of both edges giver the outer edge.
//
//	int edge_found = 0;
//	int r_start,c_start;
//	for(int r=0; r<framed_image.NumberOfRows() && (!edge_found); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns()-1 && (!edge_found); c++)
//		{
//			//Work with image indexes, but enter positions to regions. This is because the position might be
//			// negative value, while index can not!
//
//			//Find a first edge and exit the loop
//			if(framed_image[r][c]==0 && framed_image[r][c+1]!=0) 
//			{
//				r_start = r;
//				c_start = c;
//				edge_found = 1;
//			}
//		}
//	}
//	//The first point in the internal edge has coordinates [r][c+1] and the external edge [r][c]
//
//	Geometry g, g2;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//	g2.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//
//	Point p;
//	bdList<Point> l, internal_edge;
//	l.AddToLeftEnd(p(r_start,c_start));
//
//
//	//Get the external edge
//	//external_edge.AddToLeftEnd(p(r_start,c_start));
//	framed_image[r_start][c_start] = 1;
//	int rn,cn,rnn,cnn,r_pos,c_pos;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			if(framed_image[rn][cn]==0)
//			{				
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==255) 
//					{
//						framed_image[rn][cn] = 1;
//						//external_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Get the internal edge
//	l.AddToLeftEnd(p(r_start,c_start+1));
//	internal_edge.AddToLeftEnd(p(r_start,c_start+1));
//	framed_image[r_start][c_start+1] = 2;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			//If the pixel might be the internal edge
//			if(framed_image[rn][cn]==255)
//			{
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==1) 
//					{
//						framed_image[rn][cn] = 2;
//						internal_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Now grow the internal edge to the inside
//	while(!internal_edge.Empty())
//	{
//		p = internal_edge.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//		{
//			//If the pixel is not a part of internal or external edge it might be candidate for adding to region
//			if(framed_image[rn][cn]!=1 && framed_image[rn][cn]!=2)
//			{
//				//If the pixel is already in this region
//				if(framed_image[rn][cn]==255)
//				{
//					//Convert the pixel to inner edge
//					framed_image[rn][cn]=2;
//					internal_edge.AddToRightEnd(p(rn,cn));
//				}
//				//If the pixel is part of the "hole" region
//				if(framed_image[rn][cn]==0)
//				{
//					//Convert the pixel to inner edge and add it to this region
//					framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//					this->AddPointToRegion(r_pos,c_pos);
//					framed_image[rn][cn]=2;
//					internal_edge.AddToRightEnd(p(rn,cn));
//				}
//			}
//		}
//		internal_edge.DeleteLeftEnd();
//	}
//
//}
//
//
//
//
//template<class T> 
//void bdRegion2D<T>::FillInTheWholeRegionBetweenOuterInternal_4_Edge()
//{
//	assert(!Empty());
//
//	Image8U framed_image;
//	this->ConvertRegionToFramedImage(framed_image, 255);
//
//	//At this point all the region edges are NOT leaning on image edges (that is hoe the frame image is made)
//	// The idea is to start from the beginning of the image, detect the first edge and follow it internaly and
//	// externaly. The intersection of both edges giver the outer edge.
//
//	int edge_found = 0;
//	int r_start,c_start;
//	for(int r=0; r<framed_image.NumberOfRows() && (!edge_found); r++)
//	{
//		for(int c=0; c<framed_image.NumberOfColumns()-1 && (!edge_found); c++)
//		{
//			//Work with image indexes, but enter positions to regions. This is because the position might be
//			// negative value, while index can not!
//
//			//Find a first edge and exit the loop
//			if(framed_image[r][c]==0 && framed_image[r][c+1]!=0) 
//			{
//				r_start = r;
//				c_start = c;
//				edge_found = 1;
//			}
//		}
//	}
//	//The first point in the internal edge has coordinates [r][c+1] and the external edge [r][c]
//
//	Geometry g, g2;
//	g.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//	g2.SetDimensions(0,framed_image.NumberOfRows(),framed_image.NumberOfColumns());
//
//	Point p;
//	bdList<Point> l, internal_edge;
//	l.AddToLeftEnd(p(r_start,c_start));
//
//
//	//Get the external edge
//	//external_edge.AddToLeftEnd(p(r_start,c_start));
//	framed_image[r_start][c_start] = 1;
//	int rn,cn,rnn,cnn,r_pos,c_pos;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//		{
//			if(framed_image[rn][cn]==0)
//			{				
//				for(g2.ForPoint(rn,cn); g2.Get_8_Positions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==255) 
//					{
//						framed_image[rn][cn] = 1;
//						//external_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Get the internal edge
//	l.AddToLeftEnd(p(r_start,c_start+1));
//	internal_edge.AddToLeftEnd(p(r_start,c_start+1));
//	framed_image[r_start][c_start+1] = 2;
//	while(!l.Empty())
//	{
//		p = l.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_8_Positions(&rn,&cn); )
//		{
//			//If the pixel might be the internal edge
//			if(framed_image[rn][cn]==255)
//			{
//				for(g2.ForPoint(rn,cn); g2.Get_4_ClosestPositions(&rnn,&cnn); )
//				{
//					if(framed_image[rnn][cnn]==1) 
//					{
//						framed_image[rn][cn] = 2;
//						internal_edge.AddToRightEnd(p(rn,cn));
//						l.AddToRightEnd(p(rn,cn));
//						g2.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//		}
//		l.DeleteLeftEnd();
//	}
//
//	//Now grow the internal edge to the inside
//	while(!internal_edge.Empty())
//	{
//		p = internal_edge.GetLeftEnd();
//		for(g.ForPoint(p.r,p.c); g.Get_4_ClosestPositions(&rn,&cn); )
//		{
//			//If the pixel is not a part of internal or external edge it might be candidate for adding to region
//			if(framed_image[rn][cn]!=1 && framed_image[rn][cn]!=2)
//			{
//				//If the pixel is already in this region
//				if(framed_image[rn][cn]==255)
//				{
//					//Convert the pixel to inner edge
//					framed_image[rn][cn]=2;
//					internal_edge.AddToRightEnd(p(rn,cn));
//				}
//				//If the pixel is part of the "hole" region
//				if(framed_image[rn][cn]==0)
//				{
//					//Convert the pixel to inner edge and add it to this region
//					framed_image.ConvertIndexToPosition(rn,cn,&r_pos,&c_pos);
//					this->AddPointToRegion(r_pos,c_pos);
//					framed_image[rn][cn]=2;
//					internal_edge.AddToRightEnd(p(rn,cn));
//				}
//			}
//		}
//		internal_edge.DeleteLeftEnd();
//	}
//
//}


int bdRegion2D::BreakIntoDisjointRegionsOf_4_Neighborhood(bdList<bdRegion2D> &list_of_regions)
{
	list_of_regions.Reset();
	
	bdGeometry g;
	g.SetDimensions(0,m_original_image_size_CR[1],m_original_image_size_CR[0]);

	bdRegion2D original_region;
	original_region = (*this);

	//Initialize the disjoint region
	bdRegion2D disjoint_region;
	disjoint_region.m_original_image_size_CR[1] = m_original_image_size_CR[1];
	disjoint_region.m_original_image_size_CR[0] = m_original_image_size_CR[0];
	disjoint_region.m_region_value = m_region_value;

	bdList<bdDiscreteCoordinates3D> list;

	bdDiscreteCoordinates3D point;

	int start_r, start_c, rn, cn;

	while(!original_region.IsEmpty())
	{
		disjoint_region.RemoveAllPoints();
		start_r = original_region.m_row_list.GetLeftEnd().m_row_index;
		start_c = original_region.m_row_list.GetLeftEnd().m_column_list.GetLeftEnd().m_column_index;

		list.AddToRightEnd(point(start_r,start_c));
		disjoint_region.AddPoint(point.R(),point.C(),this->GetValue(point.R(),point.C()));

		while(!list.IsEmpty())
		{
			point = list.GetLeftEnd();

			for(g.ForCoordinates_4_Neighbors(point.R(),point.C()); g.Get_4_Neighbors(rn,cn); )
			{
				if(original_region.IsPointInRegion(rn,cn) && !disjoint_region.IsPointInRegion(rn,cn))
				{
					list.AddToRightEnd(point(rn,cn));
					disjoint_region.AddPoint(rn,cn,this->GetValue(rn,cn));
				}
			}
			list.DeleteLeftEnd();
		}

		//After this part we have a disjoint region
		
		//Delete the points from original_region that exist in the disjoint region
		bdRegion2DIterator it;
		for(it.SetBegin(&disjoint_region); it.IsValid(); it.MoveToNext())
		{
			original_region.RemovePoint(it.GetIndexRow(),it.GetIndexColumn());
		}

		list_of_regions.AddToRightEnd(disjoint_region);
	}

	if(list_of_regions.GetNumberOfElements()<=1) return 0;//there are no disjoint region
	else return list_of_regions.GetNumberOfElements(); //there are disjoint regions
}


int bdRegion2D::BreakIntoDisjointRegionsOf_8_Neighborhood(bdList<bdRegion2D> &list_of_regions)
{
	list_of_regions.Reset();
	
	bdGeometry g;
	g.SetDimensions(0,m_original_image_size_CR[1],m_original_image_size_CR[0]);

	bdRegion2D original_region;
	original_region = (*this);

	//Initialize the disjoint region
	bdRegion2D disjoint_region;
	disjoint_region.m_original_image_size_CR[1] = m_original_image_size_CR[1];
	disjoint_region.m_original_image_size_CR[0] = m_original_image_size_CR[0];
	disjoint_region.m_region_value = m_region_value;

	bdList<bdDiscreteCoordinates3D> list;

	bdDiscreteCoordinates3D point;

	int start_r, start_c, rn, cn;

	while(!original_region.IsEmpty())
	{
		disjoint_region.RemoveAllPoints();
		start_r = original_region.m_row_list.GetLeftEnd().m_row_index;
		start_c = original_region.m_row_list.GetLeftEnd().m_column_list.GetLeftEnd().m_column_index;

		list.AddToRightEnd(point(start_r,start_c));
		disjoint_region.AddPoint(point.R(),point.C(),this->GetValue(point.R(),point.C()));

		while(!list.IsEmpty())
		{
			point = list.GetLeftEnd();

			for(g.ForCoordinates_8_Neighbors(point.R(),point.C()); g.Get_8_Neighbors(rn,cn); )
			{
				if(original_region.IsPointInRegion(rn,cn) && !disjoint_region.IsPointInRegion(rn,cn))
				{
					list.AddToRightEnd(point(rn,cn));
					disjoint_region.AddPoint(rn,cn,this->GetValue(rn,cn));
				}
			}
			list.DeleteLeftEnd();
		}

		//After this part we have a disjoint region
		
		//Delete the points from original_region that exist in the disjoint region
		bdRegion2DIterator it;
		for(it.SetBegin(&disjoint_region); it.IsValid(); it.MoveToNext())
		{
			original_region.RemovePoint(it.GetIndexRow(),it.GetIndexColumn());
		}

		list_of_regions.AddToRightEnd(disjoint_region);
	}

	if(list_of_regions.GetNumberOfElements()<=1) return 0;//there are no disjoint region
	else return list_of_regions.GetNumberOfElements(); //there are disjoint regions
}






//template<class T>
//void bdRegion2D<T>::CreateRegionFromTheListOfDisjointRegions(bdList<bdRegion2D<T>> &list_of_regions)
//{
//	assert(!list_of_regions.Empty());
//
//	Delete();
//
//	m_original_image_size_CR[1] = list_of_regions.GetLeftEnd().m_original_image_size_CR[1];
//	m_original_image_size_CR[0] = list_of_regions.GetLeftEnd().m_original_image_size_CR[0];
//	m_region_value = list_of_regions.GetLeftEnd().m_region_value;
//
//	T el;
//	int r,c;
//	bdRegion2D<T> *p_temp_region; 
//
//	for(list_of_regions.For_StartFromLeftEnd(); list_of_regions.GetAllPointersToElemetsByMovingToRight(&p_temp_region); )
//	{
//		for(p_temp_region->For_TheWholeRegion(); p_temp_region->GetAllPointsAndElements(&r,&c,&el); )
//		{
//			AddPointToRegion(r,c,el);
//		}
//	}
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::IntersectionOfTwoRegions(bdRegion2D<T> &region1, bdRegion2D<T> &region2)
//{
//	if(region1.Empty() || region2.Empty()) 
//	{
//		DeleteAllNodes();
//		return;
//	}
//
//	int r,c;
//
//	//If region1 is smaller or equal to region2
//	if(region1.NumberOfElements()<=region2.NumberOfElements())
//	{
//		for(region1.For_TheWholeRegion(); region1.GetAllPoints(&r,&c); )
//		{
//			if(region2.PointIsInTheRegion(r,c)) AddPointToRegion(r,c);
//		}
//	}
//	//If region2 is smaller than region1
//	else
//	{
//		for(region2.For_TheWholeRegion(); region2.GetAllPoints(&r,&c); )
//		{
//			if(region1.PointIsInTheRegion(r,c)) AddPointToRegion(r,c);
//		}
//	}
//}
//
//
//
//template <class T>
//void bdRegion2D<T>::Erode_4()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfInternal_4_EdgePoints(list);
//	
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		DeletePoint(point.r, point.c);
//	}
//}
//
//
//template <class T>
//void bdRegion2D<T>::Erode_4_NotAppliedToImageEdges()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfInternal_4_EdgePoints(list);
//	
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		if(point.r!=0 && point.c!=0 && point.r!=this->NumberOfRowsOfOriginalImage()-1 && point.c!=this->NumberOfColumnsOfOriginalImage()-1)
//			DeletePoint(point.r, point.c);
//	}
//}
//
//
//
//template <class T>
//void bdRegion2D<T>::Erode_8()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfInternal_8_EdgePoints(list);
//	
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		DeletePoint(point.r, point.c);
//	}
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Dilate_4()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfExternal_4_EdgePoints(list);
//
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		AddPointToRegion(point.r,point.c);
//	}
//}
//
//
//
//template <class T>
//void bdRegion2D<T>::Dilate_4_NotAppliedToImageEdges()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfExternal_4_EdgePoints(list);
//
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		if(point.r>1 && point.c>1 && point.r<this->NumberOfRowsOfOriginalImage()-2 && point.c<this->NumberOfColumnsOfOriginalImage()-2)
//			AddPointToRegion(point.r,point.c);
//	}
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Dilate_8()
//{
//	bdList<Point> list;
//	Point point;
//
//	ListOfExternal_8_EdgePoints(list);
//
//	for(list.For_StartFromLeftEnd(); list.GetAllElemetsByMovingToRight(&point); )
//	{
//		AddPointToRegion(point.r,point.c);
//	}
//}
//
//
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Opening_4(int radius)
//{
//	assert(radius>0);
//
//	for(int i=0; i<radius; i++) Erode_4();
//
//	for(int i=0; i<radius; i++) Dilate_4();
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Opening_8(int radius)
//{
//	assert(radius>0);
//
//	for(int i=0; i<radius; i++) Erode_8();
//
//	for(int i=0; i<radius; i++) Dilate_8();
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Closing_4(int radius)
//{
//	assert(radius>0);
//
//	for(int i=0; i<radius; i++) Dilate_4();
//
//	for(int i=0; i<radius; i++) Erode_4();
//}
//
//
//
//
//template <class T>
//void bdRegion2D<T>::Closing_8(int radius)
//{
//	assert(radius>0);
//
//	for(int i=0; i<radius; i++) Dilate_8();
//
//	for(int i=0; i<radius; i++) Erode_8();
//}
//
//
//
//template <class T>
//void bdRegion2D<T>::FrameCoordinatesForRegion(int *r_start_included, int *r_end_included, int *c_start_included, int *c_end_included)
//{
//	assert(!Empty());
//
//	int rn,cn;
//
//	int is_first_point = 1;
//
//	for(this->For_TheWholeRegion(); this->GetAllPoints(&rn,&cn); )
//	{
//		if(is_first_point)
//		{
//			*r_start_included = rn;
//			*r_end_included = rn;
//			*c_start_included = cn;
//			*c_end_included = cn;
//			is_first_point = 0;
//		}
//		else
//		{
//			if(rn<*r_start_included) *r_start_included = rn;
//			if(rn>*r_end_included) *r_end_included = rn;
//			if(cn<*c_start_included) *c_start_included = cn;
//			if(cn>*c_end_included) *c_end_included = cn;
//		}
//	}
//}
//
//
//
//
//template <class T>
//ostream& operator << <>(ostream &o, bdRegion2D<T> &r)
//{ 
//	int rn,cn;
//	T *pel;
//
//	for(r.For_TheWholeRegion(); r.GetAllPointsAndPointersToElements(&rn,&cn,&pel); )
//	{
//		o<<"("<<rn<<","<<cn<<")"<<(*pel)<<" ";
//	}
//    return o;
//}



template class BD_REGION_2D_API bdList<bdRegionColumnNode>;
template class BD_REGION_2D_API bdListNode<bdRegionRowNode>;
template class BD_REGION_2D_API bdListNode<bdRegionColumnNode>;
template class BD_REGION_2D_API bdList<bdRegionRowNode>;
template class BD_REGION_2D_API bdRegularGridT<unsigned short>;//><bdRegionRowNode>;


