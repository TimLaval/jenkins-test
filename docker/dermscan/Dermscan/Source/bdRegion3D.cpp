/**************************************************************************
 Class of 3D region 

 Author: Danilo Babin
 File name: "bdRegion3D.cpp"
***************************************************************************/



#include "bdRegion3D.h"




bdRegion3DNode::bdRegion3DNode()
{
	m_node_of_slice_node = NULL;
	//m_region_2D_node = NULL;
}


bdRegion3DNode::~bdRegion3DNode()
{
	m_node_of_slice_node = NULL;
	//m_region_2D_node = NULL;
}


int bdRegion3DNode::IsNonZero()
{
	if((!m_node_of_slice_node) || (!m_region_2D_node.IsNonZero())) return 0;
	return 1;
}


int& bdRegion3DNode::GetElement()
{
	return m_region_2D_node.GetElement();
}


int bdRegion3DNode::GetIndexSlice()
{
	return m_node_of_slice_node->GetElementPointer()->m_slice_index;
}


int bdRegion3DNode::GetIndexRow()
{
	return m_region_2D_node.GetIndexRow();
}
	

int bdRegion3DNode::GetIndexColumn()
{
	return m_region_2D_node.GetIndexColumn();
}




//-------------------------------------------------------------------------------------------------------------------------------------------------------------





bdRegion3DIterator::bdRegion3DIterator()
{
	m_node_of_slice_node = NULL;
	//m_node_of_row_node = NULL;
}


bdRegion3DIterator::~bdRegion3DIterator()
{
	m_node_of_slice_node = NULL;
	//m_node_of_row_node = NULL;
}


void bdRegion3DIterator::SetBegin(bdRegion3D *region)
{
	if(!region)
	{
		m_node_of_slice_node = NULL;
		//m_node_of_row_node = NULL;
		return;
	}
	m_node_of_slice_node = region->GetSliceList()->GetLeftEndNodePointer();
	m_region_2D_node = m_node_of_slice_node->GetElementPointer()->m_slice_region_2D.GetNodeForFirstPoint();
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


int bdRegion3DIterator::IsValid()
{
	if(m_node_of_slice_node && m_region_2D_node.IsNonZero()) return 1;
	return 0;
}


void bdRegion3DIterator::MoveToNext()
{
	m_region_2D_node.MoveToNext();
	if(!m_region_2D_node.IsNonZero())
	{
		m_node_of_slice_node = m_node_of_slice_node->GetRight();
		if(m_node_of_slice_node)
		{
			m_region_2D_node = m_node_of_slice_node->GetElementPointer()->m_slice_region_2D.GetNodeForFirstPoint();
		}
	}
}





//-------------------------------------------------------------------------------------------------------------------------------------------------------------




bdRegion3D::bdRegion3D()
{
	m_original_image_size_CRS[2] = m_original_image_size_CRS[1] = m_original_image_size_CRS[0] = 0;
	m_region_value = 0;
	m_number_of_elements = 0;
}


bdRegion3D::~bdRegion3D()
{
	this->Reset();
}


int bdRegion3D::IsEmpty()
{
	if(m_number_of_elements!=0) return 0;
	else return 1;
}


int bdRegion3D::Bounds(unsigned int &output_s_min, unsigned int &output_s_max, unsigned int &output_r_min, unsigned int &output_r_max, unsigned int &output_c_min, unsigned int &output_c_max)
{
    if(this->IsEmpty()) return 0;
    
    bdRegion3DIterator it;
    it.SetBegin(this);
    output_s_min = output_s_max = it.GetIndexSlice();
    output_r_min = output_r_max = it.GetIndexRow();
    output_c_min = output_c_max = it.GetIndexColumn();
    for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
    {
        if(it.GetIndexSlice()<output_s_min) output_s_min = it.GetIndexSlice();
        if(it.GetIndexSlice()>output_s_max) output_s_max = it.GetIndexSlice();
        if(it.GetIndexRow()<output_r_min) output_r_min = it.GetIndexRow();
        if(it.GetIndexRow()>output_r_max) output_r_max = it.GetIndexRow();
        if(it.GetIndexColumn()<output_c_min) output_c_min = it.GetIndexColumn();
        if(it.GetIndexColumn()>output_c_max) output_c_max = it.GetIndexColumn();
    }
    
    return 1;
}


int bdRegion3D::BoundingPoints(bdDiscreteCoordinates3D &output_s_min, bdDiscreteCoordinates3D &output_s_max, bdDiscreteCoordinates3D &output_r_min, bdDiscreteCoordinates3D &output_r_max, bdDiscreteCoordinates3D &output_c_min, bdDiscreteCoordinates3D &output_c_max)
{
    if(this->IsEmpty()) return 0;
    
    bdRegion3DIterator it;
    it.SetBegin(this);
    output_s_min.S() = output_s_max.S() = output_r_min.S() = output_r_max.S() = output_c_min.S() = output_c_max.S() = it.GetIndexSlice();
    output_s_min.R() = output_s_max.R() = output_r_min.R() = output_r_max.R() = output_c_min.R() = output_c_max.R() = it.GetIndexRow();
    output_s_min.C() = output_s_max.C() = output_r_min.C() = output_r_max.C() = output_c_min.C() = output_c_max.C() = it.GetIndexColumn();

    for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
    {
        if(it.GetIndexSlice()<output_s_min.S())
        {
            output_s_min.S() = it.GetIndexSlice();
            output_s_min.R() = it.GetIndexRow();
            output_s_min.C() = it.GetIndexColumn();
        }
        if(it.GetIndexSlice()>output_s_max.S())
        {
            output_s_max.S() = it.GetIndexSlice();
            output_s_max.R() = it.GetIndexRow();
            output_s_max.C() = it.GetIndexColumn();
        }
        if(it.GetIndexRow()<output_r_min.R())
        {
            output_r_min.S() = it.GetIndexSlice();
            output_r_min.R() = it.GetIndexRow();
            output_r_min.C() = it.GetIndexColumn();
        }
        if(it.GetIndexRow()>output_r_max.R())
        {
            output_r_max.S() = it.GetIndexSlice();
            output_r_max.R() = it.GetIndexRow();
            output_r_max.C() = it.GetIndexColumn();
        }
        if(it.GetIndexColumn()<output_c_min.C())
        {
            output_c_min.S() = it.GetIndexSlice();
            output_c_min.R() = it.GetIndexRow();
            output_c_min.C() = it.GetIndexColumn();
        }
        if(it.GetIndexColumn()>output_c_max.C())
        {
            output_c_max.S() = it.GetIndexSlice();
            output_c_max.R() = it.GetIndexRow();
            output_c_max.C() = it.GetIndexColumn();
        }
    }
    
    return 1;
}


bdListNode<bdRegionSliceNode>* bdRegion3D::GetNodeOfSliceNode(unsigned int s_index)
{
	bdListIterator<bdRegionSliceNode> it;
	for(it.SetLeftEnd(m_slice_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_slice_index == s_index)
			return it.GetNodeAddress();
		if(it.GetElementPointer()->m_slice_index > s_index) return NULL;
	}
	return NULL;
}


bdRegion3DNode bdRegion3D::GetNode(unsigned int s_index, unsigned int r_index, unsigned int c_index)
{
	bdRegion3DNode node;
	node.m_node_of_slice_node = NULL;
	//node.m_node_of_column_node = NULL;

	//If empty, local index not found
	if(this->IsEmpty()) return node;

	//bdListNode<bdRegionSliceNode> *slice_node = this->GetNodeOfSliceNode(s_index);
	//if(!slice_node) return node;

	node.m_node_of_slice_node = this->GetNodeOfSliceNode(s_index);
	if(!node.m_node_of_slice_node) return node;

	node.m_region_2D_node = node.m_node_of_slice_node->GetElementPointer()->m_slice_region_2D.GetNode(r_index,c_index);

	return node;
}




void bdRegion3D::Reset()
{
	m_slice_list.Reset();
	m_original_image_size_CRS[2] = m_original_image_size_CRS[1] = m_original_image_size_CRS[0] = 0;
	m_region_value = 0;
	m_number_of_elements = 0;
}


void bdRegion3D::RemoveAllPoints()
{
	//It changes only slice list and does not change other member variables.
	m_slice_list.Reset();
	m_number_of_elements = 0;
}


int bdRegion3D::RemovePoint(int s_index, int r_index, int c_index)
{
	bdRegion3DNode node = this->GetNode(s_index, r_index, c_index);
	return this->RemovePoint(node);
}


int bdRegion3D::RemovePoint(bdRegion3DNode &node)
{
	if(this->IsEmpty()) return 1;

	if(!node.IsNonZero()) return 0;

	node.m_node_of_slice_node->GetElementPointer()->m_slice_region_2D.RemovePoint(node.m_region_2D_node);
	m_number_of_elements--;

	//If a region list had only one element and is now empty, remove it from the slice list
	if(node.m_node_of_slice_node->GetElementPointer()->m_slice_region_2D.IsEmpty())
		{ this->m_slice_list.DeleteNode(node.m_node_of_slice_node); }

	return 1;
}


int bdRegion3D::AddPoint(unsigned int s, unsigned int r, unsigned int c, int &element)
{
	if(this->IsEmpty()) 
	{
		bdRegionSliceNode *sn = this->m_slice_list.AddNewToRightEnd();
		sn->m_slice_index = s;
		sn->m_slice_region_2D.AddPoint(r,c,element);
		m_number_of_elements++;
		return 1;
	}

	//If the region is not empty
	//cout<<"--> find pointer location for new coordinates"<<endl;

	// Find the node of slice_node for which the index is the last one smaller or equal to 's'.
	bdListIterator<bdRegionSliceNode> it;
	bdListNode<bdRegionSliceNode> *recorded_node_of_slice_node = m_slice_list.GetLeftEndNodePointer();
	for(it.SetLeftEnd(m_slice_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_slice_index<=s) recorded_node_of_slice_node = it.GetNodeAddress();
		else { break; }
	}

	// There are 3 possible cases:

	//--- 1. The exact slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index == s)
	{
		//... add a new element to the slice list.
		recorded_node_of_slice_node->GetElementPointer()->m_slice_region_2D.AddPoint(r,c,element);
		m_number_of_elements++;
		return 1;
	}
	//--- 2. Lower slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index < s)
	{
		//... make a new slice node with given element and add to right of the recorded_node_of_slice_node.
		bdRegionSliceNode sn;
		sn.m_slice_index = s;
		sn.m_slice_region_2D.AddPoint(r,c,element);
		m_slice_list.AddToRightOfNode(recorded_node_of_slice_node,sn);
		m_number_of_elements++;
		return 1;
	}
	//--- 3. Higher slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index > s)
	{
		//... make a new slice node with given element and add to left of the recorded_node_of_slice_node.
		bdRegionSliceNode sn;
		sn.m_slice_index = s;
		sn.m_slice_region_2D.AddPoint(r,c,element);
		m_slice_list.AddToLeftOfNode(recorded_node_of_slice_node,sn);
		m_number_of_elements++;
		return 1;
	}

	return 1;//This means that adding succeeded
}


int bdRegion3D::AddPoint(unsigned int s, unsigned int r, unsigned int c, int *p_element)
{
	if(!p_element) 
	{
		return(this->AddPoint(s,r,c));
		//return 0;
	}
	return (this->AddPoint(s,r,c,*p_element));
}


int bdRegion3D::AddPoint(unsigned int s, unsigned int r, unsigned int c)
{
	int el = 0;
	return (this->AddPoint(s,r,c,el));
}


void bdRegion3D::AddRegion(bdRegion2D &region, unsigned int s)
{
	if(this->IsEmpty()) 
	{
		bdRegionSliceNode *sn = this->m_slice_list.AddNewToRightEnd();
		sn->m_slice_index = s;
		sn->m_slice_region_2D = region;
		m_number_of_elements = region.GetNumberOfElements();
	}

	//If the region is not empty
	//cout<<"--> find pointer location for new coordinates"<<endl;

	// Find the node of slice_node for which the index is the last one smaller or equal to 's'.
	bdListIterator<bdRegionSliceNode> it;
	bdListNode<bdRegionSliceNode> *recorded_node_of_slice_node = m_slice_list.GetLeftEndNodePointer();
	for(it.SetLeftEnd(m_slice_list); it.IsValid(); it.MoveRight())
	{
		if(it.GetElementPointer()->m_slice_index<=s) recorded_node_of_slice_node = it.GetNodeAddress();
		else { break; }
	}

	// There are 3 possible cases:

	//--- 1. The exact slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index == s)
	{
		//... add a new element to the slice list.
		recorded_node_of_slice_node->GetElementPointer()->m_slice_region_2D.AddRegion(region);
		m_number_of_elements += region.GetNumberOfElements();
	}
	//--- 2. Lower slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index < s)
	{
		//... make a new slice node with given element and add to right of the recorded_node_of_slice_node.
		bdRegionSliceNode sn;
		sn.m_slice_index = s;
		sn.m_slice_region_2D = region;
		m_slice_list.AddToRightOfNode(recorded_node_of_slice_node,sn);
		m_number_of_elements += region.GetNumberOfElements();
	}
	//--- 3. Higher slice index was found ---
	if(recorded_node_of_slice_node->GetElementPointer()->m_slice_index > s)
	{
		//... make a new slice node with given element and add to left of the recorded_node_of_slice_node.
		bdRegionSliceNode sn;
		sn.m_slice_index = s;
		sn.m_slice_region_2D = region;
		m_slice_list.AddToLeftOfNode(recorded_node_of_slice_node,sn);
		m_number_of_elements += region.GetNumberOfElements();
	}
}


void bdRegion3D::AddRegion(bdRegion3D &region)
{
	if(this->IsEmpty()) 
	{
		(*this) = region;
		return;
	}

	if(region.IsEmpty()) return;

	bdRegion3DIterator it;
	for(it.SetBegin(&region); it.IsValid(); it.MoveToNext())
	{
		this->AddPoint(it.GetIndexRow(),it.GetIndexColumn(),it.GetElement());
	}
}


int bdRegion3D::IsPointInRegion(unsigned int s_index, unsigned int r_index, unsigned int c_index)
{
	//If empty, local index not found
	if(this->IsEmpty()) return 0;

	bdListNode<bdRegionSliceNode> *slice_node = this->GetNodeOfSliceNode(s_index);
	if(!slice_node) return 0;

	return (slice_node->GetElementPointer()->m_slice_region_2D.IsPointInRegion(r_index,c_index));
}


int bdRegion3D::IsOverlappingWithRegion(bdRegion3D &region)
{
	bdRegion3DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		if(region.IsPointInRegion(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn())) { return 1; }
	}
	return 0;
}


int* bdRegion3D::GetValue(unsigned int s, unsigned int r, unsigned int c)
{
	bdRegion3DNode node = this->GetNode(s,r,c);
	if(node.IsNonZero()) return (&(node.GetElement()));
	return NULL;
}


bdRegion3D& bdRegion3D::operator =(bdRegion3D &region)
{
	if (&region==this) return *this;

	this->Reset();

	m_slice_list = region.m_slice_list;
	m_region_value = region.m_region_value;
	m_number_of_elements = region.m_number_of_elements;
	m_original_image_size_CRS[2] = region.m_original_image_size_CRS[2];
	m_original_image_size_CRS[1] = region.m_original_image_size_CRS[1];
	m_original_image_size_CRS[0] = region.m_original_image_size_CRS[0];

	return *this;
}


void bdRegion3D::SetSizeOfOriginalImage(unsigned int slices, unsigned int rows, unsigned int columns)
{	
	m_original_image_size_CRS[2] = slices;
	m_original_image_size_CRS[1] = rows;
	m_original_image_size_CRS[0] = columns;
}


int bdRegion3D::CreateRegionFromImage(bdImage16U &image, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CRS[2] = image.GetNumberOfSlices();
	m_original_image_size_CRS[1] = image.GetNumberOfRows();
	m_original_image_size_CRS[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(255);

	for(unsigned int s=0; s<image.GetNumberOfSlices(); s++)
	{
		for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
			{
				if(image(t,s,r,c)!=0)
				{
					this->AddPoint(s,r,c,(int&)image(t,s,r,c));
				}
			}
		}
	}
	return 1;
}


int bdRegion3D::CreateRegionFromImage(bdImage16U &image, unsigned short value_for_growing, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CRS[2] = image.GetNumberOfSlices();
	m_original_image_size_CRS[1] = image.GetNumberOfRows();
	m_original_image_size_CRS[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(value_for_growing);

	for(unsigned int s=0; s<image.GetNumberOfSlices(); s++)
	{
		for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
			{
				if(image(t,s,r,c)==value_for_growing)
				{
					this->AddPoint(s,r,c,(int&)value_for_growing);
				}
			}
		}
	}
	return 1;
}


int bdRegion3D::CreateRegion_26_FromSeed(bdImage16U &image, unsigned int c, unsigned int r, unsigned int s, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(c>=image.GetNumberOfColumns() || r>=image.GetNumberOfRows() || s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CRS[2] = image.GetNumberOfSlices();
	m_original_image_size_CRS[1] = image.GetNumberOfRows();
	m_original_image_size_CRS[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(image(t,s,r,c));

	//'temp_image' contains pixels that are in the region
	bdImage temp_image;
	temp_image.SetSize(1,image.GetNumberOfSlices(),image.GetNumberOfRows(),image.GetNumberOfColumns());
	temp_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(image.GetNumberOfSlices(), image.GetNumberOfRows(), image.GetNumberOfColumns());

	bdDiscreteCoordinates3D p(s,r,c);
	bdList<bdDiscreteCoordinates3D> list_of_points;
	list_of_points.AddToRightEnd(p);

	this->AddPoint(p.S(),p.R(),p.C(),(int&)image(t,s,r,c));
	temp_image(t,p.S(),p.R(),p.C()) = 1;

	int sn,rn,cn;
	while(!list_of_points.IsEmpty())
	{
		p = list_of_points.GetLeftEnd();
		for(g.ForCoordinates_26_Neighbors(p.S(),p.R(),p.C()); g.Get_26_Neighbors(sn,rn,cn); )
		{
			if(image(t,sn,rn,cn)!=0)
			{
				//If the point is not in the region
				if(temp_image(t,sn,rn,cn)==0) 
				{
					list_of_points.AddToRightEnd(p(sn,rn,cn));
					this->AddPoint(sn,rn,cn,(int&)image(t,sn,rn,cn));
					temp_image(t,sn,rn,cn) = 1;
				}
			}
		}
		list_of_points.DeleteLeftEnd();
	}

	return 1;
}




int bdRegion3D::CreateRegionFrom_27_Neighborhood(bdImage16U &image, unsigned int c, unsigned int r, unsigned int s, unsigned int t)
{
	if(image.IsEmpty()) return 0;
	if(c>=image.GetNumberOfColumns() || r>=image.GetNumberOfRows() || s>=image.GetNumberOfSlices() || t>=image.GetNumberOfTimeSeries()) return 0;

	this->Reset();

	m_original_image_size_CRS[2] = image.GetNumberOfSlices();
	m_original_image_size_CRS[1] = image.GetNumberOfRows();
	m_original_image_size_CRS[0] = image.GetNumberOfColumns();
	
	this->SetRegionValue(image(t,s,r,c));

	bdGeometry g;
	g.SetDimensions(image.GetNumberOfSlices(),image.GetNumberOfRows(),image.GetNumberOfColumns());

	int sn,rn,cn;
	for(g.ForCoordinates_27_Neighbors(s,r,c); g.Get_27_Neighbors(sn,rn,cn); )
	{
		if(image(t,sn,rn,cn) != 0)
		{
			this->AddPoint(sn,rn,cn,(int&)image(t,sn,rn,cn));
		}
	}
	return 1;
}


int bdRegion3D::CreateRegionFrom_27_Neighborhood(bdRegion3D &region, unsigned int c, unsigned int r, unsigned int s)
{
	if(region.IsEmpty()) return 0;
	if(c>=region.GetNumberOfColumnsOfOriginalImage() || r>=region.GetNumberOfRowsOfOriginalImage() || s>=region.GetNumberOfSlicesOfOriginalImage()) return 0;

	this->Reset();

	m_original_image_size_CRS[2] = region.GetNumberOfSlicesOfOriginalImage();
	m_original_image_size_CRS[1] = region.GetNumberOfRowsOfOriginalImage();
	m_original_image_size_CRS[0] = region.GetNumberOfColumnsOfOriginalImage();
	
	this->SetRegionValue(*(region.GetValue(s,r,c)));

	bdGeometry g;
	g.SetDimensions(region.GetNumberOfSlicesOfOriginalImage(),region.GetNumberOfRowsOfOriginalImage(),region.GetNumberOfColumnsOfOriginalImage());

	int sn,rn,cn;
	for(g.ForCoordinates_27_Neighbors(s,r,c); g.Get_27_Neighbors(sn,rn,cn); )
	{
		if(region.IsPointInRegion(sn,rn,cn))
		{
			this->AddPoint(sn,rn,cn,(*(region.GetValue(s,r,c))));
		}
	}
	return 1;
}



int bdRegion3D::CenterOfMass(unsigned int &output_s, unsigned int &output_r, unsigned int &output_c)
{
    if(this->IsEmpty()) return 0;
    
    unsigned int n = 0;
    output_s = output_r = output_c = 0;
    
    bdRegion3DIterator it;
    for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
    {
        output_s += it.GetIndexSlice();
        output_r += it.GetIndexRow();
        output_c += it.GetIndexColumn();
        n++;
    }
    if(n==0) return 0;
    
    output_s = output_s / n;
    output_r = output_r / n;
    output_c = output_c / n;
    
    return 1;
}



//
//template<class T>
//void bdRegion3D<T>::PasteRegionToImage(bdImage8U &image)
//{
//	assert(!Empty());
//	assert(!image.Empty());
//	assert(number_of_slices_of_original_image==image.NumberOfSlices() && number_of_rows_of_original_image==image.NumberOfRows() && number_of_columns_of_original_image==image.NumberOfColumns());
//
//	int s,r,c;
//
//	for(For_TheWholeRegion(); GetAllPoints(&s,&r,&c); )
//	{
//		image[s][r][c] = region_value;
//	}
//}
//


int bdRegion3D::Internal_26_EdgePointsList(bdList<bdDiscreteCoordinates3D> &list_of_egde_points)
{
	if(this->IsEmpty()) return 0;

	list_of_egde_points.Reset();

	bdGeometry g;
	g.SetDimensions(this->GetNumberOfSlicesOfOriginalImage(),this->GetNumberOfRowsOfOriginalImage(),this->GetNumberOfColumnsOfOriginalImage());

	bdRegion3DIterator it;
	for(it.SetBegin(this); it.IsValid(); it.MoveToNext())
	{
		int sn, rn, cn;
		int n_of_neighbors = 0;//we need this to detect if the edge of a region is on the edge of an image!
		bdDiscreteCoordinates3D p;
		for(g.ForCoordinates_26_Neighbors(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn()); g.Get_26_Neighbors(sn,rn,cn); )
		{
			n_of_neighbors++;

			//If the 26-neighborhood point doesn't belong to the region, this point is an edge point. 
			if(!(this->IsPointInRegion(sn,rn,cn)))
			{
				//cout<<"point  is an edge point: r,c = "<<r<<","<<c<<endl;
				list_of_egde_points.AddToLeftEnd(p(0,it.GetIndexRow(),it.GetIndexColumn()));
				n_of_neighbors = 26; //set the indicator, so that the next IF outside the FOR loop is not entered
				break;
			}
		}

		//If the number of neighbors for the point was less than 8, this point is on the edge of the image, so it IS an internal edge point!
		if(n_of_neighbors<26) list_of_egde_points.AddToLeftEnd(p(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn()));
	}

	return 1;
}


int bdRegion3D::BreakIntoDisjointRegionsOf_26_Neighborhood(bdList<bdRegion3D> &list_of_regions)
{
	list_of_regions.Reset();
	
	bdGeometry g;
	g.SetDimensions(m_original_image_size_CRS[2],m_original_image_size_CRS[1],m_original_image_size_CRS[0]);

	bdRegion3D original_region;
	original_region = (*this);

	//Initialize the disjoint region
	bdRegion3D disjoint_region;
	disjoint_region.m_original_image_size_CRS[2] = m_original_image_size_CRS[2];
	disjoint_region.m_original_image_size_CRS[1] = m_original_image_size_CRS[1];
	disjoint_region.m_original_image_size_CRS[0] = m_original_image_size_CRS[0];
	disjoint_region.m_region_value = m_region_value;

	bdList<bdDiscreteCoordinates3D> list;

	bdDiscreteCoordinates3D point;

	int start_s, start_r, start_c, sn, rn, cn;

	while(!original_region.IsEmpty())
	{
		cout<<original_region.GetNumberOfElements()<<" ";

		disjoint_region.RemoveAllPoints();
		start_s = original_region.m_slice_list.GetLeftEnd().m_slice_index;
		start_r = original_region.m_slice_list.GetLeftEnd().m_slice_region_2D.GetRowList()->GetLeftEnd().m_row_index;
		start_c = original_region.m_slice_list.GetLeftEnd().m_slice_region_2D.GetRowList()->GetLeftEnd().m_column_list.GetLeftEnd().m_column_index;

		point(start_s,start_r,start_c);
		list.AddToRightEnd(point);
		disjoint_region.AddPoint(point.S(),point.R(),point.C(),this->GetValue(point.S(),point.R(),point.C()));

		while(!list.IsEmpty())
		{
			//cout<<list.GetNumberOfElements()<<" ";
			point = list.GetLeftEnd();

			for(g.ForCoordinates_26_Neighbors(point.S(),point.R(),point.C()); g.Get_26_Neighbors(sn,rn,cn); )
			{
				if(original_region.IsPointInRegion(sn,rn,cn) && !disjoint_region.IsPointInRegion(sn,rn,cn))
				{
					//cout<<"("<<sn<<","<<rn<<","<<cn<<")";
					bdDiscreteCoordinates3D p; p(sn,rn,cn);
					list.AddToRightEnd(p);
					disjoint_region.AddPoint(sn,rn,cn,this->GetValue(sn,rn,cn));
					original_region.RemovePoint(sn,rn,cn);
				}
			}
			list.DeleteLeftEnd();
		}

		//After this part we have a disjoint region
		
		//Delete the points from original_region that exist in the disjoint region
		bdRegion3DIterator it;
		for(it.SetBegin(&disjoint_region); it.IsValid(); it.MoveToNext())
		{
			original_region.RemovePoint(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn());
			//if(!original_region.RemovePoint(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn())) cout<<"NO ";
			//if(original_region.IsPointInRegion(it.GetIndexSlice(),it.GetIndexRow(),it.GetIndexColumn())) cout<<"YES ";
		}

		list_of_regions.AddToRightEnd(disjoint_region);
	}

	if(list_of_regions.GetNumberOfElements()<=1) return 0;//there are no disjoint region
	else return list_of_regions.GetNumberOfElements(); //there are disjoint regions
}



template class BD_REGION_3D_API bdListNode<bdRegionSliceNode>;
template class BD_REGION_3D_API bdList<bdRegionSliceNode>;





