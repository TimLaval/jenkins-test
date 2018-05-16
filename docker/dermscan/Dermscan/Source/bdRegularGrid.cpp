/**************************************************************************
 Template bdRegularGrid for grids up to 4D data

 Author: Danilo Babin
 File name: "bdRegularGrid.cpp"
***************************************************************************/



#include "bdRegularGrid.h"




template<class T>
bdRegularGridT<T>::bdRegularGridT()
{
	//m_current_index_ST[0] = m_current_index_ST[1] = 0;
	m_data = NULL;
	m_is_data_created_internally = 0;
	m_dimensions_CRST[0] = m_dimensions_CRST[1] = m_dimensions_CRST[2] = m_dimensions_CRST[3] = 0;
//	m_orientation_CxCyCzRxRyRzSxSySz[0] = 0; m_orientation_CxCyCzRxRyRzSxSySz[1] = 1; m_orientation_CxCyCzRxRyRzSxSySz[2] = 0;
//	m_orientation_CxCyCzRxRyRzSxSySz[3] = 1; m_orientation_CxCyCzRxRyRzSxSySz[4] = 0; m_orientation_CxCyCzRxRyRzSxSySz[5] = 0;
//    m_orientation_CxCyCzRxRyRzSxSySz[6] = 0; m_orientation_CxCyCzRxRyRzSxSySz[7] = 0; m_orientation_CxCyCzRxRyRzSxSySz[8] = 1;

    m_orientation_CxCyCzRxRyRzSxSySz[0] = 1; m_orientation_CxCyCzRxRyRzSxSySz[1] = 0; m_orientation_CxCyCzRxRyRzSxSySz[2] = 0;
    m_orientation_CxCyCzRxRyRzSxSySz[3] = 0; m_orientation_CxCyCzRxRyRzSxSySz[4] = 1; m_orientation_CxCyCzRxRyRzSxSySz[5] = 0;
    m_orientation_CxCyCzRxRyRzSxSySz[6] = 0; m_orientation_CxCyCzRxRyRzSxSySz[7] = 0; m_orientation_CxCyCzRxRyRzSxSySz[8] = 1;

    
    m_origin_CRST[0] = m_origin_CRST[1] = m_origin_CRST[2] = m_origin_CRST[3] = 0;
	m_spacing_CRST[0] = m_spacing_CRST[1] = m_spacing_CRST[2] = m_spacing_CRST[3] = 1;
}

template<class T>
bdRegularGridT<T>::~bdRegularGridT()
{
	this->Reset();
}


template<class T>
int bdRegularGridT<T>::GetDataFromAddress(void *data_address, unsigned int n_of_columns, unsigned int n_of_rows, unsigned int n_of_slices, unsigned int n_of_time_series)
{
	if(n_of_columns<=0 || n_of_rows<=0 || n_of_slices<=0 || n_of_time_series<=0) return 0;//Check regularity of input data
	if(data_address==NULL) return 0;

//!!!!!!!!!!	////----- 2. If the grid is not empty, delete it, set the properties to default values
	//if(!this->IsEmpty()) delete [] _data;

	//Set dimensions
	m_dimensions_CRST[0] = n_of_columns;
	m_dimensions_CRST[1] = n_of_rows;
	m_dimensions_CRST[2] = n_of_slices;
	m_dimensions_CRST[3] = n_of_time_series;

	m_data = (T*) data_address;//Set the pointer to data

	this->m_pR.Set(n_of_rows);
	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions_CRST[0]; }

	this->m_pS.Set(n_of_slices);
	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]; }

	this->m_pT.Set(n_of_time_series);
	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	{ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[0]*m_dimensions_CRST[2]; }

	return 1;
}


template<class T>
void bdRegularGridT<T>::CopyFrom(bdRegularGridT<T> &m)
{
	if(m.IsEmpty()) return;
    
    bdDataObject::CopyFrom(&m);
    
	this->SetSize(m.GetNumberOfTimeSeries(),m.GetNumberOfSlices(),m.GetNumberOfRows(),m.GetNumberOfColumns());
	int n_of_elements = this->GetNumberOfDataElements();
	for(int i=0; i<n_of_elements; i++){ m_data[i] = m.m_data[i]; }
	//m_current_index_ST[0] = m.m_current_index_ST[0];
	//m_current_index_ST[1] = m.m_current_index_ST[1];
	m_orientation_CxCyCzRxRyRzSxSySz[0] = m.m_orientation_CxCyCzRxRyRzSxSySz[0]; m_orientation_CxCyCzRxRyRzSxSySz[1] = m.m_orientation_CxCyCzRxRyRzSxSySz[1];
	m_orientation_CxCyCzRxRyRzSxSySz[2] = m.m_orientation_CxCyCzRxRyRzSxSySz[2]; m_orientation_CxCyCzRxRyRzSxSySz[3] = m.m_orientation_CxCyCzRxRyRzSxSySz[3];
	m_orientation_CxCyCzRxRyRzSxSySz[4] = m.m_orientation_CxCyCzRxRyRzSxSySz[4]; m_orientation_CxCyCzRxRyRzSxSySz[5] = m.m_orientation_CxCyCzRxRyRzSxSySz[5];
    m_orientation_CxCyCzRxRyRzSxSySz[6] = m.m_orientation_CxCyCzRxRyRzSxSySz[6]; m_orientation_CxCyCzRxRyRzSxSySz[7] = m.m_orientation_CxCyCzRxRyRzSxSySz[7];
    m_orientation_CxCyCzRxRyRzSxSySz[8] = m.m_orientation_CxCyCzRxRyRzSxSySz[8];

    
    m_origin_CRST[0] = m.m_origin_CRST[0]; m_origin_CRST[1] = m.m_origin_CRST[1];
	m_origin_CRST[2] = m.m_origin_CRST[2]; m_origin_CRST[3] = m.m_origin_CRST[3];
	m_spacing_CRST[0] = m.m_spacing_CRST[0]; m_spacing_CRST[1] = m.m_spacing_CRST[1];
	m_spacing_CRST[2] = m.m_spacing_CRST[2]; m_spacing_CRST[3] = m.m_spacing_CRST[3];
}


template<class T>
bdRegularGridT<T>& bdRegularGridT<T>::operator =(bdRegularGridT<T> &m)
{
	this->CopyFrom(m);
	return(*this);
}


template<class T>
int bdRegularGridT<T>::IsEmpty()
{
	if(m_dimensions_CRST[0]==0 || m_dimensions_CRST[1]==0 || m_dimensions_CRST[2]==0 || m_dimensions_CRST[3]==0)
		return 1;
	else return 0;
}


template<class T>
int bdRegularGridT<T>::IsEqualSizeAs(bdRegularGridT<T> &grid)
{
	if(this->GetNumberOfTimeSeries()==grid.GetNumberOfTimeSeries()) return this->IsEqualSizeAs_3D(grid);
	else return 0;
}


template<class T>
int bdRegularGridT<T>::IsEqualSizeAs_3D(bdRegularGridT<T> &grid)
{
	if(this->GetNumberOfSlices()==grid.GetNumberOfSlices()) return this->IsEqualSizeAs_2D(grid);
	else return 0;
}


template<class T>
int bdRegularGridT<T>::IsEqualSizeAs_2D(bdRegularGridT<T> &grid)
{
	if(this->GetNumberOfRows()==grid.GetNumberOfRows() && this->GetNumberOfColumns()==grid.GetNumberOfColumns()) return 1;
	else return 0;
}


template<class T>
T* bdRegularGridT<T>::GetData()
{
	return m_data;
}


template<class T>
unsigned int* bdRegularGridT<T>::GetDimensions()
{
	return m_dimensions_CRST;
}

template<class T>
unsigned int bdRegularGridT<T>::GetNumberOfTimeSeries()
{
	return m_dimensions_CRST[3];
}


template<class T>
unsigned int bdRegularGridT<T>::GetNumberOfSlices()
{
	return m_dimensions_CRST[2];
}


template<class T>
unsigned int bdRegularGridT<T>::GetNumberOfRows()
{
	return m_dimensions_CRST[1];
}


template<class T>
unsigned int bdRegularGridT<T>::GetNumberOfColumns()
{
	return m_dimensions_CRST[0];
}


template<class T>
unsigned int bdRegularGridT<T>::GetNumberOfDataElements()
{
	return m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[2]*m_dimensions_CRST[3];
}


//template<class T>
//int bdRegularGridT<T>::GetNumberOfSlicesOrTimeSeries()
//{
//	if(this->GetNumberOfSlices()>1) return this->GetNumberOfSlices();
//	else return this->GetNumberOfTimeSeries();
//}

template<class T>
double* bdRegularGridT<T>::GetOrigin(){return m_origin_CRST;}

template<class T>
double bdRegularGridT<T>::GetOrigin_T(){return m_origin_CRST[3];}

template<class T>
double bdRegularGridT<T>::GetOrigin_S(){return m_origin_CRST[2];}

template<class T>
double bdRegularGridT<T>::GetOrigin_R(){return m_origin_CRST[1];}

template<class T>
double bdRegularGridT<T>::GetOrigin_C(){return m_origin_CRST[0];}

template<class T>
void bdRegularGridT<T>::SetOrigin(double t, double s, double r, double c)
{
	m_origin_CRST[0]=c; m_origin_CRST[1]=r; m_origin_CRST[2]=s; m_origin_CRST[3]=t;
}

template<class T>
double* bdRegularGridT<T>::GetSpacing(){return m_spacing_CRST;}

template<class T>
double bdRegularGridT<T>::GetSpacing_T(){return m_spacing_CRST[3];}

template<class T>
double bdRegularGridT<T>::GetSpacing_S(){return m_spacing_CRST[2];}

template<class T>
double bdRegularGridT<T>::GetSpacing_R(){return m_spacing_CRST[1];}

template<class T>
double bdRegularGridT<T>::GetSpacing_C(){return m_spacing_CRST[0];}

template<class T>
void bdRegularGridT<T>::SetSpacing(double t, double s, double r, double c)
{
	m_spacing_CRST[0]=c; m_spacing_CRST[1]=r; m_spacing_CRST[2]=s; m_spacing_CRST[3]=t;
}

template<class T>
double* bdRegularGridT<T>::GetOrientation(){return m_orientation_CxCyCzRxRyRzSxSySz;}
	
template<class T>
double bdRegularGridT<T>::GetOrientation_Cx(){return m_orientation_CxCyCzRxRyRzSxSySz[0];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Cy(){return m_orientation_CxCyCzRxRyRzSxSySz[1];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Cz(){return m_orientation_CxCyCzRxRyRzSxSySz[2];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Rx(){return m_orientation_CxCyCzRxRyRzSxSySz[3];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Ry(){return m_orientation_CxCyCzRxRyRzSxSySz[4];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Rz(){return m_orientation_CxCyCzRxRyRzSxSySz[5];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Sx(){return m_orientation_CxCyCzRxRyRzSxSySz[6];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Sy(){return m_orientation_CxCyCzRxRyRzSxSySz[7];}

template<class T>
double bdRegularGridT<T>::GetOrientation_Sz(){return m_orientation_CxCyCzRxRyRzSxSySz[8];}

template<class T>
void bdRegularGridT<T>::SetOrientationCxCyCzRxRyRz(double cx, double cy, double cz, double rx, double ry, double rz)
{
	m_orientation_CxCyCzRxRyRzSxSySz[0]=cx; m_orientation_CxCyCzRxRyRzSxSySz[1]=cy; m_orientation_CxCyCzRxRyRzSxSySz[2]=cz;
	m_orientation_CxCyCzRxRyRzSxSySz[3]=rx; m_orientation_CxCyCzRxRyRzSxSySz[4]=ry; m_orientation_CxCyCzRxRyRzSxSySz[5]=rz;
    
    
    m_orientation_CxCyCzRxRyRzSxSySz[6] = (ry * cz) - (rz * cy);
    m_orientation_CxCyCzRxRyRzSxSySz[7] = (rz * cx) - (rx * cz);
    m_orientation_CxCyCzRxRyRzSxSySz[8] = (rx * cy) - (ry * cx);
    
//    m_orientation_CxCyCzRxRyRzSxSySz[6] = - (ry * cz) + (rz * cy);
//    m_orientation_CxCyCzRxRyRzSxSySz[7] = - (rz * cx) + (rx * cz);
//    m_orientation_CxCyCzRxRyRzSxSySz[8] = - (rx * cy) + (ry * cx);
}


template<class T>
int bdRegularGridT<T>::CenterOfMass(double &out_world_x, double &out_world_y, double &out_world_z)
{
    double x_min, x_max, y_min, y_max, z_min, z_max;
    
    if(!this->Bounds_3D_WorldCoordinates(x_min, x_max, y_min, y_max, z_min, z_max)) return 0;
    
    out_world_x = (x_min + x_max ) / 2.0;
    out_world_y = (y_min + y_max ) / 2.0;
    out_world_z = (z_min + z_max ) / 2.0;
    
    return 1;
}


template<class T>
int bdRegularGridT<T>::Bounds_3D_WorldCoordinates(double &x_min, double &x_max, double &y_min, double &y_max, double &z_min, double &z_max)
{
    if(this->IsEmpty()) return 0;
    
    this->GetWorldCoordinatesForIndexes(0,0,0,z_min,y_min,x_min);
    this->GetWorldCoordinatesForIndexes(this->GetNumberOfSlices()-1, this->GetNumberOfRows()-1, this->GetNumberOfColumns()-1, z_max,y_max,x_max);
    
    // Due to various orientation possibilities, we have to check if the x_min is really smaller than x_max and if not, invert them (do also for y,z,)
    if(x_max<x_min){ double x_temp = x_max; x_max = x_min; x_min = x_temp; }
    if(y_max<y_min){ double y_temp = y_max; y_max = y_min; y_min = y_temp; }
    if(z_max<z_min){ double z_temp = z_max; z_max = z_min; z_min = z_temp; }
    
    return 1;
}


template<class T>
int bdRegularGridT<T>::SetSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(c<=0 || r<=0 || s<=0 || t<=0) return 0;//Check regularity of input data

	//----- 1. If the size is already the same, just return -----
	if((this->m_dimensions_CRST[0]==c) && (this->m_dimensions_CRST[1]==r) && (this->m_dimensions_CRST[2]==s)
		&& (this->m_dimensions_CRST[3]==t)) return 1;

//!!!!!!!!!!	////----- 2. If the grid is not empty, delete it, set the properties to default values
	//if(!this->IsEmpty()) delete [] _data;

	//Set dimensions
	m_dimensions_CRST[0] = c;
	m_dimensions_CRST[1] = r;
	m_dimensions_CRST[2] = s;
	m_dimensions_CRST[3] = t;


	m_data = new T[(c*r*s*t)] ;//Allocate new data

	m_is_data_created_internally = 1;//Set indicator that the data is internally created!

	this->m_pR.Set(r);
	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions_CRST[0]; }

	this->m_pS.Set(s);
	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]; }

	this->m_pT.Set(t);
	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	{ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[2]; }

	return 1;
}

template<class T>
int bdRegularGridT<T>::SetSizeOf3DGrid(unsigned int s, unsigned int r, unsigned int c){return SetSize(1,s,r,c);}

template<class T>
int bdRegularGridT<T>::SetSizeOf2DTimeSeries(unsigned int t, unsigned int r, unsigned int c){return SetSize(t,1,r,c);}

template<class T>
int bdRegularGridT<T>::SetSizeOf2DGrid(unsigned int r, unsigned int c){return SetSize(1,1,r,c);}

template<class T>
void bdRegularGridT<T>::SetSizeOnlyAs(bdRegularGridT<T> &grid)
{
	this->SetSize(grid.GetNumberOfTimeSeries(), grid.GetNumberOfSlices(), grid.GetNumberOfRows(), grid.GetNumberOfColumns());
}

template<class T>
void bdRegularGridT<T>::SetSizeAndPropertiesAs(bdRegularGridT<T> &grid)
{
	this->SetSize(grid.GetNumberOfTimeSeries(), grid.GetNumberOfSlices(), grid.GetNumberOfRows(), grid.GetNumberOfColumns());
	this->SetSpacing(grid.GetSpacing_T(),grid.GetSpacing_S(),grid.GetSpacing_R(),grid.GetSpacing_C());
	this->SetOrigin(grid.GetOrigin_T(),grid.GetOrigin_S(),grid.GetOrigin_R(),grid.GetOrigin_C());
	this->SetOrientationCxCyCzRxRyRz(grid.GetOrientation_Cx(), grid.GetOrientation_Cy(), grid.GetOrientation_Cz(), 
		grid.GetOrientation_Rx(), grid.GetOrientation_Ry(), grid.GetOrientation_Rz());
}

template<class T>
void bdRegularGridT<T>::SetVisualizationPropertiesToMatchInput(bdRegularGridT<T> &grid)
{
	this->SetSpacing(grid.GetSpacing_T(),grid.GetSpacing_S(),grid.GetSpacing_R(),grid.GetSpacing_C());
	this->SetOrigin(grid.GetOrigin_T(),grid.GetOrigin_S(),grid.GetOrigin_R(),grid.GetOrigin_C());
	this->SetOrientationCxCyCzRxRyRz(grid.GetOrientation_Cx(), grid.GetOrientation_Cy(), grid.GetOrientation_Cz(), 
		grid.GetOrientation_Rx(), grid.GetOrientation_Ry(), grid.GetOrientation_Rz());
}

template<class T>
T& bdRegularGridT<T>::RC(unsigned int r, unsigned int c)
{
	return TSRC(0,0,r,c);
	//return TSRC(m_current_index_ST[1],m_current_index_ST[0],r,c);
}


template<class T>
T& bdRegularGridT<T>::SRC(unsigned int s, unsigned int r, unsigned int c)
{
	return TSRC(0,s,r,c);
	//return _data[(s*_dimensions_CRST[0]*_dimensions_CRST[1] + r*_dimensions_CRST[0] + c)];
}

template<class T>
T& bdRegularGridT<T>::TRC(unsigned int t, unsigned int r, unsigned int c)
{
	return TSRC(t,0,r,c);
	//return _data[(s*_dimensions_CRST[0]*_dimensions_CRST[1] + r*_dimensions_CRST[0] + c)];
}

template<class T>
T& bdRegularGridT<T>::TSRC(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	return m_data[(m_pT[t]+m_pS[s]+m_pR[r]+c)];
}


//template<class T>
//bdRegularGridT<T>& bdRegularGridT<T>::SetS(int s)
//{
//	if(s>=0 && s<m_dimensions_CRST[2]) m_current_index_ST[0] = s;
//	//else _current_index_ST[0] = 0;
//	return (*this);
//}
//
//
//template<class T>
//bdRegularGridT<T>& bdRegularGridT<T>::SetT(int t)
//{
//	if(t>=0 && t<m_dimensions_CRST[3]) m_current_index_ST[1] = t;
//	//else _current_index_ST[1] = 0;
//	return (*this);
//}

template<class T>
T& bdRegularGridT<T>::operator()(unsigned int r, unsigned int c){return (RC(r,c));}

template<class T>
T& bdRegularGridT<T>::operator()(unsigned int s, unsigned int r, unsigned int c){return (SRC(s,r,c));}

template<class T>
T& bdRegularGridT<T>::operator()(unsigned int t, unsigned int s, unsigned int r, unsigned int c){return (TSRC(t,s,r,c));}


template<class T>
int bdRegularGridT<T>::GetIndexesForWorldCoordinates(double w_t, double w_z, double w_y, double w_x, int &out_t, int &out_s, int &out_r, int &out_c)
{
    // world_coordinates = origin  + Orientation_matrix * (index * spacing).
    
    //Orientation is calculated as follows:
    //   [ xr yr zr 1 ]
    // M=[ xc yc zc 1 ]
    //   [ xs ys zs 1 ]
    //   [ 0  0  0  1 ]
    // xs = (yr * zc) - (zr * yc)
    // ys = (zr * xc) - (xr * zc)
    // zs = (xr * yc) - (yr * xc)

    double xc = m_orientation_CxCyCzRxRyRzSxSySz[0];
    double yc = m_orientation_CxCyCzRxRyRzSxSySz[1];
    double zc = m_orientation_CxCyCzRxRyRzSxSySz[2];
    double xr = m_orientation_CxCyCzRxRyRzSxSySz[3];
    double yr = m_orientation_CxCyCzRxRyRzSxSySz[4];
    double zr = m_orientation_CxCyCzRxRyRzSxSySz[5];
    double xs = m_orientation_CxCyCzRxRyRzSxSySz[6];
    double ys = m_orientation_CxCyCzRxRyRzSxSySz[7];
    double zs = m_orientation_CxCyCzRxRyRzSxSySz[8];
    

    // Since M is orthogonal matrix, to get the inverse rotation, we just need to transpose it:
    //       [ xr xc xs 1 ]
    // M^T = [ yr yc ys 1 ]
    //       [ zr zc zs 1 ]
    //       [ 0  0  0  1 ]
    //
    //
    //       [ lx ]   [ xr*lx + xc*ly + xs*lz ]
    // M^T * [ ly ] = [ yr*lx + yc*ly + ys*lz ]
    //       [ lz ]   [ zr*lx + zc*ly + zs*lz ]
    //       [  1 ]   [            1          ]
    // where lx = w_x - origin_x, ly = w_y - origin_y, lz = w_z - origin_z.

    double lx = w_x - m_origin_CRST[0];
    double ly = w_y - m_origin_CRST[1];
    double lz = w_z - m_origin_CRST[2];
    
    double si_x = xr*lx + xc*ly + xs*lz;
    double si_y = yr*lx + yc*ly + ys*lz;
    double si_z = zr*lx + zc*ly + zs*lz;
    
    // Finally, the image index is calculated from si values:
    out_c = (int) (si_x / m_spacing_CRST[0]);
    out_r = (int) (si_y / m_spacing_CRST[1]);
    out_s = (int) (si_z / m_spacing_CRST[2]);
	out_t = (int) ((w_t - m_origin_CRST[3]) / m_spacing_CRST[3]);

    //cout<<"("<<w_z<<", "<<w_y<<", "<<w_x<<") -> ["<<out_s<<", "<<out_r<<", "<<out_c<<"]"<<endl;
    
    // Make absolute values out of the calculated ones - it could be that visualization planes have different positive/negative directions.
    if(out_s<0) out_s = -out_s;
    if(out_r<0) out_r = -out_r;
    if(out_c<0) out_c = -out_c;
    
//    cout<<out_t<<","<<out_s<<","<<out_r<<","<<out_c<<endl;
    
    if(out_t>= ((int)this->GetNumberOfTimeSeries())) return 0;
    if(out_s>= ((int)this->GetNumberOfSlices())) return 0;
    if(out_r>= ((int)this->GetNumberOfRows())) return 0;
    if(out_c>= ((int)this->GetNumberOfColumns())) return 0;

    //Without orientation taken into account, the equations are:
//  out_s = (int)((w_s - m_origin_CRST[2]) / m_spacing_CRST[2]);
//  out_r = (int)((w_r - m_origin_CRST[1]) / m_spacing_CRST[1]);
//  out_c = (int)((w_c - m_origin_CRST[0]) / m_spacing_CRST[0]);
//	if(out_t<0 || out_t>= ((int)this->GetNumberOfTimeSeries())) return 0;
//	if(out_s<0 || out_s>= ((int)this->GetNumberOfSlices())) return 0;
//	if(out_r<0 || out_r>= ((int)this->GetNumberOfRows())) return 0;
//	if(out_c<0 || out_c>= ((int)this->GetNumberOfColumns())) return 0;

	return 1;
}


template<class T>
int bdRegularGridT<T>::GetIndexesForWorldCoordinates(double w_s, double w_r, double w_c, int &out_s, int &out_r, int &out_c)
{
	int t;
	return this->GetIndexesForWorldCoordinates(m_origin_CRST[3], w_s, w_r, w_c, t, out_s, out_r, out_c);
}


template<class T>
int bdRegularGridT<T>::GetWorldCoordinatesForIndexes(unsigned int t, unsigned int s, unsigned int r, unsigned int c, double &out_w_t, double &out_w_z, double &out_w_y, double &out_w_x)
{
    // world_coordinates = origin  + Orientation_matrix * (index * spacing).
    //Orientation is calculated as follows:
    //   [ xr yr zr Tx ]
    // M=[ xc yc zc Ty ]
    //   [ xs ys zs Tz ]
    //   [ 0  0  0  1  ]
    // xs = (yr * zc) - (zr * yc)
    // ys = (zr * xc) - (xr * zc)
    // zs = (xr * yc) - (yr * xc)
    
    if(t>= (this->GetNumberOfTimeSeries())) return 0;
    if(s>= (this->GetNumberOfSlices())) return 0;
    if(r>= (this->GetNumberOfRows())) return 0;
    if(c>= (this->GetNumberOfColumns())) return 0;

    
    double xc = m_orientation_CxCyCzRxRyRzSxSySz[0];
    double yc = m_orientation_CxCyCzRxRyRzSxSySz[1];
    double zc = m_orientation_CxCyCzRxRyRzSxSySz[2];
    double xr = m_orientation_CxCyCzRxRyRzSxSySz[3];
    double yr = m_orientation_CxCyCzRxRyRzSxSySz[4];
    double zr = m_orientation_CxCyCzRxRyRzSxSySz[5];
    double xs = m_orientation_CxCyCzRxRyRzSxSySz[6];
    double ys = m_orientation_CxCyCzRxRyRzSxSySz[7];
    double zs = m_orientation_CxCyCzRxRyRzSxSySz[8];
    
    
    // World Position:
    //       [ xr yr zr 0 ]   [lx]   [ xr*lx + yr*ly + zr*lz ]
    // R   = [ xc yc zc 0 ] * [ly] = [ xc*lx + yc*ly + zc*lz ]
    //       [ xs ys zs 0 ]   [lz]   [ xs*lx + ys*ly + zs*lz ]
    //       [ 0  0  0  1 ]   [ 0]   [          0            ]
    //
    // where lx = c * spacing_x, ly = r * spacing_y, lz = s * spacing_z.
    
    double lx = c * m_spacing_CRST[0];
    double ly = r * m_spacing_CRST[1];
    double lz = s * m_spacing_CRST[2];
    
    
    out_w_x = m_origin_CRST[0] + xr*lx + yr*ly + zr*lz;
    out_w_y = m_origin_CRST[1] + xc*lx + yc*ly + zc*lz;
    out_w_z = m_origin_CRST[2] + xs*lx + ys*ly + zs*lz;
    out_w_t = m_origin_CRST[3] + t * m_spacing_CRST[3];
    
    
    
    cout<<"["<<s<<", "<<r<<", "<<c<<"] -> ("<<out_w_z<<", "<<out_w_y<<", "<<out_w_x<<")"<<endl;


    return 1;
}


template<class T>
int bdRegularGridT<T>::GetWorldCoordinatesForIndexes(unsigned int s, unsigned int r, unsigned int c, double &out_w_z, double &out_w_y, double &out_w_x)
{
    double w_t;
    return this->GetWorldCoordinatesForIndexes(0, s, r, c, w_t, out_w_z, out_w_y, out_w_x);
}



template<class T>
void bdRegularGridT<T>::Reset()
{
	if(this->IsEmpty()) return;
	//m_current_index_ST[0] = m_current_index_ST[1] = 0;
	if(m_is_data_created_internally) delete [] m_data;
	else m_data = NULL;
	m_is_data_created_internally = 0;
	m_pR.Reset();
	m_pS.Reset();
	m_pT.Reset();
	m_dimensions_CRST[0] = m_dimensions_CRST[1] = m_dimensions_CRST[2] = m_dimensions_CRST[3] = 0;
//	m_orientation_CxCyCzRxRyRzSxSySz[0] = 0; m_orientation_CxCyCzRxRyRzSxSySz[1] = 1; m_orientation_CxCyCzRxRyRzSxSySz[2] = 0;
//	m_orientation_CxCyCzRxRyRzSxSySz[3] = 1; m_orientation_CxCyCzRxRyRzSxSySz[4] = 0; m_orientation_CxCyCzRxRyRzSxSySz[5] = 0;
//    m_orientation_CxCyCzRxRyRzSxSySz[6] = 0; m_orientation_CxCyCzRxRyRzSxSySz[7] = 0; m_orientation_CxCyCzRxRyRzSxSySz[8] = 1;

    m_orientation_CxCyCzRxRyRzSxSySz[0] = 1; m_orientation_CxCyCzRxRyRzSxSySz[1] = 0; m_orientation_CxCyCzRxRyRzSxSySz[2] = 0;
    m_orientation_CxCyCzRxRyRzSxSySz[3] = 0; m_orientation_CxCyCzRxRyRzSxSySz[4] = 1; m_orientation_CxCyCzRxRyRzSxSySz[5] = 0;
    m_orientation_CxCyCzRxRyRzSxSySz[6] = 0; m_orientation_CxCyCzRxRyRzSxSySz[7] = 0; m_orientation_CxCyCzRxRyRzSxSySz[8] = 1;

    
    m_origin_CRST[0] = m_origin_CRST[1] = m_origin_CRST[2] = m_origin_CRST[3] = 0;
	m_spacing_CRST[0] = m_spacing_CRST[1] = m_spacing_CRST[2] = m_spacing_CRST[3] = 1;
}



template<class T>
int bdRegularGridT<T>::ConvertToSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	//Check if the input parameters are OK.
	if(t*s*r*c != this->GetNumberOfDataElements()) return 0;
	if(t<=0 || s<=0 || r<=0 || c<=0) return 0;

	//If the size is already the same, just return success.
	if((this->m_dimensions_CRST[0]==c) && (this->m_dimensions_CRST[1]==r) && (this->m_dimensions_CRST[2]==s)
		&& (this->m_dimensions_CRST[3]==t)) return 1;

	//Reset the arrays of indexes to data.
	m_pR.Reset();
	m_pS.Reset();
	m_pT.Reset();

	//Set new dimensions
	m_dimensions_CRST[0] = c;
	m_dimensions_CRST[1] = r;
	m_dimensions_CRST[2] = s;
	m_dimensions_CRST[3] = t;

	//if(is_S0Trange_SnTrange)
	//{
	//	//Set the pointers to data.
	//	this->m_pR.Set(r);
	//	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions_CRST[0]; }

	//	this->m_pS.Set(s);
	//	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++)
	//		{ m_pS[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[3]; }
	//	//for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]; }

	//	this->m_pT.Set(t);
	//	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++){ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]; }
	//	//for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	//	//{ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[0]*m_dimensions_CRST[2]; }
	//}
	//else
	//{
		//Set the pointers to data.
		this->m_pR.Set(r);
		for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions_CRST[0]; }

		this->m_pS.Set(s);
		for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]; }

		this->m_pT.Set(t);
		for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
		//{ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[0]*m_dimensions_CRST[2]; }
		{ m_pT[i] = i*m_dimensions_CRST[0]*m_dimensions_CRST[1]*m_dimensions_CRST[2]; }
	//}

	return 1;
}


template<class T>
void bdRegularGridT<T>::FillInWith(T value)
{
	int n_of_elements = this->GetNumberOfDataElements();
	for(int i=0; i<n_of_elements; i++){ m_data[i] = value; }
}


template<class T>
void bdRegularGridT<T>::FillInSliceWith(unsigned int t, unsigned int s, T value)
{
	for(unsigned int r=0; r<this->GetNumberOfRows(); r++)
	{
		for(unsigned int c=0; c<this->GetNumberOfColumns(); c++)
		{
			(*this)(t,s,r,c) = value;
		}
	}
}


template<class T>
bdString bdRegularGridT<T>::PrintSpacing()
{
	bdString s, s2;
	s2.NumberToString(this->GetSpacing_C(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetSpacing_R(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetSpacing_S(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetSpacing_T(),2); s.Append(s2);
	return s;
}


template<class T>
bdString bdRegularGridT<T>::PrintOrigin()
{
	bdString s, s2;
	s2.NumberToString(this->GetOrigin_C(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetOrigin_R(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetOrigin_S(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetOrigin_T(),2); s.Append(s2); 
	return s;
}


template<class T>
bdString bdRegularGridT<T>::PrintDimensions()
{
	bdString s, s2;
	s2.NumberToString(this->GetNumberOfColumns(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetNumberOfRows(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetNumberOfSlices(),2); s.Append(s2); s.Append(", ");
	s2.NumberToString(this->GetNumberOfTimeSeries(),2); s.Append(s2);
	return s;
}


template<class T>
bdString bdRegularGridT<T>::PrintOrientation()
{
	bdString s, s2;
    
    
    s2.NumberToString(this->GetOrientation_Rx(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Ry(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Rz(),15); s.Append(s2); s.Append(" / ");
    
    s2.NumberToString(this->GetOrientation_Cx(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Cy(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Cz(),15); s.Append(s2); s.Append(" (");
    
    s2.NumberToString(this->GetOrientation_Sx(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Sy(),15); s.Append(s2); s.Append(", ");
    s2.NumberToString(this->GetOrientation_Sz(),15); s.Append(s2); s.Append(")");

    
	
//    if(this->GetOrientation_Rx()<0 && this->GetOrientation_Rx()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Rx()>0 && this->GetOrientation_Rx()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Rx(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Ry()<0 && this->GetOrientation_Ry()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Ry()>0 && this->GetOrientation_Ry()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Ry(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Rz()<0 && this->GetOrientation_Rz()>-0.01) { s.Append("-0.00 / "); }
//    else
//    {
//        if(this->GetOrientation_Rz()>0 && this->GetOrientation_Rz()<0.01) { s.Append("+0.00 / "); }
//        else { s2.NumberToString(this->GetOrientation_Rz(),2); s.Append(s2); s.Append(" / "); }
//    }
//    if(this->GetOrientation_Cx()<0 && this->GetOrientation_Cx()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Cx()>0 && this->GetOrientation_Cx()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Cx(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Cy()<0 && this->GetOrientation_Cy()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Cy()>0 && this->GetOrientation_Cy()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Cy(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Cz()<0 && this->GetOrientation_Cz()>-0.01) { s.Append("-0.00 ("); }
//    else
//    {
//        if(this->GetOrientation_Cz()>0 && this->GetOrientation_Cz()<0.01) { s.Append("+0.00 ("); }
//        else { s2.NumberToString(this->GetOrientation_Cz(),2); s.Append(s2); s.Append(" ("); }
//    }
//    
//    if(this->GetOrientation_Sx()<0 && this->GetOrientation_Sx()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Sx()>0 && this->GetOrientation_Sx()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Sx(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Sy()<0 && this->GetOrientation_Sy()>-0.01) { s.Append("-0.00, "); }
//    else
//    {
//        if(this->GetOrientation_Sy()>0 && this->GetOrientation_Sy()<0.01) { s.Append("+0.00, "); }
//        else { s2.NumberToString(this->GetOrientation_Sy(),2); s.Append(s2); s.Append(", "); }
//    }
//    if(this->GetOrientation_Sz()<0 && this->GetOrientation_Sz()>-0.01) { s.Append("-0.00)"); }
//    else
//    {
//        if(this->GetOrientation_Sz()>0 && this->GetOrientation_Sz()<0.01) { s.Append("+0.00)"); }
//        else { s2.NumberToString(this->GetOrientation_Sz(),2); s.Append(s2); s.Append(")"); }
//    }

 
	return s;
}




template class BD_REGULAR_GRID_API bdRegularGridT<unsigned char>;
template class BD_REGULAR_GRID_API bdRegularGridT<unsigned short>;






