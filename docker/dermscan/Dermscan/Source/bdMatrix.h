/**************************************************************************
 Template bdMatrix for matrices up to 4D data

 Author: Danilo Babin
 File name: "bdMatrix.h"
***************************************************************************/




#ifndef BD_MATRIXT_DEF
	#define BD_MATRIXT_DEF




#include "bdArray.h"



#pragma warning(disable: 4661)




template<class T>
class bdMatrix
{
protected:

	/// Number of Columns, Rows, Slices and Time sequences.
	unsigned int m_dimensions[4];

	/// We assume that any 4D data is read/written and a single continuous memory space. 
	/// 'm_pT' is an array of pointers to different 3D images in the 4D matrix. 'm_pS' is an array of RELATIVE pointers to different
	/// 2D images in the 3D matrix. 'm_pR' is an array of RELATIVE pointers to different matrix rows in the 2D matrix.
	/// Address of the data element with index (t,s,r,c) is: address = m_pT[t] + m_pS[s] + m_pR[r] + c.
	bdArray<int> m_pT, m_pS, m_pR;

	/// Pointer to data
	T* m_data;

	/// Indicates if the _data was internaly created with 'new' function (e.g. in the SetSize method). If so, the data will 
	/// have to be released with 'delete' function.
	int m_is_data_created_internally;

	/// Take the existing data and attach it to this object. Be sure that the data stored at data_address was allocated with enough size using 'new' function.
	int GetDataFromAddress(void *data_address, int n_of_columns, int n_of_rows, int n_of_slices, int n_of_time_series);

	/// Detach the array structures from the data and returns the pointer to data
	T* DetachAndGetPointerToData();


public:
	
	/// Constructor / Destructor.
	bdMatrix();
	~bdMatrix();

	/// Copy entire matrix from input matrix.
	virtual void CopyFrom(bdMatrix<T> &m);

	/// Print type and header info.
	//ostream& PrintType(ostream &o){o<<typeid(*this).name(); return o;};
	//ostream& PrintInfo(ostream &o);

	/// Check if the matrix is empty.
	int IsEmpty();

	/// Check if the matrix is equal in size as the input matrix.
	int IsEqualSizeAs(bdMatrix<T> &matrix);
	/// Check if the matrix is equal in size as the input matrix. Works also for 2D+time matrices.
	int IsEqualSizeAs_3D(bdMatrix<T> &matrix);
	int IsEqualSizeAs_2D(bdMatrix<T> &matrix);

	/// Get pointer to data elements.
	T* GetData();

	/// Get matrix dimensions.
	unsigned int* GetDimensions();
	unsigned int GetNumberOfTimeSeries();
	unsigned int GetNumberOfSlices();
	unsigned int GetNumberOfRows();
	unsigned int GetNumberOfColumns();
	unsigned int GetNumberOfDataElements();
	

	T& operator()(int r, int c);
	T& operator()(int s, int r, int c);
	T& operator()(int t, int s, int r, int c);

	/// Set the size of the matrix.
	virtual int SetSize(int t, int s, int r, int c);
	virtual int SetSize(int s, int r, int c);
	virtual int SetSize(int r, int c);

	///Set size of the matrix to match input matrix size.
	void SetSizeAs(bdMatrix<T> &matrix);

	/// Assignment operator. 
	bdMatrix<T>& operator =(bdMatrix<T> &m);

	/// Reset the object to its intial state (just after construction).
	virtual void Reset();

	/// For an existing matrix, this method will re-arrange the pointers to get a desired number of time_series 't', slices 's',
	/// rows 'r' and columns 'c'. ATTENTION! The number of total data elements MUST REMAIN THE SAME: t*s*r*c must be
	/// equal to existing total number of elements, if not the conversion will abort and the method will return FAIL (0).
	/// Use this method to convert the externally loaded 3D data to desired type of matrix.
	virtual int ConvertToSize(int t, int s, int r, int c);

	/// Fill in the whole matrix with input value.
	void FillInWith(T value);


};




template<class T>
bdMatrix<T>::bdMatrix()
{
	m_data = NULL;
	m_is_data_created_internally = 0;
	m_dimensions[0] = m_dimensions[1] = m_dimensions[2] = m_dimensions[3] = 0;
}

template<class T>
bdMatrix<T>::~bdMatrix()
{
	this->Reset();
}


template<class T>
int bdMatrix<T>::GetDataFromAddress(void *data_address, int n_of_columns, int n_of_rows, int n_of_slices, int n_of_time_series)
{
	if(n_of_columns<=0 || n_of_rows<=0 || n_of_slices<=0 || n_of_time_series<=0) return 0;//Check regularity of input data
	if(data_address==NULL) return 0;

//!!!!!!!!!!	////----- 2. If the matrix is not empty, delete it, set the properties to default values
	//if(!this->IsEmpty()) delete [] _data;

	//Set dimensions
	m_dimensions[0] = n_of_columns;
	m_dimensions[1] = n_of_rows;
	m_dimensions[2] = n_of_slices;
	m_dimensions[3] = n_of_time_series;

	m_data = (T*) data_address;//Set the pointer to data

	this->m_pR.Set(n_of_rows);
	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions[0]; }

	this->m_pS.Set(n_of_slices);
	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions[0]*m_dimensions[1]; }

	this->m_pT.Set(n_of_time_series);
	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	{ m_pT[i] = i*m_dimensions[0]*m_dimensions[1]*m_dimensions[2]; }

	return 1;
}


template<class T>
void bdMatrix<T>::CopyFrom(bdMatrix<T> &m)
{
    if(m.IsEmpty()) return;
    this->SetSize(m.GetNumberOfTimeSeries(),m.GetNumberOfSlices(),m.GetNumberOfRows(),m.GetNumberOfColumns());
    int n_of_elements = this->GetNumberOfDataElements();
    for(int i=0; i<n_of_elements; i++){ m_data[i] = m.m_data[i]; }
}


template<class T>
int bdMatrix<T>::IsEmpty()
{
	if(m_dimensions[0]==0 || m_dimensions[1]==0 || m_dimensions[2]==0 || m_dimensions[3]==0)
		return 1;
	else return 0;
}


template<class T>
int bdMatrix<T>::IsEqualSizeAs(bdMatrix<T> &matrix)
{
	if(this->IsEqualSizeAs_3D(matrix) && this->GetNumberOfTimeSeries()==matrix.GetNumberOfTimeSeries()) return 1;
	else return 0;
}


template<class T>
int bdMatrix<T>::IsEqualSizeAs_3D(bdMatrix<T> &matrix)
{
	if(this->GetNumberOfSlices()==matrix.GetNumberOfSlices() &&
		this->IsEqualSizeAs_2D(matrix)) return 1;
	else return 0;
}


template<class T>
int bdMatrix<T>::IsEqualSizeAs_2D(bdMatrix<T> &matrix)
{
	if(this->GetNumberOfRows()==matrix.GetNumberOfRows() && this->GetNumberOfColumns()==matrix.GetNumberOfColumns()) return 1;
	else return 0;
}


template<class T>
T* bdMatrix<T>::GetData()
{
	return m_data;
}


template<class T>
unsigned int* bdMatrix<T>::GetDimensions()
{
	return m_dimensions;
}

template<class T>
unsigned int bdMatrix<T>::GetNumberOfTimeSeries()
{
	return m_dimensions[3];
}


template<class T>
unsigned int bdMatrix<T>::GetNumberOfSlices()
{
	return m_dimensions[2];
}


template<class T>
unsigned int bdMatrix<T>::GetNumberOfRows()
{
	return m_dimensions[1];
}


template<class T>
unsigned int bdMatrix<T>::GetNumberOfColumns()
{
	return m_dimensions[0];
}


template<class T>
unsigned int bdMatrix<T>::GetNumberOfDataElements()
{
	return m_dimensions[0]*m_dimensions[1]*m_dimensions[2]*m_dimensions[3];
}


template<class T>
int bdMatrix<T>::SetSize(int r, int c)
{
	return this->SetSize(1,1,r,c);
}


template<class T>
int bdMatrix<T>::SetSize(int s, int r, int c)
{
	return this->SetSize(1,s,r,c);
}


template<class T>
int bdMatrix<T>::SetSize(int t, int s, int r, int c)
{
	if(c<=0 || r<=0 || s<=0 || t<=0) return 0;//Check regularity of input data

	//----- 1. If the size is already the same, just return -----
	if((this->m_dimensions[0]==c) && (this->m_dimensions[1]==r) && (this->m_dimensions[2]==s)
		&& (this->m_dimensions[3]==t)) return 1;

	this->Reset();

//!!!!!!!!!!	////----- 2. If the matrix is not empty, delete it, set the properties to default values
	//if(!this->IsEmpty()) delete [] _data;

	//Set dimensions
	m_dimensions[0] = c;
	m_dimensions[1] = r;
	m_dimensions[2] = s;
	m_dimensions[3] = t;


	m_data = new T[(c*r*s*t)] ;//Allocate new data

	m_is_data_created_internally = 1;//Set indicator that the data is internally created!

	this->m_pR.Set(r);
	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions[0]; }

	this->m_pS.Set(s);
	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions[0]*m_dimensions[1]; }

	this->m_pT.Set(t);
	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	{ m_pT[i] = i*m_dimensions[0]*m_dimensions[1]*m_dimensions[2]; }

	return 1;
}


template<class T>
void bdMatrix<T>::SetSizeAs(bdMatrix<T> &matrix)
{
	this->SetSize(matrix.GetNumberOfTimeSeries(), matrix.GetNumberOfSlices(), matrix.GetNumberOfRows(), matrix.GetNumberOfColumns());
}


template<class T>
bdMatrix<T>& bdMatrix<T>::operator =(bdMatrix<T> &m)
{
	if (&m==this) return *this;
	this->SetSizeAs(m);
	for(unsigned int i=0; i<this->GetNumberOfDataElements(); i++)
	{
		m_data[i] = m.m_data[i];
	}
	return *this;
}


template<class T>
T& bdMatrix<T>::operator()(int r, int c)
{
	return ((*this)(0,0,r,c));
}

template<class T>
T& bdMatrix<T>::operator()(int s, int r, int c)
{
	return ((*this)(0,s,r,c));
}

template<class T>
T& bdMatrix<T>::operator()(int t, int s, int r, int c)
{
	return m_data[(m_pT[t]+m_pS[s]+m_pR[r]+c)];
}


template<class T>
void bdMatrix<T>::Reset()
{
	if(this->IsEmpty()) return;
	if(m_is_data_created_internally) delete [] m_data;
	else m_data = NULL;
	m_is_data_created_internally = 0;
	m_pR.Reset();
	m_pS.Reset();
	m_pT.Reset();
	m_dimensions[0] = m_dimensions[1] = m_dimensions[2] = m_dimensions[3] = 0;
}


template<class T>
int bdMatrix<T>::ConvertToSize(int t, int s, int r, int c)
{
	//Check if the input parameters are OK.
	if(t*s*r*c != this->GetNumberOfDataElements()) return 0;
	if(t<=0 || s<=0 || r<=0 || c<=0) return 0;

	//If the size is already the same, just return success.
	if((this->m_dimensions[0]==c) && (this->m_dimensions[1]==r) && (this->m_dimensions[2]==s)
		&& (this->m_dimensions[3]==t)) return 1;

	//Reset the arrays of indexes to data.
	m_pR.Reset();
	m_pS.Reset();
	m_pT.Reset();

	//Set new dimensions
	m_dimensions[0] = c;
	m_dimensions[1] = r;
	m_dimensions[2] = s;
	m_dimensions[3] = t;

	//Set the pointers to data.
	this->m_pR.Set(r);
	for(unsigned int i=0; i<m_pR.GetNumberOfElements(); i++){ m_pR[i] = i*m_dimensions[0]; }

	this->m_pS.Set(s);
	for(unsigned int i=0; i<m_pS.GetNumberOfElements(); i++){ m_pS[i] = i*m_dimensions[0]*m_dimensions[1]; }

	this->m_pT.Set(t);
	for(unsigned int i=0; i<m_pT.GetNumberOfElements(); i++)
	{ m_pT[i] = i*m_dimensions[0]*m_dimensions[1]*m_dimensions[2]; }

	return 1;
}


template<class T>
void bdMatrix<T>::FillInWith(T value)
{
	int n_of_elements = this->GetNumberOfDataElements();
	for(int i=0; i<n_of_elements; i++){ m_data[i] = value; }
}





//template<class T>
//int bdMatrix<T>::GetVoxelMinimumAndMaximumValue(int *min, int *max)
//{
//	if(this->IsEmpty()) return 0;
//
//	*min = *max = m_data[0];
//	int n_of_elements = this->GetNumberOfDataElements();
//	for(int i=0; i<n_of_elements; i++)
//	{
//		if(*min > m_data[i]) *min = m_data[i];
//		if(*max < m_data[i]) *max = m_data[i];
//	}
//
//	return 1;
//}


#endif
