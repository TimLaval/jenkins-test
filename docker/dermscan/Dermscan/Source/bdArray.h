/**************************************************************************
 bdArray Template

 Author: Danilo Babin
 File name: "bdArray.h"
***************************************************************************/
/**************************************************************************
 Copyright (c) Danilo Babin
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even the 
 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 ***************************************************************************/



#ifndef BD_ARRAY_DEF
	#define BD_ARRAY_DEF


#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include "bdStandard.h"

using namespace std;


template <typename T> 
class bdArray
{
protected:
	unsigned int m_noel;/// Number of elements.
	T *m_pointer;/// Pointer to elements.

public:

	/// Constructor/Destructor.
	bdArray(int r, T *i);
	bdArray();
	bdArray(int r);
	~bdArray();

	int IsEmpty();/// Is the array empty.
	unsigned int GetNumberOfElements() {return m_noel;}
	T* GetPointerToData(){return m_pointer;};
	
	/// Set the size of the array (with possible initialization from existing array). 
	/// Success 1, fail 0.
	int Set(unsigned int r);
	int Set(unsigned int r, T *i);

	/// Be sure that the data stored at data_address was allocated with the right data size!!!
	int GetDataFromAddress(T *data_address, unsigned int number_of_data_elements);

	/// Sets the m_pointer to data to NULL and m_noel to zero without releasing memory
	void DetachFromData(){m_pointer=NULL; m_noel=0;}

	/// Fills in an array of set size with the input value.
	void FillInWith(T element);

	/// Fills in an array of set size with the input argument.
	void InitializeWith(T &initial_value);

	/// Reset the array to the empty state (state just after creation).
	void Reset();

	/// Indexing operator.
	T& operator [](unsigned int r) const;

	/// Assignment operator. 
	bdArray<T>& operator =(bdArray<T> &m);
	
	friend std::ostream& operator <<(ostream &o, const bdArray<T> &a);

	/// For a set this array, we copy the elements of the input array starting from
	/// the specified index in this array.
	int CopyArrayIntoThisArrayStartingFromIndex(bdArray<T> &array_to_be_copied, unsigned int start_index = 0);

	//void SaveToMatlab_M_File(char *file_name_root);

	///// Loading is not foolproof. The input array line should look like this:
	///// array_name = [1, 4, 132, 5]
	///// Don't put space before the array name!
	///// It can load an aray from a file with multiplle arrays by finding the correct array_name.
	//int LoadFromFile_IntegerValues(char *file_name, char *array_name);

	///// Loading is not foolproof. The input array line should look like this:
	///// array_name = [1.2, 4.434, 132, 5.4021]
	///// Don't put space before the array name!
	///// It can load an aray from a file with multiplle arrays by finding the correct array_name.
	//int LoadFromFile_DoubleValues(char *file_name, char *array_name);

	///// This function does not need array name, it just searches for '[' and ']' and reads the values between them.
	//int LoadFromString_DoubleValues(char *array_string);

	///// This function does not need array name, it just searches for '[' and ']' and reads the values between them.
	///// CAUTION! The array size must be set to exact size BEFORE calling this function. This method is used for loading rows
	///// of a matrix.
	//int LoadFromString_DoubleValues_ForPreSetArraySize(char *array_string);
};



/// Warning! This will work only if the elements are comparable (operators ==, <, >, <=, >= must be defined)!

template <typename T>
class bdArraySortAndSearch
{
public:
	/// Sort the elements in the ascending order.
    void SortAscending(bdArray<T> &a);

	/// Sort the elements in the ascending order. Index change is recorded to 'index_change_order'.
    void SortAscending(bdArray<T> &a, bdArray<int> &index_change_order);
	
	/// Sort the elements in the descending order. 
    void SortDescending(bdArray<T> &a);
	
	/// Sort the elements in the descending order. Index change is recorded to 'index_change_order'.
    void SortDescending(bdArray<T> &a, bdArray<int> &index_change_order);

	/// Get minimum and maximum values in the array and their indexes. 
    int MinAndMax(bdArray<T> &a, T &min, T &max, int &index_min, int &index_max);

	/// Get the index of the given input element. If none found, return -1.
    int IndexOfElement(bdArray<T> &a, T &element);

};




//template<class T> class bdArray;  // pre-declare the template class itself
//template<class T> bdArray<T> operator+ (const bdArray<T> &a1, const bdArray<T> &a2);
template<class T> std::ostream& operator<< (ostream &o, const bdArray<T> &a); 



template <class T>
T& bdArray<T>::operator [](unsigned int r) const
{ 
	assert(r<m_noel);
    return m_pointer[r];
}


template <class T>
int bdArray<T>::IsEmpty()
{
    if(m_noel==0) return 1;
    else return 0;
}


template <class T>
int bdArray<T>::Set(unsigned int r)
{
	if(r<=0) return 0;
	if(!IsEmpty()) Reset();
	m_noel = r;
	m_pointer = new T [m_noel];
	return 1;
}


template <class T>
int bdArray<T>::Set(unsigned int r, T *i)
{
	if(r<=0) return 0;
	if(i==NULL) return 0;
	if(!IsEmpty()) Reset();
	m_noel = r;
	m_pointer = new T [m_noel];
	for(unsigned int j=0; j<m_noel; j++) m_pointer[j] = i[j];
	return 1;
}


template <class T>
int bdArray<T>::GetDataFromAddress(T *data_address, unsigned int number_of_data_elements)
{
	if(data_address == NULL) return 0;
	if(!IsEmpty()) this->Reset();
	m_noel = number_of_data_elements;
	m_pointer = data_address;
	return 1;
}


template <class T>
void bdArray<T>::FillInWith(T element)
{
	for(unsigned int i=0; i<m_noel; i++) m_pointer[i] = element;
}


template <class T>
void bdArray<T>::InitializeWith(T &initial_value)
{
	for(unsigned int i=0; i<m_noel; i++) m_pointer[i] = initial_value;
}


template <class T>
bdArray<T>::bdArray()
{
    m_noel=0;
    m_pointer=NULL;
}


template <class T>
bdArray<T>::bdArray(int r)
{
	this->Set(r);
}


template <class T>
bdArray<T>::bdArray(int r, T *i)
{
	this->Set(r,i);
}


//template <class T>
//bdArray<T>::bdArray(bdArray<T> &a)
//{
//    m_noel = a.m_noel;
//
//    m_pointer = new T [m_noel];
//
//	{
//    for(int i=0; i<m_noel; i++) m_pointer[i] = a.m_pointer[i];
//	}
//
//}


template <class T>
bdArray<T>::~bdArray()
{
    this->Reset();
}


template <class T>
void bdArray<T>::Reset()
{
	delete [] m_pointer;
	m_noel = 0;
	m_pointer = NULL;
}


//template<class T>
//T& bdArray<T>::Max()
//{
//	if(this->IsEmpty());
//
//	T* pmax;
//
//	pmax = &m_pointer[0];
//
//	for(int i=1; i<m_noel; i++)
//	{
//		if(m_pointer[i]>(*pmax)) pmax = &m_pointer[i];
//	}
//
//	return (*pmax);
//}


template <class T>
bdArray<T>& bdArray<T>::operator =(bdArray<T> &a)
{
    if (&a==this) return *this;
	this->Reset();
    m_noel = a.m_noel;
    m_pointer = new T [m_noel];
	for(unsigned int i=0; i<m_noel; i++) m_pointer[i] = a.m_pointer[i];
    return *this;
}


template <class T>
ostream& operator <<(ostream &o, const bdArray<T> &a)
{ 
	for(unsigned int i=0; i<a.m_noel; i++) o<<"["<<i<<"]="<<a.m_pointer[i]<<"  ";
    return o;
}


template <class T>
int bdArray<T>::CopyArrayIntoThisArrayStartingFromIndex(bdArray<T> &array_to_be_copied, unsigned int start_index)
{
	if(this->IsEmpty()) return 0;
	if(start_index<0 || start_index>=this->GetNumberOfElements()) return 0;
	if(array_to_be_copied.IsEmpty()) return 1;
	for(unsigned int i=start_index; i<this->GetNumberOfElements() && (i-start_index)<array_to_be_copied.GetNumberOfElements(); i++)
	{
 		m_pointer[i] = array_to_be_copied[i-start_index];
	}
	return 1;
}



//----------------------------------------------------------------------------------------------------------------------




template<class T>
int bdArraySortAndSearch<T>::MinAndMax(bdArray<T> &a, T &min, T &max, int &index_min, int &index_max)
{
	if(a.IsEmpty()) return 0;
	min = a[0];
	max = a[0];
	index_min = 0;
	index_max = 0;
	for(unsigned int i=1; i<a.GetNumberOfElements(); i++)
	{
		if(a[i]>(max)) 
		{ max = a[i]; index_max = i; }
		if(a[i]<(min)) 
		{ min = a[i]; index_min = i; }
	}
	return 1;
}


template<class T>
int bdArraySortAndSearch<T>::IndexOfElement(bdArray<T> &a, T &element)
{
	if(a.IsEmpty()) return -1;
	for(unsigned int i=0; i<a.GetNumberOfElements(); i++)
	{
		if(a[i]==element) return i;
	}
	return -1;
}


template <class T>
void bdArraySortAndSearch<T>::SortAscending(bdArray<T> &a)
{
	if(a.IsEmpty()) return;

	//indicator if any changes were made in the while loop
	int changes_made = 1;
	T temp;
	//Loop while there are changes in order
	while(changes_made)
	{
		//reset indicator, initialize pointers
		changes_made = 0;
		for(unsigned int i=0; i<a.GetNumberOfElements()-1; i++)
		{
			if( a[i] > a[i+1])
			{
				changes_made = 1;
				temp = a[i+1];
				a[i+1] = a[i];
				a[i] = temp;
			}
		}
	}
}



template <class T>
void bdArraySortAndSearch<T>::SortAscending(bdArray<T> &a, bdArray<int> &index_change_order)
{
	if(a.IsEmpty()) return;

	//Set the index array size
	index_change_order.Set(a.GetNumberOfElements());
	//Set the index values as the elements of the array
	for(unsigned int i=0; i<index_change_order.GetNumberOfElements(); i++)
	{
		index_change_order[i] = i;
	}
	//indicator if any changes were made in the while loop
	int changes_made = 1;
	T temp;
	int temp_index;
	//Loop while there are changes in order
	while(changes_made)
	{
		//reset indicator, initialize pointers
		changes_made = 0;
		for(unsigned int i=0; i<a.GetNumberOfElements()-1; i++)
		{
			if( a[i] > a[i+1])
			{
				changes_made = 1;
				temp = a[i+1];
				a[i+1] = a[i];
				a[i] = temp;
				temp_index = index_change_order[i+1];
				index_change_order[i+1] = index_change_order[i];
				index_change_order[i] = temp_index;
			}
		}
	}
}



template <class T>
void bdArraySortAndSearch<T>::SortDescending(bdArray<T> &a)
{
	if(a.IsEmpty()) return;

	//indicator if any changes were made in the while loop
	int changes_made = 1;
	T temp;
	//Loop while there are changes in order
	while(changes_made)
	{
		//reset indicator, initialize pointers
		changes_made = 0;
		for(unsigned int i=0; i<a.GetNumberOfElements()-1; i++)
		{
			if( a[i] < a[i+1])
			{
				changes_made = 1;
				temp = a[i+1];
				a[i+1] = a[i];
				a[i] = temp;
			}
		}
	}
}



template <class T>
void bdArraySortAndSearch<T>::SortDescending(bdArray<T> &a, bdArray<int> &index_change_order)
{
	if(a.IsEmpty()) return;

	//Set the index array size
	index_change_order.Set(a.GetNumberOfElements());
	//Set the index values as the elements of the array
	for(unsigned int i=0; i<index_change_order.GetNumberOfElements(); i++)
	{
		index_change_order[i] = i;
	}
	//indicator if any changes were made in the while loop
	int changes_made = 1;
	T temp;
	int temp_index;
	//Loop while there are changes in order
	while(changes_made)
	{
		//reset indicator, initialize pointers
		changes_made = 0;
		for(unsigned int i=0; i<a.GetNumberOfElements()-1; i++)
		{
			if( a[i] < a[i+1])
			{
				changes_made = 1;
				temp = a[i+1];
				a[i+1] = a[i];
				a[i] = temp;
				temp_index = index_change_order[i+1];
				index_change_order[i+1] = index_change_order[i];
				index_change_order[i] = temp_index;
			}
		}
	}
}



//template <class T>
//void bdArray<T>::SaveToMatlab_M_File(char *file_name_root)
//{
//	stringstream file_name;
//	file_name<<file_name_root<<".m"<<ends;
//	
//	char *name;
//	name = file_name.str();
//
//	ofstream output_file;
//	output_file.open(name,ios::binary);
//
//	
//	output_file<<"a = [";
//	
//	for(unsigned int i=0; i<GetNumberOfElements(); i++)
//	{
//		//output_file<<((int)(m_pointer[i]));
//		output_file<<(m_pointer[i]);
//		
//		//If the condition will be fullfiled
//		if(i==GetNumberOfElements()-1)
//		{
//			output_file<<"];"<<endl;
//		}
//		else
//		{
//			output_file<<", ";
//		}
//	}
//	output_file.close();
//}
//
//
////template <class T>
////int bdArray<T>::LoadFromFile_IntegerValues(char *file_name, char *array_name)
//template<>
//int bdArray<int>::LoadFromFile_IntegerValues(char *file_name, char *array_name)
//{
//	int line_length = 100000;
//	char *line;
//	line = new char [line_length];
//
//	for(int p=0; p<line_length-1; p++) line[p] = 'a';
//
//	//Open the file
//	ifstream input_file;
//	input_file.open(file_name,ios::binary);
//    if(!input_file) 
//	{
//		cout<<"Unable to open bdArray file: "<<file_name<<endl;
//        return 0;
//    }
//	
//	//If file has multiple arrays stored in it, find the one that matches the array name
//	int array_found = 0;
//	int i; //line character counter
//	while(!array_found)
//	{	
//		//Get a line from the file
//		if(!(input_file.getline(line, line_length-1)))
//		{
//			cout<<"Error in reading bdArray file!"<<endl;
//			input_file.close();
//			return 0;
//		}
//
//		//Check the taken line
//		//int name_detected = 1; 
//		array_found = 1;
//		for(i=0; array_name[i]!=0 && array_found; i++)
//		{
//			//If the names don't match
//			if(line[i]!=array_name[i]) array_found = 0;
//		}
//	}
//
//	//Search for bracket: "["
//	while(line[i]!='[') 
//	{
//		i++;
//		if(i==line_length-1)
//		{
//			cout<<"Error in reading bdArray file! No '[' character found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//	i++;
//
//	//Search for closing bracket: "]"
//	int i2 = i;
//	while(line[i2]!=']') 
//	{
//		i2++;
//		if(i2==line_length-1)
//		{
//			cout<<"Error in reading bdArray file! No ']' character found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//	i2--;
//
//	//Count the number of elements in the array - this is done by calculating the number of commas ","
//	//Also checks the correctness of characters in between two brackets: "[" and "]"
//	int number_of_elements_of_array = 1;
//	for(int c=i; c<=i2; c++)
//	{
//		if(line[c]==',') number_of_elements_of_array++;
//		//Check correctness of characters
//		if(line[c]!=' ' && line[c]!=',' && line[c]<'0' && line[c]>'9')
//		{
//			cout<<"Error in reading bdArray file! Incorrect characters found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//
//	this->Reset();
//	this->Set(number_of_elements_of_array);
//
//	//Fill in the whole array
//	int dec = 1;
//	int value = 0;
//	int number_of_added_elements = 0;
//	for(int r=i2 ; r>=i-1; r--)
//	{
//		if(line[r]==',' || line[r]=='[')
//		{
//			//Write the calculated value in the array
//			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
//			number_of_added_elements++;
//
//			//Reset 'dec' and 'value' for a new number 
//			dec = 1;
//			value = 0;
//		}
//		//If a number is found
//		if(line[r]>='0' && line[r]<='9')
//		{
//			value += ((int)(line[r] - '0'))*dec;
//			dec = dec * 10;
//		}
//	}
//
//	delete line;
//	input_file.close();
//	return 1;
//}
//
//
//template<>
//int bdArray<double>::LoadFromFile_DoubleValues(char *file_name, char *array_name)
//{
//	int line_length = 100000;
//	char *line;
//	line = new char [line_length];
//
//	for(int p=0; p<line_length-1; p++) line[p] = 'a';
//
//	//Open the file
//	ifstream input_file;
//	input_file.open(file_name,ios::binary);
//    if(!input_file) 
//	{
//		cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Unable to open bdArray file: "<<file_name<<endl;
//        return 0;
//    }
//	
//	//If file has multiple arrays stored in it, find the one that matches the array name
//	int array_found = 0;
//	int i; //line character counter
//	while(!array_found)
//	{	
//		//Get a line from the file
//		if(!(input_file.getline(line, line_length-1)))
//		{
//			cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Error in reading bdArray file!"<<endl;
//			input_file.close();
//			return 0;
//		}
//
//		//Check the taken line
//		//int name_detected = 1; 
//		array_found = 1;
//		for(i=0; array_name[i]!=0 && array_found; i++)
//		{
//			//If the names don't match
//			if(line[i]!=array_name[i]) array_found = 0;
//		}
//	}
//
//	//Search for bracket: "["
//	while(line[i]!='[') 
//	{
//		i++;
//		if(i==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Error in reading bdArray file! No '[' character found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//	i++;
//
//	//Search for closing bracket: "]"
//	int i2 = i;
//	while(line[i2]!=']') 
//	{
//		i2++;
//		if(i2==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Error in reading bdArray file! No ']' character found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//	i2--;
//
//	//Count the number of elements in the array - this is done by calculating the number of commas ","
//	//Also checks the correctness of characters in between two brackets: "[" and "]"
//	int number_of_elements_of_array = 1;
//	for(int c=i; c<=i2; c++)
//	{
//		if(line[c]==',') number_of_elements_of_array++;
//		//Check correctness of characters
//		if(line[c]!=' ' && line[c]!=',' && line[c]!='.' && line[c]<'0' && line[c]>'9')
//		{
//			cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Error in reading bdArray file! Incorrect characters found!"<<endl;
//			input_file.close();
//			return 0;
//		}
//	}
//
//	this->Reset();
//	this->Set(number_of_elements_of_array);
//
//	//Fill in the whole array
//	double dec = 1;
//	double value = 0;
//	int number_of_added_elements = 0;
//	int is_dot_found = 0;
//	for(int r=i2 ; r>=i-1; r--)
//	{
//		if(line[r]==',' || line[r]=='[')
//		{
//			//Write the calculated value in the array
//			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
//			number_of_added_elements++;
//
//			//Reset 'dec' and 'value' for a new number 
//			dec = 1;
//			value = 0;
//
//			//Set the dot indicator to 0
//			is_dot_found = 0;
//		}
//		//If a number is found
//		if(line[r]>='0' && line[r]<='9')
//		{
//			value += ((int)(line[r] - '0'))*dec;
//			dec = dec * 10;
//		}
//		if(line[r]=='.')
//		{
//			//If this is the first dot
//			if(!is_dot_found)
//			{
//				value = value/dec;
//				is_dot_found = 1;
//				dec = 1;
//			}
//			//If this is the second dot, the error ocurred, exit with 0
//			else
//			{
//				cout<<"bdArray<T>::LoadFromFile_DoubleValues(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
//				input_file.close();
//				return 0;
//			}
//		}
//	}
//
//	delete line;
//	input_file.close();
//	return 1;
//}
//
//
//
//template<>
//int bdArray<double>::LoadFromString_DoubleValues(char *array_string)
//{
//	int line_length = 100000;
//
//	//Search for bracket: "["
//	int i=0;
//	for(; array_string[i]!='[' && i<line_length; i++)
//	{
//		if(i==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No '[' character found!"<<endl;
//			return 0;
//		}
//	}
//	i++;
//
//	//Search for closing bracket: "]"
//	int i2 = i;
//	for(; array_string[i2]!=']' && i2<line_length; i2++)
//	{
//		if(i2==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No ']' character found!"<<endl;
//			return 0;
//		}
//	}
//	i2--;
//
//	//Count the number of elements in the array - this is done by calculating the number of commas ","
//	//Also checks the correctness of characters in between two brackets: "[" and "]"
//	int number_of_elements_of_array = 1;
//	for(int c=i; c<=i2; c++)
//	{
//		if(array_string[c]==',') number_of_elements_of_array++;
//		//Check correctness of characters
//		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && array_string[c]!='-' && array_string[c]!='+' && (!(array_string[c]>='0' && array_string[c]<='9')))
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect characters found!"<<endl;
//			return 0;
//		}
//	}
//
//	this->Reset();
//	this->Set(number_of_elements_of_array);
//
//
//	//Fill in the whole array
//	double dec = 1;
//	double value = 0;
//	int number_of_added_elements = 0;
//	int is_dot_found = 0;
//	int is_sign_found = 0;
//	for(int r=i2 ; r>=i-1; r--)
//	{
//		if(array_string[r]==',' || array_string[r]=='[')
//		{
//			//Write the calculated value in the array
//			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
//			number_of_added_elements++;
//
//			//Reset 'dec' and 'value' for a new number 
//			dec = 1;
//			value = 0;
//
//			//Set the dot indicator to 0
//			is_dot_found = 0;
//
//			//Set the sign indicator to 0
//			is_sign_found = 0;
//		}
//		//Since we search backwards, the numbers and dot can be taken into account only if the sign has not yet been found!
//		if(!is_sign_found)
//		{
//			//If MINUS sign is found, make negative value, set the indicator of found sign
//			if(array_string[r]=='-')
//			{
//				value = -value;
//				is_sign_found = 1;
//			}
//			//If PLUS sign is found, just set the indicator of found sign
//			if(array_string[r]=='+') is_sign_found = 1;
//
//			//If a number is found
//			if(array_string[r]>='0' && array_string[r]<='9')
//			{
//				value += ((int)(array_string[r] - '0'))*dec;
//				dec = dec * 10;
//			}
//			if(array_string[r]=='.')
//			{
//				//If this is the first dot
//				if(!is_dot_found)
//				{
//					value = value/dec;
//					is_dot_found = 1;
//					dec = 1;
//				}
//				//If this is the second dot, the error ocurred, exit with 0
//				else
//				{
//					cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
//					return 0;
//				}
//			}
//		}
//		else
//		{
//			//If the sign is found and the current character is not ',' or '[', report an error
//			if(array_string[r]!=',' && array_string[r]!='[' && array_string[r]!=' ')
//			{
//				cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect character found before the sign!"<<endl;
//				return 0;
//			}
//		}
//	}
//
//	return 1;
//}
//
//
//
//
//
////template<class T>
////int bdArray<T>::LoadFromString_DoubleValues(char *array_string)
////{
////	int line_length = 100000;
////
////	//Search for bracket: "["
////	int i=0;
////	for(; array_string[i]!='[' && i<line_length; i++)
////	{
////		if(i==line_length-1)
////		{
////			cout<<"bdArray<T>::LoadFromString_DoubleValues(): Error in reading bdArray file! No '[' character found!"<<endl;
////			return 0;
////		}
////	}
////	i++;
////
////	//Search for closing bracket: "]"
////	int i2 = i;
////	for(; array_string[i2]!=']' && i2<line_length; i2++)
////	{
////		if(i2==line_length-1)
////		{
////			cout<<"bdArray<T>::LoadFromString_DoubleValues(): Error in reading bdArray file! No ']' character found!"<<endl;
////			return 0;
////		}
////	}
////	i2--;
////
////	//Count the number of elements in the array - this is done by calculating the number of commas ","
////	//Also checks the correctness of characters in between two brackets: "[" and "]"
////	int number_of_elements_of_array = 1;
////	for(int c=i; c<=i2; c++)
////	{
////		if(array_string[c]==',') number_of_elements_of_array++;
////		//Check correctness of characters
////		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && (!(array_string[c]>='0' && array_string[c]<='9')))
////		{
////			cout<<"bdArray<T>::LoadFromString_DoubleValues(): Error in reading bdArray file! Incorrect characters found!"<<endl;
////			return 0;
////		}
////	}
////
////	this->Reset();
////	this->Set(number_of_elements_of_array);
////
////	//Fill in the whole array
////	double dec = 1;
////	double value = 0;
////	int number_of_added_elements = 0;
////	int is_dot_found = 0;
////	for(int r=i2 ; r>=i-1; r--)
////	{
////		if(array_string[r]==',' || array_string[r]=='[')
////		{
////			//Write the calculated value in the array
////			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
////			number_of_added_elements++;
////
////			//Reset 'dec' and 'value' for a new number 
////			dec = 1;
////			value = 0;
////
////			//Set the dot indicator to 0
////			is_dot_found = 0;
////		}
////		//If a number is found
////		if(array_string[r]>='0' && array_string[r]<='9')
////		{
////			value += ((int)(array_string[r] - '0'))*dec;
////			dec = dec * 10;
////		}
////		if(array_string[r]=='.')
////		{
////			//If this is the first dot
////			if(!is_dot_found)
////			{
////				value = value/dec;
////				is_dot_found = 1;
////				dec = 1;
////			}
////			//If this is the second dot, the error ocurred, exit with 0
////			else
////			{
////				cout<<"bdArray<T>::LoadFromString_DoubleValues(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
////				return 0;
////			}
////		}
////	}
////
////	return 1;
////}
//
//
//
//template<>
//int bdArray<double>::LoadFromString_DoubleValues_ForPreSetArraySize(char *array_string)
//{
//	int line_length = 100000;
//
//	//Search for bracket: "["
//	int i=0;
//	for(; array_string[i]!='[' && i<line_length; i++)
//	{
//		if(i==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No '[' character found!"<<endl;
//			return 0;
//		}
//	}
//	i++;
//
//	//Search for closing bracket: "]"
//	int i2 = i;
//	for(; array_string[i2]!=']' && i2<line_length; i2++)
//	{
//		if(i2==line_length-1)
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No ']' character found!"<<endl;
//			return 0;
//		}
//	}
//	i2--;
//
//	//Count the number of elements in the array - this is done by calculating the number of commas ","
//	//Also checks the correctness of characters in between two brackets: "[" and "]"
//	int number_of_elements_of_array = 1;
//	for(int c=i; c<=i2; c++)
//	{
//		if(array_string[c]==',') number_of_elements_of_array++;
//		//Check correctness of characters
//		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && array_string[c]!='-' && array_string[c]!='+' && (!(array_string[c]>='0' && array_string[c]<='9')))
//		{
//			cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect characters found!"<<endl;
//			return 0;
//		}
//	}
//
//	//If the number of numbers in the line is not the same as the size of the array, return 0.
//	if(this->GetNumberOfElements()!=number_of_elements_of_array) return 0;
//
//
//	//Fill in the whole array
//	double dec = 1;
//	double value = 0;
//	int number_of_added_elements = 0;
//	int is_dot_found = 0;
//	int is_sign_found = 0;
//	for(int r=i2 ; r>=i-1; r--)
//	{
//		if(array_string[r]==',' || array_string[r]=='[')
//		{
//			//Write the calculated value in the array
//			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
//			number_of_added_elements++;
//
//			//Reset 'dec' and 'value' for a new number 
//			dec = 1;
//			value = 0;
//
//			//Set the dot indicator to 0
//			is_dot_found = 0;
//
//			//Set the sign indicator to 0
//			is_sign_found = 0;
//		}
//		//Since we search backwards, the numbers and dot can be taken into account only if the sign has not yet been found!
//		if(!is_sign_found)
//		{
//			//If MINUS sign is found, make negative value, set the indicator of found sign
//			if(array_string[r]=='-')
//			{
//				value = -value;
//				is_sign_found = 1;
//			}
//			//If PLUS sign is found, just set the indicator of found sign
//			if(array_string[r]=='+') is_sign_found = 1;
//
//			//If a number is found
//			if(array_string[r]>='0' && array_string[r]<='9')
//			{
//				value += ((int)(array_string[r] - '0'))*dec;
//				dec = dec * 10;
//			}
//			if(array_string[r]=='.')
//			{
//				//If this is the first dot
//				if(!is_dot_found)
//				{
//					value = value/dec;
//					is_dot_found = 1;
//					dec = 1;
//				}
//				//If this is the second dot, the error ocurred, exit with 0
//				else
//				{
//					cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
//					return 0;
//				}
//			}
//		}
//		else
//		{
//			//If the sign is found and the current character is not ',' or '[', report an error
//			if(array_string[r]!=',' && array_string[r]!='[' && array_string[r]!=' ')
//			{
//				cout<<"bdArray<T>::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect character found before the sign!"<<endl;
//				return 0;
//			}
//		}
//	}
//
//	return 1;
//}



#endif