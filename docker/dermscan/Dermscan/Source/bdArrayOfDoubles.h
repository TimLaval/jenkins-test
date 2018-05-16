/**************************************************************************
 Array of double values

 Author: Danilo Babin
 File name: "bdArrayOfDoubles.h"
***************************************************************************/
/**************************************************************************
 Copyright (c) Danilo Babin
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even the 
 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 ***************************************************************************/



//To build as DLL, add:" /D "BD_ARRAY_OF_DOUBLES_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_ARRAY_OF_DOUBLES_EXPORTS
		#define BD_ARRAY_OF_DOUBLES_API __declspec(dllexport) 
	#else
		#define BD_ARRAY_OF_DOUBLES_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_ARRAY_OF_DOUBLES_EXPORTS
		#define BD_ARRAY_OF_DOUBLES_API __attribute__((visibility("default")))
	#else
		#define BD_ARRAY_OF_DOUBLES_API
	#endif
#endif



#ifndef BD_ARRAY_OF_DOUBLES_DEF
	#define BD_ARRAY_OF_DOUBLES_DEF


#include "bdArray.h"


/// Array of double values.

class BD_ARRAY_OF_DOUBLES_API bdArrayOfDoubles : public bdArray<double>
{
public:

	/// Constructor/Destructor.
	bdArrayOfDoubles();
	~bdArrayOfDoubles();

	/// Sort the array. 
	void SortAscending();
	void SortAscending(bdArray<int> &index_change_order);
	void SortDescending();
	void SortDescending(bdArray<int> &index_change_order);

	/// Get minimum and maximum values in the array and their indexes. 
	int MinAndMax(double &min, double &max, int &index_min, int &index_max);

	/// Get the index of the first occurance of input element in the array.
	/// If the element is not found, returns -1.
	int IndexOfElement(double &element);

	/// Save the array to .m file.
	void SaveToMatlab_M_File(char *file_name_root);

	/// Loading is not foolproof. The input array line should look like this:
	/// array_name = [1.2, 4.434, 132, 5.4021]
	/// Don't put space before the array name!
	/// It can load an aray from a file with multiplle arrays by finding the correct array_name.
	int LoadFromFile_DoubleValues(char *file_name, char *array_name);

	/// This function does not need array name, it just searches for '[' and ']' and reads the values between them.
	int LoadFromString_DoubleValues(char *array_string);

	/// This function does not need array name, it just searches for '[' and ']' and reads the values between them.
	/// CAUTION! The array size must be set to exact size BEFORE calling this function. This method is used for loading rows
	/// of a matrix.
	int LoadFromString_DoubleValues_ForPreSetArraySize(char *array_string);

};




#endif