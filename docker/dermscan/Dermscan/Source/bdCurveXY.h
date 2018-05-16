/**************************************************************************
 XY Curve for representing 2D signals/functions.

 Author: Danilo Babin
 File name: "bdCurveXY.h"
***************************************************************************/
/**************************************************************************
 Copyright (c) Danilo Babin
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even the 
 implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 ***************************************************************************/




//To build as DLL, add:" /D "BD_CURVE_XY_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_CURVE_XY_EXPORTS
		#define BD_CURVE_XY_API __declspec(dllexport) 
	#else
		#define BD_CURVE_XY_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_CURVE_XY_EXPORTS
		#define BD_CURVE_XY_API __attribute__((visibility("default")))
	#else
		#define BD_CURVE_XY_API
	#endif
#endif





#ifndef BD_CURVE_XY_DEF
	#define BD_CURVE_XY_DEF


#include "bdObject.h"
#include "bdArrayOfDoubles.h"


//template class BD_CURVE_XY_API bdArray<double>;


/// Used to represent curves in 2D.

class BD_CURVE_XY_API bdCurveXY : public bdDataObject
{
protected:

	/// Arrays of x and y values.
    bdArrayOfDoubles m_x, m_y;
    
    /// Labels for x and y axis.
    bdString m_x_axis_label;
    bdString m_y_axis_label;

public:

	/// Constructor.
	bdCurveXY();

	/// Destructor.
	~bdCurveXY();

	/// Check if the array is empty.
	int IsEmpty();

	/// Reset the object to its intial state (just like after construction).
	virtual void Reset();

	/// Get number of elements in the array.
	unsigned int GetNumberOfElements();
    
    /// Get Axis X label.
    bdString GetAxis_X_Label();
    
    /// Get Axis Y label.
    bdString GetAxis_Y_Label();

    /// Set Axis X label.
    void SetAxis_X_Label(bdString &label);
    void SetAxis_X_Label(const char *label);
    
    /// Set Axis Y label.
    void SetAxis_Y_Label(bdString &label);
    void SetAxis_Y_Label(const char *label);

    
	/// Get pointer to x array data.
	bdArrayOfDoubles* Get_X_Array();

	/// Get pointer to y array data.
	bdArrayOfDoubles* Get_Y_Array();

	/// Get pointer to x array data.
	double* GetPointerTo_X_Data();

	/// Get pointer to y array data.
	double* GetPointerTo_Y_Data();
    
    /// Get bound values for x axis. Bounds are [x_min,x_max].
    int X_Bounds(double &x_min, double &x_max);

    /// Get bound values for y axis. Bounds are [y_min,y_max].
    int Y_Bounds(double &y_min, double &y_max);
    
	/// Set the size of the xy array.
	virtual void Set(unsigned int r);
			
	/// For already set size of the X array, fills it in with numbers starting with 'start_value' and increasing by 'increment' (default: 0,1,2,3,...)
	void MakeStandard_X_Values(double start_value = 0, double increment = 1);
	
	/// Make an exact copy of the input xy array.
    int CopyFrom(bdCurveXY *a); // virtual int CopyFrom(bdCurveXY &a);
    int CopyFrom(bdCurveXY &a);
    
	/// Assignment operator.
	bdCurveXY& operator =(bdCurveXY &m);
    
    /// Invert Y values for the curve ('y' will become '-y').
    void Negative();

	/// Set x and y values for given input index.
	void SetValues(unsigned int index, double x_value, double y_value);

	/// Get the y value for given index.
	double& Y(unsigned int index);

    /// Get the y value for given x value.
    double Y_FromX(double x_value);
    
    /// Get the index in the array of values of given input x value (or first larger value). Assumes that x array is ASCENDING values array.
    /// If exact value is not found, returns the index of the first larger value or 0 if there is no lower value.
    unsigned int IndexOfX(double x_value);
    
	/// Get the x value for given index.
	double& X(unsigned int index);
    
    /// The whole signal is moved in y axis by the given constant value (the given value is added to each sample of the curve).
    void AddConstantValue(double y_const);
    
    /// Find the maximum Y value. Alternatively, the index of the maximum value is also given as output.
    int Max_Y_Value(double &max_y_value, unsigned int *index_of_max_y_value = NULL);
    
    /// Find the minimum Y value. Alternatively, the index of the minimum value is also given as output.
    int Min_Y_Value(double &min_y_value, unsigned int *index_of_min_y_value = NULL);

    /// Find the maximum Y convex value (means that the samples at index-1 and index+1 are smaller or equal to max). Alternatively, the index of the maximum value is also given as output.
    int Max_Y_ConvexValue(double &max_y_value, unsigned int *index_of_max_y_value = NULL);
    
    /// Find the minimum Y convex value (means that the samples at index-1 and index+1 are larger or equal to min). Alternatively, the index of the minimum value is also given as output.
    int Min_Y_ConvexValue(double &min_y_value, unsigned int *index_of_min_y_value = NULL);

    
    /// Sorts x values (corresponding y values are also sorted accordingly).
    void Sort_X();
    
	
	int RegroupXValuesToGivenNumberOfRanges(unsigned int number_of_ranges, bdCurveXY &output);
	
	/// Shift elements to rigth. The right end element is placed to left end.
	int ShiftElementsToRightCircular();

	/// Shift elements to left. The left end element is placed to right end.
	int ShiftElementsToLeftCircular();

//	/// Save to MATLAB file.
//	void SaveToMatlab_M_File(const char *file_name);//, char *file_name_root);

    /// Save to a default file type.
    void SaveToDefaultFileType(const char *file_path);
    
	/// Save to SINGLE ARRAY MATLAB file.
	void SaveToSingleArrayMatlab_M_File(const char *file_name);
    int LoadSingleArrayFile_v3(const char *file_name);
	int LoadSingleArrayFile_v2(const char *file_name);
    int LoadSingleArrayFile_v1(const char *file_name);
    
    /// General save method.
    int LoadSingleArrayFile(const char *file_name);

    
	/// Loading is not foolproof. There should be no space before the array name!
	/// It can load an aray from a file with multiplle arrays by finding the correct array_name.
	/// The input array line should look like this: array_name = [1.2, 4.434, 132, 5.4021]
	//If the array name for X is not supplied (array_nameX), x is made as standard scale starting from 0 with increment 1.
	int LoadFromFile_DoubleValues(char *file_name, char *array_nameY, char *array_nameX = NULL);

	//friend bdCurveXY operator + <>(const bdCurveXY &a1, const bdCurveXY &a2);
	//friend std::ostream& operator << <>(ostream &o, const Array<T> &a);
};


#endif