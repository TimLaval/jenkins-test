/**************************************************************************
 Class of 2D Image (slice) tags

 Author: Danilo Babin
 File name: "bdImageSliceTag.h"
***************************************************************************/


// To build as DLL, add:" /D "BD_IMAGE_SLICE_TAG_EXPORTS" "
// in command line build options of the project.


#if defined(_MSC_VER) //  Microsoft 
	#ifdef BD_IMAGE_SLICE_TAG_EXPORTS
		#define BD_IMAGE_SLICE_TAG_API __declspec(dllexport) 
	#else
		#define BD_IMAGE_SLICE_TAG_API __declspec(dllimport) 
	#endif
#else // GCC
	#ifdef BD_IMAGE_SLICE_TAG_EXPORTS
		#define BD_IMAGE_SLICE_TAG_API __attribute__((visibility("default")))
	#else
		#define BD_IMAGE_SLICE_TAG_API
	#endif
#endif



#ifndef BD_IMAGE_SLICE_TAG_DEF
	#define BD_IMAGE_SLICE_TAG_DEF



#include <iostream>
#include "bdArray.h"
#include "bdString.h"


template class BD_IMAGE_SLICE_TAG_API bdArray<bdString>;
template class BD_IMAGE_SLICE_TAG_API bdArray<double>;


class BD_IMAGE_SLICE_TAG_API bdImageSliceTagNumber
{
public:
    bdString m_tag_number1;// 4 character hexadecimal number e.g. "0008"
    bdString m_tag_number2;// 4 character hexadecimal number e.g. "0013"
    bdString m_tag_description;//e.g. "Creation Time"
};


class BD_IMAGE_SLICE_TAG_API bdImageSliceTagValue
{
public:
    bdString m_value;// actual value of the tag.
    

    /// Get value as a double number. If element is NAN return fail 0, else success 1.
    int GetValueAsDouble(double &output);

    /// Get value as an array of double numbers. If no numbers are found return 0, else success returns number of added values.
    int GetValuesAsArrayOfDoubles(bdArray<double> &output_values);
};



/// Class of Image Slice Tags. Contains regular and proprietary DICOM fields.

class BD_IMAGE_SLICE_TAG_API bdImageSliceTag
{
public:
    
//	/// Arrays containing names of tags.
//	static bdArray<bdString> m_dicom_tag_names;
//	static bdArray<bdString> m_siemens_tag_names;

    
    /// Arrays containing names of tags.
    static bdArray<bdImageSliceTagNumber> m_dicom_tag_names;

    
//	/// Arrays containing tag values
//	bdArray<double> m_dicom_tag_values;
//	bdArray<double> m_siemens_tag_values;


    /// Arrays containing tag values
    bdArray<bdImageSliceTagValue> m_dicom_tag_values;
    
    
	/// Indicator to say if the loaded header is Siemens dicom.
	int m_is_siemens_dicom;


	/// Constructor.
	bdImageSliceTag();

	/// Destructor.
	~bdImageSliceTag();


	/// For input tag name (any: dicom, siemens...) return the pointer to value. If none found, return NULL. 
	bdImageSliceTagValue* Value(bdString &tag_number1, bdString &tag_number2);
	bdImageSliceTagValue* Value(const char *tag_number1, const char *tag_number2);

    /// For input tag name (any: dicom, siemens...) return the string value. If none found, return empty string.
    bdString* ValueString(bdString &tag_number1, bdString &tag_number2);
    bdString* ValueString(const char *tag_number1, const char *tag_number2);

    /// For input tag name (any: dicom, siemens...) return an array of containing double values. If none found, return empty array.
    int ValueDoubles(bdString &tag_number1, bdString &tag_number2, bdArray<double> &output_values);
    int ValueDoubles(const char *tag_number1, const char *tag_number2, bdArray<double> &output_values);
    
    
//	/// Get DICOM tag value for given index;
//	double& ValueDicom(int index);

	/// For input DICOM tag name, returns its index. If not found, return -1.
	int IndexOfDicomTagNumber(bdString &tag_number1, bdString &tag_number2);
	int IndexOfDicomTagNumber(const char *tag_number1, const char *tag_number2);

//	/// For input DICOM hexadecimal tag number (e.g. "0008,0001"), returns its index. If not found, return -1.
//	int IndexOfDicomTagHexadecimal(bdString &tag_hexadecimal);
//	int IndexOfDicomTagHexadecimal(const char *tag_hexadecimal);
//
//	/// Get SIEMENS tag value for given index;
//	double& ValueSiemens(int index);
//
//	/// For input SIEMENS tag name, returns its index. If not found, return -1.
//	int IndexOfSiemensTagName(bdString &tag_name);
//	int IndexOfSiemensTagName(const char *tag_name);

	/// Print out.
	void Print(bdString &s);

	/// Saving (printing) to a string to be used when saveing to a file.
	void SaveToString(bdString &s);

	/// Load saved tags from a string. Each tag should be bordered with '<' and '>'.
	/// E.g.: <(0018,1060)[Trigger Time]{0}><(0020,0013)[Image Number]{5}> .
    /// Warning: Loading is done based ONLY on TAG NUMBERS (description is NOT used), e.g. string <(0020,0013){5}> will load correctly.
	int LoadFromSavedString(bdString &s);

	/// Friend stream operators.
	BD_IMAGE_SLICE_TAG_API friend std::ostream& operator <<(std::ostream &o, const bdImageSliceTag &t);
	BD_IMAGE_SLICE_TAG_API friend std::stringstream& operator <<(std::stringstream &o, const bdImageSliceTag &t);
};







#endif