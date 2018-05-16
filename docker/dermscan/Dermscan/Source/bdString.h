/**************************************************************************
 String class

 Author: Danilo Babin
 File name: "bdString.h"
***************************************************************************/




//To build as DLL, add:" /D "BD_STRING_EXPORTS" "
// in command line build options of the project.



#if defined(_MSC_VER) //  Microsoft
	#ifdef BD_STRING_EXPORTS
		#define BD_STRING_API __declspec(dllexport)
	#else
		#define BD_STRING_API __declspec(dllimport)
	#endif
#else // consider GCC
    #ifdef BD_STRING_EXPORTS
		#define BD_STRING_API __attribute__((visibility("default")))
	#else
		#define BD_STRING_API
	#endif
#endif


#ifndef BD_STRING_DEF
	#define BD_STRING_DEF


#include <sstream>
#include <list>
#include "bdList.h"
#include "bdStandard.h"



template class BD_STRING_API std::allocator<char>;
template class BD_STRING_API std::basic_string<char>;


class BD_STRING_API bdString
{
public:

	/// Data is stored as a standard string, which is accessible publicly.
	std::string m_str;

public:

	/// Constructor.
	bdString(){};

	/// Destructor.
	~bdString(){};

	/// Get std::string.
	std::string& GetStdString(){return m_str;};
	
	/// Clear the string (make it empty).
	void Clear(){this->m_str.clear();};

	/// Check if string is empty.
	int IsEmpty(){return this->m_str.empty();};

	/// Get number of characters in the string.
	unsigned int NumberOfCharacters() {return ((unsigned int)this->m_str.length());};

    /// Assign string/characters to existing string.
    void Assign(bdString &s){this->m_str.assign(s.m_str);};
    void Assign(const char c){this->m_str.assign(&c,0,1);};
    void Assign(const char *s){this->m_str.assign(s);};
    void Assign(const std::string &s){this->m_str.assign(s);};
    
	/// Append string/characters to existing string.
	void Append(bdString &s){this->m_str.append(s.m_str);};
	void Append(const char c){this->m_str.append(&c,0,1);};
	void Append(const char *s){this->m_str.append(s);};
	void Append(const std::string &s){this->m_str.append(s);};
    void AppendNumber(int number);
    void AppendNumber(double number);
    void AppendDirSeparator();
    
    
	/// Convert number int/double to string.
	void NumberToString(int number);
	void NumberToString(double number, int number_of_decimals = 2);
	//void NumberToString(double number, unsigned int number_of_decimal_places = 3);

	/// Check if the string has a certain pattern. To do so, we use '#' to indicate any number in the original string or '%' to indicate any 
	/// non-number character in the original string. E.g. pattern "data_###%" covers strings "data_001g", "data_569x", "data_042h".
	int IsOfPattern(bdString &pattern_string);

	/// Convert string to number int/double. Success 1, fail 0.
	int StringToInt(int &output_number);
	int StringToDouble(double &output_number);
    

    /// Convert string to position coordinates. Success 1, fail 0.
    int StringToPosition(double &out_coord1, double &out_coord2, double &out_coord3);

	/// From the string all int numbers are extracted and stored to list, other characters are discarded.
	void ExtractAllIntNumbersToList(bdList<int> &list_of_extracted_int_numbers);
	void ExtractAllIntNumbersToList(std::list<int> &list_of_extracted_int_numbers);

    /// Extract strings representing hexadecimal numbers.
    void ExtractHexadecimalIntNumberStringsToList(bdList<bdString> &list_hexadecimal_number_strings);
    
    /// Extract strings representing int numbers with dots (signes '+' and '-' are NOT extracted) - use to extract serial numbers, e.g.:
    /// "26 1.2.840 UI 10008.5.1.4.1.1.4" will be extracted as: "26", "1.2.840", "10008.5.1.4.1.1.4".
    void ExtractSeriesOfIntNumberStringsWithDots(bdList<bdString> &list_number_strings_with_dots);

	/// Extracts the last int number from the string and records the rest of the string to 'remaining_string'. 
	/// E.g. for string ''data_5_00003_', remaining_string = 'data_5_' and number=3.
	void ExtractLastIntNumber(bdString &remaining_string, int &number);

	/// From a string containing integer numbers and plain words (strings), extracts the first non-number string.
	/// E.g. 'data00005_04' results in 'data'. 
	void ExtractFirstNonIntNumberString(bdString &extracted_string);

	/// From the string all int numbers are extracted and stored to list, other characters are discarded.
	void ExtractAllDoubleNumbersToList(std::list<double> &list_of_extracted_doubles);
	void ExtractAllDoubleNumbersToList(bdList<double> &list_of_extracted_doubles);
    
    /// Check if the string begins with the given input string.
    int HasPrefix(const char *prefix);
	
	/// Returns the contents of the string as a const C-style string. A null terminator is appended
	/// The C-style string is owned by this->m_str and should not be deleted.
	const char* C_String(){return this->m_str.c_str();};

	/// Replace one character with another.
	void ReplaceCharacter(char character_to_replace, char replacement_character);

	/// Convert existing path to Windows path type by replacing '/' with '\'.
	void PathToWindowsPath();

	/// Convert existing path to Linux path type by replacing '\' with '/'.
	void PathToLinuxPath();

	/// Gets the string between given characters. Useful for reading number arrays, e.g. for '[1, -5, 2]' string,
	/// start character '[' and end character ']' the output is '1, -5, 2'. Success 1, fail 0.
	int ExtractStringBetweenCharacters(char start_character, char end_character, bdString &output);

    /// Gets the string that is common for this and input string.
    /// e.g. for "data_0324" and "data_1224", the output is "data_".
    int ExtractFirstCommonStringPartWithInputString(bdString &input, bdString &output);

	/// It is assumed that the directory is of pattern: 'x:\dir1\dir2\' or '\dir1\dir2\' ('\' can be replaced with '/').
	int ExtractDirectory(bdString &output_dir);

	/// For file path of pattern: 'x:\dir\file.png' extracts: 'x:\dir', 'file' and 'png'.
	int ExtractFileNamePathParts(bdString &out_dir, bdString &out_file_name, bdString &out_extension);

	/// Inserts a number at places marked with #. E.g. n=3 for 'name###_set' converts to 'name003_set'.
	int InsertNumberToString(int n);
    
    /// Inverts letters (small become capital and capital become small).
    void InvertCapitalization();

	/// Add the tag to the beginning of the string between two 'tag_character's. 
	/// E.g. for tag_character='#', string "#my_name" and m_tag "my_tag001", the output tagged string is "#my_tag001#my_name".
	void AddTag(char tag_character, bdString &input_tag);

	/// Extracts the tag from the beginning of the string which is found between two 'tag_character's. 
	/// E.g. for tag_character='#' and string "#my_tag001#my_name", the output m_tag string is "my_tag001".
	int ExtractTag(char tag_character, bdString &output_tag);

	/// Removes the tag from the beginning of the string which is found between two 'tag_character's. 
	/// E.g. for tag_character='#' and string "#my_tag001#my_name", the output string is "my_name".
	int RemoveTag(char tag_character, bdString &output_string);
    int RemoveTag(char tag_character1, char tag_character2, bdString &output_string);

	/// The string has to be of type (tag_character='*'): "*option1**option2**option3*". The extracted string will be:
	/// "option1", "option2" and "option3".
	int ExtractOptions(char tag_character, bdList<bdString> &output_options_list);

    /// The string has to be of type (tag_character1='[', tag_character2=']'): "[option1][option2][option3]". The extracted string will be:
    /// "option1", "option2" and "option3".
    int ExtractOptions(char tag_character1, char tag_character2, bdList<bdString> &output_options_list);
    
	/// Removes the part of the string that is (consequitively) equal in both strings and gives the remainder from 'this' string
	/// as an output, e.g.: this:"str_first", input_string:"str_second" -> output:"first".
	int DifferenceFromInputString(bdString &input_string, bdString &output_difference_string);

	/// Removes the part of the string/path that is (consequitively) equal in both strings/paths and gives the remainder from 'this' string/path
	/// as an output, e.g.: this:"c:\one\two", input_string:"c:/one/three" -> output:"two".
	int DifferenceFromFilePath(bdString &input_path, bdString &output_difference_path);

	char& operator [](int r){return this->m_str[r];};
	bdString& operator ()(char *s);
    bdString& operator ()(const char *s);
	
	bdString& operator =(bdString &s);
	bdString& operator =(char *s);
	bdString& operator =(const char *s);
	
	int operator ==(bdString &s);
	int operator ==(char *s);
	int operator ==(const char *s);

	int operator !=(bdString &s);
	int operator !=(char *s);
    int operator !=(const char *s);

	int operator <(bdString &s);
	int operator >(bdString &s);

	BD_STRING_API friend std::ostream& operator << (std::ostream &o, bdString &s);
	BD_STRING_API friend std::stringstream& operator << (std::stringstream &o, bdString &s);
};


#endif