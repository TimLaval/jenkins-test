/**************************************************************************
 String class implementation

 Author: Danilo Babin
 File name: "bdString.cpp"
***************************************************************************/



#include "bdString.h"



//bdString::bdString(char *s)
//{
//	for(unsigned int i=0; i<bdString_MAX_STRING_SIZE && s[i]!=0; i++) list.AddToRightEnd(s[i]);
//}
//
//
//bdString::bdString(const char *s)
//{
//	for(unsigned int i=0; i<bdString_MAX_STRING_SIZE && s[i]!=0; i++) list.AddToRightEnd(((char&)(s[i])));
//}



//void bdString::Append(char *s, unsigned int number_of_characters)
//{
//	if(s==NULL) return;
//	unsigned int n;
//	if(number_of_characters>bdString_MAX_STRING_SIZE) n = bdString_MAX_STRING_SIZE;
//	else n = number_of_characters;
//
//	for(unsigned int i=0; i<n && s[i]!=0; i++) list.AddToRightEnd(s[i]);
//}


//void bdString::Append(const std::string &s)
//{
//	// Fill string up with random lower case characters
//	char c;
//	for (int i=0; i<((int)s.length()); i++) 
//	{
//		c = (char) s[i];
//		this->list.AddToRightEnd(c);
//	}
//}
//
//
//void bdString::NumberToString(int number)
//{
//	this->Delete();
//
//	if(number==0)
//	{
//		char ch = '0';
//		this->list.AddToLeftEnd(ch);
//		return;
//	}
//
//	int d = number;
//	if(d<0) d = -d;
//	while(d!=0)
//	{
//		int c = d%10;
//		char ch = '0' + ((char) c);
//		this->list.AddToLeftEnd(ch);
//		d = d/10;
//	}
//
//	if(number<0)
//	{
//		char ch = '-';
//		this->list.AddToLeftEnd(ch);
//	}
//}
//
//
//void bdString::NumberToString(double number, unsigned int number_of_decimal_places)
//{
//	this->Delete();
//	unsigned int n;
//	if(number_of_decimal_places>5) n = 5;
//	else n = number_of_decimal_places;
//	//Determine the base starting value
//	int b = 1;
//	for(unsigned int i=0; i<n; i++) b = b*10;
//
//	double d = number;
//	int before_dot, after_dot;
//	before_dot = (int)d;
//	after_dot = (int)( (d-((double)before_dot)) * ((double) b) );
//	if(after_dot<0) after_dot = -after_dot;
//	bdString int_string;
//	int_string.NumberToString(before_dot);
//	this->Append(int_string);
//	this->Append(".");
//	int_string.NumberToString(after_dot);
//	this->Append(int_string);
//
//	//Check if the last characters are '0' and if so, delete them.
//	int is_non_zero_char_found = 0;
//	while(!is_non_zero_char_found)
//	{
//		char ch;
//		ch = this->list.GetRightEnd();
//		if(ch=='0' || ch=='.') this->list.DeleteRightEnd();
//		else is_non_zero_char_found = 1;
//	}
//}



void bdString::AppendNumber(int number)
{
    bdString bds; bds.NumberToString(number);
    this->Append(bds);
}


void bdString::AppendNumber(double number)
{
    bdString bds; bds.NumberToString(number);
    this->Append(bds);    
}


void bdString::AppendDirSeparator()
{
#if defined(_MSC_VER) //  Microsoft
    this->Append("\\");
#else // consider GCC
    this->Append("/");
#endif

}


void bdString::NumberToString(int number)
{
	this->m_str.clear();
	std::ostringstream os;
	os<<number;
	this->m_str = os.str();
}


void bdString::NumberToString(double number, int number_of_decimals)
{
    this->Clear();

    // Record the sign and make input positive.
    if(number<0)
    {
        this->Append("-");
        number = -number;
    }
    
    // Convert double to int and from it to string.
    bdString int_string_to_add;
    for(unsigned int i=0; i<number_of_decimals+1; i++)
    {
        int n = ((int) number);
        
        bdString int_string;
        int_string.NumberToString(n);
        int_string_to_add.Append(int_string);
        
        if(n>0)
        {
            this->Append(int_string_to_add);
            int_string_to_add.Clear();
        }
        
        number = (number - ((double)n)) * 10.0;
        
        if(i==0)
        {
            this->Append(int_string_to_add);
            int_string_to_add.Assign(".");
        }
    }
    
    
    
    
    
//	this->m_str.clear();
//	std::ostringstream os;
//	os<<number;
//	std::string s = os.str();
//	int j=0;
//	int is_dot_found = 0;
//	for(unsigned int i=0; i<s.size() && j<number_of_decimals; i++)
//	{
//		this->m_str.push_back(s[i]);
//		if(!is_dot_found)
//		{
//			if(s[i]=='.') { is_dot_found = 1; }
//		}
//		else j++;
//	}
    
}




int bdString::IsOfPattern(bdString &pattern_string)
{
	if(this->IsEmpty()) return 0;
	//If the number of characters don't match, the pattern will not match eather.
	if(this->NumberOfCharacters()!=pattern_string.NumberOfCharacters()) return 0;

	char ch, ch2;
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		ch = (*this)[i];
		ch2 = pattern_string[i];

		if(ch2=='#') { if(!(ch>='0' && ch<='9')) return 0; }
		else
		{
			if(ch2=='%') { if(ch>='0' && ch<='9') return 0; }
			else { if(ch!=ch2) return 0; }
		}
	}
	return 1;
}


//int bdString::StringToInt(int &output_number)
//{
//	//cout<<"StringToInt string:"<<(*this)<<endl;
//
//	if(this->IsEmpty()) return 0;
//	if(this->NumberOfCharacters()==1)//If there is only 1 character and it's not a number, return 0 (this is in case the character is '-' or '+')
//	{
//		if(!(this->list.GetLeftEnd()>='0' && this->list.GetLeftEnd()<='9')) return 0;
//	}
//
//	char ch;
//	int m = 1;
//	int basis = 1;
//	int is_sign_found = 0;//After the first character, this indicator is set to 1.
//	output_number = 0;
//	int is_non_number_character_found = 0;
//	for(this->list.For_StartFromRightEnd(); this->list.GetAllElemetsByMovingToLeft(&ch); )
//	{
//		if(!is_sign_found)//This is ecexuted just in the first iteration.
//		{
//			if(ch=='+');
//			else
//			{
//				if(ch=='-') m = -1;
//				else
//				{
//					if(ch>='0' && ch<='9')
//					{
//						int i = ch - '0';
//						output_number += i*basis;
//						basis = basis * 10;
//					}
//					else
//					{
//						is_non_number_character_found = 1;
//						this->list.ExitLoopAtNextFunctionCall();
//					}
//				}
//			}
//			is_sign_found = 1;
//		}
//		else
//		{
//			if(ch>='0' && ch<='9')
//			{
//				int i = ch - '0';
//				output_number += i*basis;
//				basis = basis * 10;
//			}
//			else
//			{
//				is_non_number_character_found = 1;
//				this->list.ExitLoopAtNextFunctionCall();
//			}
//		}
//	}
//
//	output_number = output_number*m;
//
//	//If there was a non-number character in the string, return 0.
//	if(is_non_number_character_found) return 0;
//	
//	return 1;
//}
//
//
//
//int bdString::StringToDouble(double &output_number)
//{//THIS IS NOT WORKING - THERE IS A BUG SOMEWHERE!!!
//	//Double value can be represented as converting two int values (one before the dot and one after it)
//
//	bdString s1, s2, *ps;//'ps' is the pointer to the string we currently append to, we initialize it with s1 address
//	ps = &s1;
//	//Find the dot and generate the strings s1 and s2:
//	char ch;
//	for(this->list.For_StartFromLeftEnd(); this->list.GetAllElemetsByMovingToRight(&ch); )
//	{
//		if(ch=='.') ps = &s2;
//		else ps->Append(ch);
//	}
//
//	//At this point s1 should contain the first int and s2 the second (if exists)
//
//	int i1=0, i2=0;
//	if(!s1.StringToInt(i1)) return 0;//If error is found, return fail.
//
//	if(!s2.IsEmpty())//If the second string is not empty
//	{
//		if(!s2.StringToInt(i2)) return 0;//If error is found, return fail
//		if(s2[0]<'0' || s2[0]>'9') return 0;//If non-zero character is found as the first in the string, return fail.
//
//		output_number = (double) i2;
//		double base = 1;
//		for(int j=0; j<s2.NumberOfCharacters(); j++)
//		{
//			base = base * 10;
//		}
//		output_number = output_number / base;
//		output_number = output_number + ((double)i1);
//	}
//	else output_number = (double) i1;
//
//	return 1;
//}


int bdString::StringToInt(int &output_number)
{
	if(this->IsEmpty()) return 0;
	std::istringstream is(this->m_str);
	return (is>>output_number) ? true : false;
}


int bdString::StringToDouble(double &output_number)
{
	if(this->IsEmpty()) return 0;
	std::istringstream is(this->m_str);
	return (is>>output_number) ? true : false;
}


int bdString::StringToPosition(double &out_coord1, double &out_coord2, double &out_coord3)
{
    bdList<double> extracted_position;
    this->ExtractAllDoubleNumbersToList(extracted_position);
    if(extracted_position.GetNumberOfElements()>=3)
    {
        out_coord1 = extracted_position[0];
        out_coord2 = extracted_position[1];
        out_coord3 = extracted_position[2];
        return 1;
    }
    return 0;
}


void bdString::ExtractAllIntNumbersToList(bdList<int> &list_of_extracted_int_numbers)
{
	//cout<<"ExtractAllIntNumbersToList string:"<<(*this)<<endl;

	list_of_extracted_int_numbers.Reset();

	if(this->IsEmpty()) return;

	bdList<bdString> l;//We extract all strings that represent numbers to this list

	int is_sign_found = 0;
	char ch;
	l.AddNewToRightEnd();
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		ch = (*this)[i];

		if(!is_sign_found)
		{
			if(ch=='+' || ch=='-' || (ch>='0' && ch<='9'))
			{
				l.GetRightEnd().Append(ch);
				is_sign_found = 1;
			}
		}
		else
		{
			if(ch>='0' && ch<='9') l.GetRightEnd().Append(ch);
			else
			{
				l.AddNewToRightEnd();
				is_sign_found = 0;
			}
		}
	}

	//Now we convert all the strings in list l to int values and if conversion is valid, we enter values to output list.
	while(!l.IsEmpty())
	{
		int n;
		if(l.GetLeftEnd().StringToInt(n)) list_of_extracted_int_numbers.AddToRightEnd(n);
		l.DeleteLeftEnd();
	}
}


void bdString::ExtractAllIntNumbersToList(std::list<int> &list_of_extracted_int_numbers)
{
	list_of_extracted_int_numbers.clear();
	if(this->IsEmpty()) return;

	bdString temp_string;
	for(unsigned int i=0; i<m_str.length(); i++)
	{
		char ch = (*this)[i];

		//If the current character is one of the possible characters for an int value...
		if(ch=='+' || ch=='-' || (ch>='0' && ch<='9'))
		{
			//...copy it into the temp_string 
			temp_string.Append(ch);
		}
		//If we have found the character that is not part of the number... 
		else
		{
			//...try to convert the temp_string into an int and if success, write it to output list.
			int n;
			if(temp_string.StringToInt(n)) { list_of_extracted_int_numbers.push_back(n); }
			temp_string.Clear();
		}
	}

	// If the string ended with a number...
	if(!temp_string.IsEmpty())
	{
		//...try to convert the temp_string into an int and if success, write it to output list.
		int n;
		if(temp_string.StringToInt(n)) { list_of_extracted_int_numbers.push_back(n); }
		temp_string.Clear();
	}
}


void bdString::ExtractHexadecimalIntNumberStringsToList(bdList<bdString> &list_hexadecimal_number_strings)
{
    list_hexadecimal_number_strings.Reset();
    
    if(this->IsEmpty()) return;
    
    // We extract all strings that represent numbers to input/output list
    
    int is_sign_found = 0;
    char ch;
    bdString *ps = list_hexadecimal_number_strings.AddNewToRightEnd();
    for(unsigned int i=0; i<this->m_str.length(); i++)
    {
        ch = (*this)[i];
        
        if(!is_sign_found)
        {
            if(ch=='+' || ch=='-' || (ch>='0' && ch<='9') || (ch>='a' && ch<='f') || (ch>='A' && ch<='F'))
            {
                ps->Append(ch);
                is_sign_found = 1;
            }
        }
        else
        {
            if((ch>='0' && ch<='9') || (ch>='a' && ch<='f') || (ch>='A' && ch<='F')) ps->Append(ch);
            else
            {
                ps = list_hexadecimal_number_strings.AddNewToRightEnd();
                is_sign_found = 0;
            }
        }
    }
    
    // Check if the last string in the list is empty (it can happen if the original string ended with a non-hexadecimal character)
    if(ps->IsEmpty()) list_hexadecimal_number_strings.DeleteRightEnd();
}


void bdString::ExtractSeriesOfIntNumberStringsWithDots(bdList<bdString> &list_number_strings_with_dots)
{
    list_number_strings_with_dots.Reset();
    
    if(this->IsEmpty()) return;
    
    // We extract all strings that represent numbers to input/output list
    
    char ch;
    int is_number_found = 0;
    bdString *ps = list_number_strings_with_dots.AddNewToRightEnd();
    for(unsigned int i=0; i<this->m_str.length(); i++)
    {
        ch = (*this)[i];
        if((ch>='0' && ch<='9') || ch=='.')
        {
            ps->Append(ch);
            is_number_found = 1;
        }
        else
        {
            if(is_number_found)
            {
                ps = list_number_strings_with_dots.AddNewToRightEnd();
                is_number_found = 0;
            }
        }
    
    }
    
    // Check if the last string in the list is empty (it can happen if the original string ended with a non-serial-number character)
    if(ps->IsEmpty()) list_number_strings_with_dots.DeleteRightEnd();
    
}


void bdString::ExtractLastIntNumber(bdString &remaining_string, int &number)
{
	if(this->IsEmpty()) return;

	bdList<char> list_of_chars;

	int is_number_found = 0;
	int i = (int) (m_str.length()-1);
	for( ; i>=0; i--)
	{
		char ch = (*this)[i];

		//If the current character is one of the possible characters for an int value...
		if(ch=='+' || ch=='-' || (ch>='0' && ch<='9'))
		{
			//...copy it into the list 
			list_of_chars.AddToLeftEnd(ch);
			if(ch>='0' && ch<='9') is_number_found = 1;
		}
		//If we have found the character that is not part of the number... 
		else
		{
			if(is_number_found)
			{
				//...try to convert the list of chars into an int.
				bdString temp_string;
				bdListIterator<char> it;
				for(it.SetLeftEnd(list_of_chars); it.IsValid(); it.MoveRight()) { temp_string.Append(it.GetElement()); }
				temp_string.StringToInt(number);
				break;
			}
		}
	}

	//copy the rest of the string to output (from index 0 till index 'i').
	if(is_number_found)
	{
		for(int j=0; j<=i; j++) { remaining_string.Append((*this)[j]); }
	}
	else
	{
		remaining_string = (*this);
	}
}


void bdString::ExtractFirstNonIntNumberString(bdString &extracted_string)
{
	if(this->IsEmpty()) return;
	extracted_string.Clear();
	
	int is_first_non_number_character_found = 0;

	for(unsigned int i=0; i<m_str.length(); i++)
	{
		char ch = (*this)[i];

		// If the first non-number character has not yet been found...
		if(!is_first_non_number_character_found) 
		{
			//If the current character is not one of the possible characters for an int value...
			if(!(ch=='+' || ch=='-' || (ch>='0' && ch<='9')))
			{
				is_first_non_number_character_found = 1;
				//...copy it into the temp_string 
				extracted_string.Append(ch);
			}
		}
		// If the first non-number character has been found (so output string is not empty)...
		else
		{
			//If the current character is not one of the possible characters for an int value...
			if(!(ch=='+' || ch=='-' || (ch>='0' && ch<='9')))
			{
				//...copy it into the temp_string 
				extracted_string.Append(ch);
			}
			// ...otherwise just exit function
			else return;
		}
	}
}


void bdString::ExtractAllDoubleNumbersToList(std::list<double> &list_of_extracted_doubles)
{
	list_of_extracted_doubles.clear();
	if(this->IsEmpty()) return;

	bdString temp_string;
	for(unsigned int i=0; i<m_str.length(); i++)
	{
		char ch = (*this)[i];

		//If the current character is one of the possible characters for a double value...
		if(ch=='+' || ch=='-' || (ch>='0' && ch<='9') || ch=='.' || ch=='e')
		{
			//...copy it into the temp_string 
			temp_string.Append(ch);
		}
		//If we have found the character that is not part of the number... 
		else
		{
			//...try to convert the temp_string into a double and if success, write it to output list.
			double d;
			if(temp_string.StringToDouble(d)) { list_of_extracted_doubles.push_back(d); }
			temp_string.Clear();
		}
	}

	// If the string ended with a number...
	if(!temp_string.IsEmpty())
	{
		//...try to convert the temp_string into a double and if success, write it to output list.
		double d;
		if(temp_string.StringToDouble(d)) { list_of_extracted_doubles.push_back(d); }
		temp_string.Clear();
	}
}


void bdString::ExtractAllDoubleNumbersToList(bdList<double> &list_of_extracted_doubles)
{
	list_of_extracted_doubles.Reset();
	if(this->IsEmpty()) return;

	bdString temp_string;
	for(unsigned int i=0; i<m_str.length(); i++)
	{
		char ch = (*this)[i];

		//If the current character is one of the possible characters for a double value...
		if(ch=='+' || ch=='-' || (ch>='0' && ch<='9') || ch=='.' || ch=='e')
		{
			//...copy it into the temp_string 
			temp_string.Append(ch);
		}
		//If we have found the character that is not part of the number... 
		else
		{
			//...try to convert the temp_string into a double and if success, write it to output list.
			double d;
			if(temp_string.StringToDouble(d)) { list_of_extracted_doubles.AddToRightEnd(d); }
			temp_string.Clear();
		}
	}

	// If the string ended with a number...
	if(!temp_string.IsEmpty())
	{
		//...try to convert the temp_string into a double and if success, write it to output list.
		double d;
		if(temp_string.StringToDouble(d)) { list_of_extracted_doubles.AddToRightEnd(d); }
		temp_string.Clear();
	}
}


int bdString::HasPrefix(const char *prefix)
{
    bdString s;
    s.Assign(prefix);
    if(this->NumberOfCharacters() < s.NumberOfCharacters()) return 0;
    for(unsigned int i=0; i<s.NumberOfCharacters(); i++)
    {
        if((*this)[i]!=s[i]) return 0;
    }
    return 1;
}



//int bdString::ToCharString(char *p_output_string)
//{
//	if(p_output_string==NULL) return 0;
//
//	int i=0;
//	for(this->list.For_StartFromLeftEnd(); this->list.GetAllElemetsByMovingToRight(&(p_output_string[i])); i++);
//	
//	p_output_string[i] = 0;
//	return 1;
//}


void bdString::ReplaceCharacter(char character_to_replace, char replacement_character)
{
	//m_str.replace( m_str.begin(), m_str.end(), character_to_replace, replacement_character);
    for(unsigned int i=0; i<this->m_str.length(); i++)
    {
        if(this->m_str[i]==character_to_replace) { this->m_str[i] = replacement_character; }
    }
}


void bdString::PathToWindowsPath()
{
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		if(this->m_str[i]=='/') { this->m_str[i] = '\\'; }
	}
}


void bdString::PathToLinuxPath()
{
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		if(this->m_str[i]=='\\') { this->m_str[i] = '/'; }
	}
}



int bdString::ExtractStringBetweenCharacters(char start_character, char end_character, bdString &output)
{
	output.Clear();
	int is_start_character_found = 0;
	int is_end_character_found = 0;
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		char ch = this->m_str[i];
		if(!is_start_character_found)
		{ 
			if(ch == start_character) is_start_character_found = 1; 
		}
		else
		{
			if(ch == end_character)
			{
				is_end_character_found = 1;
				break;
			}
			else output.Append(ch);
		}
	}	

	if(is_end_character_found) return 1;
	else
	{
		output.Clear();
		return 0;
	}
}


int bdString::ExtractFirstCommonStringPartWithInputString(bdString &input, bdString &output)
{
    if(this->IsEmpty() || input.IsEmpty()) return 0;
    output.Clear();
    for(unsigned int i=0; i<this->m_str.length() && i<input.m_str.length(); i++)
    {
        if(this->m_str[i]==input.m_str[i]) { output.Append(this->m_str[i]); }
        else { break; }
    }
    
    return 1;
}


int bdString::ExtractDirectory(bdString &output_dir)
{
	if(this->IsEmpty()) return 0;
	if(!( (((*this)[0]>='a' && (*this)[0]<='z')||((*this)[0]>='A' && (*this)[0]<='Z')) || (*this)[0]=='\\' || (*this)[0]=='/') ) return 0;

	//Check the double dot if the path starts with a letter
	if( ((*this)[0]>='a' && (*this)[0]<='z') || ((*this)[0]>='A' && (*this)[0]<='Z') )
	{
		if((*this)[1]!=':') return 0;
	}

	char ch;
	int index_of_last_bracket = -1;
	int number_of_brackets = 0; 
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		ch = (*this)[i];
		if(ch=='\\' || ch=='/')
		{
			index_of_last_bracket = i;
			number_of_brackets++;
		}
	}

	if(number_of_brackets==0) return 0;

	//Just copy the string to the last bracket...
	output_dir.Clear();
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		ch = (*this)[i];
		//if(i==index_of_last_bracket) return 1;
		if(i==index_of_last_bracket+1) return 1;//Changed this to return the slash as the last character!
		else output_dir.Append(ch);
	}

	return 1;
}


int bdString::ExtractFileNamePathParts(bdString &out_dir, bdString &out_file_name, bdString &out_extension)
{
	out_dir.Clear();
	out_extension.Clear();
	out_file_name.Clear();

	this->ExtractDirectory(out_dir);
	bdString str_file_name_with_ext;
	this->DifferenceFromFilePath(out_dir, str_file_name_with_ext);
	
	int index_of_dot = -1;
	for(int i=str_file_name_with_ext.NumberOfCharacters()-1; i>=0; i--) 
	{ 
		if(str_file_name_with_ext[i]=='.') 
		{
			index_of_dot = i;
			break;
		}
	}
	if(index_of_dot>=((int)str_file_name_with_ext.NumberOfCharacters())) return 0;
    if(index_of_dot<0)
    {
        out_file_name = str_file_name_with_ext;
        out_extension.Append("");
        return 1;
    }
    for(int i=0; ((unsigned int)i)<str_file_name_with_ext.NumberOfCharacters(); i++)
	{
		if(i<index_of_dot) out_file_name.Append(str_file_name_with_ext[i]);
		else
		{
			if(i>index_of_dot) out_extension.Append(str_file_name_with_ext[i]);
		}
	}
	return 1;
}


int bdString::InsertNumberToString(int n)
{
	//char *pch;
	//int i=0;
	int num_pos = 0;
	int is_num_sign_found = 0;
	int number_of_num_characters = 0;
	//for(this->list.For_StartFromLeftEnd(); this->list.GetAllPointersToElemetsByMovingToRight(&pch); )
	for(unsigned int i=0; i<this->m_str.length(); i++ )
	{
		if(!is_num_sign_found)
		{
			//if((*pch)=='#')
			if((*this)[i]=='#')
			{
				is_num_sign_found = 1;
				number_of_num_characters++;
				num_pos = i;
				//*pch = '0';
				(*this)[i] = '0';
			}
		}
		else
		{
			//if((*pch)=='#')
			if((*this)[i]=='#')
			{
				number_of_num_characters++;
				//*pch = '0';
				(*this)[i] = '0';

			}
			//else this->list.ExitLoopAtNextFunctionCall();
			else break;
		}
		//i++;
	}

	if(!is_num_sign_found) return 0;

	bdString bds;
	bds.NumberToString(n);

	for(int j=0; j<number_of_num_characters && j<((int)bds.NumberOfCharacters()); j++)
	{
		(*this)[num_pos+number_of_num_characters-1-j] = bds[bds.NumberOfCharacters()-1-j];
	}

	return 1;
}


void bdString::InvertCapitalization()
{
    for(unsigned int i=0; i<this->m_str.length(); i++ )
    {
        if(islower((*this)[i])) (*this)[i] = toupper((*this)[i]);
        else (*this)[i] = tolower((*this)[i]);
    }
}


void bdString::AddTag(char tag_character, bdString &input_tag)
{
	std::string temp = this->m_str;

	this->m_str.clear();
	this->m_str.assign(tag_character + input_tag.m_str + tag_character + temp);
}


int bdString::ExtractTag(char tag_character, bdString &output_tag)
{
	output_tag.Clear();

	int is_first_tag_sign_found = 0;
	int is_second_tag_sign_found = 0;
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		char ch = this->m_str[i];
		if(!is_first_tag_sign_found)
		{
			if(ch==tag_character) is_first_tag_sign_found = 1;
		}
		else
		{
			if(ch==tag_character)
			{
				is_second_tag_sign_found = 1;
				break;
			}
			else output_tag.Append(ch);
		}
	}	

	if(is_second_tag_sign_found) return 1;
	else
	{
		output_tag.Clear();
		return 0;
	}
}


int bdString::RemoveTag(char tag_character, bdString &output_string)
{
	output_string.Clear();

	char ch;
	int is_first_tag_sign_found = 0;
	int is_second_tag_sign_found = 0;
	int index_of_first_non_tag_sharacter = -1;
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		ch = this->m_str[i];
		if(!is_first_tag_sign_found)
		{
			if(ch==tag_character) is_first_tag_sign_found = 1;
		}
		else
		{
			if(ch==tag_character) 
			{
				is_second_tag_sign_found = 1;
				index_of_first_non_tag_sharacter = i+1;
				break;
			}
		}
	}	

	if(is_second_tag_sign_found && index_of_first_non_tag_sharacter<((int)this->m_str.length())) 
	{
		output_string.m_str.append(this->m_str.c_str(),index_of_first_non_tag_sharacter,this->m_str.length()-index_of_first_non_tag_sharacter);
		return 1;
	}
	else return 0;
}


int bdString::RemoveTag(char tag_character1, char tag_character2, bdString &output_string)
{
    output_string.Clear();
    
    char ch;
    int is_first_tag_sign_found = 0;
    int is_second_tag_sign_found = 0;
    int index_of_first_non_tag_sharacter = -1;
    for(unsigned int i=0; i<this->m_str.length(); i++)
    {
        ch = this->m_str[i];
        if(!is_first_tag_sign_found)
        {
            if(ch==tag_character1) is_first_tag_sign_found = 1;
        }
        else
        {
            if(ch==tag_character2)
            {
                is_second_tag_sign_found = 1;
                index_of_first_non_tag_sharacter = i+1;
                break;
            }
        }
    }
    
    if(is_second_tag_sign_found && index_of_first_non_tag_sharacter<((int)this->m_str.length()))
    {
        output_string.m_str.append(this->m_str.c_str(),index_of_first_non_tag_sharacter,this->m_str.length()-index_of_first_non_tag_sharacter);
        return 1;
    }
    else return 0;
}



int bdString::ExtractOptions(char tag_character, bdList<bdString> &output_options_list)
{
	if(this->IsEmpty()) return 0;

	//!!!! IT SHOULD CHECK HERE IF IT WILL BE POSSIBLE TO EXTRACT AT LEAST ONE OPTION FROM THE STRING!!!!
	//I.E. CHECK IF IT CONTAINS AT LEAST 2 tag_characters!!!!!

	bdString t1, t2, *p1, *p2, *p_temp, tag_string;
	t1 = *this;
	p1 = &t1; p2 = &t2;
	while(p1->ExtractTag(tag_character,tag_string))
	{
		output_options_list.AddToRightEnd(tag_string);
		p1->RemoveTag(tag_character,*p2);
		p_temp = p1;
		p1 = p2;
		p2 = p_temp;
	}
	return 1;
}


int bdString::ExtractOptions(char tag_character1, char tag_character2, bdList<bdString> &output_options_list)
{
    if(this->IsEmpty()) return 0;
    
    //!!!! IT SHOULD CHECK HERE IF IT WILL BE POSSIBLE TO EXTRACT AT LEAST ONE OPTION FROM THE STRING!!!!
    //I.E. CHECK IF IT CONTAINS AT LEAST 2 tag_characters!!!!!
    
    bdString t1, t2, *p1, *p2, *p_temp, tag_string;
    t1 = *this;
    p1 = &t1; p2 = &t2;
    while(p1->ExtractStringBetweenCharacters(tag_character1, tag_character2, tag_string))
    {
        output_options_list.AddToRightEnd(tag_string);
        p1->RemoveTag(tag_character1,tag_character2,*p2);
        p_temp = p1;
        p1 = p2;
        p2 = p_temp;
    }
    return 1;
}


int bdString::DifferenceFromInputString(bdString &input_string, bdString &output_difference_string)
{
	if(input_string.IsEmpty() || this->IsEmpty()) return 0;
	output_difference_string.Clear();

	//Find the index of the first different character
	int index_of_first_different_char = -1;
	for(unsigned int i=0; i<this->m_str.length(); i++)
	{
		if(input_string[i]!=this->m_str[i])
		{
			index_of_first_different_char = i;
			break;
		}
	}

	//COPY THE STRING this->m_str FROM INDEX index_of_first_different_char TO output_difference_string.
	if(index_of_first_different_char>=0) output_difference_string.m_str.append(this->m_str.c_str(),index_of_first_different_char,this->m_str.length()-index_of_first_different_char);

	return 1;
}


int bdString::DifferenceFromFilePath(bdString &input_path, bdString &output_difference_path)
{
	if(input_path.IsEmpty() || this->IsEmpty()) return 0;
	output_difference_path.Clear();

	char c1, c2;
	int index_from_which_strings_different = -1;
	for(unsigned int i=0; i<this->m_str.length() && index_from_which_strings_different<0; i++)
	{
		c1 = this->m_str[i];
		c2 = input_path.m_str[i];
		if(!((c1=='\\' || c1=='/') && (c2=='\\' || c2=='/')) && c1!=c2) index_from_which_strings_different = i;
	}

	//If strings are not different, return 1.
	if(index_from_which_strings_different<0) return 1;

	//If strings are different, copy the different part
	output_difference_path.m_str.append(this->m_str.c_str(),index_from_which_strings_different,this->m_str.length()-index_from_which_strings_different);

	return 1;
}


bdString& bdString::operator ()(char *s)
{
	if(s==NULL) return *this;
	this->Clear();
	this->Append(s);
   
    return *this;
}


bdString& bdString::operator ()(const char *s)
{
    if(s==NULL) return *this;
    this->Clear();
    this->Append(s);
    
    return *this;
}


bdString& bdString::operator =(bdString &s)
{
    if (&s==this) return *this;

	this->Clear();
	this->Append(s);
   
    return *this;
}


bdString& bdString::operator =(char *s)
{
	if(s==NULL) return *this;

	this->Clear();
	this->Append(s);  
    return *this;
}


bdString& bdString::operator =(const char *s)
{
	if(s==NULL) return *this;

	this->Clear();
	this->Append(s);  
    return *this;
}



int bdString::operator ==(bdString &s)
{
	if(&s==this) return 1;

    if(this->IsEmpty() && s.IsEmpty()) return 1;
    else
    {
        if(this->IsEmpty() || s.IsEmpty()) return 1;
    }
    
    if(this->NumberOfCharacters()!=s.NumberOfCharacters()) return 0;
       
       
	//The 'compare' method returns 0 if the strings are equal!
	if(!this->m_str.compare(s.m_str)) return 1;
	else return 0;

 //   if(&s==this) return 1;
	//if(this->NumberOfCharacters()!=s.NumberOfCharacters()) return 0;

	//this->list.CurrentNodeSetToLeftEnd();
	//s.list.CurrentNodeSetToLeftEnd();
	//int is_checking = 1;
	//while(is_checking)
	//{
	//	char c1, c2;
	//	this->list.CurrentNodeGetElement(&c1);
	//	s.list.CurrentNodeGetElement(&c2);
	//	if(c1!=c2) return 0;

	//	if(this->list.CurrentNodeMoveToRight())
	//	{
	//		s.list.CurrentNodeMoveToRight();
	//	}
	//	else is_checking = 0;
	//}

	////If it came to here the strings are identical
	//return 1;
}



int bdString::operator ==(char *s)
{
	if(s==NULL) return 0;

	//The 'compare' method returns 0 if the strings are equal!
	if(!this->m_str.compare(s)) return 1;
	else return 0;

	//if(this->IsEmpty()) return 0;

	//this->list.CurrentNodeSetToLeftEnd();
	//int is_checking = 1;
	//int i = 0;
	//while(is_checking)
	//{
	//	if(s[i]==0) return 0;

	//	char c1;
	//	this->list.CurrentNodeGetElement(&c1);
	//	if(s[i]!=c1) return 0;

	//	if(!this->list.CurrentNodeMoveToRight()) is_checking = 0;
	//	i++;
	//}

	////If it came to here the strings are identical
	//return 1;
}


int bdString::operator ==(const char *s)
{
	if(s==NULL) return 0;

	//The 'compare' method returns 0 if the strings are equal!
	if(!this->m_str.compare(s)) return 1;
	else return 0;
}



int bdString::operator !=(bdString &s)
{
	//if((*this)==s) return 0;
	//else return 1;
    return (!((*this)==s));
}


int bdString::operator !=(char *s)
{
	//if((*this)==s) return 0;
	//else return 1;
    return (!((*this)==s));
}


int bdString::operator !=(const char *s)
{
    //if((*this)==s) return 0;
    //else return 1;
    return (!((*this)==s));
}


int bdString::operator <(bdString &s)
{
	if(&s==this) return 0;

	unsigned int length1, length2, length;
	length1 = (unsigned int) m_str.length();
	length2 = (unsigned int) s.m_str.length();
	if(length1<length2) length = length1;
	else length = length2;

	char c1, c2;
	for(unsigned int i=0; i<length; i++)
	{
		c1 = m_str[i];
		c2 = s.m_str[i];
		if(c1>=c2) return 0;
	}
	//if(length1>length2) return 0;

	return 1;
}


int bdString::operator >(bdString &s)
{
	if(&s==this) return 0;

	unsigned int length1, length2, length;
	length1 = (unsigned int) m_str.length();
	length2 = (unsigned int) s.m_str.length();
	if(length1<length2) length = length1;
	else length = length2;

	char c1, c2;
	for(unsigned int i=0; i<length; i++)
	{
		c1 = m_str[i];
		c2 = s.m_str[i];
		if(c1<=c2) return 0;
	}
	//if(length1>length2) return 0;

	return 1;
}


std::ostream& operator << (std::ostream &o, bdString &s)
{ 
	o<<s.m_str;
    return o;
}


std::stringstream& operator << (std::stringstream &o, bdString &s)
{ 
	o<<s.m_str;
    return o;
}