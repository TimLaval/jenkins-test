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


#include "bdArrayOfDoubles.h"



bdArrayOfDoubles::bdArrayOfDoubles()
{
}


bdArrayOfDoubles::~bdArrayOfDoubles()
{
}


void bdArrayOfDoubles::SortAscending()
{
    bdArraySortAndSearch<double> ss;
	ss.SortAscending(*this);
}


void bdArrayOfDoubles::SortAscending(bdArray<int> &index_change_order)
{
    bdArraySortAndSearch<double> ss;
    ss.SortAscending((*this),index_change_order);
}


void bdArrayOfDoubles::SortDescending()
{
    bdArraySortAndSearch<double> ss;
	ss.SortDescending((*this));
}


void bdArrayOfDoubles::SortDescending(bdArray<int> &index_change_order)
{
    bdArraySortAndSearch<double> ss;
	ss.SortDescending((*this),index_change_order);
}


int bdArrayOfDoubles::MinAndMax(double &min, double &max, int &index_min, int &index_max)
{
    bdArraySortAndSearch<double> ss;
	return (ss.MinAndMax((*this), min, max, index_min, index_max));
}


int bdArrayOfDoubles::IndexOfElement(double &element)
{
    bdArraySortAndSearch<double> ss;
	return (ss.IndexOfElement((*this), element));
}


int bdArrayOfDoubles::LoadFromFile_DoubleValues(char *file_name, char *array_name)
{
	int line_length = 100000;
	char *line;
	line = new char [line_length];

	for(int p=0; p<line_length-1; p++) line[p] = 'a';

	//Open the file
	ifstream input_file;
	input_file.open(file_name,ios::binary);
    if(!input_file) 
	{
		cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Unable to open bdArray file: "<<file_name<<endl;
        return 0;
    }
	
	//If file has multiple arrays stored in it, find the one that matches the array name
	int array_found = 0;
	int i=0; //line character counter
	while(!array_found)
	{	
		//Get a line from the file
		if(!(input_file.getline(line, line_length-1)))
		{
			cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Error in reading bdArray file!"<<endl;
			input_file.close();
			return 0;
		}

		//Check the taken line
		//int name_detected = 1; 
		array_found = 1;
		for(i=0; array_name[i]!=0 && array_found; i++)
		{
			//If the names don't match
			if(line[i]!=array_name[i]) array_found = 0;
		}
	}

	//Search for bracket: "["
	while(line[i]!='[')
	{
		i++;
		if(i==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Error in reading bdArray file! No '[' character found!"<<endl;
			input_file.close();
			return 0;
		}
	}
	i++;

	//Search for closing bracket: "]"
	int i2 = i;
	while(line[i2]!=']') 
	{
		i2++;
		if(i2==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Error in reading bdArray file! No ']' character found!"<<endl;
			input_file.close();
			return 0;
		}
	}
	i2--;

	//Count the number of elements in the array - this is done by calculating the number of commas ","
	//Also checks the correctness of characters in between two brackets: "[" and "]"
	int number_of_elements_of_array = 1;
	for(int c=i; c<=i2; c++)
	{
		if(line[c]==',') number_of_elements_of_array++;
		//Check correctness of characters
		if(line[c]!=' ' && line[c]!=',' && line[c]!='.' && line[c]<'0' && line[c]>'9')
		{
			cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Error in reading bdArray file! Incorrect characters found!"<<endl;
			input_file.close();
			return 0;
		}
	}

	this->Reset();
	this->Set(number_of_elements_of_array);

	//Fill in the whole array
	double dec = 1;
	double value = 0;
	int number_of_added_elements = 0;
	int is_dot_found = 0;
	for(int r=i2 ; r>=i-1; r--)
	{
		if(line[r]==',' || line[r]=='[')
		{
			//Write the calculated value in the array
			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
			number_of_added_elements++;

			//Reset 'dec' and 'value' for a new number 
			dec = 1;
			value = 0;

			//Set the dot indicator to 0
			is_dot_found = 0;
		}
		//If a number is found
		if(line[r]>='0' && line[r]<='9')
		{
			value += ((int)(line[r] - '0'))*dec;
			dec = dec * 10;
		}
		if(line[r]=='.')
		{
			//If this is the first dot
			if(!is_dot_found)
			{
				value = value/dec;
				is_dot_found = 1;
				dec = 1;
			}
			//If this is the second dot, the error ocurred, exit with 0
			else
			{
				cout<<"bdArrayOfDoubles::LoadFromFile_DoubleValues(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
				input_file.close();
				return 0;
			}
		}
	}

	delete line;
	input_file.close();
	return 1;
}




int bdArrayOfDoubles::LoadFromString_DoubleValues(char *array_string)
{
	int line_length = 100000;

	//Search for bracket: "["
	int i=0;
	for(; array_string[i]!='[' && i<line_length; i++)
	{
		if(i==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No '[' character found!"<<endl;
			return 0;
		}
	}
	i++;

	//Search for closing bracket: "]"
	int i2 = i;
	for(; array_string[i2]!=']' && i2<line_length; i2++)
	{
		if(i2==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No ']' character found!"<<endl;
			return 0;
		}
	}
	i2--;

	//Count the number of elements in the array - this is done by calculating the number of commas ","
	//Also checks the correctness of characters in between two brackets: "[" and "]"
	int number_of_elements_of_array = 1;
	for(int c=i; c<=i2; c++)
	{
		if(array_string[c]==',') number_of_elements_of_array++;
		//Check correctness of characters
		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && array_string[c]!='-' && array_string[c]!='+' && (!(array_string[c]>='0' && array_string[c]<='9')))
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect characters found!"<<endl;
			return 0;
		}
	}

	this->Reset();
	this->Set(number_of_elements_of_array);


	//Fill in the whole array
	double dec = 1;
	double value = 0;
	int number_of_added_elements = 0;
	int is_dot_found = 0;
	int is_sign_found = 0;
	for(int r=i2 ; r>=i-1; r--)
	{
		if(array_string[r]==',' || array_string[r]=='[')
		{
			//Write the calculated value in the array
			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
			number_of_added_elements++;

			//Reset 'dec' and 'value' for a new number 
			dec = 1;
			value = 0;

			//Set the dot indicator to 0
			is_dot_found = 0;

			//Set the sign indicator to 0
			is_sign_found = 0;
		}
		//Since we search backwards, the numbers and dot can be taken into account only if the sign has not yet been found!
		if(!is_sign_found)
		{
			//If MINUS sign is found, make negative value, set the indicator of found sign
			if(array_string[r]=='-')
			{
				value = -value;
				is_sign_found = 1;
			}
			//If PLUS sign is found, just set the indicator of found sign
			if(array_string[r]=='+') is_sign_found = 1;

			//If a number is found
			if(array_string[r]>='0' && array_string[r]<='9')
			{
				value += ((int)(array_string[r] - '0'))*dec;
				dec = dec * 10;
			}
			if(array_string[r]=='.')
			{
				//If this is the first dot
				if(!is_dot_found)
				{
					value = value/dec;
					is_dot_found = 1;
					dec = 1;
				}
				//If this is the second dot, the error ocurred, exit with 0
				else
				{
					cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
					return 0;
				}
			}
		}
		else
		{
			//If the sign is found and the current character is not ',' or '[', report an error
			if(array_string[r]!=',' && array_string[r]!='[' && array_string[r]!=' ')
			{
				cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect character found before the sign!"<<endl;
				return 0;
			}
		}
	}

	return 1;
}





//int bdArrayOfDoubles::LoadFromString_DoubleValues(char *array_string)
//{
//	int line_length = 100000;
//
//	//Search for bracket: "["
//	int i=0;
//	for(; array_string[i]!='[' && i<line_length; i++)
//	{
//		if(i==line_length-1)
//		{
//			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues(): Error in reading bdArray file! No '[' character found!"<<endl;
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
//			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues(): Error in reading bdArray file! No ']' character found!"<<endl;
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
//		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && (!(array_string[c]>='0' && array_string[c]<='9')))
//		{
//			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues(): Error in reading bdArray file! Incorrect characters found!"<<endl;
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
//		}
//		//If a number is found
//		if(array_string[r]>='0' && array_string[r]<='9')
//		{
//			value += ((int)(array_string[r] - '0'))*dec;
//			dec = dec * 10;
//		}
//		if(array_string[r]=='.')
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
//				cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
//				return 0;
//			}
//		}
//	}
//
//	return 1;
//}




int bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(char *array_string)
{
	int line_length = 100000;

	//Search for bracket: "["
	int i=0;
	for(; array_string[i]!='[' && i<line_length; i++)
	{
		if(i==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No '[' character found!"<<endl;
			return 0;
		}
	}
	i++;

	//Search for closing bracket: "]"
	int i2 = i;
	for(; array_string[i2]!=']' && i2<line_length; i2++)
	{
		if(i2==line_length-1)
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! No ']' character found!"<<endl;
			return 0;
		}
	}
	i2--;

	//Count the number of elements in the array - this is done by calculating the number of commas ","
	//Also checks the correctness of characters in between two brackets: "[" and "]"
	int number_of_elements_of_array = 1;
	for(int c=i; c<=i2; c++)
	{
		if(array_string[c]==',') number_of_elements_of_array++;
		//Check correctness of characters
		if(array_string[c]!=' ' && array_string[c]!=',' && array_string[c]!='.' && array_string[c]!='-' && array_string[c]!='+' && (!(array_string[c]>='0' && array_string[c]<='9')))
		{
			cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect characters found!"<<endl;
			return 0;
		}
	}

	//If the number of numbers in the line is not the same as the size of the array, return 0.
	if(this->GetNumberOfElements()!=number_of_elements_of_array) return 0;


	//Fill in the whole array
	double dec = 1;
	double value = 0;
	int number_of_added_elements = 0;
	int is_dot_found = 0;
	int is_sign_found = 0;
	for(int r=i2 ; r>=i-1; r--)
	{
		if(array_string[r]==',' || array_string[r]=='[')
		{
			//Write the calculated value in the array
			this->m_pointer[this->GetNumberOfElements()-1-number_of_added_elements] = value;
			number_of_added_elements++;

			//Reset 'dec' and 'value' for a new number 
			dec = 1;
			value = 0;

			//Set the dot indicator to 0
			is_dot_found = 0;

			//Set the sign indicator to 0
			is_sign_found = 0;
		}
		//Since we search backwards, the numbers and dot can be taken into account only if the sign has not yet been found!
		if(!is_sign_found)
		{
			//If MINUS sign is found, make negative value, set the indicator of found sign
			if(array_string[r]=='-')
			{
				value = -value;
				is_sign_found = 1;
			}
			//If PLUS sign is found, just set the indicator of found sign
			if(array_string[r]=='+') is_sign_found = 1;

			//If a number is found
			if(array_string[r]>='0' && array_string[r]<='9')
			{
				value += ((int)(array_string[r] - '0'))*dec;
				dec = dec * 10;
			}
			if(array_string[r]=='.')
			{
				//If this is the first dot
				if(!is_dot_found)
				{
					value = value/dec;
					is_dot_found = 1;
					dec = 1;
				}
				//If this is the second dot, the error ocurred, exit with 0
				else
				{
					cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect '.' character found!"<<endl;
					return 0;
				}
			}
		}
		else
		{
			//If the sign is found and the current character is not ',' or '[', report an error
			if(array_string[r]!=',' && array_string[r]!='[' && array_string[r]!=' ')
			{
				cout<<"bdArrayOfDoubles::LoadFromString_DoubleValues_ForPreSetArraySize(): Error in reading bdArray file! Incorrect character found before the sign!"<<endl;
				return 0;
			}
		}
	}

	return 1;
}
