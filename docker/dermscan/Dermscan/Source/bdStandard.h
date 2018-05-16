/**************************************************************************
 Various practical (extended standard) functions

 Author: Danilo Babin
 File name: "bdStandard.h"
***************************************************************************/





#ifndef BD_STANDARD_DEF
	#define BD_STANDARD_DEF
	

/// Way to turn definition values to strings.	
#ifndef BD_STRINGIZE_DEF
	#define BD_STRINGIZE_DEF
	#define BD_STRINGIZE2(s) #s
	#define BD_STRINGIZE(s) BD_STRINGIZE2(s)
#endif
	


#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>



namespace bdStandard
{
	template <typename T>
	inline bool FromString(const std::string &s, T &t)
	{
		std::istringstream is(s);
		return (is>>t) ? true : false;
	};

	template <typename T>
	inline std::string ToString(T t)
	{
		std::ostringstream os;
		os<<t;
		return os.str();
	};

	template <typename T>
	inline void ToString(T &t, std::string &output_string)
	{
		output_string.clear();
		std::ostringstream os;
		os<<t;
		output_string = os.str();
	};

//	//Convert integer to string. Returns 0 if 'string_size' is not large enough to store given int.
//	int IntToString(int i, char *output_string, int string_size);
//
//	//Checks if the strings are equal to character 0 or to given string_size.
//	int AreStringsEqual(const char *s1, const char *s2, int string_size);
//
//	//Calculates the number of non-zero bits in the INT number (the sign should not be taken into account).
//	int NumberOfNonZeroBits(int v);
//};


//int bdStandard::IntToString(int i, char *output_string, int string_size)
inline int IntToString(int i, char *output_string, int string_size)
{
	if(string_size<3) return 0;

	int b = 10;
	int index = 0;
	int d = i;

	//Determine the base starting value
	int n;
	for(n=2; (d/b) != 0; n++)
	{
		b = b*10;
	}
	b = b/10;

	//Check if the string size can handle this number
	if(n+2>string_size)//This condition should not be the same for negative and positive numbers, but here we take a more strict condition.
	{
		return 0;
	}

	if(i<0)
	{
		output_string[index] = '-';
		index++;
	}
	
	while(b!=0 && index<string_size-1)
	{
		int c = d/b;
		output_string[index] = '0' + ((char) c);
		if(b>1) d = (d/b)%(b/10);
		b = b/10;
		index++;
	}

	output_string[index] = 0;

	return 1;
};



//int bdStandard::AreStringsEqual(const char *s1, const char *s2, int string_size)
inline int AreStringsEqual(const char *s1, const char *s2, int string_size)
{
	if(string_size<2) return 0;
	for(int i=0; (s1[i]!=0 || s2[i]!=0) && i<string_size; i++)
	{
		if(s1[i]!=s2[i]) return 0;
	}
	return 1;
};

//int bdStandard::NumberOfNonZeroBits(int v)
inline int NumberOfNonZeroBits(int v)
{
	int x = v;
	int n_of_non_zero_bits = 0;
	while(x!=0)
	{
		if(x%2==1) n_of_non_zero_bits++;
		x = x/2;
	}
	return n_of_non_zero_bits;
};


};

#endif