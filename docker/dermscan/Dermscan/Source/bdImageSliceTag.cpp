/**************************************************************************
 Class of 2D Image (slice) tags

 Author: Danilo Babin
 File name: "bdImageSliceTag.cpp"
***************************************************************************/


#include "bdImageSliceTag.h"


bdArray<bdImageSliceTagNumber> bdImageSliceTag::m_dicom_tag_names;
//bdArray<bdString> bdImageSliceTag::m_siemens_tag_names;



int bdImageSliceTagValue::GetValueAsDouble(double &output)
{
    return m_value.StringToDouble(output);
}


int bdImageSliceTagValue::GetValuesAsArrayOfDoubles(bdArray<double> &output_values)
{
    bdList<double> list_of_extracted_doubles;
    m_value.ExtractAllDoubleNumbersToList(list_of_extracted_doubles);
    output_values.Set(list_of_extracted_doubles.GetNumberOfElements());
    bdListIterator<double> itd;
    unsigned int i = 0;
    for(itd.SetLeftEnd(list_of_extracted_doubles); itd.IsValid(); itd.MoveRight())
    {
        output_values[i] = itd.GetElement();
        i++;
    }
    return output_values.GetNumberOfElements();
}



//------------------------------------------------------------------------------------------------------------------------------------------------------




bdImageSliceTag::bdImageSliceTag()
{
	if(m_dicom_tag_names.IsEmpty())
	{
		m_dicom_tag_names.Set(8);//if you change the value here, do not forget to change it bellow for values array!
		int i=0;
//		m_dicom_tag_names[i].m_tag_number1.Append("0008"); m_dicom_tag_names[i].m_tag_number2.Append("0013");
//        m_dicom_tag_names[i].m_tag_description.Append("Creation Time"); i++;
//        m_dicom_tag_names[i].m_tag_number1.Append("0008"); m_dicom_tag_names[i].m_tag_number2.Append("1140");
//        m_dicom_tag_names[i].m_tag_description.Append("Referenced Image Sequence"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0018"); m_dicom_tag_names[i].m_tag_number2.Append("1060");
        m_dicom_tag_names[i].m_tag_description.Append("Trigger Time"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0020"); m_dicom_tag_names[i].m_tag_number2.Append("0013");
        m_dicom_tag_names[i].m_tag_description.Append("Image Number"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0020"); m_dicom_tag_names[i].m_tag_number2.Append("0032");
        m_dicom_tag_names[i].m_tag_description.Append("Patient Position"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0020"); m_dicom_tag_names[i].m_tag_number2.Append("0037");
        m_dicom_tag_names[i].m_tag_description.Append("Image Orientation RC"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0020"); m_dicom_tag_names[i].m_tag_number2.Append("1041");
        m_dicom_tag_names[i].m_tag_description.Append("Slice Location"); i++;
        
        m_dicom_tag_names[i].m_tag_number1.Append("0020"); m_dicom_tag_names[i].m_tag_number2.Append("0052");
        m_dicom_tag_names[i].m_tag_description.Append("Frame of Reference UID"); i++;

        m_dicom_tag_names[i].m_tag_number1.Append("0008"); m_dicom_tag_names[i].m_tag_number2.Append("1150");
        m_dicom_tag_names[i].m_tag_description.Append("Referenced SOP Class UID"); i++;

        m_dicom_tag_names[i].m_tag_number1.Append("0008"); m_dicom_tag_names[i].m_tag_number2.Append("1155");
        m_dicom_tag_names[i].m_tag_description.Append("Referenced SOP Instance UID"); i++;

	}

//	if(m_siemens_tag_names.IsEmpty())
//	{
//		m_siemens_tag_names.Set(8);//if you change the value here, do not forget to change it bellow for values array!
//		int i=0;
//		m_siemens_tag_names[i].Append("dPhaseFOV"); i++;// Siemens MR CSA Series Header Filed Of View for rows.
//		m_siemens_tag_names[i].Append("dReadoutFOV"); i++;// Siemens MR CSA Series Header Filed Of View for columns.
//		m_siemens_tag_names[i].Append("sPosition.dSag"); i++;// Siemens MR CSA Series Header Slice position Sagittal.
//		m_siemens_tag_names[i].Append("sPosition.dCor"); i++;// Siemens MR CSA Series Header Slice position Coronal.
//		m_siemens_tag_names[i].Append("sPosition.dTra"); i++;// Siemens MR CSA Series Header Slice position Transversal.
//		m_siemens_tag_names[i].Append("sNormal.dSag"); i++;// Siemens MR CSA Series Header Slice normal Sagittal.
//		m_siemens_tag_names[i].Append("sNormal.dCor"); i++;// Siemens MR CSA Series Header Slice normal Coronal.
//		m_siemens_tag_names[i].Append("sNormal.dTra"); i++;// Siemens MR CSA Series Header Slice normal Transversal.
//	}

	// Set values for dicom
	m_dicom_tag_values.Set(m_dicom_tag_names.GetNumberOfElements());

    //	for(unsigned int i=0; i<m_dicom_tag_values.GetNumberOfElements(); i++) { m_dicom_tag_values[i] = 0; }

//	// Set values for siemens
//	m_is_siemens_dicom = 0;
//	m_siemens_tag_values.Set(8);
//	for(unsigned int i=0; i<m_siemens_tag_values.GetNumberOfElements(); i++) { m_siemens_tag_values[i] = 0; }
}


bdImageSliceTag::~bdImageSliceTag()
{
}


bdImageSliceTagValue* bdImageSliceTag::Value(bdString &tag_number1, bdString &tag_number2)
{
	int index = this->IndexOfDicomTagNumber(tag_number1,tag_number2);
    if(index>=0)
    {
        return(&(this->m_dicom_tag_values[index]));
    }
	return NULL;
}


bdImageSliceTagValue* bdImageSliceTag::Value(const char *tag_number1, const char *tag_number2)
{
    int index = this->IndexOfDicomTagNumber(tag_number1,tag_number2);
    if(index>=0)
    {
        return(&(this->m_dicom_tag_values[index]));
    }
    return NULL;
}


bdString* bdImageSliceTag::ValueString(bdString &tag_number1, bdString &tag_number2)
{
    bdImageSliceTagValue *v = this->Value(tag_number1,tag_number2);
    if(v) return &(v->m_value);
    return NULL;
}

bdString* bdImageSliceTag::ValueString(const char *tag_number1, const char *tag_number2)
{
    bdImageSliceTagValue *v = this->Value(tag_number1,tag_number2);
    if(v) return &(v->m_value);
    return NULL;
}


int bdImageSliceTag::ValueDoubles(bdString &tag_number1, bdString &tag_number2, bdArray<double> &output_values)
{
    bdImageSliceTagValue *v = this->Value(tag_number1,tag_number2);
    if(v) { return v->GetValuesAsArrayOfDoubles(output_values); }
    //cout<<"bdImageSliceTag::ValueDoubles(): Error geting bdImageSliceTagValue for tags "<<tag_number1<<","<<tag_number2<<endl;
    return 0;
}


int bdImageSliceTag::ValueDoubles(const char *tag_number1, const char *tag_number2, bdArray<double> &output_values)
{
    bdImageSliceTagValue *v = this->Value(tag_number1,tag_number2);
    if(v) { return v->GetValuesAsArrayOfDoubles(output_values); }
    return 0;
}



//double& bdImageSliceTag::ValueDicom(int index)
//{
//	return m_dicom_tag_values[index];
//}

	
int bdImageSliceTag::IndexOfDicomTagNumber(bdString &tag_number1, bdString &tag_number2)
{
    for(unsigned int i=0; i<m_dicom_tag_names.GetNumberOfElements(); i++)
    {
        if((m_dicom_tag_names[i].m_tag_number1==tag_number1) && (m_dicom_tag_names[i].m_tag_number2==tag_number2))
            return i;
    }
    return -1;
}


int bdImageSliceTag::IndexOfDicomTagNumber(const char *tag_number1, const char *tag_number2)
{
    for(unsigned int i=0; i<m_dicom_tag_names.GetNumberOfElements(); i++)
    {
        if((m_dicom_tag_names[i].m_tag_number1==tag_number1) && (m_dicom_tag_names[i].m_tag_number2==tag_number2))
            return i;
    }
    return -1;
}


//int bdImageSliceTag::IndexOfDicomTagHexadecimal(bdString &tag_hexadecimal)
//{
//	for(unsigned int i=0; i<m_dicom_tag_names.GetNumberOfElements(); i++)
//	{
//		bdString t;
//		m_dicom_tag_names[i].ExtractStringBetweenCharacters('(',')',t);
//		if(t==tag_hexadecimal) return i;
//	}
//	return -1;
//}
//
//
//int bdImageSliceTag::IndexOfDicomTagHexadecimal(const char *tag_hexadecimal)
//{
//	bdString t; t.Append(tag_hexadecimal);
//	return this->IndexOfDicomTagHexadecimal(t);
//}
//
//	
//double& bdImageSliceTag::ValueSiemens(int index)
//{
//	return m_siemens_tag_values[index];
//}
//
//
//int bdImageSliceTag::IndexOfSiemensTagName(bdString &tag_name)
//{
//    bdArraySortAndSearch<bdString> ss;
//	return (ss.IndexOfElement(m_siemens_tag_names, tag_name));
//}
//
//
//int bdImageSliceTag::IndexOfSiemensTagName(const char *tag_name)
//{
//	bdString t; t.Append(tag_name);
//	return this->IndexOfDicomTagName(t);
//}


void bdImageSliceTag::Print(bdString &s)
{
	bdString temp;
	for(unsigned int i=0; i<m_dicom_tag_names.GetNumberOfElements(); i++)
	{
		s.Append("("); s.Append(m_dicom_tag_names[i].m_tag_number1); s.Append(","); s.Append(m_dicom_tag_names[i].m_tag_number2); s.Append(") ");
        s.Append(m_dicom_tag_names[i].m_tag_description); s.Append(": ");
        s.Append(m_dicom_tag_values[i].m_value); s.Append("\n");
	}

//	if(m_is_siemens_dicom)
//	{
//		s.Append("\nSiemens CSA Series MR Header:\n");
//		for(unsigned int i=0; i<m_siemens_tag_names.GetNumberOfElements(); i++)
//		{
//			s.Append("  "); s.Append(m_siemens_tag_names[i]); s.Append(": "); temp.NumberToString(m_siemens_tag_values[i]); s.Append(temp); s.Append("\n");
//		}
//	}
}


void bdImageSliceTag::SaveToString(bdString &s)
{
	bdString temp;
	for(unsigned int i=0; i<m_dicom_tag_names.GetNumberOfElements(); i++)
	{
        s.Append("<");
        s.Append("("); s.Append(m_dicom_tag_names[i].m_tag_number1); s.Append(","); s.Append(m_dicom_tag_names[i].m_tag_number2); s.Append(")");
        s.Append("["); s.Append(m_dicom_tag_names[i].m_tag_description); s.Append("]");
        m_dicom_tag_values[i].m_value.ReplaceCharacter('\n', ' ');
        s.Append("{"); s.Append(m_dicom_tag_values[i].m_value); s.Append("}");
        s.Append(">");
        
//        s.Append("$["); s.Append(m_dicom_tag_names[i].m_tag_description);
//        s.Append(" ("); s.Append(m_dicom_tag_names[i].m_tag_number1); s.Append(","); s.Append(m_dicom_tag_names[i].m_tag_number2);
//        s.Append(")]{"); s.Append(m_dicom_tag_values[i].m_value); s.Append("}$");
////        s.Append("$["); s.Append(m_dicom_tag_names[i]); s.Append("]{"); temp.NumberToString(m_dicom_tag_values[i]); s.Append(temp); s.Append("}$");
	}

//	if(m_is_siemens_dicom)
//	{
//		//s.Append("$[Siemens]$");
//		for(unsigned int i=0; i<m_siemens_tag_names.GetNumberOfElements(); i++)
//		{
//			s.Append("$["); s.Append(m_siemens_tag_names[i]); s.Append("]{"); temp.NumberToString(m_siemens_tag_values[i]); s.Append(temp); s.Append("}$");
//		}
//	}
}


int bdImageSliceTag::LoadFromSavedString(bdString &s)
{
	if(s.IsEmpty()) return 0;

	bdList<bdString> list_of_tags;
	s.ExtractOptions('<','>',list_of_tags);
	bdListIterator<bdString> it;
	for(it.SetLeftEnd(list_of_tags); it.IsValid(); it.MoveRight())
	{
        bdString tag_numbers;
        it.GetElementPointer()->ExtractStringBetweenCharacters('(',')',tag_numbers);
        
        bdList<bdString> tag_numbers_list;
        tag_numbers.ExtractHexadecimalIntNumberStringsToList(tag_numbers_list);
        if(tag_numbers_list.GetNumberOfElements()!=2) return 0;
        
        // find the tag numbers, and if regular add value.
        int index = this->IndexOfDicomTagNumber(tag_numbers_list.GetLeftEnd(), tag_numbers_list.GetRightEnd());
        if(index>=0)// If it is a regular DICOM tag name
        {
            if(!it.GetElementPointer()->ExtractStringBetweenCharacters('{','}',m_dicom_tag_values[index].m_value)) return 0;
        }
        
//		// Find the tag name
//        bdArraySortAndSearch<bdString> ss;
//		int index = ss.IndexOfElement(m_dicom_tag_names,tag_name);
//		if(index>=0)// If it is a regular DICOM tag name
//		{
//			double number;
//			tag_value.StringToDouble(number);
//			m_dicom_tag_values[index] = number;
//		}
//		else //if not found in regular dicom fields, search the Siemens dicom fields...
//		{
//			index = ss.IndexOfElement(m_siemens_tag_names,tag_name);
//			if(index>=0)// If it is a Siemens DICOM tag name
//			{
//				m_is_siemens_dicom = 1;
//				double number;
//				tag_value.StringToDouble(number);
//				m_siemens_tag_values[index] = number;
//			}
//			//else std::cout<<"NO ";
//		}
	}

	return 1;
}


std::ostream& operator <<(std::ostream &o, const bdImageSliceTag &t)
{

    for(unsigned int i=0; i<t.m_dicom_tag_names.GetNumberOfElements(); i++)
    {
        o<<"("<<t.m_dicom_tag_names[i].m_tag_number1<<","<<t.m_dicom_tag_names[i].m_tag_number2<<") "<<t.m_dicom_tag_names[i].m_tag_description<<": "<<t.m_dicom_tag_values[i].m_value<<"\n";
    }
//	for(unsigned int i=0; i<t.m_dicom_tag_names.GetNumberOfElements(); i++)
//	{
//		o<<(t.m_dicom_tag_names[i])<<": "<<(t.m_dicom_tag_values[i])<<endl;
//	}

//	if(t.m_is_siemens_dicom)
//	{
//		o<<endl<<"Siemens CSA Series MR Header:"<<endl;
//		for(unsigned int i=0; i<t.m_siemens_tag_names.GetNumberOfElements(); i++)
//		{
//			o<<"  "<<(t.m_siemens_tag_names[i])<<": "<<(t.m_siemens_tag_values[i])<<endl;
//		}
//	}

    return o;
}


std::stringstream& operator << (std::stringstream &o, bdImageSliceTag &t)
{
    for(unsigned int i=0; i<t.m_dicom_tag_names.GetNumberOfElements(); i++)
    {
        o<<"("<<t.m_dicom_tag_names[i].m_tag_number1<<","<<t.m_dicom_tag_names[i].m_tag_number2<<") "<<t.m_dicom_tag_names[i].m_tag_description<<": "<<t.m_dicom_tag_values[i].m_value<<"\n";
    }

//	if(t.m_is_siemens_dicom)
//	{
//		o<<endl<<"Siemens CSA Series MR Header:"<<endl;
//		for(unsigned int i=0; i<t.m_siemens_tag_names.GetNumberOfElements(); i++)
//		{
//			o<<"  "<<(t.m_siemens_tag_names[i])<<": "<<(t.m_siemens_tag_values[i])<<endl;
//		}
//	}

    return o;
}