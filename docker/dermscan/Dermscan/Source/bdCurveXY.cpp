/**************************************************************************
 XY Curve for representing 2D signals/functions.

 Author: Danilo Babin
 File name: "bdCurveXY.h"
***************************************************************************/


#include "bdCurveXY.h"



bdCurveXY::bdCurveXY()
{
    m_x_axis_label.Assign("X Axis");
    m_y_axis_label.Assign("Y Axis");
}


bdCurveXY::~bdCurveXY()
{
	this->Reset();
}


int bdCurveXY::IsEmpty()
{
	if(m_x.IsEmpty() || (m_y.GetNumberOfElements()!=m_x.GetNumberOfElements())) return 1;
	else return 0;
}


void bdCurveXY::Reset()
{
	m_x.Reset();
	m_y.Reset();
}


unsigned int bdCurveXY::GetNumberOfElements() 
{
	return m_x.GetNumberOfElements();
}


bdString bdCurveXY::GetAxis_X_Label()
{
    return m_x_axis_label;
}


bdString bdCurveXY::GetAxis_Y_Label()
{
    return m_y_axis_label;
}


void bdCurveXY::SetAxis_X_Label(bdString &label)
{
    m_x_axis_label.Assign(label);
}


void bdCurveXY::SetAxis_X_Label(const char *label)
{
    m_x_axis_label.Assign(label);
}


void bdCurveXY::SetAxis_Y_Label(bdString &label)
{
    m_y_axis_label.Assign(label);
}


void bdCurveXY::SetAxis_Y_Label(const char *label)
{
    m_y_axis_label.Assign(label);    
}


bdArrayOfDoubles* bdCurveXY::Get_X_Array()
{
	return &m_x;
}


bdArrayOfDoubles* bdCurveXY::Get_Y_Array()
{
	return &m_y;
}


double* bdCurveXY::GetPointerTo_X_Data()
{
	return (m_x.GetPointerToData());
}


double* bdCurveXY::GetPointerTo_Y_Data()
{
	return (m_y.GetPointerToData());
}


int bdCurveXY::X_Bounds(double &x_min, double &x_max)
{
    if(this->IsEmpty()) return 0;
    x_min = X(0);
    x_max = x_min;
    for(int i=1; i<this->GetNumberOfElements(); i++)
    {
        if(X(i)<x_min) x_min = X(i);
        if(X(i)>x_max) x_max = X(i);
    }
    return 1;
}


int bdCurveXY::Y_Bounds(double &y_min, double &y_max)
{
    if(this->IsEmpty()) return 0;
    y_min = Y(0);
    y_max = y_min;
    for(int i=1; i<this->GetNumberOfElements(); i++)
    {
        if(Y(i)<y_min) y_min = Y(i);
        if(Y(i)>y_max) y_max = Y(i);
    }
    return 1;
}


void bdCurveXY::Set(unsigned int r)
{
	m_x.Set(r); m_y.Set(r);
}



void bdCurveXY::MakeStandard_X_Values(double start_value, double increment)
{
	//assert(!this->x.Empty());

	m_x[0] = start_value;
	for(unsigned int r=1; r<m_x.GetNumberOfElements(); r++)
	{ m_x[r] = m_x[r-1] + increment; }
}



int bdCurveXY::CopyFrom(bdCurveXY *a)
{
    if(!a) return 0;
    if(a==this) return 1;
    this->Reset();
    bdDataObject::CopyFrom(a);
    m_x = a->m_x;
    m_y = a->m_y;
    m_x_axis_label = a->m_x_axis_label;
    m_y_axis_label = a->m_y_axis_label;
    return 1;
}



int bdCurveXY::CopyFrom(bdCurveXY &a)
{
    return (this->CopyFrom(&a));
}


bdCurveXY& bdCurveXY::operator =(bdCurveXY &m)
{
	this->CopyFrom(m);
	return(*this);
}


void bdCurveXY::Negative()
{
    for(int j=0; j<this->GetNumberOfElements(); j++)
    {
        this->Y(j) = - this->Y(j);
    }
}


void bdCurveXY::SetValues(unsigned int index, double x_value, double y_value)
{
	m_x[index]=x_value; m_y[index]=y_value;
}


double& bdCurveXY::Y(unsigned int index)
{
	return m_y[index];
}


double bdCurveXY::Y_FromX(double x_value)
{
    // Find the closest samples (on x axis) in the input cuve
    int sample_index = -1;//set to impossible value to indicate if the index is found or not.
    for(int j=0; j<this->GetNumberOfElements(); j++)
    {
        if(this->X(j)>x_value)
        {
            sample_index = j;
            break;
        }
    }

    //Calculate the output sample value:
    //-- 1. sample is on the left, the signal is 0 --
    if(sample_index==0)
    {
        return 0;
    }
    //-- 2. sample is on the right, the signal is 0 --
    if(sample_index<0)
    {
        return 0;
    }
    //-- 3. sample is in between 2 existing samples, we do interpolation --
    //if(sample_index>0)
    //{
    return ( (this->Y(sample_index)- this->Y(sample_index-1)) * (x_value-this->X(sample_index-1)) / (this->X(sample_index)-this->X(sample_index-1)) + this->Y(sample_index-1) );
    //}
}


unsigned int bdCurveXY::IndexOfX(double x_value)
{
    for(int j=0; j<this->GetNumberOfElements(); j++)
    {
        if(this->X(j)>x_value) { return j; }
    }
    return 0;
}


double& bdCurveXY::X(unsigned int index)
{
	return m_x[index];
}


void bdCurveXY::AddConstantValue(double y_const)
{
    for(unsigned int i=0; i<this->GetNumberOfElements(); i++) { Y(i) += y_const; }
}


int bdCurveXY::Max_Y_Value(double &max_y_value, unsigned int *index_of_max_y_value)
{
    if(this->IsEmpty()) return 0;
    max_y_value = Y(0);
    for(unsigned int i=0; i<this->GetNumberOfElements(); i++)
    {
        if(Y(i) > max_y_value)
        {
            max_y_value = Y(i);
            if(index_of_max_y_value) *index_of_max_y_value = i;
        }
    }
    return 1;
}


int bdCurveXY::Min_Y_Value(double &min_y_value, unsigned int *index_of_min_y_value)
{
    if(this->IsEmpty()) return 0;
    min_y_value = Y(0);
    for(unsigned int i=0; i<this->GetNumberOfElements(); i++)
    {
        if(Y(i) < min_y_value)
        {
            min_y_value = Y(i);
            if(index_of_min_y_value) *index_of_min_y_value = i;
        }
    }
    return 1;
}


int bdCurveXY::Max_Y_ConvexValue(double &max_y_value, unsigned int *index_of_max_y_value)
{
    if(this->GetNumberOfElements()<3) return 0;
    
    max_y_value = 0;
    int is_max_value_candidate_found = 0;
    for(unsigned int i=1; i<this->GetNumberOfElements()-1; i++)
    {
        if((Y(i) >= Y(i-1)) && (Y(i) >= Y(i+1)))//if a convex point was found
        {
            if(Y(i)>max_y_value || (!is_max_value_candidate_found))
            {
                is_max_value_candidate_found = 1;
                max_y_value = Y(i);
                if(index_of_max_y_value) *index_of_max_y_value = i;
            }
        }
    }
    return is_max_value_candidate_found;
}


int bdCurveXY::Min_Y_ConvexValue(double &min_y_value, unsigned int *index_of_min_y_value)
{
    if(this->GetNumberOfElements()<3) return 0;
    
    min_y_value = 0;
    int is_min_value_candidate_found = 0;
    for(unsigned int i=1; i<this->GetNumberOfElements()-1; i++)
    {
        if((Y(i) <= Y(i-1)) && (Y(i) <= Y(i+1)))//if a convex point was found
        {
            if(Y(i)<min_y_value || (!is_min_value_candidate_found))
            {
                is_min_value_candidate_found = 1;
                min_y_value = Y(i);
                if(index_of_min_y_value) *index_of_min_y_value = i;
            }
        }
    }
    return is_min_value_candidate_found;
}



void bdCurveXY::Sort_X()
{
    bdArray<int> index_change_order;
    this->Get_X_Array()->SortAscending(index_change_order);
    bdArrayOfDoubles temp;
    temp = *(this->Get_Y_Array());
    for(unsigned int i=0; i<index_change_order.GetNumberOfElements(); i++)
    {
        Y(i) = temp[(index_change_order[i])];
    }
}


int bdCurveXY::RegroupXValuesToGivenNumberOfRanges(unsigned int number_of_ranges, bdCurveXY &output)
{
	if(this->IsEmpty()) return 0;
	if(this->GetNumberOfElements()<number_of_ranges) return 0;
	if(number_of_ranges<2) return 0;

	//int single_range = this->GetNumberOfElements()/number_of_ranges;
	//if(single_range<1) return 0;

	double single_range = (m_x[this->GetNumberOfElements()-1] - m_x[0]) / number_of_ranges;
	if(single_range<0) single_range = -single_range;
	if(single_range<1) return 0;



	output.Reset();
	output.Set(number_of_ranges);
	output.m_x.FillInWith(0);
	output.m_y.FillInWith(0);
	int j=0;
	for(unsigned int i=0; i<this->GetNumberOfElements(); i++)
	{
		j = (int) (m_x[i] / single_range);
		//if(output.m_x[i/single_range]<m_x[i]) output.m_x[i/single_range] = m_x[i];
		//output.m_y[i/single_range] += m_y[i];
		if(output.m_x[j]<m_x[i]) output.m_x[j] = m_x[i];
		output.m_y[j] += m_y[i];
	}
	
	return 1;
}


int bdCurveXY::ShiftElementsToRightCircular()
{
	if(this->IsEmpty()) return 0;
	if(this->GetNumberOfElements()==1) return 1;

	double dx = m_x[1] - m_x[0];
	double y_temp = m_y[(this->GetNumberOfElements()-1)];
	for(int i=this->GetNumberOfElements()-1; i>=1; i--)
	{
		m_x[i] = m_x[i-1]; m_y[i] = m_y[i-1];
	}
	m_x[0] = m_x[1] - dx;
	m_y[0] = y_temp;

	return 1;
}


int bdCurveXY::ShiftElementsToLeftCircular()
{
	if(this->IsEmpty()) return 0;
	if(this->GetNumberOfElements()==1) return 1;

	double dx = m_x[(this->GetNumberOfElements()-1)] - m_x[(this->GetNumberOfElements()-2)];
	double y_temp = m_y[0];
	for(unsigned int i=0; i<this->GetNumberOfElements()-1; i++)
	{
		m_x[i] = m_x[i+1]; m_y[i] = m_y[i+1];
	}
	m_x[(this->GetNumberOfElements()-1)] = m_x[(this->GetNumberOfElements()-2)] + dx;
	m_y[(this->GetNumberOfElements()-1)] = y_temp;

	return 1;
}


//void bdCurveXY::SaveToMatlab_M_File(const char *file_name)//, char *file_name_root)
//{
//    bdString bds_file, bds_dir, bds_file_name, bds_ext;
//    bds_file.Assign(file_name);
//    bds_file.ExtractFileNamePathParts(bds_dir, bds_file_name, bds_ext);
//    if(!(bds_ext=="m" || bds_ext=="M"))
//    {
//        bds_file.Append(".m");
//    }
//    
//	ofstream output_file;
//	output_file.open(bds_file.C_String(),ios::binary);
//
//	//Print out the X array
//	output_file<<"m_x = [";
//	
//	for(unsigned int i=0; i<m_x.GetNumberOfElements(); i++)
//	{
//		output_file<<m_x[i];
//		
//		//If the condition will be fullfiled
//		if(i==m_x.GetNumberOfElements()-1)
//		{
//			output_file<<"];"<<endl;
//		}
//		else
//		{
//			output_file<<", ";
//		}
//	}
//
//	//Print out the Y array
//    output_file<<"m_y = [";
//	
//	for(unsigned int i=0; i<m_y.GetNumberOfElements(); i++)
//	{
//		output_file<<m_y[i];
//		
//		//If the condition will be fullfiled
//		if(i==m_y.GetNumberOfElements()-1)
//		{
//			output_file<<"];"<<endl;
//		}
//		else
//		{
//			output_file<<", ";
//		}
//	}
//
//	output_file.close();
//}


void bdCurveXY::SaveToDefaultFileType(const char *file_path)
{
    this->SaveToSingleArrayMatlab_M_File(file_path);
}


void bdCurveXY::SaveToSingleArrayMatlab_M_File(const char *file_name)
{
    bdString bds_file, bds_dir, bds_file_name, bds_ext;
    bds_file.Assign(file_name);
    bds_file.ExtractFileNamePathParts(bds_dir, bds_file_name, bds_ext);
    if(!(bds_ext=="m" || bds_ext=="M"))
    {
        bds_file.Append(".m");
    }
    
	ofstream output_file;
	output_file.open(bds_file.C_String(),ios::binary);

	output_file<<"%Single Array File: "<<endl;
    
    output_file<<"%v"<<endl;// version is recorded
    output_file<<"%2.0"<<endl;

    output_file<<"%l"<<endl;// axes labels are recorded
    output_file<<"%"<<this->GetAxis_X_Label().C_String()<<endl;
    output_file<<"%"<<this->GetAxis_Y_Label().C_String()<<endl;
    
    output_file<<"%t"<<endl;// tag is recorded
    bdString bds_tag;
    this->SaveTags(bds_tag);
    output_file<<"%"<<bds_tag<<endl;


	//Print out the X array
	output_file<<"m_x = [";
	
	for(unsigned int i=0; i<m_x.GetNumberOfElements(); i++)
	{
		output_file<<m_x[i];
		
		//If the condition will be fullfiled
		if(i==m_x.GetNumberOfElements()-1) output_file<<"];"<<endl;
		else output_file<<", ";
	}

	//Print out the Y array
	output_file<<"m_y = [";
	
	for(unsigned int i=0; i<m_y.GetNumberOfElements(); i++)
	{
		output_file<<m_y[i];
		
		//If the condition will be fullfiled
		if(i==m_y.GetNumberOfElements()-1) output_file<<"];"<<endl;
		else output_file<<", ";
	}

	output_file.close();
}


int bdCurveXY::LoadSingleArrayFile_v3(const char *file_name)
{
    ifstream input_file;
    input_file.open(file_name,ios::binary);
    if(!input_file)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v3() : Unable to open single array file: "<<file_name<<endl;
        return 0;
    }
    
    
    char text_line[100000];//If you change this value, you must also change the amount of loaded characters, see below.
    unsigned int amount_of_characters_to_load = 100000;
    
    if(!input_file.getline(text_line, amount_of_characters_to_load)) return 0;
    if(text_line[0]!='%' || text_line[1]!='S' || text_line[2]!='i' || text_line[3]!='n' || text_line[4]!='g' || text_line[5]!='l' || text_line[6]!='e' ||
       text_line[7]!=' ' || text_line[8]!='A' || text_line[9]!='r' || text_line[10]!='r' || text_line[11]!='a' || text_line[12]!='y' ||
       text_line[13]!=' ' || text_line[14]!='F' || text_line[15]!='i' || text_line[16]!='l' || text_line[17]!='e' ||  text_line[18]!=':')
    { input_file.close(); return 0; }
    
    // load version
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%' && text_line[1]!='v') { input_file.close(); return 0; }
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    bdString bds;
    bds.Append(&(text_line[1]));
    bdList<int> list;
    bds.ExtractAllIntNumbersToList(list);
    if(list.IsEmpty())
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v3(): No version value found!"<<endl;
        input_file.close();
        return 0;
    }
    if(list.GetLeftEnd()>3)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v3(): Version higher than supported, version: "<<list.GetLeftEnd()<<endl;
        input_file.close();
        return 0;
    }
    
    // read axes labels.
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%' && text_line[1]!='l') { input_file.close(); return 0; }
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    this->SetAxis_X_Label(&(text_line[1]));
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    this->SetAxis_Y_Label(&(text_line[1]));
    
    // read object tag.
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%' && text_line[1]!='t') { input_file.close(); return 0; }
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    this->LoadTags(&(text_line[1]));
    
    //Load X array
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(!m_x.LoadFromString_DoubleValues(text_line))
    {
        //cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading X array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
        cout<<"bdCurveXY::LoadSingleArrayFile_v3(): error loading X array! File name: "<<file_name<<endl;
        cout<<"text_line = "<<text_line<<endl;
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    //Load Y array
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(!m_y.LoadFromString_DoubleValues(text_line))
    {
        //cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading Y array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
        cout<<"bdCurveXY::LoadSingleArrayFile_v3(): error loading Y array! File name: "<<file_name<<endl;
        cout<<"text_line = "<<text_line<<endl;
        //Delete both arrays
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    if(m_x.GetNumberOfElements()!=m_y.GetNumberOfElements())
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v3(): loaded X and Y arrays are not of the same size: X= "<<m_x.GetNumberOfElements()<<", Y= "<<m_y.GetNumberOfElements()<<endl;
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    if(m_x.GetNumberOfElements()==0)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v3()_v2: no elemetns loaded, array is empty!"<<endl;
        input_file.close();
        return 0;
    }
    
    input_file.close();
    return 1;
}



int bdCurveXY::LoadSingleArrayFile_v2(const char *file_name)
{
    ifstream input_file;
    input_file.open(file_name,ios::binary);
    if(!input_file)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v2() : Unable to open single array file: "<<file_name<<endl;
        return 0;
    }
    
    
    char text_line[100000];//If you change this value, you must also change the amount of loaded characters, see below.
    unsigned int amount_of_characters_to_load = 100000;
    
    if(!input_file.getline(text_line, amount_of_characters_to_load)) return 0;
    if(text_line[0]!='%' || text_line[1]!='S' || text_line[2]!='i' || text_line[3]!='n' || text_line[4]!='g' || text_line[5]!='l' || text_line[6]!='e' ||
       text_line[7]!=' ' || text_line[8]!='A' || text_line[9]!='r' || text_line[10]!='r' || text_line[11]!='a' || text_line[12]!='y' ||
       text_line[13]!=' ' || text_line[14]!='F' || text_line[15]!='i' || text_line[16]!='l' || text_line[17]!='e' ||  text_line[18]!=':')
        { input_file.close(); return 0; }
    
    // load version
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%' && text_line[1]!='v') { input_file.close(); return 0; }
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    bdString bds;
    bds.Append(&(text_line[1]));
    bdList<int> list;
    bds.ExtractAllIntNumbersToList(list);
    if(list.IsEmpty())
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v2(): No version value found!"<<endl;
        input_file.close();
        return 0;
    }
    if(list.GetLeftEnd()>2)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v2(): Version higher than supported, version: "<<list.GetLeftEnd()<<endl;
        input_file.close();
        return 0;
    }
    
    // read object tag.
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%' && text_line[1]!='t') { input_file.close(); return 0; }
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(text_line[0]!='%') { input_file.close(); return 0; }
    this->LoadTags(&(text_line[1]));
    
    //Load X array
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(!m_x.LoadFromString_DoubleValues(text_line))
    {
        //cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading X array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
        cout<<"bdCurveXY::LoadSingleArrayFile_v2(): error loading X array! File name: "<<file_name<<endl;
        cout<<"text_line = "<<text_line<<endl;
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    //Load Y array
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if(!m_y.LoadFromString_DoubleValues(text_line))
    {
        //cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading Y array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
        cout<<"bdCurveXY::LoadSingleArrayFile_v2(): error loading Y array! File name: "<<file_name<<endl;
        cout<<"text_line = "<<text_line<<endl;
        //Delete both arrays
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    if(m_x.GetNumberOfElements()!=m_y.GetNumberOfElements())
    {
        cout<<"bdCurveXY::LoadSingleArrayFile_v2(): loaded X and Y arrays are not of the same size: X= "<<m_x.GetNumberOfElements()<<", Y= "<<m_y.GetNumberOfElements()<<endl;
        m_x.Reset();
        m_y.Reset();
        input_file.close();
        return 0;
    }
    
    if(m_x.GetNumberOfElements()==0)
    {
        cout<<"bdCurveXY::LoadSingleArrayFile()_v2: no elemetns loaded, array is empty!"<<endl;
        input_file.close();
        return 0;
    }
    
    input_file.close();
    return 1;
}



int bdCurveXY::LoadSingleArrayFile_v1(const char *file_name)
{
	ifstream input_file;
	input_file.open(file_name,ios::binary);
    if(!input_file)
	{
		cout<<"bdCurveXY::LoadSingleArrayFile_v1() : Unable to open single array file: "<<file_name<<endl;
		return 0;
    }

	char text_line[100000];//If you change this value, you must also change the amount of loaded characters, see below.

	if(!input_file.getline(text_line, 99999)) return 0;
	if(text_line[0]!='%' || text_line[1]!='S' || text_line[2]!='i' || text_line[3]!='n' || text_line[4]!='g' || text_line[5]!='l' || text_line[6]!='e' ||
		text_line[7]!=' ' || text_line[8]!='A' || text_line[9]!='r' || text_line[10]!='r' || text_line[11]!='a' || text_line[12]!='y' ||
		text_line[13]!=' ' || text_line[14]!='F' || text_line[15]!='i' || text_line[16]!='l' || text_line[17]!='e' ||  text_line[18]!=':') return 0;

	//Load X array
	if(!input_file.getline(text_line, 99999)) return 0;
	if(!m_x.LoadFromString_DoubleValues(text_line))
	{
		//cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading X array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
		cout<<"bdCurveXY::LoadSingleArrayFile_v1(): error loading X array! File name: "<<file_name<<endl;
		cout<<"text_line = "<<text_line<<endl;
		m_x.Reset();
		m_y.Reset();
		return 0;
	}
	
	//Load Y array
	if(!input_file.getline(text_line, 99999)) return 0;
	if(!m_y.LoadFromString_DoubleValues(text_line))
	{
		//cout<<"Array2D<Tx,Ty>::LoadSingleArrayFile(): error loading Y array! File name: "<<file_name<<", array_name: "<<array_name<<endl;
		cout<<"bdCurveXY::LoadSingleArrayFile_v1(): error loading Y array! File name: "<<file_name<<endl;
		cout<<"text_line = "<<text_line<<endl;
		//Delete both arrays
		m_x.Reset();
		m_y.Reset();
		return 0;
	}

	if(m_x.GetNumberOfElements()!=m_y.GetNumberOfElements())
	{
		cout<<"bdCurveXY::LoadSingleArrayFile_v1(): loaded X and Y arrays are not of the same size: X= "<<m_x.GetNumberOfElements()<<", Y= "<<m_y.GetNumberOfElements()<<endl;
		m_x.Reset();
		m_y.Reset();
		return 0;
	}

	if(m_x.GetNumberOfElements()==0)
	{
		cout<<"bdCurveXY::LoadSingleArrayFile_v1(): no elemetns loaded, array is empty!"<<endl;
		return 0;
	}

	input_file.close();
	return 1;
}


int bdCurveXY::LoadSingleArrayFile(const char *file_name)
{
    if(this->LoadSingleArrayFile_v3(file_name)) return 1;
    if(this->LoadSingleArrayFile_v2(file_name)) return 1;
    if(this->LoadSingleArrayFile_v1(file_name)) return 1;
    return 0;
}


int bdCurveXY::LoadFromFile_DoubleValues(char *file_name, char *array_nameY, char *array_nameX)
{
	if(!m_y.LoadFromFile_DoubleValues(file_name, array_nameY))
	{
		cout<<"int bdCurveXY::LoadFromFile_DoubleValues(): error loading Y array! File name: "<<file_name<<", array_nameY: "<<array_nameY<<endl;
		return 0;
	}

	if(array_nameX==NULL)
	{
		m_x.Set(m_y.GetNumberOfElements());
		this->MakeStandard_X_Values();
	}
	else
	{
		if(!m_x.LoadFromFile_DoubleValues(file_name, array_nameX))
		{
			cout<<"int bdCurveXY::LoadFromFile_DoubleValues(): error loading X array! File name: "<<file_name<<", array_nameX: "<<array_nameX<<endl;
			return 0;
		}
	}
	return 1;
}
