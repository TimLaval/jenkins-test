/**************************************************************************
 Template bdImage for images up to 4D data

 Author: Danilo Babin
 File name: "bdImage.cpp"
***************************************************************************/



#include "bdImage.h"


template<class T>
bdImageT<T>::bdImageT()
{
	//m_current_index_ST[0] = m_current_index_ST[1] = 0;
	this->m_data = NULL;
	this->m_is_data_created_internally = 0;
	this->m_dimensions_CRST[0] = this->m_dimensions_CRST[1] = this->m_dimensions_CRST[2] = this->m_dimensions_CRST[3] = 0;
    this->m_orientation_CxCyCzRxRyRzSxSySz[0] = 1; this->m_orientation_CxCyCzRxRyRzSxSySz[1] = 0; this->m_orientation_CxCyCzRxRyRzSxSySz[2] = 0;
    this->m_orientation_CxCyCzRxRyRzSxSySz[3] = 0; this->m_orientation_CxCyCzRxRyRzSxSySz[4] = 1; this->m_orientation_CxCyCzRxRyRzSxSySz[5] = 0;
    this->m_orientation_CxCyCzRxRyRzSxSySz[6] = 0; this->m_orientation_CxCyCzRxRyRzSxSySz[7] = 0; this->m_orientation_CxCyCzRxRyRzSxSySz[8] = 1;
    this->m_origin_CRST[0] = this->m_origin_CRST[1] = this->m_origin_CRST[2] = this->m_origin_CRST[3] = 0;
	this->m_spacing_CRST[0] = this->m_spacing_CRST[1] = this->m_spacing_CRST[2] = this->m_spacing_CRST[3] = 1;
}


template<class T>
bdImageT<T>::~bdImageT()
{
	this->Reset();
}


template<class T>
void bdImageT<T>::CopyFrom(bdImageT<T> &m)
{
	if(m.IsEmpty()) return;
	bdRegularGridT<T>::CopyFrom(m);
	//this->SetSize(m.GetNumberOfTimeSeries(),m.GetNumberOfSlices(),m.GetNumberOfRows(),m.GetNumberOfColumns());
	//this->SetVisualizationPropertiesAs(m);
	for(unsigned int i=0; i<this->GetNumberOfDataElements(); i++){ this->m_data[i] = m.m_data[i]; }
	m_tags = m.m_tags;
}


template<class T>
bdImageT<T>& bdImageT<T>::operator =(bdImageT<T> &m)
{
	this->CopyFrom(m);
	return(*this);
}


template<>
void bdImageT<unsigned char>::GetVoxelValueFullRange(unsigned int *full_range_min, unsigned int *full_range_max)
{
	*full_range_min = 0; *full_range_max = 255;
}


template<>
void bdImageT<unsigned short>::GetVoxelValueFullRange(unsigned int *full_range_min, unsigned int *full_range_max)
{
	*full_range_min = 0; *full_range_max = 65535;
}


template<class T>
int bdImageT<T>::GetVoxelMinimumAndMaximumValue(unsigned int *min, unsigned int *max)
{
	if(this->IsEmpty()) return 0;

	*min = *max = this->m_data[0];
	int n_of_elements = this->GetNumberOfDataElements();
	for(int i=0; i<n_of_elements; i++)
	{
		if(*min > this->m_data[i]) *min = this->m_data[i];
		if(*max < this->m_data[i]) *max = this->m_data[i];
	}

	return 1;
}


template<class T>
int bdImageT<T>::GetVoxelValueForWorldCoordinates(double w_t, double w_s, double w_r, double w_c, T &out_value)
{
	int index_t = (int)((w_t - this->m_origin_CRST[3]) / this->m_spacing_CRST[3]);
	int index_s = (int)((w_s - this->m_origin_CRST[2]) / this->m_spacing_CRST[2]);
	int index_r = (int)((w_r - this->m_origin_CRST[1]) / this->m_spacing_CRST[1]);
	int index_c = (int)((w_c - this->m_origin_CRST[0]) / this->m_spacing_CRST[0]);

	if(!this->GetIndexesForWorldCoordinates(w_t, w_s, w_r, w_c, index_t, index_s, index_r, index_c)) return 0;
	//if(index_t<0 || index_t>= this->GetNumberOfTimeSeries()) return 0;
	//if(index_s<0 || index_s>= this->GetNumberOfSlices()) return 0;
	//if(index_r<0 || index_r>= this->GetNumberOfRows()) return 0;
	//if(index_c<0 || index_c>= this->GetNumberOfColumns()) return 0;

	out_value = (*this)(index_t,index_s,index_r,index_c);
	return 1;
}


template<class T>
int bdImageT<T>::GetVoxelValueForWorldCoordinates(double w_s, double w_r, double w_c, T &out_value)
{
	return this->GetVoxelValueForWorldCoordinates(this->m_origin_CRST[3], w_s, w_r, w_c, out_value);
}


template<class T>
bdImageSliceTag& bdImageT<T>::Tag(unsigned int t, unsigned int s)
{
    if((t*this->GetNumberOfSlices() + s)>=m_tags.GetNumberOfElements()) return m_tag_default;
	return m_tags[(t*this->GetNumberOfSlices() + s)];
}


template<class T>
bdImageSliceTag& bdImageT<T>::Tag(unsigned int s)
{
	return m_tags[s];
}


template<class T>
int bdImageT<T>::SetSize(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(!bdRegularGridT<T>::SetSize(t,s,r,c)) return 0;
	m_tags.Set(t*s);
	return 1;
}


template<class T>
void bdImageT<T>::SetSizeAndPropertiesAs(bdImageT<T> &image)
{
	bdRegularGridT<T>::SetSizeAndPropertiesAs(image);
	//m_tags.Set(image.GetNumberOfSlices()*image.GetNumberOfTimeSeries());
    m_tags = image.m_tags;
}


template<class T>
int bdImageT<T>::SetSizeOf2DTimeSeries(unsigned int t, unsigned int r, unsigned int c)
{
	if(!bdRegularGridT<T>::SetSizeOf2DTimeSeries(t,r,c)) return 0;
	m_tags.Set(t);
	return 1;
}


template<class T>
int bdImageT<T>::SetSizeOf3DGrid(unsigned int s, unsigned int r, unsigned int c)
{
	if(!bdRegularGridT<T>::SetSizeOf3DGrid(s,r,c)) return 0;
	m_tags.Set(s);
	return 1;
}


template<class T>
void bdImageT<T>::SetSizeOnlyAs(bdImageT<T> &image)
{
	bdRegularGridT<T>::SetSizeOnlyAs(image);
	m_tags.Set(image.GetNumberOfSlices()*image.GetNumberOfTimeSeries());
}


template<class T>
void bdImageT<T>::Reset()
{
	bdRegularGridT<T>::Reset();
	m_tags.Reset();
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_8_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(this->IsEmpty()) return 0;
	int n = 0;

	if(r>=1)//if(r-1>=0)
	{
		if(c>=1) {if((*this)(t,s,r-1,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r-1,c+1)) n++;}//if(c+1<this->GetNumberOfColumns()) ...
		if((*this)(t,s,r-1,c)) n++;
	}
	if(r<this->GetNumberOfRows()-1)
	{
		if(c>=1) {if((*this)(t,s,r+1,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r+1,c+1)) n++;}
		if((*this)(t,s,r+1,c)) n++;
	}
	if(c>=1) {if((*this)(t,s,r,c-1)) n++;}
	if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r,c+1)) n++;}

	return n;
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_8_Neighborhood(unsigned int s, unsigned int r, unsigned int c)
{
	return this->NumberOfNonZeroVexelsIn_8_Neighborhood(0,s,r,c);
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_9_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	int n = this->NumberOfNonZeroVexelsIn_8_Neighborhood(t,s,r,c);
	if((*this)(t,s,r,c)) n++;
	return n;
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_9_Neighborhood(unsigned int s, unsigned int r, unsigned int c)
{
	return this->NumberOfNonZeroVexelsIn_9_Neighborhood(0,s,r,c);
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_26_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	if(this->IsEmpty()) return 0;
	int n = 0;

	if(s>=1)
	{
		if(r>=1)
		{
			if(c>=1) {if((*this)(t,s-1,r-1,c-1)) n++;}
			if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s-1,r-1,c+1)) n++;}
			if((*this)(t,s-1,r-1,c)) n++;
		}
		if(r<this->GetNumberOfRows()-1)
		{
			if(c>=1) {if((*this)(t,s-1,r+1,c-1)) n++;}
			if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s-1,r+1,c+1)) n++;}
			if((*this)(t,s-1,r+1,c)) n++;
		}
		if(c>=1) {if((*this)(t,s-1,r,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s-1,r,c+1)) n++;}
		if((*this)(t,s-1,r,c)) n++;
	}
	if(s<this->GetNumberOfSlices()-1)
	{
		if(r>=1)
		{
			if(c>=1) {if((*this)(t,s+1,r-1,c-1)) n++;}
			if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s+1,r-1,c+1)) n++;}
			if((*this)(t,s+1,r-1,c)) n++;
		}
		if(r<this->GetNumberOfRows()-1)
		{
			if(c>=1) {if((*this)(t,s+1,r+1,c-1)) n++;}
			if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s+1,r+1,c+1)) n++;}
			if((*this)(t,s+1,r+1,c)) n++;
		}
		if(c>=1) {if((*this)(t,s+1,r,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s+1,r,c+1)) n++;}
		if((*this)(t,s+1,r,c)) n++;
	}
	if(r>=1)
	{
		if(c>=1) {if((*this)(t,s,r-1,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r-1,c+1)) n++;}
		if((*this)(t,s,r-1,c)) n++;
	}
	if(r<this->GetNumberOfRows()-1)
	{
		if(c>=1) {if((*this)(t,s,r+1,c-1)) n++;}
		if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r+1,c+1)) n++;}
		if((*this)(t,s,r+1,c)) n++;
	}
	if(c>=1) {if((*this)(t,s,r,c-1)) n++;}
	if(c<this->GetNumberOfColumns()-1) {if((*this)(t,s,r,c+1)) n++;}

	return n;
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_26_Neighborhood(unsigned int s, unsigned int r, unsigned int c)
{
	//return this->NumberOfNonZeroVexelsIn_26_Neighborhood(m_current_index_ST[1],s,r,c);
	return this->NumberOfNonZeroVexelsIn_26_Neighborhood(0,s,r,c);
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_27_Neighborhood(unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	int n = this->NumberOfNonZeroVexelsIn_26_Neighborhood(t,s,r,c);
	if((*this)(t,s,r,c)) n++;
	return n;
}


template<class T>
int bdImageT<T>::NumberOfNonZeroVexelsIn_27_Neighborhood(unsigned int s, unsigned int r, unsigned int c)
{
	//return this->NumberOfNonZeroVexelsIn_27_Neighborhood(m_current_index_ST[1],s,r,c);
	return this->NumberOfNonZeroVexelsIn_27_Neighborhood(0,s,r,c);
}


template<class T>
void bdImageT<T>::SortSlices(const char *criteria1_tag_number1, const char *criteria1_tag_number2, const char *criteria2_tag_number1, const char *criteria2_tag_number2, const char *criteria3_tag_number1, const char *criteria3_tag_number2)
{
	//This array will indicate how the indexes changed during the sort.
	bdArray<unsigned int> index_array;
	index_array.Set(m_tags.GetNumberOfElements());
	for(unsigned int i=0; i<index_array.GetNumberOfElements(); i++) index_array[i] = i;

	// Do the bubble sort based on the given criteria.
	int is_change_made = 1;
	while(is_change_made)
	{
		is_change_made = 0;
		for(unsigned int i=0; i<m_tags.GetNumberOfElements()-1; i++)
		{
            bdArray<double> value_array1;
            if(!m_tags[i].ValueDoubles(criteria1_tag_number1, criteria1_tag_number2, value_array1))
            {
                cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria1_tag_number1<<","<<criteria1_tag_number2<<endl;
                return;
            }
            bdArray<double> value_array2;
            //double v1 = value_array[0]; //double v1 = *(m_tags[i].Value(tag_criteria1));
            if(!m_tags[i+1].ValueDoubles(criteria1_tag_number1, criteria1_tag_number2, value_array2))
            {
                cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria1_tag_number1<<","<<criteria1_tag_number2<<endl;
                return;
            }
            //double v2 = value_array[0]; //double v2 = *(m_tags[i+1].Value(tag_criteria1));

			//cout<<" v1="<<v1<<" v2="<<v2<<" ";
            double v1 = 0, v2 = 0;
            for(unsigned int j=0; j<value_array1.GetNumberOfElements() && j<value_array2.GetNumberOfElements(); j++)
            {
                if(value_array1[j]+0.00005<value_array2[j] || value_array1[j]-0.00005>value_array2[j])
                {
                    v1 = value_array1[j]; v2 = value_array2[j];
                    break;
                }
            }

			//Sorting based on the criteria1
			if(v1>v2) 
			{
				//Invert the tag order
				bdImageSliceTag temp_tag;
				temp_tag = m_tags[i+1];
				m_tags[i+1] = m_tags[i];
				m_tags[i] = temp_tag;

				//Invert the index array elements order
				unsigned int temp_i;
				temp_i = index_array[i+1];
				index_array[i+1] = index_array[i];
				index_array[i] = temp_i;

				is_change_made = 1;
			}
			else
			{
				if(v1==v2 && criteria2_tag_number1 && criteria2_tag_number2)
				{
//					v1 = *(m_tags[i].Value(tag_criteria2));
//					v2 = *(m_tags[i+1].Value(tag_criteria2));
                    
                    if(!m_tags[i].ValueDoubles(criteria2_tag_number1, criteria2_tag_number2, value_array1))
                    {
                        cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria2_tag_number1<<","<<criteria2_tag_number2<<endl;
                        return;
                    }
                    //v1 = value_array[0];
                    if(!m_tags[i+1].ValueDoubles(criteria2_tag_number1, criteria2_tag_number2, value_array2))
                    {
                        cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria2_tag_number1<<","<<criteria2_tag_number2<<endl;
                        return;
                    }
                    //v2 = value_array[0];

                    v1 = 0, v2 = 0;
                    for(unsigned int j=0; j<value_array1.GetNumberOfElements() && j<value_array2.GetNumberOfElements(); j++)
                    {
                        if(value_array1[j]+0.00005<value_array2[j] || value_array1[j]-0.00005>value_array2[j])
                        {
                            v1 = value_array1[j]; v2 = value_array2[j];
                            break;
                        }
                    }

					//Sorting based on the criteria2
					if(v1>v2) 
					{
						//Invert the tag order
						bdImageSliceTag temp_tag;
						temp_tag = m_tags[i+1];
						m_tags[i+1] = m_tags[i];
						m_tags[i] = temp_tag;

						//Invert the index array elements order
						unsigned int temp_i;
						temp_i = index_array[i+1];
						index_array[i+1] = index_array[i];
						index_array[i] = temp_i;

						is_change_made = 1;
					}
					else
					{
						if(v1==v2 && criteria3_tag_number1 && criteria3_tag_number2)
						{
//							v1 = *(m_tags[i].Value(tag_criteria3));
//							v2 = *(m_tags[i+1].Value(tag_criteria3));
                            
                            if(!m_tags[i].ValueDoubles(criteria3_tag_number1, criteria3_tag_number2, value_array1))
                            {
                                cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria3_tag_number1<<","<<criteria3_tag_number2<<endl;
                                return;
                            }
                            //v1 = value_array[0];
                            if(!m_tags[i+1].ValueDoubles(criteria3_tag_number1, criteria3_tag_number2, value_array2))
                            {
                                cout<<"bdImageT<T>::SortSlices(): Error getting values from tags:"<<criteria3_tag_number1<<","<<criteria3_tag_number2<<endl;
                                return;
                            }
                            //v2 = value_array[0];
                            
                            v1 = 0, v2 = 0;
                            for(unsigned int j=0; j<value_array1.GetNumberOfElements() && j<value_array2.GetNumberOfElements(); j++)
                            {
                                if(value_array1[j]+0.00005<value_array2[j] || value_array1[j]-0.00005>value_array2[j])
                                {
                                    v1 = value_array1[j]; v2 = value_array2[j];
                                    break;
                                }
                            }


							//Sorting based on the criteria3
							if(v1>v2) 
							{
								//Invert the tag order
								bdImageSliceTag temp_tag;
								temp_tag = m_tags[i+1];
								m_tags[i+1] = m_tags[i];
								m_tags[i] = temp_tag;

								//Invert the index array elements order
								unsigned int temp_i;
								temp_i = index_array[i+1];
								index_array[i+1] = index_array[i];
								index_array[i] = temp_i;

								is_change_made = 1;
							}
						}
					}
				}
			}
		}
	}

	//for(unsigned int i=0; i<index_array.GetNumberOfElements(); i++) cout<<" "<<index_array[i]<<" ";

	//Tags are sorted, the index change order is in index_array, so we now have to sort slices.
	bdImageT<T> temp_img;//temp_img.SetSize(1,1,this->GetNumberOfRows(),this->GetNumberOfColumns());
	temp_img.CopyFrom(*this);
	for(unsigned int i=0; i<index_array.GetNumberOfElements(); i++)
	{
		int old_index = index_array[i];// new index is 'i'.		
		int s_old = old_index % this->GetNumberOfSlices();
		int t_old = old_index / this->GetNumberOfSlices();
		int s_new = i % this->GetNumberOfSlices();
		int t_new = i / this->GetNumberOfSlices();

		for(unsigned int r=0; r<this->GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<this->GetNumberOfColumns(); c++)
			{
				(*this)(t_new,s_new,r,c) = temp_img(t_old,s_old,r,c);
			}
		}
	}
}


template<class T>
int bdImageT<T>::LoadInfoAndTags(const char *file_name, int is_loading_image_properties)
{
	if(this->IsEmpty()) return 0;

	ifstream input_file;
	input_file.open(file_name,ios::binary);
    if(!input_file)
	{
		cout<<"bdImageT<T>::LoadSliceTags() : Unable to open file: "<<file_name<<endl;
		return 0;
    }

	char text_line[100000];//If you change this value, you must also change the amount of loaded characters, see below.
	int line_size = 99999;

	if(!input_file.getline(text_line, line_size)) return 0;
	if(text_line[0]!='%' || text_line[1]!='b' || text_line[2]!='d' || text_line[3]!='I' || text_line[4]!='m' || text_line[5]!='a' || text_line[6]!='g' ||
		text_line[7]!='e' || text_line[8]!=' ' || text_line[9]!='I' || text_line[10]!='n' || text_line[11]!='f' || text_line[12]!='o' ||
		text_line[13]!=':') { input_file.close(); return 0; }

	if(!input_file.getline(text_line, line_size)) return 0;
	if(text_line[0]!='v') return 0;
	if(!input_file.getline(text_line, line_size)) return 0;
	bdString bds;
	bds.Append(text_line);
	bdList<int> list;
	bds.ExtractAllIntNumbersToList(list);
	if(list.IsEmpty())
	{
		cout<<"bdImageT<T>::LoadFromTextFileOfPositions(): No version value found!"<<endl;
		return 0;
	}
	if(list.GetLeftEnd()>1)
	{
		cout<<"bdImageT<T>::LoadFromTextFileOfPositions(): Version higher than supported, version: "<<list.GetLeftEnd()<<endl;
		return 0;
	}

	if(!input_file.getline(text_line, line_size)) return 0;
	if(text_line[0]!='h') return 0;
	this->m_history.Clear();
	if(!input_file.getline(text_line, line_size)) return 0;
	while(text_line[0]=='%')
	{
		this->m_history.Append(&(text_line[1]));
		this->m_history.Append("\n");
		if(!input_file.getline(text_line, line_size)) return 0;
	}
	if(text_line[0]!='t') return 0;

	for( unsigned int index_of_tag = 0; input_file.getline(text_line, line_size) && index_of_tag<m_tags.GetNumberOfElements(); index_of_tag++)
	{
		bds.Clear();
		bds.Append(text_line);
		m_tags[index_of_tag].LoadFromSavedString(bds);
	}

	return 1;
}


template<class T>
void bdImageT<T>::SaveInfoAndTags(const char *file_name)
{
	ofstream file;
	file.open(file_name,ios::binary);

	file<<"%bdImage Info:" <<endl;
	file<<"v"<<endl;
	file<<"1.0"<<endl;
	//file<<"d"<<endl;// Indicates that from this moment on, the origin is recorded
	//file<<(this->GetNumberOfTimeSeries())<<endl;
	//file<<(this->GetNumberOfSlices())<<endl;
	//file<<(this->GetNumberOfRows())<<endl;
	//file<<(this->GetNumberOfColumns())<<endl;
	//file<<"o"<<endl;// Indicates that from this moment on, the origin is recorded
	//file<<(this->GetOrigin_T())<<endl;
	//file<<(this->GetOrigin_S())<<endl;
	//file<<(this->GetOrigin_R())<<endl;
	//file<<(this->GetOrigin_C())<<endl;
	//file<<"s"<<endl;// Indicates that from this moment on, the spacing is recorded
	//file<<(this->GetSpacing_T())<<endl;
	//file<<(this->GetSpacing_S())<<endl;
	//file<<(this->GetSpacing_R())<<endl;
	//file<<(this->GetSpacing_C())<<endl;
	//file<<"O"<<endl;// Indicates that from this moment on, the orientation is recorded
	//file<<(this->GetOrientation_Cx())<<endl;
	//file<<(this->GetOrientation_Cy())<<endl;
	//file<<(this->GetOrientation_Cz())<<endl;
	//file<<(this->GetOrientation_Rx())<<endl;
	//file<<(this->GetOrientation_Ry())<<endl;
	//file<<(this->GetOrientation_Rz())<<endl;
	file<<"h"<<endl;// Indicates that from this moment on, history is recorded.
	if(this->m_history.IsEmpty()) file<<"%"<<endl;
	else file<<"%"<<this->m_history.C_String();//we do not add 'endl' because each entry in history should end with 'endl'.
	file<<"t"<<endl;// Indicates that from this moment on, slice tags are recorded.
	// Record tags into a separate file.
	for(unsigned int t=0; t<this->GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<this->GetNumberOfSlices(); s++)
		{
			bdString bds;
			this->Tag(t,s).SaveToString(bds);
			file<<bds.C_String()<<endl;
		}
	}

	file.close();
}


template<class T>
int bdImageT<T>::SaveToRawVTK16UFile(const char *file_name)
{
//    ofstream output_file(file_name,ios::binary);
//    if(!output_file)
//    {
//        cout<<"Error: bdImage::SaveToRaw16UFile(): Unable to open '"<<file_name<<"'!"<<endl;
//        return 0;
//    }
//    
//    output_file<<"%bdImage Raw File: "<<endl;
//    
//    //----- Write version number -----
//    output_file<<"version"<<endl;//Letter 'v' indicates that the version number is recorded
//    output_file<<"1.0"<<endl;
//    
//    //----- Write dimensions -----
//    output_file<<"dimensions"<<endl;
//    output_file<<(this->GetNumberOfTimeSeries())<<endl;
//    output_file<<(this->GetNumberOfSlices())<<endl;
//    output_file<<(this->GetNumberOfRows())<<endl;
//    output_file<<(this->GetNumberOfColumns())<<endl;
//    
//    //----- Write origin -----
//    output_file<<"origin"<<endl;
//    output_file<<(this->GetOrigin_T())<<endl;
//    output_file<<(this->GetOrigin_S())<<endl;
//    output_file<<(this->GetOrigin_R())<<endl;
//    output_file<<(this->GetOrigin_C())<<endl;
//    
//    //----- Write spacing -----
//    output_file<<"spacing"<<endl;
//    output_file<<(this->GetSpacing_T())<<endl;
//    output_file<<(this->GetSpacing_S())<<endl;
//    output_file<<(this->GetSpacing_R())<<endl;
//    output_file<<(this->GetSpacing_C())<<endl;
//    
//    //----- Write orientation (CxCyCzRxRyRz) -----
//    output_file<<"orientation"<<endl;
//    output_file<<(this->GetOrientation_Cx())<<endl;
//    output_file<<(this->GetOrientation_Cy())<<endl;
//    output_file<<(this->GetOrientation_Cz())<<endl;
//    output_file<<(this->GetOrientation_Rx())<<endl;
//    output_file<<(this->GetOrientation_Ry())<<endl;
//    output_file<<(this->GetOrientation_Rz())<<endl;
//
//    
//    //----- Write all the values -----
//    output_file<<"values"<<endl;
//    
//    //Save all the points consisting of 3 index values and m_vector_dimensions[2]*m_vector_dimensions[1]*m_vector_dimensions[0] float values
//    for(unsigned int t=0; t<this->GetNumberOfTimeSeries(); t++)
//    {
//        for(unsigned int s=0; s<this->GetNumberOfSlices(); s++)
//        {
//            for(unsigned int r=0; r<this->GetNumberOfRows(); r++)
//            {
//                for(unsigned int c=0; c<this->GetNumberOfColumns(); c++)
//                {
//                    output_file.write((char*) &((*this)(t,s,r,c)), sizeof(unsigned short));
//                }
//            }
//        }
//    }
//    output_file.close();
//    return 1;
    

    
    ofstream output_file(file_name,ios::binary);
    if(!output_file)
    {
        cout<<"Error: bdImage::SaveToRaw16UFile(): Unable to open '"<<file_name<<"'!"<<endl;
        return 0;
    }
    output_file<<"# vtk DataFile Version 3.0"<<endl;
    output_file<<"T"<<this->GetNumberOfTimeSeries()<<"t";
    output_file<<"S"<<this->GetNumberOfSlices()<<"s";
    output_file<<"O"<<this->GetOrientation_Cx()<<" "<<this->GetOrientation_Cy()<<" "<<this->GetOrientation_Cz()<<" "<<this->GetOrientation_Rx()<<" "<<this->GetOrientation_Ry()<<" "<<this->GetOrientation_Rz()<<"o";
    output_file<<"B"<<this->GetOrigin_T()<<"b";
    output_file<<"I"<<this->GetSpacing_T()<<"i"<<endl;
    output_file<<"BINARY"<<endl;//or: "ASCII\n"
    output_file<<"DATASET STRUCTURED_POINTS"<<endl;
    output_file<<"DIMENSIONS "<<this->GetNumberOfColumns()<<" "<<this->GetNumberOfRows()<<" "<<(this->GetNumberOfTimeSeries()*this->GetNumberOfSlices())<<endl;
    output_file<<"SPACING "<<this->GetSpacing_C()<<" "<<this->GetSpacing_R()<<" "<<(this->GetSpacing_T()*this->GetSpacing_S())<<endl;
    output_file<<"ORIGIN "<<this->GetOrigin_C()<<" "<<this->GetOrigin_R()<<" "<<this->GetOrigin_S()<<endl;
    output_file<<"POINT_DATA "<<(this->GetNumberOfColumns()*this->GetNumberOfRows()*this->GetNumberOfTimeSeries()*this->GetNumberOfSlices())<<endl;
    output_file<<"SCALARS scalars unsigned_short"<<endl;
    output_file<<"LOOKUP_TABLE default\n";
    
    //Save all the points consisting of 3 index values and m_vector_dimensions[2]*m_vector_dimensions[1]*m_vector_dimensions[0] float values
    for(unsigned int t=0; t<this->GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<this->GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<this->GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<this->GetNumberOfColumns(); c++)
                {
                    output_file.write((char*) &((*this)(t,s,r,c)), sizeof(unsigned short));
                }
            }
        }
    }
    output_file.close();
    return 1;
}



template<class T>
int bdImageT<T>::LoadRawVTK16UFile(const char *file_name)
{
    this->Reset();
    
    char text_line[100000];//If you change this value, you must also change the amount of loaded characters, see below.
    int amount_of_characters_to_load = 99999;
    
    
    ifstream input_file;
    input_file.open(file_name,ios::binary);
    if(!input_file)
    {
        cout<<"bdImageT<T>::LoadRawVTK16UFile(): Unable to open file: "<<file_name<<endl;
        return 0;
    }
    
    //cout<<"1"<<endl;

    // Load the first line, but skip its processing.
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }

    // The second line reads additional fields
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    bdString bds, bds2;
    bds.Append(text_line);
    int t_size = 1, s_size = 1;
    double t_origin = 0, t_spacing = 0;
    int is_loading_time_series = 1;
    if(!bds.IsEmpty())
    {
        if(!bds.ExtractStringBetweenCharacters('T','t',bds2)) is_loading_time_series = 0;
        else bds2.StringToInt(t_size);
        if(!bds.ExtractStringBetweenCharacters('S','s',bds2)) is_loading_time_series = 0;
        else bds2.StringToInt(s_size);
        if(bds.ExtractStringBetweenCharacters('O','o',bds2))
        {
            bdList<double> l;
            bds2.ExtractAllDoubleNumbersToList(l);
            if(l.GetNumberOfElements()==6) this->SetOrientationCxCyCzRxRyRz(l[0],l[1],l[2],l[3],l[4],l[5]);
        }
        if(!bds.ExtractStringBetweenCharacters('B','b',bds2)) is_loading_time_series = 0;
        else bds2.StringToDouble(t_origin);
        if(!bds.ExtractStringBetweenCharacters('I','i',bds2)) is_loading_time_series = 0;
        else bds2.StringToDouble(t_spacing);
    }
    
    // This line must state: "BINARY"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='B' && text_line[1]=='I' && text_line[2]=='N' && text_line[3]=='A' && text_line[4]=='R' && text_line[5]=='Y') ){ input_file.close(); return 0; }

    // This line must state: "DATASET STRUCTURED_POINTS"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='D' && text_line[1]=='A' && text_line[2]=='T' && text_line[3]=='A' && text_line[4]=='S' && text_line[5]=='E' && text_line[6]=='T' && text_line[7]==' '
          && text_line[8]=='S' && text_line[9]=='T' && text_line[10]=='R' && text_line[11]=='U' && text_line[12]=='C' && text_line[13]=='T' && text_line[14]=='U' && text_line[15]=='R' && text_line[16]=='E'
          && text_line[17]=='D' && text_line[18]=='_' && text_line[19]=='P' && text_line[20]=='O' && text_line[21]=='I' && text_line[22]=='N' && text_line[23]=='T' && text_line[24]=='S') ){ input_file.close(); return 0; }

    // This line must state: "DIMENSIONS"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='D' && text_line[1]=='I' && text_line[2]=='M' && text_line[3]=='E' && text_line[4]=='N' && text_line[5]=='S' && text_line[6]=='I' && text_line[7]=='O' && text_line[8]=='N' && text_line[9]=='S') )
        { input_file.close(); return 0; }
    else
    {
        bds.Clear();
        bds.Append(text_line);
        bdList<int> list_of_ints;
        bds.ExtractAllIntNumbersToList(list_of_ints);
        if(is_loading_time_series) { this->SetSize(t_size, s_size, list_of_ints[1],list_of_ints[0]); }
        else { this->SetSize(1,list_of_ints[2],list_of_ints[1],list_of_ints[0]); }
    }
    
    // This line must state: "SPACING"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='S' && text_line[1]=='P' && text_line[2]=='A' && text_line[3]=='C' && text_line[4]=='I' && text_line[5]=='N' && text_line[6]=='G') )
    { input_file.close(); return 0; }
    else
    {
        bds.Clear();
        bds.Append(text_line);
        bdList<double> list_of_doubles;
        bds.ExtractAllDoubleNumbersToList(list_of_doubles);
        if(is_loading_time_series) { this->SetSpacing(t_spacing,list_of_doubles[2],list_of_doubles[1],list_of_doubles[0]); }
        else { this->SetSpacing(1,list_of_doubles[2],list_of_doubles[1],list_of_doubles[0]); }
    }
    
    // This line must state: "ORIGIN"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='O' && text_line[1]=='R' && text_line[2]=='I' && text_line[3]=='G' && text_line[4]=='I' && text_line[5]=='N') )
    { input_file.close(); return 0; }
    else
    {
        bds.Clear();
        bds.Append(text_line);
        bdList<double> list_of_doubles;
        bds.ExtractAllDoubleNumbersToList(list_of_doubles);
        if(is_loading_time_series) { this->SetOrigin(t_origin,list_of_doubles[2],list_of_doubles[1],list_of_doubles[0]); }
        else { this->SetOrigin(0,list_of_doubles[2],list_of_doubles[1],list_of_doubles[0]); }
    }

    // This line must state: "POINT_DATA"
    int total_number_of_data_points;
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='P' && text_line[1]=='O' && text_line[2]=='I' && text_line[3]=='N' && text_line[4]=='T' && text_line[5]=='_' && text_line[6]=='D' && text_line[7]=='A' && text_line[8]=='T' && text_line[9]=='A') )
    { input_file.close(); return 0; }
    else
    {
        bds.Clear();
        bds.Append(text_line);
        bdList<int> list_of_ints;
        bds.ExtractAllIntNumbersToList(list_of_ints);
        total_number_of_data_points = list_of_ints[0];
        if(total_number_of_data_points!=this->GetNumberOfColumns()*this->GetNumberOfRows()*this->GetNumberOfTimeSeries()*this->GetNumberOfSlices()){ input_file.close(); return 0; }
    }
    
    //cout<<"3"<<endl;
    
    // This line must state: "SCALARS"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='S' && text_line[1]=='C' && text_line[2]=='A' && text_line[3]=='L' && text_line[4]=='A' && text_line[5]=='R' && text_line[6]=='S') )
        { input_file.close(); return 0; }

    // This line must state: "LOOKUP_TABLE"
    if(!input_file.getline(text_line, amount_of_characters_to_load)) { input_file.close(); return 0; }
    if( !(text_line[0]=='L' && text_line[1]=='O' && text_line[2]=='O' && text_line[3]=='K' && text_line[4]=='U' && text_line[5]=='P' && text_line[6]=='_'
          && text_line[7]=='T' && text_line[8]=='A' && text_line[9]=='B' && text_line[10]=='L' && text_line[11]=='E') )
    { input_file.close(); return 0; }
    

    //Reading the data
    for(unsigned int t=0; t<this->GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<this->GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<this->GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<this->GetNumberOfColumns(); c++)
                {
                    unsigned short value;
                    input_file.read((char*) &value, sizeof(unsigned short));
                    (*this)(t,s,r,c) = value;
                }
            }
        }
    }
   
    input_file.close();
    
    
    return 1;
}



template class BD_IMAGE_API bdImageT<unsigned char>;
template class BD_IMAGE_API bdImageT<unsigned short>;







	//void SaveRawFile(char *file_name_root);
	//void LoadRawFile(char *file_name_root);


	////Sorts slices by values in Image Number Tag.
	//void SortSlicesBy_ImageNumber_Tag();
	//void SortSlicesBy_ImageNumber_Tag(Image3D<T> &output_image);
	////Sorts slices by values in Slice Location Tag.
	//void SortSlicesBy_SliceLocation_Tag();
	//void SortSlicesBy_SliceLocation_Tag(Image3D<T> &output_image);
	////Sorts the slices by giving priority to 'created time' tag and uses 'image number' tag as the second criterion.
	//void SortSlicesBy_CreationTime_Tag_And_ImageNumber_Tag();
	//void SortSlicesBy_CreationTime_Tag_And_ImageNumber_Tag(Image3D<T> &output_image);
	////Finds the lowest 'S' value of Image Patient Position for all slices and sets the Origin of the whole data set by origin (S,R,C) of the found slice.
	//void SetOriginFrom_ImagePositionPatient_Tag();
	////Finds the lowest value of Slice Location and sets the only the 'S' Origin value of the whole data set to the found value.
	//void Set_S_OriginFrom_SliceLocation_Tag();
	//int Is_SliceLocationTag_TheSameForEverySlice();
	//int Set_S_SpacingAsAverageTriggerTimeInterval();


	////Inverts the order of the data stored in memory (for loading PNGs).
	//void _InvertMemoryOrder();

	////For a smaller 2D image finds its position in the 3D image that contains the 2D image. If position found returns 1, if not returns 0.
	//int FindPositionOf2DImage(Image<T> &input_image, Point &output_position_of_2D_image);
	//void InitializeWith(T &initial_value);
	//void Invert();
	//int Invert3D();//Invert according to the min_max range of the whole 3D image. 
	//void Normalize();
	//int Normalize3D();//Normalize according to the min_max range of the whole 3D image. 
	//int RescaleToRangeOfValues(double lower_value, double upper_value, Image3D<T> &output_image);
	//void Threshold(unsigned int threshold);
	//void Threshold(int threshold, Image3D<T> &output_image);
	////Interpolates slices of this image, number_of_surrounding_slices determines how many slices around the current one are taken into account,
	//// distance_coefficient is the weight difference of each one (the weight decreases as we go away from the current slice).
	//// Default: distance_coefficient = 0, means that all slices will be treated equaly, 
	//// otherwise: new = a(0) + distance_coefficient/n*a(n) + distance_coefficient/n*a(-n), n=[1,...,number_of_surrounding_slices];
	//void InterpolateSlices(Image3D<T> &output_image, int number_of_surrounding_slices, double distance_coefficient = 0);
	////Use this to arrange slice order in images where the slices were not well ordered (1 or 2 parts of the image should be at the different
	//// place in the image). 'slice1' indicates the starting slice and the part copied after it are the slices with higher number up to the 
	//// border defined either by NumberOfSlices() of the image or values 'slice2' or 'slice3'. Values 'slice1' and 'slice2' must be different
	//// and in [0,NumberOfSlices()-1] range. 'slice3' must be different than 'slice1' or 'slice2' and is -1 by default, when it is not used.
	//void ReArrangeSliceOrder(Image3D<T> &output, int slice1, int slice2, int slice3=-1);

	////Performs masking of this image with mask_image and stores result to output_image 
	//int Mask(Image3D<T> &mask_image, Image3D<T> &output_image);
	//int MinNonZeroValue();//Finds the minimum non-zero value in the image. If image empty or all are zeros, returns 0.
	//int MaxNonZeroValue();//Finds the maximum non-zero value in the image. If image empty or all are zeros, returns 0.
	//void Add(Image3D<unsigned short> &input_image1, Image3D<unsigned short> &input_image2);//Add values from input image and store to output image
	////Shifts the image by given shift values for slices, rows and columns. Positions that are not copied are filled with fill_in_value
	//int Shift(Image3D<T> &output_image, int slice_shift, int row_shift, int column_shift, T fill_in_value);
	//int Flip_S_Axis(Image3D<T> &output_image);
	//int Flip_R_Axis(Image3D<T> &output_image);
	//int Flip_C_Axis(Image3D<T> &output_image);
	//int RotateAround_C_Axis_Minus90Degrees(Image3D<T> &output);//'S' becomes '-R', 'R' becoomes 'S', 'C' stays 'C'.
	//int RotateAround_C_Axis_Plus90Degrees(Image3D<T> &output);//'S' becomes 'R', 'R' becoomes '-S', 'C' stays 'C'.
	//int RotateAround_R_Axis_Plus90Degrees(Image3D<T> &output);//'S' becomes '-R', 'R' stays 'R', 'C' becomes 'S'.
	//int RotateAround_R_Axis_Minus90Degrees(Image3D<T> &output);//'S' becomes 'R', 'R' stays 'R', 'C' becomes '-S'.
	//int RotateAround_S_Axis_Plus90Degrees(Image3D<T> &output);//'S' stays 'S', 'R' becomes 'C', 'C' becomes '-R'.
	//int RotateAround_S_Axis_Minus90Degrees(Image3D<T> &output);//'S' stays 'S', 'R' becomes '-C', 'C' becomes 'R'.
	//void Rotate_SRC_To_RSC();
	//void Rotate_SRC_To_CRS();
	//void Rotate_SRC_To_SCR();
	//int CutOutPartOf3DImage(Image3D<T> &output_image, int s_start_included, int s_end_included, int r_start_included, int r_end_included, int c_start_included, int c_end_included);
	//int binary_NumberOfPixelsIn_26_Neighborhood(int s, int r, int c);
	//int binary_NumberOfPixelsIn_27_Neighborhood(int s, int r, int c);
	//void binary_GetObjectBounds(int *s_start, int *s_end, int *r_start, int *r_end, int *c_start, int *c_end);

	////Fills in this images with values from input images, where first image has priority (if both 
	//// input images are non-zero in the same voxel the value copied is from the first input image).
	//int MakeCombinationOfImages(Image3D<T> &input_image_with_priority, Image3D<T> &input_image, int second_input_image_label = -1);

	//Image<T>& operator [](int s);
	//T& operator ()(int s, int r, int c);
	//Image3D<T>& operator =(Image3D<T> &m);

	////int LoadFromSlices(const char *file_name_root, int number_of_slices, char *extension);
	////int LoadFromSlices(const char *file_name_root, int start_slice_index_included, int end_slice_index_included, char *extension);
	//int LoadFrom_PNG_Slices(const char *file_name_root, int _number_of_slices);
	//int LoadFrom_PNG_Slices(const char *file_name_root, int slice_start_index, int slice_end_index);
	//int LoadVTKFileStructured(char *file_name);//OVO JE ZA HASTE!!!!!!!
	//int Load_VTK_FileStructuredPoints(char *file_name);
	//int Load_VTK_FileStructuredPoints(const char *file_name);
	//int Load_VTK_v3_FileStructuredPoints(char *file_name);
	//void SaveToSlices(char *file_name_root, char *extension);
	//void SaveToSlices(char *file_name_root, int sequence_number, char *extension);
	//void SaveToVTKFile(const char *file_name);
	

	////In the skeletonized image with paths that are not connected, ceates connected paths
	//int binary_ConnectPathsInSkeletonizedImage(Image3D<T> &output);

	////In the current binary image finds the largest 26-connected region of values larger than 'threshold' and stores it into 'output_image'. LEGACY FUNCTION.
	//void binary_ExtractLargest_26_ConnectedRegion(Image3D<T> &output_image, int threshold=0);

	////In the current binary image finds the 26-connected region of values larger than 'threshold' from seed point and stores it into 'output_image'.
	//int binary_Extract_26_ConnectedRegionFromSeedPoint(Image3D<T> &output_image, Point &seed_point_indexes, int threshold=0);

	////In the current binary image finds the largest 6-connected region and stores it into 'output_image'. LEGACY FUNCTION.
	//void binary_ExtractLargest_6_ConnectedRegion(Image3D<T> &output_image);
