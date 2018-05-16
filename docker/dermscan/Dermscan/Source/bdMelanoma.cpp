/**************************************************************************
 Methods for ABCDE features of melanoma.

 Author: Danilo Babin
 File name: "bdMelanoma.cpp"
***************************************************************************/



#include "bdMelanoma.h"




int bdMelanoma::IsMasked(bdImage &mask, int r, int c)
{
    if( (mask(0,0,r,c) && mask(1,0,r,c)) || (mask(0,0,r,c) && mask(2,0,r,c)) || (mask(1,0,r,c) && mask(2,0,r,c)) ) return 1;
    return 0;
}



int bdMelanoma::CreateMaskImage(bdImage &segmented_image, bdImage &mask_image)
{
    if(segmented_image.IsEmpty()) return 0;
    
    mask_image.SetSize(1,1,segmented_image.GetNumberOfRows(),segmented_image.GetNumberOfColumns());
    mask_image.SetVisualizationPropertiesToMatchInput(segmented_image);
    mask_image.FillInWith(0);
    
    for(unsigned int r=0; r<segmented_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_image.GetNumberOfColumns(); c++)
        {
            if(this->IsMasked(segmented_image,r,c))
            {
                mask_image(r,c) = 255;
            }
        }
    }
    
    return 1;
    
}


int bdMelanoma::ProjectPositionToLine(double *line_pos1_src, double *line_pos2_src, double *position_to_project_src, double *output_projected_position_src, double &n)
{
    if(line_pos1_src[0]==line_pos2_src[0] && line_pos1_src[1]==line_pos2_src[1] && line_pos1_src[2]==line_pos2_src[2]) return 0;
    
    //'n' is the coefficient stating where the position falls on the vector defined by 'vector_position1' and 'vector_position2'
    
    //The derived formula is:
    // n(x1^2 - 2x1x2 + x2^2) + x1x2- x1^2 + x1x3 - x2x3 + n(y1^2 - 2y1y2 + y2^2) + y1y2- y1^2 + y1y3 - y2y3 + n(z1^2 - 2z1z2 + z2^2) + z1z2- z1^2 + z1z3 - z2z3 = 0
    
    //'p1' will be the part NOT multiplied by 'n': x1x2- x1^2 + x1x3 - x2x3 + y1y2- y1^2 + y1y3 - y2y3 + z1z2- z1^2 + z1z3 - z2z3.
    double p1 = line_pos1_src[2]*line_pos2_src[2] - line_pos1_src[2]*line_pos1_src[2] + line_pos1_src[2]*position_to_project_src[2] - line_pos2_src[2]*position_to_project_src[2] +
    line_pos1_src[1]*line_pos2_src[1] - line_pos1_src[1]*line_pos1_src[1] + line_pos1_src[1]*position_to_project_src[1] - line_pos2_src[1]*position_to_project_src[1] +
    line_pos1_src[0]*line_pos2_src[0] - line_pos1_src[0]*line_pos1_src[0] + line_pos1_src[0]*position_to_project_src[0] - line_pos2_src[0]*position_to_project_src[0];
    
    //'p2' will be the part multiplied by 'n': (x1^2 - 2x1x2 + x2^2) + (y1^2 - 2y1y2 + y2^2) + (z1^2 - 2z1z2 + z2^2).
    double p2 = line_pos1_src[2]*line_pos1_src[2] - 2*line_pos1_src[2]*line_pos2_src[2] + line_pos2_src[2]*line_pos2_src[2] +
    line_pos1_src[1]*line_pos1_src[1] - 2*line_pos1_src[1]*line_pos2_src[1] + line_pos2_src[1]*line_pos2_src[1] +
    line_pos1_src[0]*line_pos1_src[0] - 2*line_pos1_src[0]*line_pos2_src[0] + line_pos2_src[0]*line_pos2_src[0];
    
    //Now the equation becomes: n*p2 + p1 = 0  -->  n = -p1/p2.
    
    //Obviously, 'n' can't be calculated if p2=0, (there is a zero vector or the position is such that no singe projection can be made), so return 0.
    if(p2==0) return 0;
    
    n = -p1/p2;
    
    //If 'n' is in range [0,1] the projection position is ON the vector, if greater than 1 it is outside vector in the direction of the vector,
    // if smaller than 0 it is outside the vector in the opposite direction than the vector direction.
    
    if(n==0)
    {
        output_projected_position_src[0] = line_pos1_src[0];
        output_projected_position_src[1] = line_pos1_src[1];
        output_projected_position_src[2] = line_pos1_src[2];
    }
    else
    {
        double dc = line_pos2_src[2] - line_pos1_src[2];
        double dr = line_pos2_src[1] - line_pos1_src[1];
        double ds = line_pos2_src[0] - line_pos1_src[0];
        
        output_projected_position_src[2] = line_pos1_src[2] + n*dc;
        output_projected_position_src[1] = line_pos1_src[1] + n*dr;
        output_projected_position_src[0] = line_pos1_src[0] + n*ds;
    }
    
    return 1;
}





int bdMelanoma::EqualizeBrightnessOfChannels(bdImage &input, bdImage &output)
{
//    // Rescale the slice values to fall within 0,255.
//    unsigned int lower_value = 0, upper_value = 255;
//    bdGIP;
//    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    output.SetSizeAndPropertiesAs(input);
    
    
    bdArray< bdArray<int> > array_of_histograms;
    bdArray< int > array_of_indexes_of_maximum_values;

    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    array_of_indexes_of_maximum_values.Set(input.GetNumberOfTimeSeries());
    
    unsigned int lower_value = 0, upper_value = 255;//65535;
    input.GetVoxelMinimumAndMaximumValue(&lower_value, &upper_value);
    
    // For each slice calculate the histogram and the maxima of the histogram.
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0, upper_value, t, 0);
        
        // Go through the histogram and find the max value and its index.
        unsigned int max = 0;
        array_of_indexes_of_maximum_values[t] = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements()-1; i++)
        {
            if(array_of_histograms[t][i]>max)
            {
                max = array_of_histograms[t][i];
                array_of_indexes_of_maximum_values[t] = i;
            }
        }
    }
    
    //cout<<" 3 ";
    
    
    // Go through the maxima array and find the index of each peak (maximum alue in the histogram) and record these indexes together wth the max index value.
    
 
    
    // Find the max value and record its index...
    unsigned int max_recorded_index = 0;
    for(unsigned int t=0; t<array_of_indexes_of_maximum_values.GetNumberOfElements(); t++)
    {
        if(array_of_indexes_of_maximum_values[t]>max_recorded_index)
        {
            max_recorded_index = array_of_indexes_of_maximum_values[t];
        }
        
        cout<<" ["<<t<<"]="<<array_of_indexes_of_maximum_values[t]<<"  ";
    }
    
    cout<<" max_recorded_index="<<max_recorded_index<<"  ";
    
    
    //output.CopyFrom(input);
    
    // Go through each of the image bands and to each pixel gray value add the value of the difference between the index of the maximum found peak index and the index of maximum value for that channel.
    
    for(unsigned int t=0; t<array_of_indexes_of_maximum_values.GetNumberOfElements(); t++)
    {
        unsigned short diff = max_recorded_index - array_of_indexes_of_maximum_values[t];
        
        cout<<" ["<<t<<"]_diff="<<diff<<"  ";
        
        double factor = ((double)max_recorded_index) / ((double)(array_of_indexes_of_maximum_values[t]+1));
        
        cout<<" ["<<t<<"]_factor="<<factor<<"  ";
        
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
//                int v = input(t,0,r,c) + diff;
//                if(v>upper_value) v = upper_value;
//                output(t,0,r,c) = v;
                
                unsigned short v = ((double)input(t,0,r,c))  * factor;
                if(v>upper_value) v = upper_value;
                output(t,0,r,c) = v;

            }
        }
    }
    
    return 1;
}


int bdMelanoma::DilateBinaryImageWithCircleSE(bdImage &input, bdImage &output, unsigned int squared_radius)
{
    if(input.IsEmpty()) return 0;
    
    unsigned int s = 0;
    
    output.CopyFrom(input);
    
    bdGeometry g;
    g.SetDimensions(1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    
    unsigned int t = 0;
    //for(unsigned int t=0; t<input_image.GetNumberOfTimeSeries(); t++)
    //{
        //for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
        //{
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    if(input(t,s,r,c)!=0)
                    {
                        //check if the pixel is an edge pixel
                        int rn,cn;
                        int is_edge_pixel = 0;
                        for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn);)
                        {
                            if(input(rn,cn)==0)
                            {
                                is_edge_pixel = 1;
                                break;
                            }
                        }
                        
                        // if the pixel is an edge pixel, perform dilation
                        if(is_edge_pixel)
                        {
                            for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn);)
                            {
                                output(rn,cn) = input(r,c);
                            }
                        }
                    }
                }
            }
        //}
    //}
    return 1;
}



int bdMelanoma::Histogram(bdImage &input, bdArray<int> &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t, unsigned int s)
{
    if(input.IsEmpty()) return 0;
    
    //unsigned int full_range_min, full_range_max;
    //input.GetVoxelValueFullRange(&full_range_min, &full_range_max);

    output_histogram.Set(range_max+1);
    output_histogram.FillInWith(0);
    //for(unsigned int t=0; t<input_image.GetNumberOfTimeSeries(); t++)
    //{
        //for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
        //{
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    if(input(t,s,r,c)>=range_min && input(t,s,r,c)<=range_max)
                    {
                        output_histogram[(input(t,s,r,c))] += 1;
                    }
                }
            }
        //}
    //}
    
    return 1;
}


int bdMelanoma::Histogram(bdImage &input, bdCurveXY &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t, unsigned int s)
{
    if(range_max<= range_min) return 0;
    
    bdArray<int> histogram;
    
    if(!this->Histogram(input, histogram, range_min, range_max, t, s)) return 0;
    
    output_histogram.Set(range_max-range_min+1);
    
    int i=0;
    int a=0;
    for(i=output_histogram.GetNumberOfElements()-1, a=histogram.GetNumberOfElements()-1; i>0 && a>0; i--, a--)
    {
        output_histogram.SetValues(i, i, histogram[a]);
    }
    
    return 1;
}



int bdMelanoma::Histogram(bdImage &input, bdImage &mask, bdArray<int> &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask)
{
    if(input.IsEmpty()) return 0;
    
    //unsigned int full_range_min, full_range_max;
    //input.GetVoxelValueFullRange(&full_range_min, &full_range_max);
    
    output_histogram.Set(range_max+1);
    output_histogram.FillInWith(0);
    //for(unsigned int t=0; t<input_image.GetNumberOfTimeSeries(); t++)
    //{
    //for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
    //{
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            if(input(t_input,s_input,r,c)>=range_min && input(t_input,s_input,r,c)<=range_max)
            {
                if(mask(t_mask,s_mask,r,c)!=0)
                {
                    output_histogram[(input(t_input,s_input,r,c))] += 1;
                }
            }
        }
    }
    //}
    //}
    
    return 1;
}


int bdMelanoma::Histogram(bdImage &input, bdImage &mask, bdCurveXY &output_histogram, unsigned int range_min, unsigned int range_max, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask)
{
    if(range_max<=range_min) return 0;
    
    bdArray<int> histogram;
    
    if(!this->Histogram(input, mask, histogram, range_min, range_max, t_input, s_input, t_mask, s_mask)) return 0;
    
    output_histogram.Set(range_max-range_min+1);
    
    int i=0;
    int a=0;
    for(i=output_histogram.GetNumberOfElements()-1, a=histogram.GetNumberOfElements()-1; i>0 && a>0; i--, a--)
    {
        output_histogram.SetValues(i, i, histogram[a]);
    }
    
    return 1;
}



int bdMelanoma::HistogramMaxima(bdArray<int> &input_histogram, bdArray<int> &output_maxima)
{
    if(input_histogram.IsEmpty()) return 0;
    
    output_maxima.Set(input_histogram.GetNumberOfElements());
    output_maxima.FillInWith(0);
    
    
    for(unsigned int i=1; i<input_histogram.GetNumberOfElements()-1; i++)
    {
        // search if the current value is the maxima to is left
        int is_maximum_to_left = -1;
        int counter_to_left = 0;
        int is_broken_to_left = 0;
        for(unsigned int j=1; j<=i; j++)
        {
            if(is_maximum_to_left<0)// have we determined if this is a maximum?
            {
                if(input_histogram[i]<input_histogram[i-j])
                {
                    is_maximum_to_left = 0;
                }
                if(input_histogram[i]>input_histogram[i-j])
                {
                    is_maximum_to_left = 1;
                }
            }
            counter_to_left++;
            if(input_histogram[i]<input_histogram[i-j])
            {
                is_broken_to_left = 1;
                break;
            }
        }
        
        if(!is_broken_to_left) counter_to_left = 255;
        
        // If the current value is a maximum after searching to its left...
        if(is_maximum_to_left)
        {
            //... search to its right
            int is_maximum_to_right = 1;
            int counter_to_right = 0;
            int is_broken_to_right = 0;
            for(unsigned int k=1; k<input_histogram.GetNumberOfElements()-1-i; k++)// we deliberately do not check till the end of the range, because there is often cutoff peak.
            {
                if(is_maximum_to_right<0)// have we determined if this is a maximum? (if is_maximum_to_right<0 we still haven't)
                {
                    if(input_histogram[i]<input_histogram[i+k])
                    {
                        is_maximum_to_right = 0;
                    }
                    if(input_histogram[i]>input_histogram[i+k])
                    {
                        is_maximum_to_right = 1;
                    }
                }
                counter_to_right++;
                if(input_histogram[i]<input_histogram[i+k])
                {
                    is_broken_to_right = 1;
                    break;
                }
            }
            
            if(!is_broken_to_right) counter_to_right = 255;
            
            // If the current value is a maximum, assign it the smaller counter value
            if(is_maximum_to_left && is_maximum_to_right)
            {
                if(counter_to_right<counter_to_left)
                {
                    output_maxima[i] = counter_to_right;
                }
                else
                {
                    output_maxima[i] = counter_to_left;
                }
            }
            else
            {
                output_maxima[i] = 0;
            }
        }
    }
    
    return 1;
}


int bdMelanoma::EnvelopeOfHistogram(bdArray<int> &input_histogram, unsigned int start_index, unsigned int end_index, unsigned int window_size, bdArray<int> &output_envelope)
{
    if(input_histogram.IsEmpty()) return 0;
    
    output_envelope.Set(input_histogram.GetNumberOfElements());
    output_envelope.CopyArrayIntoThisArrayStartingFromIndex(input_histogram);
    
    for(unsigned int i=0; i<input_histogram.GetNumberOfElements(); i++)
    {
        if(i>=window_size)// this is to say that: i-window_size>=0
        {
            for(unsigned int j=1; j<window_size; j++)
            {
                int v = (input_histogram[i] * (window_size - j)) / j;
                if(v>output_envelope[i-j]) output_envelope[i-j] = v;
            }
        }
        if(i+window_size<input_histogram.GetNumberOfElements())
        {
            for(unsigned int j=1; j<window_size; j++)
            {
                int v = (input_histogram[i+window_size] * (i + j)) / (i+window_size);
                if(v>output_envelope[i+j]) output_envelope[i+j] = v;
            }
        }
    }
    
    return 1;
}


int bdMelanoma::ConnectedComponents_2D_WithSeedIndexes(bdImage &input_image, int threshold, unsigned int t, bdDiscreteCoordinates3D &seed_indexes, int &output_number_of_components, int &output_region_size)
{
    if(input_image.IsEmpty()) return 0;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(input_image, threshold, thresholded_image);
    
    
    //Initialize output image.
    bdImage output_image;
    output_image.SetSize(thresholded_image.GetNumberOfTimeSeries(),thresholded_image.GetNumberOfSlices(),thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
    output_image.FillInWith(0);
    
    unsigned short n_of_regions = 1;
    for(unsigned int s=0; s<thresholded_image.GetNumberOfSlices(); s++)
    {
        for(unsigned int r=0; r<thresholded_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<thresholded_image.GetNumberOfColumns(); c++)
            {
                //Find the 3D region
                if(thresholded_image(t,s,r,c)>threshold && output_image(t,s,r,c)==0)
                {
                    bdRegion3D reg;
                    reg.CreateRegion_26_FromSeed(thresholded_image,c,r,s,t);
                    
                    bdRegion3DIterator itr;
                    for(itr.SetBegin(&reg); itr.IsValid(); itr.MoveToNext())
                    {
                        output_image(itr.GetIndexSlice(),itr.GetIndexRow(),itr.GetIndexColumn()) = n_of_regions;
                    }
                    n_of_regions++;
                    //cout<<" n:"<<n_of_regions<<" ";
                }
            }
        }
    }
    
    output_region_size = 0;
    
    bdGeometry g;
    g.SetDimensions(thresholded_image.GetNumberOfSlices(),thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
    
    int rn, cn;
    int is_region_found = 0;
    for(int squared_radius=0; squared_radius<g.GetNumberOfSphereElements() && !is_region_found; squared_radius++)
    {
        for(g.ForCoordinates_Circle(seed_indexes.R(), seed_indexes.C(), squared_radius); g.Get_Circle(squared_radius+1, rn, cn); )
        {
            if(output_image(t,0,rn,cn)!=0)
            {
                bdRegion3D reg;
                reg.CreateRegion_26_FromSeed(thresholded_image,cn,rn,0,t);
                output_region_size = reg.GetNumberOfElements();
                is_region_found = 1;
                break;
            }
        }
    }
    
    
    
    output_number_of_components = n_of_regions-1;
    
    return 1;
}



int bdMelanoma::LevelComponentsOfRGBImage(bdImage &input, unsigned int level_value, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    
    output.CopyFrom(input);
    
    unsigned int voxel_value_min, voxel_value_max;
    input.GetVoxelValueFullRange(&voxel_value_min, &voxel_value_max);
    if(level_value>voxel_value_max) level_value = voxel_value_max;
    if(level_value<voxel_value_min) level_value = voxel_value_min;
    
    unsigned int s = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
                if(input(t,s,r,c) > level_value) output(t,s,r,c) = level_value;
            }
        }
    }
    
    return 1;
}


int bdMelanoma::MeanOperatorPerRGBComponent(bdImage &input, bdImage &mask, unsigned int squared_circular_SE_radius, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    
    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    output.FillInWith(0);

    bdGeometry g;
    g.SetDimensions(1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    
    unsigned int s = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
                if(mask(r,c)!=0)
                {
                    unsigned int sum = 0;
                    unsigned int number_of_voxels = 0;
                    int rn,cn;
                    for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_circular_SE_radius, rn,cn);)
                    {
                        sum += input(t,s,rn,cn);
                        number_of_voxels++;
                    }
                    
                    output(t,s,r,c) = (unsigned short) (sum/number_of_voxels);
                    
                }
            }
        }
    }
    
    
    return 1;
}



int bdMelanoma::NormalizeComponentsOfRGBImage(bdImage &input, unsigned int lower_value, unsigned int upper_value, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    
    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    output.FillInWith(0);
    
    unsigned int voxel_value_min, voxel_value_max;
    input.GetVoxelValueFullRange(&voxel_value_min, &voxel_value_max);
    if(upper_value>voxel_value_max) upper_value = voxel_value_max;
    if(lower_value<voxel_value_min) lower_value = voxel_value_min;
    
    unsigned int s = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        //for(unsigned int s=0; s<input.GetNumberOfSlices(); s++)
        //{
            // Find min and max value for this slice
            unsigned int min = voxel_value_max, max = 0;
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    if(input(t,s,r,c)<min) min = input(t,s,r,c);
                    if(input(t,s,r,c)>max) max = input(t,s,r,c);
                }
            }
            
            unsigned int range = max-min;
            
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    //cout<<"["<<t<<","<<s<<","<<r<<","<<c<<"]="<<input_image(t,s,r,c)<<"  ";
                    output(t,s,r,c) = input(t,s,r,c) - min;
                    //cout<<"("<<t<<","<<s<<","<<r<<","<<c<<")="<<output_image(t,s,r,c)<<"  ";
                    output(t,s,r,c) = (unsigned short)( lower_value + ((((double)(output(t,s,r,c))) * (upper_value-lower_value))) / ((double)range) );
                }
            }
            
        //}
    }

    
    
    return 1;
}



int bdMelanoma::NormalizeWholeRGBRangeImage(bdImage &input, unsigned int lower_value, unsigned int upper_value, bdImage &output)
{
    if(input.IsEmpty()) return 0;
 
    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    output.FillInWith(0);
    
    unsigned int voxel_value_min, voxel_value_max;
    input.GetVoxelValueFullRange(&voxel_value_min, &voxel_value_max);
    if(upper_value>voxel_value_max) upper_value = voxel_value_max;
    if(lower_value<voxel_value_min) lower_value = voxel_value_min;
    
    unsigned int s = 0;
    
    
    unsigned int min = voxel_value_max, max = 0;
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
                if(input(t,s,r,c)<min) min = input(t,s,r,c);
                if(input(t,s,r,c)>max) max = input(t,s,r,c);
            }
        }
    }
    
    unsigned int range = max-min;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
                //cout<<"["<<t<<","<<s<<","<<r<<","<<c<<"]="<<input_image(t,s,r,c)<<"  ";
                output(t,s,r,c) = input(t,s,r,c) - min;
                //cout<<"("<<t<<","<<s<<","<<r<<","<<c<<")="<<output_image(t,s,r,c)<<"  ";
                output(t,s,r,c) = (unsigned short)( lower_value + ((((double)(output(t,s,r,c))) * (upper_value-lower_value))) / ((double)range) );
            }
        }
    }
    
    
    return 1;
}



double bdMelanoma::SumOfVoxelValuesForAllRGBComponents(bdImage &input, bdImage &mask)
{
    if(input.IsEmpty()) return 0;
    
    double sum =0;
    
    unsigned int s = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
            {
                if(mask(r,c)!=0)
                {
                    sum += input(t,s,r,c);
                }
            }
        }
    }
    
    return sum;
}




//int bdMelanoma::SegmentMelanoma_Thresholding(bdImage &input_image, bdImage &output)
//{
//    cout<<" 1 ";
//    
//	if(input_image.IsEmpty()) return 0;
//    
//    
//    bdGIP gip;
//    bdImage downsampled;
//    gip.DownsampleByCoefficient(input_image, 1, 10, 10, downsampled);
//    
//    
//    bdImage input;
//    gip.Negative(downsampled, input);
//    
//    //gip.
//    
//	//if(seed_points.GetNumberOfElements()!=2) return 0;
//
//	//bdImage img;
//	//bdGIP gip;
//	//gip.Extract3DImage(input,0,img);
//
//	if(this->IsAbortRequested()) return 0;
//	this->SetProgressCounterRelativeValue(5);
//
//
//    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
//	output.SetVisualizationPropertiesToMatchInput(input);
//	output.FillInWith(0);
//    
//    
//    bdArray< bdArray<int> > array_of_histograms;
//    array_of_histograms.Set(input.GetNumberOfTimeSeries());
//
//    // For each slice calculate the histogram
//    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
//    {
//        this->Histogram(input, array_of_histograms[t], t, 0);
//    }
//    
//        cout<<" 2 ";
//    
//    // For each histogram calculate variance and find the index of the max variance (actually MINIMUM variance because we inverted the images)
//    double max_variance = -1;//initial non-valid value.
//    unsigned int max_variance_index = 0;
//    
//    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
//    {
//        // calculate mean
//        double mean = 0, sum = 0, n=0;
//        for(unsigned int i=0; i<array_of_histograms[t].GetNumberOfElements(); i++)
//        {
//            sum += array_of_histograms[t][i] * i;
//            n += i;
//        }
//        mean = sum / n;
//        
//        // calculate variance
//        double variance = 0;
//        for(unsigned int i=0; i<array_of_histograms[t].GetNumberOfElements(); i++)
//        {
//            variance += (array_of_histograms[t][i] - mean) * (array_of_histograms[t][i] - mean) * i;
//        }
//        variance = variance / n;
//        
//        cout<<"variance["<<t<<"]="<<variance<<"  ";
//
//
//        // WE USE INVERTED IMAGES, SO WE SHOULD SEARCH FOR THE SMALLEST VARIANCE!
//        if(variance<max_variance || max_variance<0)
//        {
//            max_variance_index = t;
//            max_variance = variance;
//        }
//    }
//    
//        cout<<" 3 ";
//    
//    // We continue to work on the slice that has the highest variance, so extract the slice
//    bdImage slice;
//    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
//    slice.SetVisualizationPropertiesToMatchInput(input);
//    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
//        {
//            slice(0,0,r,c) = input(max_variance_index,0,r,c);
//        }
//    }
//    
//        cout<<" 4 ";
//    
//    // Rescale the slice values to fall within 0,255.
//    unsigned int lower_threshold = 0, upper_threshold = 255;
//    
//    bdImage rescaled_slice;
//    gip.RescaleWholeRange(slice, lower_threshold, upper_threshold, rescaled_slice);
//    
//        cout<<" 5 ";
//    
//    // Calculate the number of connected components for each threshold value in order to determine the optimal threshold
//    bdBIP bip;
//    bdArray<int> histogram_of_CCs;
//    histogram_of_CCs.Set(upper_threshold+1);
//    histogram_of_CCs.FillInWith(0);
//    bdImage connected_components;
//    //for(unsigned int threshold = lower_threshold+1; threshold<upper_threshold; threshold++)
//    for(unsigned int threshold = lower_threshold+1; threshold<upper_threshold; threshold++)
//    {
//        int number_of_CCs = 0;
//        int size_of_CC = 0;
//        
//        //bip.ConnectedComponents_3D(rescaled_slice, threshold, connected_components, 0, &number_of_CCs);
//        bdDiscreteCoordinates3D seed_indexes;
//        seed_indexes.S() = 0;
//        seed_indexes.R() = rescaled_slice.GetNumberOfRows()/2;
//        seed_indexes.C() = rescaled_slice.GetNumberOfColumns()/2;
//        
//        this->ConnectedComponents_2D_WithSeedIndexes(rescaled_slice, threshold, 0, seed_indexes, number_of_CCs, size_of_CC);
//        
//        histogram_of_CCs[threshold] = number_of_CCs;
//        if(number_of_CCs==0) number_of_CCs = 1;
//        //cout<<"["<<threshold<<"]:"<<size_of_CC<<"/"<<number_of_CCs<<"="<<(size_of_CC/number_of_CCs)<<"  ";
//        cout<<"["<<threshold<<"]:"<<(number_of_CCs)<<"  ";
//    }
//    
//    
//    // temporarily!!!!!!!!!!!!!!!!!
//    output = slice;
//    
//    
//    
//	if(this->IsAbortRequested()) return 0;
//	this->SetProgressCounterRelativeValue(100);
//
//	return 1;
//}




int bdMelanoma::Downsample(bdImage &input, bdImage &output, int size_of_new_pixel)
{
    if(input.IsEmpty()) return 0;
    
    int dimension_r_of_output_image = input.GetNumberOfRows() / size_of_new_pixel;
    int dimension_c_of_output_image = input.GetNumberOfColumns() / size_of_new_pixel;
    
    int half_of_new_pixel_area = (size_of_new_pixel*size_of_new_pixel)*8/10;//*2/3;
    
    output.SetSize(1,1, dimension_r_of_output_image, dimension_c_of_output_image);
    output.SetVisualizationPropertiesToMatchInput(input);
    
    for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
        {
            int n_of_non_zero_pixels = 0;
            for(unsigned int rn=0; rn<size_of_new_pixel; rn++)
            {
                for(unsigned int cn=0; cn<size_of_new_pixel; cn++)
                {
                    if(input( (r*size_of_new_pixel+rn),(c*size_of_new_pixel+cn) )!=0)
                    {
                        n_of_non_zero_pixels++;
                    }
                }
            }
            
            if(n_of_non_zero_pixels > half_of_new_pixel_area)
            {
                output(r,c) = 255;
            }
            else
            {
                output(r,c) = 0;
            }
        }
    }
    
    return 1;
}


int bdMelanoma::Downsample_Mean(bdImage &input, bdImage &output, int size_of_new_pixel)
{
    if(input.IsEmpty()) return 0;
    
    int dimension_r_of_output_image = input.GetNumberOfRows() / size_of_new_pixel;
    int dimension_c_of_output_image = input.GetNumberOfColumns() / size_of_new_pixel;
    
    output.SetSize(input.GetNumberOfTimeSeries(),1, dimension_r_of_output_image, dimension_c_of_output_image);
    output.SetVisualizationPropertiesToMatchInput(input);
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                int n_of_pixels = 0;
                int sum = 0;
                for(unsigned int rn=0; rn<size_of_new_pixel; rn++)
                {
                    for(unsigned int cn=0; cn<size_of_new_pixel; cn++)
                    {
                        sum += input(t, 0, (r*size_of_new_pixel+rn), (c*size_of_new_pixel+cn));
                        n_of_pixels++;
                    }
                }
                
                if(n_of_pixels==0)
                {
                    output(t,0,r,c) = 0;
                }
                else
                {
                    output(t,0,r,c) = sum / n_of_pixels;
                }
            }
        }
    }
    
    return 1;
}


int bdMelanoma::Downsample_Max(bdImage &input, bdImage &output, int size_of_new_pixel)
{
    if(input.IsEmpty()) return 0;
    
    int dimension_r_of_output_image = input.GetNumberOfRows() / size_of_new_pixel;
    int dimension_c_of_output_image = input.GetNumberOfColumns() / size_of_new_pixel;
    
    output.SetSize(input.GetNumberOfTimeSeries(),1, dimension_r_of_output_image, dimension_c_of_output_image);
    output.SetVisualizationPropertiesToMatchInput(input);
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                int max = 0;
                for(unsigned int rn=0; rn<size_of_new_pixel; rn++)
                {
                    for(unsigned int cn=0; cn<size_of_new_pixel; cn++)
                    {
                        if(input(t, 0, (r*size_of_new_pixel+rn), (c*size_of_new_pixel+cn)) > max) max = input(t, 0, (r*size_of_new_pixel+rn), (c*size_of_new_pixel+cn));
                    }
                }
                
                output(t,0,r,c) = max;
            }
        }
    }
    
    return 1;
}


int bdMelanoma::MaskWithDownsampledImage(bdImage &input, bdImage &downsampled_mask, bdImage &output, int size_of_downsampled_pixel)
{
    if(input.IsEmpty()) return 0;
    if(downsampled_mask.IsEmpty()) return 0;
    
    output.SetSizeAndPropertiesAs(input);
    output.FillInWith(0);
    
    for(unsigned int r=0; r<downsampled_mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<downsampled_mask.GetNumberOfColumns(); c++)
        {
            if(downsampled_mask(r,c)!=0)
            {
                for(unsigned int rn=0; rn<size_of_downsampled_pixel; rn++)
                {
                    for(unsigned int cn=0; cn<size_of_downsampled_pixel; cn++)
                    {
                        output( (r*size_of_downsampled_pixel+rn),(c*size_of_downsampled_pixel+cn) ) = input( (r*size_of_downsampled_pixel+rn),(c*size_of_downsampled_pixel+cn) );
                    }
                }
            }
        }
    }
    
    return 1;
}


int bdMelanoma::FillInSegmentedMelanoma(bdImage &original_image, bdImage &segmented_image, bdImage &output)
{
    if(segmented_image.IsEmpty()) return 0;
    
    bdGeometry g;
    g.SetDimensions(1,segmented_image.GetNumberOfRows(),segmented_image.GetNumberOfColumns());
    
    cout<<" 1 ";
    
    // Find ALL borders of the melanoma and store them to an image
    bdImage border_image;
    border_image.SetSize(1,1,segmented_image.GetNumberOfRows(),segmented_image.GetNumberOfColumns());
    border_image.SetVisualizationPropertiesToMatchInput(segmented_image);
    border_image.FillInWith(0);
    
    for(unsigned int r=0; r<segmented_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_image.GetNumberOfColumns(); c++)
        {
            // if pixels belongs to the mask...
            if(segmented_image(r,c)!=0 )
            {
                //... check if it is a border
                int rn,cn;
                for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
                {
                    if(segmented_image(rn,cn)==0)
                    {
                        border_image(r,c) = 255;
                        break;
                    }
                }
            }
        }
    }
    
    cout<<" 2 ";
    
    // Extract the largest border image.
    bdBIP bip;
    bdImage largest_border_image;
    bip.ExtractLargest_26_ConnectedComponent(border_image, 1, largest_border_image);
    
    
    cout<<" 3 ";
    
    // Go through the borders images and process those borders that are not the largest...
    output.CopyFrom(segmented_image);
    for(unsigned int r=0; r<segmented_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_image.GetNumberOfColumns(); c++)
        {
            // if the pixels belongs to a border which is not the largest...
            if(border_image(r,c)!=0 && largest_border_image(r,c)==0)
            {
                //... find the background pixel...
                int rn,cn;
                for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
                {
                    if(output(rn,cn)==0)
                    {
                        //... and fill in the background region
                        bdVoxel first_voxel;
                        first_voxel.R() = rn;
                        first_voxel.C() = cn;
                        
                        bdList< bdVoxel > temp_list;
                        //bdImage temp_image;
                        //temp_image.SetSizeAndPropertiesAs(segmented_image);
                        //temp_image.FillInWith(0);
                        temp_list.AddToLeftEnd(first_voxel);
                        //temp_image(first_voxel.R(),first_voxel.C()) = 255;
                        output(first_voxel.R(),first_voxel.C()) = original_image(first_voxel.R(),first_voxel.C());
                        while(!temp_list.IsEmpty())
                        {
                            bdVoxel v = temp_list.GetLeftEnd();
                            
                            // Find the next background pixel and enter it to the list
                            {
                                int rnn,cnn;
                                for(g.ForCoordinates_4_Neighbors(v.R(),v.C()); g.Get_4_Neighbors(rnn,cnn); )
                                {
                                    //if(largest_border_image(rnn,cnn)!=0 && temp_image(rnn,cnn)==0)
                                    if(output(rnn,cnn)==0)
                                    {
                                        output(rnn,cnn)=original_image(rnn,cnn);
                                        if(output(rnn,cnn)==0) output(rnn,cnn) = 1;
                                        bdVoxel v1;
                                        temp_list.AddToRightEnd(v1(rnn,cnn));
                                        //break;
                                    }
                                }
                            }
                            
                            temp_list.DeleteLeftEnd();
                        }
                        
                        break;
                    }
                }
            }
        }
    }
    
    //output.CopyFrom(largest_border_image);
    
    return 1;
    
}



int bdMelanoma::SegmentMelanoma_Thresholding(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    // Downsample the image to work faster.
    bdGIP gip;
//    bdImage downsampled;
//    gip.DownsampleByCoefficient(input_image, 1, 10, 10, downsampled);
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input); //gip.RescaleWholeRange(downsampled, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    

    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
    //cout<<" 4 ";
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
 
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    // The threshold value is the index in between index_of_max_value1 and index_of_max_value2
    unsigned int threshold = (index_of_max_value1 + index_of_max_value2) / 2;
    
    cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"("<<index_of_max_value1<<","<<index_of_max_value2<<") ";
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    

    bdBIP bip;
    bdImage extracted_image;
    bip.ExtractLargest_26_ConnectedComponent(inverted_image, inverted_threshold, output);
    
    
//    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
//    output.SetVisualizationPropertiesToMatchInput(input);
//    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
//        {
//            slice(0,0,r,c) = input(max_variance_index,0,r,c);
//        }
//    }
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}


int bdMelanoma::SegmentMelanoma_Thresholding2(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
//    //cout<<" 4 ";
//    
//    // Calculate maxima in the histogram of the chosen slice.
//    bdArray<int> maxima;
//    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
    
    
    // Go through the histogram of chosen slice to determine the maximum peak in the histogram.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<array_of_histograms[max_variance_index].GetNumberOfElements(); i++)
    {
        if(array_of_histograms[max_variance_index][i]>max_value1)
        {
            max_value1 = array_of_histograms[max_variance_index][i];
            index_of_max_value1 = i;
        }
    }
    
    // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
    // in array of maxima that preceed the maximum peak.
    unsigned int max_difference = 0;
    unsigned int index_of_max_difference = 0;
    for(unsigned int i=index_of_max_value1-1; i>0; i--)
    {
        int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
        if(d<0) d = -d;
        if(d>max_difference)
        {
            max_difference = d;
            index_of_max_difference = i;
        }
    }
    
    // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
    // y = kx + n ;
    // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
    // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
    
    double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
    double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
    
    // For y=0: x = -n/k;
    double x_cross = - n/k;
    cout<<"  x_cross="<<x_cross<<" ";
    
    
    // The threshold value is the index in of the crossing with x axis
    unsigned int threshold = (unsigned int) x_cross;
    
    cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max_diff="<<index_of_max_difference<<") ";
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    
    
    // Find a non-zero pixel that is the closest to the image center
    bdGeometry g;
    g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
    unsigned int center_r = input_image.GetNumberOfRows() / 2;
    unsigned int center_c = input_image.GetNumberOfColumns() / 2;
    int r,c;
    int found_r = 0, found_c = 0;
    for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, r, c); )
    {
        if(inverted_image(r,c)>inverted_threshold)
        {
            found_r = r;
            found_c = c;
            break;
        }
    }
    
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    
    bdBIP bip;
    bdImage extracted_image;
    //bip.ExtractLargest_26_ConnectedComponent(inverted_image, inverted_threshold, output);
    bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, output, 0, found_r, found_c);
    
    //ExtractLargest_26_ConnectedComponent(inverted_image, inverted_threshold, output);

    
    
    //    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
    //    output.SetVisualizationPropertiesToMatchInput(input);
    //    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    //    {
    //        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
    //        {
    //            slice(0,0,r,c) = input(max_variance_index,0,r,c);
    //        }
    //    }
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}



int bdMelanoma::SegmentMelanoma_Thresholding3(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
    //    //cout<<" 4 ";
    //
    //    // Calculate maxima in the histogram of the chosen slice.
    //    bdArray<int> maxima;
    //    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
    
    
    // Go through the histogram of chosen slice to determine the maximum peak in the histogram.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<array_of_histograms[max_variance_index].GetNumberOfElements(); i++)
    {
        if(array_of_histograms[max_variance_index][i]>max_value1)
        {
            max_value1 = array_of_histograms[max_variance_index][i];
            index_of_max_value1 = i;
        }
    }
    
    // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
    // in array of maxima that preceed the maximum peak.
    unsigned int max_difference = 0;
    unsigned int index_of_max_difference = 0;
    for(unsigned int i=index_of_max_value1-1; i>0; i--)
    {
        int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
        if(d<0) d = -d;
        if(d>max_difference)
        {
            max_difference = d;
            index_of_max_difference = i;
        }
    }
    
    // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
    // y = kx + n ;
    // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
    // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
    
    double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
    double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
    
    // For y=0: x = -n/k;
    double x_cross = - n/k;
    cout<<"  x_cross="<<x_cross<<" ";
    
    
    // The threshold value is the index in of the crossing with x axis
    unsigned int threshold = (unsigned int) x_cross;
    
    cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max_diff="<<index_of_max_difference<<") ";
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    
    
    // Find a non-zero pixel that is the closest to the image center
    bdGeometry g;
    g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
    unsigned int center_r = input_image.GetNumberOfRows() / 2;
    unsigned int center_c = input_image.GetNumberOfColumns() / 2;
    int rn,cn;
    int found_r = 0, found_c = 0;
    for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
    {
        if(inverted_image(rn,cn)>inverted_threshold)
        {
            found_r = rn;
            found_c = cn;
            break;
        }
    }
    
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    
    bdBIP bip;
    bdImage extracted_image;
    bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
    
    
    bdMelanoma melanoma;
    bdImage downsampled_image;
    melanoma.Downsample(extracted_image, downsampled_image, 50);
    
    
    bdImage extracted_downsampled_image;
    bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
    
    // Dilate
    bdImage dilated_extracted_downsampled_image;
    dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
    dilated_extracted_downsampled_image.FillInWith(0);
    bdGeometry g2;
    g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
    for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
        {
            if(extracted_downsampled_image(r,c)!=0)
            {
                for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                {
                    dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                }
            }
        }
    }
    
    
    
    bdImage masked_image;
    melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
    
    


    output.CopyFrom(masked_image);

    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}



int bdMelanoma::SegmentMelanoma_Thresholding4(bdImage &input_image, bdImage &output, bdImage &output2)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
//    //    //cout<<" 4 ";
//    //
//    //    // Calculate maxima in the histogram of the chosen slice.
//    //    bdArray<int> maxima;
//    //    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
//    
//    
//    // Go through the histogram of chosen slice to determine the maximum peak in the histogram.
//    // Find the max value...
//    unsigned int max_value1 = 0;
//    unsigned int index_of_max_value1 = 0;
//    for(unsigned int i=0; i<array_of_histograms[max_variance_index].GetNumberOfElements(); i++)
//    {
//        if(array_of_histograms[max_variance_index][i]>max_value1)
//        {
//            max_value1 = array_of_histograms[max_variance_index][i];
//            index_of_max_value1 = i;
//        }
//    }
    
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }

    
    
    //----- Determine threshold value 1 -----
    unsigned int threshold;
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that preceed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold = (unsigned int) x_cross;
        
        cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------
    
    //----- Determine threshold value 2 -----
    unsigned int threshold2;
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that succeed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
        {
            int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold2 = (unsigned int) x_cross;
        
        cout<<" ["<<max_variance_index<<"]threshold2="<<threshold2<<"(index_of_max="<<index_of_max_value1<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------

    
    //----- Perform extraction 1 -----
    {
        bdImage inverted_image;
        gip.Negative(slice, inverted_image);
        int inverted_threshold = 255 - threshold;
        
        
        // Find a non-zero pixel that is the closest to the image center
        bdGeometry g;
        g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
        unsigned int center_r = input_image.GetNumberOfRows() / 2;
        unsigned int center_c = input_image.GetNumberOfColumns() / 2;
        int rn,cn;
        int found_r = 0, found_c = 0;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            if(inverted_image(rn,cn)>inverted_threshold)
            {
                found_r = rn;
                found_c = cn;
                break;
            }
        }
        
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
        
        bdBIP bip;
        bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        bdMelanoma melanoma;
        bdImage downsampled_image;
        melanoma.Downsample(extracted_image, downsampled_image, 50);
        
        
        bdImage extracted_downsampled_image;
        bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
        
        // Dilate
        bdImage dilated_extracted_downsampled_image;
        dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
        dilated_extracted_downsampled_image.FillInWith(0);
        bdGeometry g2;
        g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
        for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
            {
                if(extracted_downsampled_image(r,c)!=0)
                {
                    for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                    {
                        dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                    }
                }
            }
        }
        
        bdImage masked_image;
        melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
        
        output.CopyFrom(masked_image);
    }
    //------------------------------
    
    
    //----- Perform extraction 2 -----
    {
        bdImage inverted_image;
        gip.Negative(slice, inverted_image);
        int inverted_threshold = 255 - threshold2;
        
        
        // Find a non-zero pixel that is the closest to the image center
        bdGeometry g;
        g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
        unsigned int center_r = input_image.GetNumberOfRows() / 2;
        unsigned int center_c = input_image.GetNumberOfColumns() / 2;
        int rn,cn;
        int found_r = 0, found_c = 0;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            if(inverted_image(rn,cn)>inverted_threshold)
            {
                found_r = rn;
                found_c = cn;
                break;
            }
        }
        
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
        
        bdBIP bip;
        bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        bdMelanoma melanoma;
        bdImage downsampled_image;
        melanoma.Downsample(extracted_image, downsampled_image, 50);
        
        
        bdImage extracted_downsampled_image;
        bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
        
        // Dilate
        bdImage dilated_extracted_downsampled_image;
        dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
        dilated_extracted_downsampled_image.FillInWith(0);
        bdGeometry g2;
        g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
        for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
            {
                if(extracted_downsampled_image(r,c)!=0)
                {
                    for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                    {
                        dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                    }
                }
            }
        }
        
        bdImage masked_image;
        melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
        
        output2.CopyFrom(masked_image);
    }
    //------------------------------
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}


int bdMelanoma::SegmentMelanoma_Thresholding5(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
    
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    
    unsigned int threshold;
    //----- Determine threshold value 1 -----
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that preceed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold = (unsigned int) x_cross;
        
        cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------
    
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
    
    
    // Find a non-zero pixel that is the closest to the image center
    int found_r = -1, found_c = -1; // indexes are set to an error value.
    {
        bdGeometry g;
        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
        int rn,cn;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            //if(inverted_image(rn,cn)>inverted_threshold)
            if(thresholded_image(rn,cn)>0)
            {
                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
                found_r = rn;
                found_c = cn;
                int rnn, cnn;
                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
                {
                    thresholded_image(rnn,cnn) = 255;
                }
                
                break;
            }
        }
    }
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    //----- Perform extraction -----
    bdImage extracted_image;
    
    // Check if a non-zero pixel was found. If so, perform extraction
    if(found_r>=0 && found_c>=0 )
    {
        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
    }
    
    
    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
    
    
    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
    if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
    {
        //----- Determine threshold value 2 -----
        {
            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
            // in array of maxima that succeed the maximum peak.
            unsigned int max_difference = 0;
            unsigned int index_of_max_difference = 0;
            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
            {
                int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
                if(d<0) d = -d;
                if(d>max_difference)
                {
                    max_difference = d;
                    index_of_max_difference = i;
                }
            }
            
            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
            // y = kx + n ;
            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
            
            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
            
            // For y=0: x = -n/k;
            double x_cross = - n/k;
            //cout<<"  x_cross="<<x_cross<<" ";
            
            
            // The threshold value is the index in of the crossing with x axis
            threshold = (unsigned int) x_cross;
            
            cout<<" ["<<max_variance_index<<"]threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
        }
        //------------------------------
        
        
        inverted_threshold = 255 - threshold;
        
        // Find a non-zero pixel that is the closest to the image center
        found_r = -1; found_c = -1; // indexes are set to an error value.
        {
            bdGeometry g;
            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
            unsigned int center_r = input_image.GetNumberOfRows() / 2;
            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
            int rn,cn;
            
            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
            {
                if(inverted_image(rn,cn)>inverted_threshold)
                {
                    found_r = rn;
                    found_c = cn;
                    break;
                }
            }
        }
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
    }
    
    
    //----- Perform extraction -----
    {
        //bdBIP bip;
        //bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        // Check if we need to perform masking by looking if the segmented region reaches the image edge
        int is_segmented_lesion_touching_image_edge = 0;
        //---- Check image edges -----
        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
        {
            if(extracted_image(r,0)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
        }
        if(!is_segmented_lesion_touching_image_edge)
        {
            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
            {
                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
            }
        }
        //---------------------------
        
        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
        
        // If the lesion is not touching image edge, we consider it a correct segmentation.
        if(is_segmented_lesion_touching_image_edge)
        {
            bdMelanoma melanoma;
            bdImage downsampled_image;
            melanoma.Downsample(extracted_image, downsampled_image, 50);
            
            
            bdImage extracted_downsampled_image;
            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
            
            // Dilate
            bdImage dilated_extracted_downsampled_image;
            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
            dilated_extracted_downsampled_image.FillInWith(0);
            bdGeometry g2;
            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
                {
                    if(extracted_downsampled_image(r,c)!=0)
                    {
                        int rn,cn;
                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                        {
                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                        }
                    }
                }
            }
            
            bdImage masked_image;
            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
            
            output.CopyFrom(masked_image);
        }
        else
        {
            output.CopyFrom(extracted_image);
            return 1;
        }
        
    }
    //------------------------------
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}





int bdMelanoma::SegmentMelanoma_Thresholding6(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    
    
    bdArray< bdArray<int> > array_of_histograms;
    array_of_histograms.Set(input.GetNumberOfTimeSeries());
    
    // For each slice calculate the histogram
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(input, array_of_histograms[t], 0,255, t, 0);
    }
    
    //cout<<" 2 ";
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum += array_of_histograms[t][i] * i;
            n += array_of_histograms[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<array_of_histograms[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * array_of_histograms[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
    
    //cout<<" 3 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            slice(0,0,r,c) = input(max_variance_index,0,r,c);
        }
    }
    
    
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(array_of_histograms[max_variance_index], maxima);
    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    
    unsigned int threshold;
    //----- Determine threshold value 1 -----
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that preceed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold = (unsigned int) x_cross;
        
        cout<<" ["<<max_variance_index<<"]threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------



    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);

    
    // Find a non-zero pixel that is the closest to the image center
    int found_r = -1, found_c = -1; // indexes are set to an error value.
    {
        bdGeometry g;
        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
        int rn,cn;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            //if(inverted_image(rn,cn)>inverted_threshold)
            if(thresholded_image(rn,cn)>0)
            {
                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
                found_r = rn;
                found_c = cn;
                int rnn, cnn;
                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
                {
                    thresholded_image(rnn,cnn) = 255;
                }
                
                break;
            }
        }
    }
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    //----- Perform extraction -----
    bdImage extracted_image;
    
    // Check if a non-zero pixel was found. If so, perform extraction
    if(found_r>=0 && found_c>=0 )
    {
        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
    }

    
    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
    

    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
    if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
    {
        //----- Determine threshold value 2 -----
        {
            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
            // in array of maxima that succeed the maximum peak.
            unsigned int max_difference = 0;
            unsigned int index_of_max_difference = 0;
            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
            {
                int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
                if(d<0) d = -d;
                if(d>max_difference)
                {
                    max_difference = d;
                    index_of_max_difference = i;
                }
            }
            
            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
            // y = kx + n ;
            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
            
            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
            
            // For y=0: x = -n/k;
            double x_cross = - n/k;
            //cout<<"  x_cross="<<x_cross<<" ";
            
            
            // The threshold value is the index in of the crossing with x axis
            threshold = (unsigned int) x_cross;
            
            cout<<" ["<<max_variance_index<<"]threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
        }
        //------------------------------

        
        inverted_threshold = 255 - threshold;
        
        // Find a non-zero pixel that is the closest to the image center
        found_r = -1; found_c = -1; // indexes are set to an error value.
        {
            bdGeometry g;
            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
            unsigned int center_r = input_image.GetNumberOfRows() / 2;
            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
            int rn,cn;
            
            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
            {
                if(inverted_image(rn,cn)>inverted_threshold)
                {
                    found_r = rn;
                    found_c = cn;
                    break;
                }
            }
        }
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";

    }

    
    //----- Perform extraction -----
    {
        //bdBIP bip;
        //bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        // Check if we need to perform masking by looking if the segmented region reaches the image edge
        int is_segmented_lesion_touching_image_edge = 0;
        //---- Check image edges -----
        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
        {
            if(extracted_image(r,0)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
        }
        if(!is_segmented_lesion_touching_image_edge)
        {
            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
            {
                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
            }
        }
        //---------------------------
        
        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
        
        // If the lesion is not touching image edge, we consider it a correct segmentation.
        if(!is_segmented_lesion_touching_image_edge)
        {
            bdMelanoma melanoma;
            melanoma.FillInSegmentedMelanoma(input,extracted_image,output);
            //output.CopyFrom(extracted_image);
            return 1;
        }
        
        
        bdMelanoma melanoma;
        bdImage downsampled_image;
        melanoma.Downsample(extracted_image, downsampled_image, 50);
        
        
        bdImage extracted_downsampled_image;
        bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
        
        // Dilate
        bdImage dilated_extracted_downsampled_image;
        dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
        dilated_extracted_downsampled_image.FillInWith(0);
        bdGeometry g2;
        g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
        for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
            {
                if(extracted_downsampled_image(r,c)!=0)
                {
                    int rn,cn;
                    for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                    {
                        dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                    }
                }
            }
        }
        
        bdImage masked_image;
        melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
        
        //output.CopyFrom(masked_image);

        
        bdImage extracted_upsampled_image;
        bip.ExtractLargest_26_ConnectedComponent(masked_image, 0, extracted_upsampled_image);

        
        
        melanoma.FillInSegmentedMelanoma(input,extracted_upsampled_image,output);
        
        
    }
    //------------------------------

    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}


//int bdMelanoma::SegmentMelanoma_Thresholding7(bdImage &input_image, bdImage &output)
//{
//    //cout<<" 1 ";
//    
//    if(input_image.IsEmpty()) return 0;
//    
//    
//    bdGIP gip;
//    
//    
////    // Rescale the slice values to fall within 0,255.
////    unsigned int lower_value = 0, upper_value = 255;
////    bdImage input;
////    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
//    
//    // Find average value for blue and green channel
//    unsigned int sum_green = 0, sum_blue = 0;
//    for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
//        {
//            sum_green += input_image(1,0,r,c);
//            sum_blue += input_image(2,0,r,c);
//        }
//    }
//    unsigned int avr_green = sum_green / (input_image.GetNumberOfRows() * input_image.GetNumberOfColumns());
//    unsigned int avr_blue = sum_blue / (input_image.GetNumberOfRows() * input_image.GetNumberOfColumns());
//    
//    
//    cout<<"avr(b,g)="<<avr_blue<<","<<avr_green;
//    
//    bdImage in_img;
//    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
//    in_img.SetSize(1,1, input_image.GetNumberOfRows(), input_image.GetNumberOfColumns());
//    in_img.SetVisualizationPropertiesToMatchInput(input_image);
//    for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
//        {
//            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
//
//            // MAX OF BLUE AND GREEN CORRECTED
//            //if(avr_green>avr_blue)
//            {
//                if(input_image(1,0,r,c) > input_image(2,0,r,c) +1000) in_img(0,0,r,c) = input_image(1,0,r,c);
//                else in_img(0,0,r,c) = input_image(2,0,r,c);// + avr_green-avr_blue;
//            }
//            //else
//            //{
//            //    if(input_image(1,0,r,c) > input_image(2,0,r,c) ) in_img(0,0,r,c) = (input_image(1,0,r,c)+avr_green-avr_blue);
//            //    else in_img(0,0,r,c) = input_image(2,0,r,c);
//            //}
//            
////            // MAX OF BLUE AND GREEN
////            if(input_image(1,0,r,c)>input_image(2,0,r,c)) in_img(0,0,r,c) = input_image(1,0,r,c);
////            else in_img(0,0,r,c) = input_image(2,0,r,c);
//            
//            
////            // BLUE
////            in_img(0,0,r,c) = input_image(2,0,r,c);
//            
////            // MAX OF ALL CHANNELS
////            if(input_image(1,0,r,c)>input_image(2,0,r,c)) in_img(0,0,r,c) = input_image(1,0,r,c);
////            else in_img(0,0,r,c) = input_image(2,0,r,c);
////            if(input_image(0,0,r,c)>in_img(0,0,r,c)) in_img(0,0,r,c) = input_image(0,0,r,c);
////            
////            // MIN OF ALL CHANNELS
////            if(input_image(1,0,r,c)<input_image(2,0,r,c)) in_img(0,0,r,c) = input_image(1,0,r,c);
////            else in_img(0,0,r,c) = input_image(2,0,r,c);
////            if(input_image(0,0,r,c)<in_img(0,0,r,c)) in_img(0,0,r,c) = input_image(0,0,r,c);
//
//            
//            
//            //if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//            //else slice(0,0,r,c) = input(0,0,r,c);
//        }
//    }
//    
//    
////    bdImage in_img2;
////    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
////    in_img2.SetSize(1,1, input_image.GetNumberOfRows(), input_image.GetNumberOfColumns());
////    in_img2.SetVisualizationPropertiesToMatchInput(input_image);
////    for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
////    {
////        for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
////        {
////            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
////            
////            //if(input_image(1,0,r,c)>input_image(2,0,r,c)) in_img(0,0,r,c) = input_image(1,0,r,c);
////            //else in_img(0,0,r,c) = input_image(2,0,r,c);
////            
////            //            // BLUE
////            //            in_img(0,0,r,c) = input_image(2,0,r,c);
////            
////            //            // MAX OF ALL CHANNELS
////            //            if(input_image(1,0,r,c)>input_image(2,0,r,c)) in_img(0,0,r,c) = input_image(1,0,r,c);
////            //            else in_img(0,0,r,c) = input_image(2,0,r,c);
////            //            if(input_image(0,0,r,c)>in_img(0,0,r,c)) in_img(0,0,r,c) = input_image(0,0,r,c);
////            
////            // MIN OF ALL CHANNELS
////            if(input_image(1,0,r,c)<input_image(2,0,r,c)) in_img2(0,0,r,c) = input_image(1,0,r,c);
////            else in_img2(0,0,r,c) = input_image(2,0,r,c);
////            if(input_image(0,0,r,c)<in_img2(0,0,r,c)) in_img2(0,0,r,c) = input_image(0,0,r,c);
////            
////            
////            
////            //if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
////            //else slice(0,0,r,c) = input(0,0,r,c);
////        }
////    }
////    
////    bdImage res_img;
////    res_img.SetSize(1,1, input_image.GetNumberOfRows(), input_image.GetNumberOfColumns());
////    res_img.SetVisualizationPropertiesToMatchInput(input_image);
////    for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
////    {
////        for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
////        {
////            res_img(0,0,r,c) = in_img(0,0,r,c) - in_img2(0,0,r,c);
////        }
////    }
//    
//    output.CopyFrom(in_img);
//    //output.CopyFrom(res_img);
//
//
//    
//    
//    // Rescale the slice values to fall within 0,255.
//    bdImage slice;
//    unsigned int lower_value = 0, upper_value = 255;
//    bdImage input;
//    gip.RescaleWholeRange(in_img, lower_value, upper_value, slice);
//
//    
//    
//    if(this->IsAbortRequested()) return 0;
//    this->SetProgressCounterRelativeValue(5);
//    
//    
//////    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
//////    bdImage slice;
//////    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
//////    slice.SetVisualizationPropertiesToMatchInput(input);
//////    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
//////    {
//////        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
//////        {
//////            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
//////            
//////            if(input(1,0,r,c)>input(2,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//////            else slice(0,0,r,c) = input(2,0,r,c);
//////            
//////            //if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//////            //else slice(0,0,r,c) = input(0,0,r,c);
//////        }
//////    }
////    
////    //Calculate the histogram of the slice
////    bdArray<int> histogram_of_slice;
////    this->Histogram(slice, histogram_of_slice, 0,255, 0, 0);
////    
////    
////    
////    
////    
////    // Calculate maxima in the histogram of the chosen slice.
////    bdArray<int> maxima;
////    this->HistogramMaxima(histogram_of_slice, maxima);
////    
////    
////    
////    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
////    // Find the max value...
////    unsigned int max_value1 = 0;
////    unsigned int index_of_max_value1 = 0;
////    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
////    {
////        if(maxima[i]>max_value1)
////        {
////            max_value1 = maxima[i];
////            index_of_max_value1 = i;
////        }
////    }
////    
////    // Find the next max value...
////    unsigned int max_value2 = 0;
////    unsigned int index_of_max_value2 = 0;
////    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
////    {
////        if(maxima[i]>max_value2)
////        {
////            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
////            //if(i!=index_of_max_value1)
////            //{
////            if(i>index_of_max_value1)
////            {
////                if(i-index_of_max_value1>10)
////                {
////                    max_value2 = maxima[i];
////                    index_of_max_value2 = i;
////                }
////            }
////            if(i<index_of_max_value1)
////            {
////                if(index_of_max_value1-i>10)
////                {
////                    max_value2 = maxima[i];
////                    index_of_max_value2 = i;
////                }
////            }
////            //}
////        }
////    }
////    
////    
////    unsigned int threshold;
////    //----- Determine threshold value 1 -----
////    {
////        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
////        // in array of maxima that preceed the maximum peak.
////        unsigned int max_difference = 0;
////        unsigned int index_of_max_difference = 0;
////        for(unsigned int i=index_of_max_value1-1; i>0; i--)
////        {
////            //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
////            int d = histogram_of_slice[i+1] - histogram_of_slice[i];
////            if(d<0) d = -d;
////            if(d>max_difference)
////            {
////                max_difference = d;
////                index_of_max_difference = i;
////            }
////        }
////        
////        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
////        // y = kx + n ;
////        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
////        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
////        
////        //        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
////        //        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
////        double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
////        double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
////        
////        
////        // For y=0: x = -n/k;
////        double x_cross = - n/k;
////        //cout<<"  x_cross="<<x_cross<<" ";
////        
////        
////        // The threshold value is the index in of the crossing with x axis
////        threshold = (unsigned int) x_cross;
////        
////        cout<<" threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
////    }
////    //------------------------------
////    
////    
////    
////    bdImage inverted_image;
////    gip.Negative(slice, inverted_image);
////    int inverted_threshold = 255 - threshold;
////    
////    bdImage thresholded_image;
////    bdBIP bip;
////    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
////    
////    
////    // Find a non-zero pixel that is the closest to the image center
////    int found_r = -1, found_c = -1; // indexes are set to an error value.
////    {
////        bdGeometry g;
////        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
////        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
////        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
////        int rn,cn;
////        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
////        {
////            //if(inverted_image(rn,cn)>inverted_threshold)
////            if(thresholded_image(rn,cn)>0)
////            {
////                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
////                found_r = rn;
////                found_c = cn;
////                int rnn, cnn;
////                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
////                {
////                    thresholded_image(rnn,cnn) = 255;
////                }
////                
////                break;
////            }
////        }
////    }
////    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
////    
////    //----- Perform extraction -----
////    bdImage extracted_image;
////    
////    // Check if a non-zero pixel was found. If so, perform extraction
////    if(found_r>=0 && found_c>=0 )
////    {
////        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
////    }
////    
////    
////    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
////    
////    
////    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
////    if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
////    {
////        //----- Determine threshold value 2 -----
////        {
////            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
////            // in array of maxima that succeed the maximum peak.
////            unsigned int max_difference = 0;
////            unsigned int index_of_max_difference = 0;
////            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
////            {
////                //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
////                int d = histogram_of_slice[i+1] - histogram_of_slice[i];
////                
////                if(d<0) d = -d;
////                if(d>max_difference)
////                {
////                    max_difference = d;
////                    index_of_max_difference = i;
////                }
////            }
////            
////            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
////            // y = kx + n ;
////            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
////            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
////            
////            //            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
////            //            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
////            double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
////            double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
////            
////            
////            // For y=0: x = -n/k;
////            double x_cross = - n/k;
////            //cout<<"  x_cross="<<x_cross<<" ";
////            
////            
////            // The threshold value is the index in of the crossing with x axis
////            threshold = (unsigned int) x_cross;
////            
////            cout<<" threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
////        }
////        //------------------------------
////        
////        
////        inverted_threshold = 255 - threshold;
////        
////        // Find a non-zero pixel that is the closest to the image center
////        found_r = -1; found_c = -1; // indexes are set to an error value.
////        {
////            bdGeometry g;
////            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
////            unsigned int center_r = input_image.GetNumberOfRows() / 2;
////            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
////            int rn,cn;
////            
////            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
////            {
////                if(inverted_image(rn,cn)>inverted_threshold)
////                {
////                    found_r = rn;
////                    found_c = cn;
////                    break;
////                }
////            }
////        }
////        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
////        
////    }
////    
////    
////    //----- Perform extraction -----
////    {
////        //bdBIP bip;
////        //bdImage extracted_image;
////        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
////        
////        
////        // Check if we need to perform masking by looking if the segmented region reaches the image edge
////        int is_segmented_lesion_touching_image_edge = 0;
////        //---- Check image edges -----
////        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
////        {
////            if(extracted_image(r,0)!=0)
////            {
////                is_segmented_lesion_touching_image_edge = 1;
////                break;
////            }
////            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
////            {
////                is_segmented_lesion_touching_image_edge = 1;
////                break;
////            }
////        }
////        if(!is_segmented_lesion_touching_image_edge)
////        {
////            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
////            {
////                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
////                {
////                    is_segmented_lesion_touching_image_edge = 1;
////                    break;
////                }
////                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
////                {
////                    is_segmented_lesion_touching_image_edge = 1;
////                    break;
////                }
////            }
////        }
////        //---------------------------
////        
////        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
////        
////        // If the lesion is not touching image edge, we consider it a correct segmentation.
////        if(is_segmented_lesion_touching_image_edge)
////        {
////            bdMelanoma melanoma;
////            bdImage downsampled_image;
////            melanoma.Downsample(extracted_image, downsampled_image, 50);
////            
////            
////            bdImage extracted_downsampled_image;
////            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
////            
////            // Dilate
////            bdImage dilated_extracted_downsampled_image;
////            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
////            dilated_extracted_downsampled_image.FillInWith(0);
////            bdGeometry g2;
////            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
////            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
////            {
////                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
////                {
////                    if(extracted_downsampled_image(r,c)!=0)
////                    {
////                        int rn,cn;
////                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
////                        {
////                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
////                        }
////                    }
////                }
////            }
////            
////            bdImage masked_image;
////            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
////            
////            output.CopyFrom(masked_image);
////        }
////        else
////        {
////            output.CopyFrom(extracted_image);
////            return 1;
////        }
////        
////    }
////    //------------------------------
//    
//    
//    //output.CopyFrom(slice);
//    //output.CopyFrom(in_img);
//
//    
//    
//    if(this->IsAbortRequested()) return 0;
//    this->SetProgressCounterRelativeValue(100);
//    
//    return 1;
//}




int bdMelanoma::SegmentMelanoma_Thresholding7(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            // Use just the blue channel.
            //slice(0,0,r,c) = input(2,0,r,c);
            
            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;

            // combined B and G channels
            if(input(1,0,r,c)>input(2,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
            else slice(0,0,r,c) = input(2,0,r,c);

            // combined R and G channels
//            if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//            else slice(0,0,r,c) = input(0,0,r,c);
            
            
        }
    }
    
    //Calculate the histogram of the slice
    bdArray<int> histogram_of_slice;
    this->Histogram(slice, histogram_of_slice, 0,255, 0, 0);
    
    


    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(histogram_of_slice, maxima);

    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    
    unsigned int threshold;
    //----- Determine threshold value 1 -----
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that preceed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            int d = histogram_of_slice[i+1] - histogram_of_slice[i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
//        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
//        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
        double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);

        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold = (unsigned int) x_cross;
        
        cout<<" threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------
    
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
    
    
    // Find a non-zero pixel that is the closest to the image center
    int found_r = -1, found_c = -1; // indexes are set to an error value.
    {
        bdGeometry g;
        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
        int rn,cn;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            //if(inverted_image(rn,cn)>inverted_threshold)
            if(thresholded_image(rn,cn)>0)
            {
                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
                found_r = rn;
                found_c = cn;
                int rnn, cnn;
                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
                {
                    thresholded_image(rnn,cnn) = 255;
                }
                
                break;
            }
        }
    }
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    //----- Perform extraction -----
    bdImage extracted_image;
    
    // Check if a non-zero pixel was found. If so, perform extraction
    if(found_r>=0 && found_c>=0 )
    {
        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
    }
    
    
    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
    
    
    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
    if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
    {
        //----- Determine threshold value 2 -----
        {
            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
            // in array of maxima that succeed the maximum peak.
            unsigned int max_difference = 0;
            unsigned int index_of_max_difference = 0;
            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
            {
                //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
                int d = histogram_of_slice[i+1] - histogram_of_slice[i];

                if(d<0) d = -d;
                if(d>max_difference)
                {
                    max_difference = d;
                    index_of_max_difference = i;
                }
            }
            
            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
            // y = kx + n ;
            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
            
//            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
//            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
            double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
            double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);

            
            // For y=0: x = -n/k;
            double x_cross = - n/k;
            //cout<<"  x_cross="<<x_cross<<" ";
            
            
            // The threshold value is the index in of the crossing with x axis
            threshold = (unsigned int) x_cross;
            
            cout<<" threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
        }
        //------------------------------
        
        
        inverted_threshold = 255 - threshold;
        
        // Find a non-zero pixel that is the closest to the image center
        found_r = -1; found_c = -1; // indexes are set to an error value.
        {
            bdGeometry g;
            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
            unsigned int center_r = input_image.GetNumberOfRows() / 2;
            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
            int rn,cn;
            
            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
            {
                if(inverted_image(rn,cn)>inverted_threshold)
                {
                    found_r = rn;
                    found_c = cn;
                    break;
                }
            }
        }
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
    }
    
    
    //----- Perform extraction -----
    {
        //bdBIP bip;
        //bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        // Check if we need to perform masking by looking if the segmented region reaches the image edge
        int is_segmented_lesion_touching_image_edge = 0;
        //---- Check image edges -----
        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
        {
            if(extracted_image(r,0)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
        }
        if(!is_segmented_lesion_touching_image_edge)
        {
            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
            {
                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
            }
        }
        //---------------------------
        
        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
        
        // If the lesion is not touching image edge, we consider it a correct segmentation.
        if(is_segmented_lesion_touching_image_edge)
        {
            bdMelanoma melanoma;
            bdImage downsampled_image;
            melanoma.Downsample(extracted_image, downsampled_image, 50);
            
            
            bdImage extracted_downsampled_image;
            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
            
            // Dilate
            bdImage dilated_extracted_downsampled_image;
            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
            dilated_extracted_downsampled_image.FillInWith(0);
            bdGeometry g2;
            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
                {
                    if(extracted_downsampled_image(r,c)!=0)
                    {
                        int rn,cn;
                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                        {
                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                        }
                    }
                }
            }
            
            bdImage masked_image;
            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
            
            output.CopyFrom(masked_image);
        }
        else
        {
            output.CopyFrom(extracted_image);
            return 1;
        }
        
    }
    //------------------------------
    
    
    //output.CopyFrom(slice);
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}




int bdMelanoma::SegmentMelanoma_Thresholding8(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            // Use just the blue channel.
            slice(0,0,r,c) = input(2,0,r,c);
            
            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
            
            // combined B and G channels
            //            if(input(1,0,r,c)>input(2,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
            //            else slice(0,0,r,c) = input(2,0,r,c);
            
            // combined R and G channels
            //            if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
            //            else slice(0,0,r,c) = input(0,0,r,c);
            
            
        }
    }
    
    //Calculate the histogram of the slice
    bdArray<int> histogram_of_slice;
    this->Histogram(slice, histogram_of_slice, 0,255, 0, 0);
    
    
    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
    //bdArray<int> histogram_of_slice_modified;
    //histogram_of_slice_modified.Set(histogram_of_slice.GetNumberOfElements());
    unsigned int last_non_zero_value = 1;
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        if(histogram_of_slice[i]==0)
        {
            histogram_of_slice[i] = last_non_zero_value;
        }
        else
        {
            last_non_zero_value = histogram_of_slice[i];
        }
    }
    
    
    
    
    
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(histogram_of_slice, maxima);
    
    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    
    unsigned int threshold;
    //----- Determine threshold value 1 -----
    {
        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
        // in array of maxima that preceed the maximum peak.
        unsigned int max_difference = 0;
        unsigned int index_of_max_difference = 0;
        for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
            int d = histogram_of_slice[i+1] - histogram_of_slice[i];
            if(d<0) d = -d;
            if(d>max_difference)
            {
                max_difference = d;
                index_of_max_difference = i;
            }
        }
        
        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
        // y = kx + n ;
        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
        
        //        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
        //        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
        double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
        double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
        
        
        // For y=0: x = -n/k;
        double x_cross = - n/k;
        //cout<<"  x_cross="<<x_cross<<" ";
        
        
        // The threshold value is the index in of the crossing with x axis
        threshold = (unsigned int) x_cross;
        
        cout<<" threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
    }
    //------------------------------
    
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
    
    
    // Find a non-zero pixel that is the closest to the image center
    int found_r = -1, found_c = -1; // indexes are set to an error value.
    {
        bdGeometry g;
        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
        int rn,cn;
        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        {
            //if(inverted_image(rn,cn)>inverted_threshold)
            if(thresholded_image(rn,cn)>0)
            {
                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
                found_r = rn;
                found_c = cn;
                int rnn, cnn;
                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
                {
                    thresholded_image(rnn,cnn) = 255;
                }
                
                break;
            }
        }
    }
    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    //----- Perform extraction -----
    bdImage extracted_image;
    
    // Check if a non-zero pixel was found. If so, perform extraction
    if(found_r>=0 && found_c>=0 )
    {
        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
    }
    
    
    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
    
    
    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
    if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
    {
        //----- Determine threshold value 2 -----
        {
            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
            // in array of maxima that succeed the maximum peak.
            unsigned int max_difference = 0;
            unsigned int index_of_max_difference = 0;
            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
            {
                //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
                int d = histogram_of_slice[i+1] - histogram_of_slice[i];
                
                if(d<0) d = -d;
                if(d>max_difference)
                {
                    max_difference = d;
                    index_of_max_difference = i;
                }
            }
            
            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
            // y = kx + n ;
            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
            
            //            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
            //            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
            double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
            double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
            
            
            // For y=0: x = -n/k;
            double x_cross = - n/k;
            //cout<<"  x_cross="<<x_cross<<" ";
            
            
            // The threshold value is the index in of the crossing with x axis
            threshold = (unsigned int) x_cross;
            
            cout<<" threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
        }
        //------------------------------
        
        
        inverted_threshold = 255 - threshold;
        
        // Find a non-zero pixel that is the closest to the image center
        found_r = -1; found_c = -1; // indexes are set to an error value.
        {
            bdGeometry g;
            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
            unsigned int center_r = input_image.GetNumberOfRows() / 2;
            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
            int rn,cn;
            
            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
            {
                if(inverted_image(rn,cn)>inverted_threshold)
                {
                    found_r = rn;
                    found_c = cn;
                    break;
                }
            }
        }
        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
    }
    
    
    //----- Perform extraction -----
    {
        //bdBIP bip;
        //bdImage extracted_image;
        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        
        
        // Check if we need to perform masking by looking if the segmented region reaches the image edge
        int is_segmented_lesion_touching_image_edge = 0;
        //---- Check image edges -----
        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
        {
            if(extracted_image(r,0)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
        }
        if(!is_segmented_lesion_touching_image_edge)
        {
            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
            {
                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
            }
        }
        //---------------------------
        
        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
        
        // If the lesion is not touching image edge, we consider it a correct segmentation.
        if(is_segmented_lesion_touching_image_edge)
        {
            bdMelanoma melanoma;
            bdImage downsampled_image;
            melanoma.Downsample(extracted_image, downsampled_image, 50);
            
            
            bdImage extracted_downsampled_image;
            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
            
            // Dilate
            bdImage dilated_extracted_downsampled_image;
            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
            dilated_extracted_downsampled_image.FillInWith(0);
            bdGeometry g2;
            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
                {
                    if(extracted_downsampled_image(r,c)!=0)
                    {
                        int rn,cn;
                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                        {
                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                        }
                    }
                }
            }
            
            bdImage masked_image;
            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
            
            output.CopyFrom(masked_image);
        }
        else
        {
            output.CopyFrom(extracted_image);
            return 1;
        }
        
    }
    //------------------------------
    
    
    //output.CopyFrom(slice);
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}



//int bdMelanoma::SegmentMelanoma_Thresholding9(bdImage &input_image, bdImage &output)
//{
//    //cout<<" 1 ";
//    
//    if(input_image.IsEmpty()) return 0;
//    
//    
//    bdGIP gip;
//    
//    
//    // Rescale the slice values to fall within 0,255.
//    unsigned int lower_value = 0, upper_value = 255;
//    bdImage input;
//    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
//    
//    
//    
//    if(this->IsAbortRequested()) return 0;
//    this->SetProgressCounterRelativeValue(5);
//    
//    
//    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
//    bdImage slice;
//    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
//    slice.SetVisualizationPropertiesToMatchInput(input);
//    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
//        {
//            // Use just the blue channel.
//            slice(0,0,r,c) = input(2,0,r,c);
//            
//            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
//            
//            // combined B and G channels
//            //            if(input(1,0,r,c)>input(2,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//            //            else slice(0,0,r,c) = input(2,0,r,c);
//            
//            // combined R and G channels
//            //            if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
//            //            else slice(0,0,r,c) = input(0,0,r,c);
//            
//            
//        }
//    }
//    
//    //Calculate the histogram of the slice
//    bdArray<int> histogram_of_slice;
//    this->Histogram(slice, histogram_of_slice, 0,255, 0, 0);
//    
//    
//    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
//    unsigned int last_non_zero_value = 1;
//    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
//    {
//        if(histogram_of_slice[i]==0)
//        {
//            histogram_of_slice[i] = last_non_zero_value;
//        }
//        else
//        {
//            last_non_zero_value = histogram_of_slice[i];
//        }
//    }
//    
//    
//    
//    
//    
//    
//    // Calculate maxima in the histogram of the chosen slice.
//    bdArray<int> maxima;
//    this->HistogramMaxima(histogram_of_slice, maxima);
//    
//    
//    
//    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
//    // Find the max value...
//    unsigned int max_value1 = 0;
//    unsigned int index_of_max_value1 = 0;
//    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
//    {
//        if(maxima[i]>max_value1)
//        {
//            max_value1 = maxima[i];
//            index_of_max_value1 = i;
//        }
//    }
//    
//    // Find the next max value...
//    unsigned int max_value2 = 0;
//    unsigned int index_of_max_value2 = 0;
//    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
//    {
//        if(maxima[i]>max_value2)
//        {
//            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
//            //if(i!=index_of_max_value1)
//            //{
//            if(i>index_of_max_value1)
//            {
//                if(i-index_of_max_value1>10)
//                {
//                    max_value2 = maxima[i];
//                    index_of_max_value2 = i;
//                }
//            }
//            if(i<index_of_max_value1)
//            {
//                if(index_of_max_value1-i>10)
//                {
//                    max_value2 = maxima[i];
//                    index_of_max_value2 = i;
//                }
//            }
//            //}
//        }
//    }
//    
//    
//    unsigned int threshold;
////    //----- Determine threshold value 1 -----
////    {
////        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
////        // in array of maxima that preceed the maximum peak.
////        unsigned int max_difference = 0;
////        unsigned int index_of_max_difference = 0;
////        for(unsigned int i=index_of_max_value1-1; i>0; i--)
////        {
////            //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
////            int d = histogram_of_slice[i+1] - histogram_of_slice[i];
////            if(d<0) d = -d;
////            if(d>max_difference)
////            {
////                max_difference = d;
////                index_of_max_difference = i;
////            }
////        }
////        
////        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
////        // y = kx + n ;
////        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
////        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
////        
////        //        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
////        //        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
////        double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
////        double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
////        
////        
////        // For y=0: x = -n/k;
////        double x_cross = - n/k;
////        //cout<<"  x_cross="<<x_cross<<" ";
////        
////        
////        // The threshold value is the index in of the crossing with x axis
////        threshold = (unsigned int) x_cross;
////        
////        cout<<" threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
////    }
////    //------------------------------
////
//    
//    //----- Determine threshold value 1 -----
//    {
//        // Now that we have the maximum peak, search for the bending (foot of the curve) between the maxima in the histogram.
//
//        double line_pos1_src[3], line_pos2_src[3];
//        line_pos1_src[0] = 0; line_pos1_src[1] = histogram_of_slice[index_of_max_value1]; line_pos1_src[2] = index_of_max_value1;
//        line_pos2_src[0] = 0; line_pos2_src[1] = histogram_of_slice[index_of_max_value2]; line_pos2_src[2] = index_of_max_value2;
//        double max_squared_distance = 0;
//        unsigned int index_of_max_distance = index_of_max_value1-1;
//        for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)// for(unsigned int i=index_of_max_value1-1; i>0; i--)
//        {
//            double position_to_project_src[3], projected_position_src[3], n;
//            position_to_project_src[0] = 0; position_to_project_src[1] = histogram_of_slice[i]; position_to_project_src[2] = i;
//            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
//            if(n>0 && n<1)
//            {
//                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
//                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
//                
//                if(d>max_squared_distance)
//                {
//                    max_squared_distance = d;
//                    index_of_max_distance = i;
//                }
//            }
//        }
//        
//        // The threshold value is the index in of the crossing with x axis
//        threshold = index_of_max_distance;
//        
//        cout<<" threshold(index_of_max_distance)="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<") ";
//    }
//    //------------------------------
//    
//
//    
//    
//    
//    bdImage inverted_image;
//    gip.Negative(slice, inverted_image);
//    int inverted_threshold = 255 - threshold;
//    
//    bdImage thresholded_image;
//    bdBIP bip;
//    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
//    
//    
//    // Find a non-zero pixel that is the closest to the image center
//    // Create an array of seeds based on the image size and distance from the image center
//    bdArray< bdDiscreteCoordinates3D > seeds;
//    unsigned int search_distance_step = 4;
//    seeds.Set(1 + search_distance_step*8);
//    int center_r = thresholded_image.GetNumberOfRows() / 2;
//    int center_c = thresholded_image.GetNumberOfColumns() / 2;
//    unsigned int seed_index = 0;
//    seeds[seed_index].R() = center_r;
//    seeds[seed_index].C() = center_c;
//    seed_index++;
//    for(unsigned int d=1; d<=search_distance_step; d++)
//    {
//        int dR = (thresholded_image.GetNumberOfRows()*d) / 20;
//        int dC = (thresholded_image.GetNumberOfColumns()*d) / 20;
//        
//        seeds[seed_index].R() = center_r - dR;
//        seeds[seed_index].C() = center_c - dC;
//        seed_index++;
//        
//        seeds[seed_index].R() = center_r - dR;
//        seeds[seed_index].C() = center_c;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r - dR;
//        seeds[seed_index].C() = center_c + dC;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r;
//        seeds[seed_index].C() = center_c - dC;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r;
//        seeds[seed_index].C() = center_c + dC;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r + dR;
//        seeds[seed_index].C() = center_c - dC;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r + dR;
//        seeds[seed_index].C() = center_c;
//        seed_index++;
//
//        seeds[seed_index].R() = center_r + dR;
//        seeds[seed_index].C() = center_c + dC;
//        seed_index++;
//    }
//    
//    
//    // Find the region from seeds that has a large enough area to qualify for segmentation result.
//    bdImage extracted_image;
//    int is_region_found = 0;
//    for(seed_index = 0; seed_index<seeds.GetNumberOfElements(); seed_index++)
//    {
//        if(thresholded_image(seeds[seed_index].R(),seeds[seed_index].C())>0)
//        {
//            cout<<"  found_r,c=("<<seeds[seed_index].R()<<","<<seeds[seed_index].C()<<") ";
//            
//            //----- Perform extraction -----
//            bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
//            int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
//            cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
//            if(area_of_extracted_region > 1000)
//            {
//                is_region_found = 1;
//                break;
//            }
//        }
//    }
//    
//    
//    
////    // Find a non-zero pixel that is the closest to the image center
////    int found_r = -1, found_c = -1; // indexes are set to an error value.
////    {
////        bdGeometry g;
////        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
////        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
////        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
////        int rn,cn;
////        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
////        {
////            //if(inverted_image(rn,cn)>inverted_threshold)
////            if(thresholded_image(rn,cn)>0)
////            {
////                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
////                found_r = rn;
////                found_c = cn;
////                int rnn, cnn;
////                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
////                {
////                    thresholded_image(rnn,cnn) = 255;
////                }
////                
////                break;
////            }
////        }
////    }
////    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
//    
////    //----- Perform extraction -----
////    //bdImage extracted_image;
////    
////    // Check if a non-zero pixel was found. If so, perform extraction
////    if(found_r>=0 && found_c>=0 )
////    {
////        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
////    }
//    
//    
//    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
//    
//    
//    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
////    int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
////    cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
////    //if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
////    if(found_r<0 || found_c<0 || area_of_extracted_region<1000)
//    
//    if(!is_region_found)
//    {
//        //----- Determine threshold value 2 -----
//        {
//            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
//            // in array of maxima that succeed the maximum peak.
//            unsigned int max_difference = 0;
//            unsigned int index_of_max_difference = 0;
//            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
//            {
//                //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
//                int d = histogram_of_slice[i+1] - histogram_of_slice[i];
//                
//                if(d<0) d = -d;
//                if(d>max_difference)
//                {
//                    max_difference = d;
//                    index_of_max_difference = i;
//                }
//            }
//            
//            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
//            // y = kx + n ;
//            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
//            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
//            
//            //            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
//            //            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
//            double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
//            double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
//            
//            
//            // For y=0: x = -n/k;
//            double x_cross = - n/k;
//            //cout<<"  x_cross="<<x_cross<<" ";
//            
//            
//            // The threshold value is the index in of the crossing with x axis
//            threshold = (unsigned int) x_cross;
//            
//            cout<<" threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
//        }
//        //------------------------------
//        
//        
//        inverted_threshold = 255 - threshold;
//        
////        // Find a non-zero pixel that is the closest to the image center
////        found_r = -1; found_c = -1; // indexes are set to an error value.
////        {
////            bdGeometry g;
////            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
////            unsigned int center_r = input_image.GetNumberOfRows() / 2;
////            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
////            int rn,cn;
////            
////            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
////            {
////                if(inverted_image(rn,cn)>inverted_threshold)
////                {
////                    found_r = rn;
////                    found_c = cn;
////                    break;
////                }
////            }
////        }
////        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
//        
//        
//        // Find the region from seeds that has a large enough area to qualify for segmentation result.
//        for(seed_index = 0; seed_index<seeds.GetNumberOfElements(); seed_index++)
//        {
//            if(inverted_image(seeds[seed_index].R(),seeds[seed_index].C())>inverted_threshold)
//            {
//                cout<<"  found_r,c=("<<seeds[seed_index].R()<<","<<seeds[seed_index].C()<<") ";
//                
//                //----- Perform extraction -----
//                //bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
//                bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
//                int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
//                cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
//                if(area_of_extracted_region > 1000)
//                {
//                    is_region_found = 1;
//                    break;
//                }
//            }
//        }
//    }
//    
//    
//    if(!is_region_found)
//    {
//        output.CopyFrom(slice);
//        return 0;
//    }
//    
//    
//    //----- Perform extraction -----
//    {
////        //bdBIP bip;
////        //bdImage extracted_image;
////        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
////        
//        
//        // Check if we need to perform masking by looking if the segmented region reaches the image edge
//        int is_segmented_lesion_touching_image_edge = 0;
//        //---- Check image edges -----
//        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
//        {
//            if(extracted_image(r,0)!=0)
//            {
//                is_segmented_lesion_touching_image_edge = 1;
//                break;
//            }
//            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
//            {
//                is_segmented_lesion_touching_image_edge = 1;
//                break;
//            }
//        }
//        if(!is_segmented_lesion_touching_image_edge)
//        {
//            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
//            {
//                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
//                {
//                    is_segmented_lesion_touching_image_edge = 1;
//                    break;
//                }
//                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
//                {
//                    is_segmented_lesion_touching_image_edge = 1;
//                    break;
//                }
//            }
//        }
//        //---------------------------
//        
//        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
//        
//        // If the lesion is not touching image edge, we consider it a correct segmentation.
//        if(is_segmented_lesion_touching_image_edge)
//        {
//            bdMelanoma melanoma;
//            bdImage downsampled_image;
//            melanoma.Downsample(extracted_image, downsampled_image, 50);
//            
//            
//            bdImage extracted_downsampled_image;
//            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
//            
//            // Dilate
//            bdImage dilated_extracted_downsampled_image;
//            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
//            dilated_extracted_downsampled_image.FillInWith(0);
//            bdGeometry g2;
//            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
//            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
//            {
//                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
//                {
//                    if(extracted_downsampled_image(r,c)!=0)
//                    {
//                        int rn,cn;
//                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
//                        {
//                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
//                        }
//                    }
//                }
//            }
//            
//            bdImage masked_image;
//            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
//            
//            output.CopyFrom(masked_image);
//        }
//        else
//        {
//            output.CopyFrom(extracted_image);
//            return 1;
//        }
//        
//    }
//    //------------------------------
//    
//    
//    //output.CopyFrom(slice);
//    
//    
//    if(this->IsAbortRequested()) return 0;
//    this->SetProgressCounterRelativeValue(100);
//    
//    return 1;
//}
//
//



int bdMelanoma::SegmentMelanoma_Thresholding9(bdImage &input_image, bdImage &output)
{
    //cout<<" 1 ";
    
    if(input_image.IsEmpty()) return 0;
    
    
    bdGIP gip;
    
    
    // Rescale the slice values to fall within 0,255.
    unsigned int lower_value = 0, upper_value = 255;
    bdImage input;
    gip.RescaleWholeRange(input_image, lower_value, upper_value, input);
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(5);
    
    
    // Create a 'combined' slice from the RGB components by looking at max values of G and B channels (to erase markers)
    bdImage slice;
    slice.SetSize(1,1, input.GetNumberOfRows(), input.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(input);
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            // Use just the blue channel.
            slice(0,0,r,c) = input(2,0,r,c);
            
            //slice(0,0,r,c) = (input(0,0,r,c) + input(1,0,r,c) + input(2,0,r,c)) / 3;
            
            // combined B and G channels
            //            if(input(1,0,r,c)>input(2,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
            //            else slice(0,0,r,c) = input(2,0,r,c);
            
            // combined R and G channels
            //            if(input(1,0,r,c)>input(0,0,r,c)) slice(0,0,r,c) = input(1,0,r,c);
            //            else slice(0,0,r,c) = input(0,0,r,c);
            
            
        }
    }
    
    //Calculate the histogram of the slice
    bdArray<int> histogram_of_slice;
    this->Histogram(slice, histogram_of_slice, 0,255, 0, 0);
    
    
    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
    unsigned int last_non_zero_value = 1;
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        if(histogram_of_slice[i]==0)
        {
            histogram_of_slice[i] = last_non_zero_value;
        }
        else
        {
            last_non_zero_value = histogram_of_slice[i];
        }
    }
    
    
    
    
    
    
    // Calculate maxima in the histogram of the chosen slice.
    bdArray<int> maxima;
    this->HistogramMaxima(histogram_of_slice, maxima);
    
    
    
    // Go through the maxima array and histogram of chosen slice to determine the 2 peaks to use for threshold calculation.
    // Find the max value...
    unsigned int max_value1 = 0;
    unsigned int index_of_max_value1 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value1)
        {
            max_value1 = maxima[i];
            index_of_max_value1 = i;
        }
    }
    
    // Find the next max value...
    unsigned int max_value2 = 0;
    unsigned int index_of_max_value2 = 0;
    for(unsigned int i=0; i<maxima.GetNumberOfElements(); i++)
    {
        if(maxima[i]>max_value2)
        {
            // max_value2 must be different than max_value1 n terms of index and given distance of indexes!
            //if(i!=index_of_max_value1)
            //{
            if(i>index_of_max_value1)
            {
                if(i-index_of_max_value1>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            if(i<index_of_max_value1)
            {
                if(index_of_max_value1-i>10)
                {
                    max_value2 = maxima[i];
                    index_of_max_value2 = i;
                }
            }
            //}
        }
    }
    
    
    unsigned int threshold;
    //    //----- Determine threshold value 1 -----
    //    {
    //        // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
    //        // in array of maxima that preceed the maximum peak.
    //        unsigned int max_difference = 0;
    //        unsigned int index_of_max_difference = 0;
    //        for(unsigned int i=index_of_max_value1-1; i>0; i--)
    //        {
    //            //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
    //            int d = histogram_of_slice[i+1] - histogram_of_slice[i];
    //            if(d<0) d = -d;
    //            if(d>max_difference)
    //            {
    //                max_difference = d;
    //                index_of_max_difference = i;
    //            }
    //        }
    //
    //        // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
    //        // y = kx + n ;
    //        // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
    //        // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
    //
    //        //        double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
    //        //        double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
    //        double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
    //        double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
    //
    //
    //        // For y=0: x = -n/k;
    //        double x_cross = - n/k;
    //        //cout<<"  x_cross="<<x_cross<<" ";
    //
    //
    //        // The threshold value is the index in of the crossing with x axis
    //        threshold = (unsigned int) x_cross;
    //
    //        cout<<" threshold="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
    //    }
    //    //------------------------------
    //
    
    //----- Determine threshold value 1 -----
    {
        // Now that we have the maximum peak, search for the bending (foot of the curve) between the maxima in the histogram.
        
        double line_pos1_src[3], line_pos2_src[3];
        line_pos1_src[0] = 0; line_pos1_src[1] = histogram_of_slice[index_of_max_value1]; line_pos1_src[2] = index_of_max_value1;
        line_pos2_src[0] = 0; line_pos2_src[1] = histogram_of_slice[index_of_max_value2]; line_pos2_src[2] = index_of_max_value2;
        double max_squared_distance = 0;
        unsigned int index_of_max_distance = index_of_max_value1-1;
        for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)// for(unsigned int i=index_of_max_value1-1; i>0; i--)
        {
            double position_to_project_src[3], projected_position_src[3], n;
            position_to_project_src[0] = 0; position_to_project_src[1] = histogram_of_slice[i]; position_to_project_src[2] = i;
            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
            if(n>0 && n<1)
            {
                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
                
                if(d>max_squared_distance)
                {
                    max_squared_distance = d;
                    index_of_max_distance = i;
                }
            }
        }
        
        // The threshold value is the index in of the crossing with x axis
        threshold = index_of_max_distance;
        
        cout<<" threshold(index_of_max_distance)="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<") ";
    }
    //------------------------------
    
    
    
    
    
    bdImage inverted_image;
    gip.Negative(slice, inverted_image);
    int inverted_threshold = 255 - threshold;
    
    bdImage thresholded_image;
    bdBIP bip;
    bip.Threshold(inverted_image, inverted_threshold, thresholded_image);
    
    
    // Find a non-zero pixel that is the closest to the image center
    // Create an array of seeds based on the image size and distance from the image center
    bdArray< bdDiscreteCoordinates3D > seeds;
    unsigned int search_distance_step = 4;
    seeds.Set(1 + search_distance_step*8);
    int center_r = thresholded_image.GetNumberOfRows() / 2;
    int center_c = thresholded_image.GetNumberOfColumns() / 2;
    unsigned int seed_index = 0;
    seeds[seed_index].R() = center_r;
    seeds[seed_index].C() = center_c;
    seed_index++;
    for(unsigned int d=1; d<=search_distance_step; d++)
    {
        int dR = (thresholded_image.GetNumberOfRows()*d) / 20;
        int dC = (thresholded_image.GetNumberOfColumns()*d) / 20;
        
        seeds[seed_index].R() = center_r - dR;
        seeds[seed_index].C() = center_c - dC;
        seed_index++;
        
        seeds[seed_index].R() = center_r - dR;
        seeds[seed_index].C() = center_c;
        seed_index++;
        
        seeds[seed_index].R() = center_r - dR;
        seeds[seed_index].C() = center_c + dC;
        seed_index++;
        
        seeds[seed_index].R() = center_r;
        seeds[seed_index].C() = center_c - dC;
        seed_index++;
        
        seeds[seed_index].R() = center_r;
        seeds[seed_index].C() = center_c + dC;
        seed_index++;
        
        seeds[seed_index].R() = center_r + dR;
        seeds[seed_index].C() = center_c - dC;
        seed_index++;
        
        seeds[seed_index].R() = center_r + dR;
        seeds[seed_index].C() = center_c;
        seed_index++;
        
        seeds[seed_index].R() = center_r + dR;
        seeds[seed_index].C() = center_c + dC;
        seed_index++;
    }
    
    
    // Find the region from seeds that has a large enough area to qualify for segmentation result.
    bdImage extracted_image;
    int is_region_found = 0;
    for(seed_index = 0; seed_index<seeds.GetNumberOfElements(); seed_index++)
    {
        if(thresholded_image(seeds[seed_index].R(),seeds[seed_index].C())>0)
        {
            cout<<"  found_r,c=("<<seeds[seed_index].R()<<","<<seeds[seed_index].C()<<") ";
            
            //----- Perform extraction -----
            bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
            int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
            cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
            if(area_of_extracted_region > 1000)
            {
                is_region_found = 1;
                break;
            }
        }
    }
    
    
    
    //    // Find a non-zero pixel that is the closest to the image center
    //    int found_r = -1, found_c = -1; // indexes are set to an error value.
    //    {
    //        bdGeometry g;
    //        g.SetDimensions(1,thresholded_image.GetNumberOfRows(),thresholded_image.GetNumberOfColumns());
    //        unsigned int center_r = thresholded_image.GetNumberOfRows() / 2;
    //        unsigned int center_c = thresholded_image.GetNumberOfColumns() / 2;
    //        int rn,cn;
    //        for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
    //        {
    //            //if(inverted_image(rn,cn)>inverted_threshold)
    //            if(thresholded_image(rn,cn)>0)
    //            {
    //                //cout<<" rn,cn="<<rn<<","<<cn<<"  ";
    //                found_r = rn;
    //                found_c = cn;
    //                int rnn, cnn;
    //                for(g.ForCoordinates_Circle(rn, cn, 350); g.Get_Circle(400, rnn, cnn); )
    //                {
    //                    thresholded_image(rnn,cnn) = 255;
    //                }
    //
    //                break;
    //            }
    //        }
    //    }
    //    cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
    
    //    //----- Perform extraction -----
    //    //bdImage extracted_image;
    //
    //    // Check if a non-zero pixel was found. If so, perform extraction
    //    if(found_r>=0 && found_c>=0 )
    //    {
    //        bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, found_r, found_c);
    //    }
    
    
    // Check if we have found a valid reagion by examining its size (should be more than a certain threshold)
    
    
    // Check if a non-zero pixel was found and if the extracted region is large enough. If not, we need to calculate a new threshold value.
    //    int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
    //    cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
    //    //if(found_r<0 || found_c<0 || bip.AreaOfObject_3D(extracted_image)<10000)
    //    if(found_r<0 || found_c<0 || area_of_extracted_region<1000)
    
    if(!is_region_found)
    {
        //----- Determine threshold value 2 -----
        {
            // Now that we have the maximum peak, search for the highest difference between the 2 neighboring values
            // in array of maxima that succeed the maximum peak.
            unsigned int max_difference = 0;
            unsigned int index_of_max_difference = 0;
            for(unsigned int i=index_of_max_value1; i<index_of_max_value2; i++)
            {
                //int d = array_of_histograms[max_variance_index][i+1] - array_of_histograms[max_variance_index][i];
                int d = histogram_of_slice[i+1] - histogram_of_slice[i];
                
                if(d<0) d = -d;
                if(d>max_difference)
                {
                    max_difference = d;
                    index_of_max_difference = i;
                }
            }
            
            // Calculate where the line (defined by the 2 neighbors with max difference value) would cross the x (index) axis.
            // y = kx + n ;
            // k = dy/dx -> k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference]
            // n = y - kx -> n = array_of_histograms[max_variance_index][index_of_max_difference] - k*index_of_max_difference
            
            //            double k = array_of_histograms[max_variance_index][index_of_max_difference+1] - array_of_histograms[max_variance_index][index_of_max_difference];
            //            double n = ((double)array_of_histograms[max_variance_index][index_of_max_difference]) - k * ((double)index_of_max_difference);
            double k = histogram_of_slice[index_of_max_difference+1] - histogram_of_slice[index_of_max_difference];
            double n = ((double)histogram_of_slice[index_of_max_difference]) - k * ((double)index_of_max_difference);
            
            
            // For y=0: x = -n/k;
            double x_cross = - n/k;
            //cout<<"  x_cross="<<x_cross<<" ";
            
            
            // The threshold value is the index in of the crossing with x axis
            threshold = (unsigned int) x_cross;
            
            cout<<" threshold2="<<threshold<<"(index_of_max="<<index_of_max_value1<<",index_of_max2="<<index_of_max_value2<<",index_of_max_diff="<<index_of_max_difference<<") ";
        }
        //------------------------------
        
        
        inverted_threshold = 255 - threshold;
        
        //        // Find a non-zero pixel that is the closest to the image center
        //        found_r = -1; found_c = -1; // indexes are set to an error value.
        //        {
        //            bdGeometry g;
        //            g.SetDimensions(1,input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
        //            unsigned int center_r = input_image.GetNumberOfRows() / 2;
        //            unsigned int center_c = input_image.GetNumberOfColumns() / 2;
        //            int rn,cn;
        //
        //            for(g.ForCoordinates_Circle(center_r, center_c, 0); g.Get_Circle(400, rn, cn); )
        //            {
        //                if(inverted_image(rn,cn)>inverted_threshold)
        //                {
        //                    found_r = rn;
        //                    found_c = cn;
        //                    break;
        //                }
        //            }
        //        }
        //        cout<<"  found_r,c=("<<found_r<<","<<found_c<<") ";
        
        
        // Find the region from seeds that has a large enough area to qualify for segmentation result.
        for(seed_index = 0; seed_index<seeds.GetNumberOfElements(); seed_index++)
        {
            if(inverted_image(seeds[seed_index].R(),seeds[seed_index].C())>inverted_threshold)
            {
                cout<<"  found_r,c=("<<seeds[seed_index].R()<<","<<seeds[seed_index].C()<<") ";
                
                //----- Perform extraction -----
                //bip.Extract_26_ConnectedComponent_3D(thresholded_image, 0, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
                bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, seeds[seed_index].R(),seeds[seed_index].C());
                int area_of_extracted_region = bip.AreaOfObject_3D(extracted_image);
                cout<<"  area_of_extracted_region="<<area_of_extracted_region<<"  ";
                if(area_of_extracted_region > 1000)
                {
                    is_region_found = 1;
                    break;
                }
            }
        }
    }
    
    
    if(!is_region_found)
    {
        output.CopyFrom(slice);
        return 0;
    }
    
    
    //----- Perform extraction -----
    {
        //        //bdBIP bip;
        //        //bdImage extracted_image;
        //        bip.Extract_26_ConnectedComponent_3D(inverted_image, inverted_threshold, extracted_image, 0, found_r, found_c);
        //
        
        // Check if we need to perform masking by looking if the segmented region reaches the image edge
        int is_segmented_lesion_touching_image_edge = 0;
        //---- Check image edges -----
        for(unsigned int r=0; r<extracted_image.GetNumberOfRows(); r++)
        {
            if(extracted_image(r,0)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
            if(extracted_image(r,extracted_image.GetNumberOfColumns()-1)!=0)
            {
                is_segmented_lesion_touching_image_edge = 1;
                break;
            }
        }
        if(!is_segmented_lesion_touching_image_edge)
        {
            for(unsigned int c=0; c<extracted_image.GetNumberOfColumns(); c++)
            {
                if(extracted_image(extracted_image.GetNumberOfRows()-1,c)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
                if(extracted_image(extracted_image.GetNumberOfRows()-1,0)!=0)
                {
                    is_segmented_lesion_touching_image_edge = 1;
                    break;
                }
            }
        }
        //---------------------------
        
        cout<<"  is_segmented_lesion_touching_image_edge="<<is_segmented_lesion_touching_image_edge<<" ";
        
        // If the lesion is not touching image edge, we consider it a correct segmentation.
        if(is_segmented_lesion_touching_image_edge)
        {
            bdMelanoma melanoma;
            bdImage downsampled_image;
            melanoma.Downsample(extracted_image, downsampled_image, 50);
            
            
            bdImage extracted_downsampled_image;
            bip.ExtractLargest_26_ConnectedComponent(downsampled_image, 0, extracted_downsampled_image);
            
            // Dilate
            bdImage dilated_extracted_downsampled_image;
            dilated_extracted_downsampled_image.SetSizeAndPropertiesAs(extracted_downsampled_image);
            dilated_extracted_downsampled_image.FillInWith(0);
            bdGeometry g2;
            g2.SetDimensions(1,extracted_downsampled_image.GetNumberOfRows(),extracted_downsampled_image.GetNumberOfColumns());
            for(unsigned int r=0; r<extracted_downsampled_image.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<extracted_downsampled_image.GetNumberOfColumns(); c++)
                {
                    if(extracted_downsampled_image(r,c)!=0)
                    {
                        int rn,cn;
                        for(g2.ForCoordinates_4_Neighbors(r, c); g2.Get_4_Neighbors(rn, cn); )
                        {
                            dilated_extracted_downsampled_image(rn,cn) = extracted_downsampled_image(r,c);
                        }
                    }
                }
            }
            
            bdImage masked_image;
            melanoma.MaskWithDownsampledImage(extracted_image, dilated_extracted_downsampled_image, masked_image, 50);
            
            output.CopyFrom(masked_image);
        }
        else
        {
            output.CopyFrom(extracted_image);
            //return 1;
        }
        
    }
    //------------------------------
    
    
    //output.CopyFrom(slice);
    
    //----- Make a RGB image by thresholding also Red and Green channels -----
    bdImage temp;
    temp.CopyFrom(output);
    output.SetSize(3,1,temp.GetNumberOfRows(),temp.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(temp);
    output.FillInWith(0);
    for(unsigned int r=0; r<temp.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<temp.GetNumberOfColumns(); c++)
        {
            if(temp(r,c)!=0)
            {
                if(input_image(0,0,r,c)<threshold) { output(0,0,r,c) = 255; }
                if(input_image(1,0,r,c)<threshold) { output(1,0,r,c) = 255; }
                if(input_image(2,0,r,c)<threshold) { output(2,0,r,c) = 255; }
                
//                if(input_image(0,0,r,c)<threshold) { output(0,0,r,c) = 255; output(1,0,r,c) = 0; output(2,0,r,c) = 0; }
//                else
//                {
//                    if(input_image(1,0,r,c)<threshold) { output(0,0,r,c) = 0; output(1,0,r,c) = 255; output(2,0,r,c) = 0; }
//                    else
//                    {
//                        output(0,0,r,c) = 0; output(1,0,r,c) = 0; output(2,0,r,c) = 255;
//                    }
//                }
            }
        }
    }
    //----------
    
    
    
    
    if(this->IsAbortRequested()) return 0;
    this->SetProgressCounterRelativeValue(100);
    
    return 1;
}








// // this version does not work
//int bdMelanoma::SegmentMelanomaInnerStructure(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
//{
//    unsigned int s = 0;
//
//    output.SetSizeAndPropertiesAs(original);
//    output.FillInWith(0);
//
//    
//    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
//    {
//        bdRegion2D segmented_region;
//        segmented_region.CreateRegionFromImage(segmented_melanoma,s,0);
//        
//        for(unsigned int threshold = 180; threshold>50; threshold-=30)
//        {
//            cout<<"thresold="<<threshold<<endl;
//            
//            bdRegion2D thresholded_region;
//            
//            bdRegion2DIterator itr;
//            for(itr.SetBegin(&segmented_region); itr.IsValid(); itr.MoveToNext())
//            {
//                if(original(t,s,itr.GetIndexRow(),itr.GetIndexColumn())>threshold)
//                {
//                    thresholded_region.AddPoint(itr.GetIndexRow(),itr.GetIndexColumn());
//                }
//            }
//            
//            bdList<bdRegion2D> list_of_regions;
//            thresholded_region.BreakIntoDisjointRegionsOf_8_Neighborhood(list_of_regions);
//            
//            bdListIterator<bdRegion2D> itl;
//            for(itl.SetLeftEnd(list_of_regions); itl.IsValid(); itl.MoveRight())
//            {
//                if(itl.GetElementPointer()->GetNumberOfElements() < segmented_region.GetNumberOfElements()/4)
//                {
//                    for(itr.SetBegin(itl.GetElementPointer()); itr.IsValid(); itr.MoveToNext())
//                    {
//                        if(output(t,s,itr.GetIndexRow(),itr.GetIndexColumn())==0)
//                        {
//                            output(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = threshold;
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    return 1;
//}



//int bdMelanoma::SegmentMelanomaInnerStructure(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
//{
//    unsigned int s = 0;
//    
//    output.SetSizeAndPropertiesAs(original);
//    output.FillInWith(0);
//    
//    unsigned int t = 0;
////    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
//    {
//        cout<<" Creating list ";
//        // put the masked voxels from the original imgae to the list
//        bdList< bdVoxel > list_of_voxels;
//        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//            {
//                if(segmented_melanoma(t,s,r,c)!=0)
//                {
//                    bdVoxel *voxel = list_of_voxels.AddNewToRightEnd();
//                    (*voxel)(r,c);
//                    voxel->SetValue(original(t,s,r,c));
//                }
//            }
//        }
//        
//        cout<<" Sorting list ";
//        
//        // bubble sort the voxels in the list in ascending order of their values
//        bdListNode< bdVoxel > *ending_node = list_of_voxels.GetRightEndNodePointer();
//        for(int is_change_made = 1; is_change_made; )
//        {
//            is_change_made = 0;
//            
//            bdListNode< bdVoxel > *node1, *node2;
//            node1 = list_of_voxels.GetLeftEndNodePointer();
//            node2 = node1->GetRight();
//            for( ; node1!=ending_node && ending_node; node1=node2, node2=node2->GetRight())
//            {
//                if(node1->GetElementPointer()->V()>node2->GetElementPointer()->V())
//                {
//                    bdVoxel v;
//                    v.CopyFrom(node2->GetElement());
//                    node2->GetElementPointer()->CopyFrom(node1->GetElement());
//                    node1->GetElementPointer()->CopyFrom(v);
//                    is_change_made = 1;
//                }
//            }
//            ending_node = ending_node->GetLeft();
//        }
//        
//        cout<<" Segmenting ";
//        
//        // Go through the sorted list and for a pixel that has no neighbors in the output image give a new
//        // label value. For a pixel that already has a neighbor with a label value give that label value.
//        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
//        unsigned int label_value = 1;
//        bdGeometry g;
//        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
//        bdListIterator< bdVoxel > it;
//        for(it.SetLeftEnd(list_of_voxels); it.IsValid(); it.MoveRight())
//        {
//            unsigned int r = it.GetElementPointer()->R();
//            unsigned int c = it.GetElementPointer()->C();
//            
//            unsigned int found_label_value = 0;
//            
//            int rn, cn;
//            for(g.ForCoordinates_9_Neighbors(r,c); g.Get_9_Neighbors(rn,cn); )
//            {
//                if(output(t,s,rn,cn)!=0)
//                {
//                    if(!found_label_value)
//                    {
//                        found_label_value = output(t,s,rn,cn);
//                        output(t,s,r,c) = found_label_value;
//                    }
//                    else
//                    {
//                        if(found_label_value!=output(t,s,rn,cn))
//                        {
//                            output(t,s,r,c) = 0;
//                            break;
//                        }
//                    }
//                }
//            }
//            
//            // if no label was found in neighborhood of the pixel, assign a new label value
//            if(!found_label_value)
//            {
//                output(t,s,r,c) = label_value;
//                label_value++;
//            }
//        }
//    }
//    
//    return 1;
//}





int bdMelanoma::SegmentMelanomaInnerStructure(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    output.SetSizeAndPropertiesAs(original);
    output.FillInWith(0);
    
    unsigned int t = 0;
    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
    {
        bdArray<int> histogram;
        unsigned int range_min = 0;
        unsigned int range_max = 65535;
        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
        
        // create array of sorted voxel values. This will be used to access voxels in order of their values.
        bdArray< bdArray< bdVoxel > > sorted_array;
        sorted_array.Set(histogram.GetNumberOfElements());
        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
        {
            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
        }
        
        cout<<" Creating list ";
        // put the masked voxels from the original image to the sorted_array
        bdArray<int> index_array;
        index_array.Set(sorted_array.GetNumberOfElements());
        index_array.FillInWith(0);
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(t,s,r,c)!=0)
                {
                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
                    (index_array[original(t,s,r,c)])++;
                }
            }
        }
        
        cout<<" Segmenting ";
        
        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
        // label value. For a pixel that already has a neighbor with a label value give that label value.
        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
        unsigned int label_value = 1;
        bdGeometry g;
        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
        
        for(unsigned int i=0; i<sorted_array.GetNumberOfElements(); i++)
        {
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                
                unsigned int found_label_value = 0;
                
                int rn, cn;
                for(g.ForCoordinates_9_Neighbors(r,c); g.Get_9_Neighbors(rn,cn); )
                {
                    if(output(t,s,rn,cn)!=0)
                    {
                        if(!found_label_value)
                        {
                            found_label_value = output(t,s,rn,cn);
                            output(t,s,r,c) = found_label_value;
                        }
                        else
                        {
//                            if(found_label_value!=output(t,s,rn,cn))
//                            {
//                                output(t,s,r,c) = 0;
//                                break;
//                            }
                            if(found_label_value!=output(t,s,rn,cn))
                            {
                                if(found_label_value<output(t,s,rn,cn))
                                {
                                    output(t,s,r,c) = found_label_value;
                                }
                                else
                                {
                                    output(t,s,r,c) = output(t,s,rn,cn);
                                }
                                break;
                            }
                            
                        }
                    }
                }
                
                // if no label was found in neighborhood of the pixel, assign a new label value
                if(!found_label_value)
                {
                    output(t,s,r,c) = label_value;
                    label_value++;
                }
                
            }
        }
    }
    
    return 1;
}



int bdMelanoma::SegmentMelanomaInnerStructure2(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    output.SetSizeAndPropertiesAs(original);
    output.FillInWith(0);
    
    unsigned int t = 0;
    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
    {
        bdArray<int> histogram;
        unsigned int range_min = 0;
        unsigned int range_max = 65535;
        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
        
        // create array of sorted voxel values. This will be used to access voxels in order of their values.
        bdArray< bdArray< bdVoxel > > sorted_array;
        sorted_array.Set(histogram.GetNumberOfElements());
        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
        {
            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
        }
        
        cout<<" Creating list ";
        // put the masked voxels from the original image to the sorted_array
        bdArray<int> index_array;
        index_array.Set(sorted_array.GetNumberOfElements());
        index_array.FillInWith(0);
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(t,s,r,c)!=0)
                {
                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
                    (index_array[original(t,s,r,c)])++;
                }
            }
        }
        
        cout<<" Segmenting ";
        
        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
        // label value. For a pixel that already has a neighbor with a label value give that label value.
        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
        bdGeometry g;
        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
        
        for(unsigned int i=0; i<sorted_array.GetNumberOfElements(); i++)// 'i' is the threshold value
        {
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                
                //unsigned int found_label_value = 0;
                
                int rn, cn;
                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
                {
                    if(output(t,s,rn,cn)==0)//if not visited yet
                    {
                        output(t,s,r,c) = 255;
                        break;
                    }
                }
            }
        }
    }
    
    return 1;
}



int bdMelanoma::SegmentMelanomaInnerStructure3(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    output.SetSizeAndPropertiesAs(original);
    output.FillInWith(0);
    
    unsigned int t = 0;
    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
    {
        bdArray<int> histogram;
        unsigned int range_min = 0;
        unsigned int range_max = 65535;
        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
        
        // create array of sorted voxel values. This will be used to access voxels in order of their values.
        bdArray< bdArray< bdVoxel > > sorted_array;
        sorted_array.Set(histogram.GetNumberOfElements());
        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
        {
            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
        }
        
        cout<<" Creating list ";
        // put the masked voxels from the original image to the sorted_array
        bdArray<int> index_array;
        index_array.Set(sorted_array.GetNumberOfElements());
        index_array.FillInWith(0);
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(t,s,r,c)!=0)
                {
                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
                    (index_array[original(t,s,r,c)])++;
                }
            }
        }
        
        cout<<" Segmenting ";
        
        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
        // label value. For a pixel that already has a neighbor with a label value give that label value.
        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
        bdGeometry g;
        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
        
        bdImage thresholded_img;
        thresholded_img.SetSizeAndPropertiesAs(original);
        thresholded_img.FillInWith(0);
        
        int step = 25;
        for(unsigned int thr=40; thr<150; thr+=step)// 'i' is the threshold value
        {
            cout<<"thr="<<thr<<endl;
            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
            {
                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
                {
                    
                    unsigned int r = sorted_array[i][j].R();
                    unsigned int c = sorted_array[i][j].C();

                    // perform trhesholdeing
                    thresholded_img(t,s,r,c) = 255;
                    
//                    int rn, cn;
//                    for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
//                    {
//                        if(output(t,s,rn,cn)==0)//if not visited yet
//                        {
//                            output(t,s,r,c) = 255;
//                            break;
//                        }
//                    }
                }
            }
            
            bdImage closed_img;
            bdBIP bip;
            bip.Closing_Sphere(thresholded_img, closed_img, 64);
            
            // Go through the current images and record places where the closed image is non-zero and thresholded image is zero (but which hasn't yet been entered into the output)
            for(unsigned int r=0; r<thresholded_img.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<thresholded_img.GetNumberOfColumns(); c++)
                {
                    if(output(r,c)==0)
                    {
                        if(thresholded_img(r,c)==0 && closed_img(r,c)!=0)
                        {
                            output(r,c) = 255;
                        }
                    }
                }
            }
        }
    }
    
    return 1;
}



int bdMelanoma::SegmentMelanomaInnerStructure4(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    output.SetSizeAndPropertiesAs(original);
    output.FillInWith(0);
    
    unsigned short edge_label = 350;//65535;
    unsigned short thresholded_label = 255;
    
    unsigned int t = 0;
    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
    {
        bdArray<int> histogram;
        unsigned int range_min = 0;
        unsigned int range_max = 65535;
        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
        
        // create array of sorted voxel values. This will be used to access voxels in order of their values.
        bdArray< bdArray< bdVoxel > > sorted_array;
        sorted_array.Set(histogram.GetNumberOfElements());
        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
        {
            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
        }
        
        //cout<<" Creating list ";
        // put the masked voxels from the original image to the sorted_array
        bdArray<int> index_array;
        index_array.Set(sorted_array.GetNumberOfElements());
        index_array.FillInWith(0);
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(t,s,r,c)!=0)
                {
                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
                    (index_array[original(t,s,r,c)])++;
                }
            }
        }
        
        //cout<<" Segmenting ";
        
        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
        // label value. For a pixel that already has a neighbor with a label value give that label value.
        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
        bdGeometry g;
        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
        
        bdImage thresholded_img;
        thresholded_img.SetSizeAndPropertiesAs(original);
        thresholded_img.FillInWith(0);
        
        int step = 6;
        unsigned int thr=30+step;
        
        // Perform initial thresholding to fill in parts that do not fall in range of examined threshold values
        for(unsigned int i=0; i<thr-step; i++)// 'i' is the threshold value
        {
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                
                // perform trhesholdeing
                thresholded_img(t,s,r,c) = thresholded_label;
            }
        }
        
        //for(unsigned int thr=40+step; thr<100; thr=thr+step)// 'i' is the threshold value
        for( ; thr<150; thr=thr+step)
        {
            cout<<" thr="<<thr<<" ";
            
            // create the thresholded image
            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
            {
                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
                {
                    
                    unsigned int r = sorted_array[i][j].R();
                    unsigned int c = sorted_array[i][j].C();
                    
                    // perform trhesholdeing
                    thresholded_img(t,s,r,c) = thresholded_label;
                }
            }
            
            // find edges of the thresholded image, add pixels to list for fast access
            bdList< bdVoxel > list_of_edge_pixels;
            for(unsigned int i=thr-step; i<thr; i++)//for(unsigned int i=40; i<thr; i++)// 'i' is the threshold value
            {
                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
                {
                    unsigned int r = sorted_array[i][j].R();
                    unsigned int c = sorted_array[i][j].C();
                    int rn, cn;
                    for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
                    {
                        if(thresholded_img(t,s,rn,cn)==0)//if edge found
                        {
                            list_of_edge_pixels.AddToRightEnd(sorted_array[i][j]);
                            thresholded_img(t,s,r,c) = edge_label;
                            break;
                        }
                    }
                    
                }
            }
            // At this point: edge pixels are in the list_of_edge_pixels and are marked with value edge_label in the thresholded_img.
            
            // Filter the list of edge pixels to maintin only a single pixel per edge.
            bdList< bdVoxel > list_of_filtered_edge_pixels;
            {
                bdListIterator< bdVoxel > it;
                for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
                {
                    unsigned int r = it.GetElementPointer()->R();
                    unsigned int c = it.GetElementPointer()->C();
                    
                    if(thresholded_img(t,s,r,c) == edge_label)
                    {
                        list_of_filtered_edge_pixels.AddToRightEnd(it.GetElement());
                        
                        bdList< bdVoxel > temp_list;
                        temp_list.AddToRightEnd(it.GetElement());
                        thresholded_img(t,s,r,c) = thresholded_label;
                        while(!temp_list.IsEmpty())
                        {
                            bdVoxel v = temp_list.GetLeftEnd();
                            
                            int rnn, cnn;
                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
                            {
                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
                                {
                                    bdVoxel v_temp;
                                    v_temp(s,rnn,cnn);
                                    temp_list.AddToRightEnd(v_temp);
                                    thresholded_img(t,s,rnn,cnn) = thresholded_label;
                                }
                            }
                            
                            temp_list.DeleteLeftEnd();
                        }
                    }
                }
            }
            // At this point: list_of_filtered_edge_pixels contains list of pixels where each pixel represents one contour (taken from the contour).
            // Edge pixels are no longer marked in the thresholded_img.
            
            // Go through the list_of_filtered_edge_pixels and check if the pixel belongs to the outer edge by growing a region
            // and checking if it touches the border of the segmented region (pixel value 0 in segmented image).
            bdListIterator< bdVoxel > it;
            for(it.SetLeftEnd(list_of_filtered_edge_pixels); it.IsValid(); it.MoveRight())
            {
                //unsigned int r = it.GetElementPointer()->R();
                //unsigned int c = it.GetElementPointer()->C();

// OVDE STAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                // list of pixels belonging to the region
                bdList< bdVoxel > region_pixels_list;
                //region_pixels_list.AddToRightEnd(it.GetElement());
                
                bdList< bdVoxel > temp_list;
                temp_list.AddToRightEnd(it.GetElement());
                //thresholded_img(t,s,r,c) = 128;//thresholded_label;
                while(!temp_list.IsEmpty())
                {
                    bdVoxel v = temp_list.GetLeftEnd();
                    
                    int rnn, cnn;
                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
                    {
                        if(thresholded_img(t,s,rnn,cnn)==0)//if edge found
                        {
                            // if the segmented image is 0 it means that we have reached the end of the segmented melanoma, so the edge is outer...
                            if(segmented_melanoma(t,s,rnn,cnn)==0)
                            {
                                //... and since the edge is outer, we need to undo the region growing...
                                while(!region_pixels_list.IsEmpty())
                                {
                                    bdVoxel v2 = region_pixels_list.GetLeftEnd();
                                    thresholded_img(t,s,v2.R(),v2.C()) = 0;
                                    output(t,s,v2.R(),v2.C()) = 0;
                                    region_pixels_list.DeleteLeftEnd();
                                }
                                
                                //region_pixels_list.Reset();
                                temp_list.Reset();
                                break;
                            }
                            // if the segmented image is non-zero we consider this edge an inner egde (untill proven otherwise)...
                            else
                            {
                                //... so perform segmentation.
                                thresholded_img(t,s,rnn,cnn) = 128;
                                output(t,s,rnn,cnn) = 255;
                                
                                bdVoxel v_temp;
                                v_temp(s,rnn,cnn);
                                temp_list.AddToRightEnd(v_temp);
                                region_pixels_list.AddToRightEnd(v_temp);
                            }
                        }
                    }
                    
                    temp_list.DeleteLeftEnd();
                }
                
                
            }
            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
        }

        
            
//            // Go through the list_of_edge_pixels and check if the pixel belongs to the outer edge
//            // the idea here is to project horizontal and vertical liines from the pixel to the edges of the image.
//            // In case there is no object pixel between the current pixel and the edge, the pixel belongs to the outer edge,
//            // and we can extract the outer edges by region growing. The remaining edges are inner edges.
//            bdListIterator< bdVoxel > it;
//            for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
//            {
//                unsigned int r = it.GetElementPointer()->R();
//                unsigned int c = it.GetElementPointer()->C();
//                
//                if(thresholded_img(t,s,r,c) == edge_label)
//                {
//                
//                    int is_object_found = 0;
//                    for(int rn = r-1; rn>=0; rn-- )
//                    {
//                        if(thresholded_img(t,s,rn,c)!=0)
//                        {
//                            is_object_found = 1;
//                            break;
//                        }
//                    }
//                    
//                    if(is_object_found)
//                    {
//                        is_object_found = 0;
//                        for(int rn = r+1; rn<thresholded_img.GetNumberOfRows(); rn++ )
//                        {
//                            if(thresholded_img(t,s,rn,c)!=0)
//                            {
//                                is_object_found = 1;
//                                break;
//                            }
//                        }
//                    }
//                    
//                    if(is_object_found)
//                    {
//                        is_object_found = 0;
//                        for(int cn = c-1; cn>=0; cn-- )
//                        {
//                            if(thresholded_img(t,s,r,cn)!=0)
//                            {
//                                is_object_found = 1;
//                                break;
//                            }
//                        }
//                    }
//                    
//                    if(is_object_found)
//                    {
//                        is_object_found = 0;
//                        for(int cn = c+1; cn<thresholded_img.GetNumberOfColumns(); cn++ )
//                        {
//                            if(thresholded_img(t,s,r,cn)!=0)
//                            {
//                                is_object_found = 1;
//                                break;
//                            }
//                        }
//                    }
//                    
//                    // if no object is found, the current pixel belongs to the outer edge, so we have to track other pixels by region growing.
//                    if(!is_object_found)
//                    {
//                        bdList< bdVoxel > temp_list;
//                        temp_list.AddToRightEnd(it.GetElement());
//                        thresholded_img(t,s,r,c) = 128;//thresholded_label;
//                        while(!temp_list.IsEmpty())
//                        {
//                            bdVoxel v = temp_list.GetLeftEnd();
//                            
//                            int rnn, cnn;
//                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//                            {
//                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
//                                {
//                                    //cout<<" first_RG ";
//                                    bdVoxel v_temp;
//                                    v_temp(s,rnn,cnn);
//                                    temp_list.AddToRightEnd(v_temp);
//                                    thresholded_img(t,s,rnn,cnn) = 128;//thresholded_label;
//                                }
//                            }
//                        
//                            temp_list.DeleteLeftEnd();
//                        }
//                    }
//                }
//            }
//            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
            
//            // Go through the segmented region and grow non-segmented regions adjacent to inner edges. Store these regions to output and thresholded_img,
//            // so that these regions are not re-visited in future iterations.
//            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
//            {
//                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
//                {
//                    unsigned int r = sorted_array[i][j].R();
//                    unsigned int c = sorted_array[i][j].C();
//                    
//                    if(thresholded_img(t,s,r,c)==edge_label)// if it's an edge pixel...
//                    {
//                        //... find the background voxel (if any)...
//                        int rn, cn;
//                        for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
//                        {
//                            if(thresholded_img(t,s,rn,cn)==0)//if background found...
//                            {
//                                //cout<<" second_RG ";
//                                //... perform region growing on bckground label.
//                                bdList< bdVoxel > temp_list;
//                                bdVoxel v;
//                                v(rn,cn);
//                                temp_list.AddToRightEnd(v);
//                                thresholded_img(t,s,rn,cn) = thresholded_label;
//                                output(t,s,rn,cn) = thresholded_label;
//                                while(!temp_list.IsEmpty())
//                                {
//                                    v = temp_list.GetLeftEnd();
//                                    
//                                    int rnn, cnn;
//                                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//                                    {
//                                        if(thresholded_img(t,s,rnn,cnn)==0)//if background pixel
//                                        {
//                                            bdVoxel v_temp;
//                                            v_temp(s,rnn,cnn);
//                                            temp_list.AddToRightEnd(v_temp);
//                                            thresholded_img(t,s,rnn,cnn) = thresholded_label;
//                                            output(t,s,rnn,cnn) = thresholded_label;
//                                        }
//                                    }
//                                    
//                                    temp_list.DeleteLeftEnd();
//                                }
//                            }
//                        }
//                        
//                    }
//                }
//            }
//        }
        //output.CopyFrom(thresholded_img);
    }
    
    return 1;
}



//int bdMelanoma::SegmentMelanomaInnerStructure5(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
//{
//    unsigned int s = 0;
//    
//    output.SetSizeAndPropertiesAs(original);
//    output.FillInWith(0);
//    
//    unsigned short edge_label = 350;//65535;
//    unsigned short thresholded_label = 255;
//    
//    unsigned int t = 0;
//    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
//    {
//        bdArray<int> histogram;
//        unsigned int range_min = 0;
//        unsigned int range_max = 65535;
//        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
//        
//        // create array of sorted voxel values. This will be used to access voxels in order of their values.
//        bdArray< bdArray< bdVoxel > > sorted_array;
//        sorted_array.Set(histogram.GetNumberOfElements());
//        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
//        {
//            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
//        }
//        
//        //cout<<" Creating list ";
//        // put the masked voxels from the original image to the sorted_array
//        bdArray<int> index_array;
//        index_array.Set(sorted_array.GetNumberOfElements());
//        index_array.FillInWith(0);
//        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//            {
//                if(segmented_melanoma(t,s,r,c)!=0)
//                {
//                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
//                    (index_array[original(t,s,r,c)])++;
//                }
//            }
//        }
//        
//        //cout<<" Segmenting ";
//        
//        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
//        // label value. For a pixel that already has a neighbor with a label value give that label value.
//        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
//        bdGeometry g;
//        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
//        
//        bdImage thresholded_img;
//        thresholded_img.SetSizeAndPropertiesAs(original);
//        thresholded_img.FillInWith(0);
//        
//        int step = 6;
//        unsigned int thr=30+step;
//        
//        // Perform initial thresholding to fill in parts that do not fall in range of examined threshold values
//        for(unsigned int i=0; i<thr-step; i++)// 'i' is the threshold value
//        {
//            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
//            {
//                
//                unsigned int r = sorted_array[i][j].R();
//                unsigned int c = sorted_array[i][j].C();
//                
//                // perform trhesholdeing
//                thresholded_img(t,s,r,c) = thresholded_label;
//            }
//        }
//        
//        //for(unsigned int thr=40+step; thr<100; thr=thr+step)// 'i' is the threshold value
//        for( ; thr<150; thr=thr+step)
//        {
//            cout<<" thr="<<thr<<" ";
//            
//            // create the thresholded image
//            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
//            {
//                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
//                {
//                    
//                    unsigned int r = sorted_array[i][j].R();
//                    unsigned int c = sorted_array[i][j].C();
//                    
//                    // perform trhesholdeing
//                    thresholded_img(t,s,r,c) = thresholded_label;
//                }
//            }
//            
//            // find edges of the thresholded image, add pixels to list for fast access
//            bdList< bdVoxel > list_of_edge_pixels;
//            for(unsigned int i=thr-step; i<thr; i++)//for(unsigned int i=40; i<thr; i++)// 'i' is the threshold value
//            {
//                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
//                {
//                    unsigned int r = sorted_array[i][j].R();
//                    unsigned int c = sorted_array[i][j].C();
//                    int rn, cn;
//                    for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
//                    {
//                        if(thresholded_img(t,s,rn,cn)==0)//if edge found
//                        {
//                            list_of_edge_pixels.AddToRightEnd(sorted_array[i][j]);
//                            thresholded_img(t,s,r,c) = edge_label;
//                            break;
//                        }
//                    }
//                    
//                }
//            }
//            // At this point: edge pixels are in the list_of_edge_pixels and are marked with value edge_label in the thresholded_img.
//            
//            // Filter the list of edge pixels to maintin only a single pixel per edge.
//            bdList< bdVoxel > list_of_filtered_edge_pixels;
//            {
//                bdListIterator< bdVoxel > it;
//                for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
//                {
//                    unsigned int r = it.GetElementPointer()->R();
//                    unsigned int c = it.GetElementPointer()->C();
//                    
//                    if(thresholded_img(t,s,r,c) == edge_label)
//                    {
//                        list_of_filtered_edge_pixels.AddToRightEnd(it.GetElement());
//                        
//                        bdList< bdVoxel > temp_list;
//                        temp_list.AddToRightEnd(it.GetElement());
//                        thresholded_img(t,s,r,c) = thresholded_label;
//                        while(!temp_list.IsEmpty())
//                        {
//                            bdVoxel v = temp_list.GetLeftEnd();
//                            
//                            int rnn, cnn;
//                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//                            {
//                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
//                                {
//                                    bdVoxel v_temp;
//                                    v_temp(s,rnn,cnn);
//                                    temp_list.AddToRightEnd(v_temp);
//                                    thresholded_img(t,s,rnn,cnn) = thresholded_label;
//                                }
//                            }
//                            
//                            temp_list.DeleteLeftEnd();
//                        }
//                    }
//                }
//            }
//            // At this point: list_of_filtered_edge_pixels contains list of pixels where each pixel represents one contour (taken from the contour).
//            // Edge pixels are no longer marked in the thresholded_img.
//            
//            // Go through the list_of_filtered_edge_pixels and check if the pixel belongs to the outer edge by growing a region
//            // and checking if it touches the border of the segmented region (pixel value 0 in segmented image).
//            bdListIterator< bdVoxel > it;
//            for(it.SetLeftEnd(list_of_filtered_edge_pixels); it.IsValid(); it.MoveRight())
//            {
//                //unsigned int r = it.GetElementPointer()->R();
//                //unsigned int c = it.GetElementPointer()->C();
//                
//                // OVDE STAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                
//                // list of pixels belonging to the region
//                bdList< bdVoxel > region_pixels_list;
//                //region_pixels_list.AddToRightEnd(it.GetElement());
//                
//                bdList< bdVoxel > temp_list;
//                temp_list.AddToRightEnd(it.GetElement());
//                //thresholded_img(t,s,r,c) = 128;//thresholded_label;
//                while(!temp_list.IsEmpty())
//                {
//                    bdVoxel v = temp_list.GetLeftEnd();
//                    
//                    int rnn, cnn;
//                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//                    {
//                        if(thresholded_img(t,s,rnn,cnn)==0)//if edge found
//                        {
//                            // if the segmented image is 0 it means that we have reached the end of the segmented melanoma, so the edge is outer...
//                            if(segmented_melanoma(t,s,rnn,cnn)==0)
//                            {
//                                //... and since the edge is outer, we need to undo the region growing...
//                                while(!region_pixels_list.IsEmpty())
//                                {
//                                    bdVoxel v2 = region_pixels_list.GetLeftEnd();
//                                    thresholded_img(t,s,v2.R(),v2.C()) = 0;
//                                    output(t,s,v2.R(),v2.C()) = 0;
//                                    region_pixels_list.DeleteLeftEnd();
//                                }
//                                
//                                //region_pixels_list.Reset();
//                                temp_list.Reset();
//                                break;
//                            }
//                            // if the segmented image is non-zero we consider this edge an inner egde (untill proven otherwise)...
//                            else
//                            {
//                                //... so perform segmentation.
//                                thresholded_img(t,s,rnn,cnn) = 128;
//                                output(t,s,rnn,cnn) = 255;
//                                
//                                bdVoxel v_temp;
//                                v_temp(s,rnn,cnn);
//                                temp_list.AddToRightEnd(v_temp);
//                                region_pixels_list.AddToRightEnd(v_temp);
//                            }
//                        }
//                    }
//                    
//                    temp_list.DeleteLeftEnd();
//                }
//                
//                
//            }
//            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
//        }
//        
//        
//        
//        //            // Go through the list_of_edge_pixels and check if the pixel belongs to the outer edge
//        //            // the idea here is to project horizontal and vertical liines from the pixel to the edges of the image.
//        //            // In case there is no object pixel between the current pixel and the edge, the pixel belongs to the outer edge,
//        //            // and we can extract the outer edges by region growing. The remaining edges are inner edges.
//        //            bdListIterator< bdVoxel > it;
//        //            for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
//        //            {
//        //                unsigned int r = it.GetElementPointer()->R();
//        //                unsigned int c = it.GetElementPointer()->C();
//        //
//        //                if(thresholded_img(t,s,r,c) == edge_label)
//        //                {
//        //
//        //                    int is_object_found = 0;
//        //                    for(int rn = r-1; rn>=0; rn-- )
//        //                    {
//        //                        if(thresholded_img(t,s,rn,c)!=0)
//        //                        {
//        //                            is_object_found = 1;
//        //                            break;
//        //                        }
//        //                    }
//        //
//        //                    if(is_object_found)
//        //                    {
//        //                        is_object_found = 0;
//        //                        for(int rn = r+1; rn<thresholded_img.GetNumberOfRows(); rn++ )
//        //                        {
//        //                            if(thresholded_img(t,s,rn,c)!=0)
//        //                            {
//        //                                is_object_found = 1;
//        //                                break;
//        //                            }
//        //                        }
//        //                    }
//        //
//        //                    if(is_object_found)
//        //                    {
//        //                        is_object_found = 0;
//        //                        for(int cn = c-1; cn>=0; cn-- )
//        //                        {
//        //                            if(thresholded_img(t,s,r,cn)!=0)
//        //                            {
//        //                                is_object_found = 1;
//        //                                break;
//        //                            }
//        //                        }
//        //                    }
//        //
//        //                    if(is_object_found)
//        //                    {
//        //                        is_object_found = 0;
//        //                        for(int cn = c+1; cn<thresholded_img.GetNumberOfColumns(); cn++ )
//        //                        {
//        //                            if(thresholded_img(t,s,r,cn)!=0)
//        //                            {
//        //                                is_object_found = 1;
//        //                                break;
//        //                            }
//        //                        }
//        //                    }
//        //
//        //                    // if no object is found, the current pixel belongs to the outer edge, so we have to track other pixels by region growing.
//        //                    if(!is_object_found)
//        //                    {
//        //                        bdList< bdVoxel > temp_list;
//        //                        temp_list.AddToRightEnd(it.GetElement());
//        //                        thresholded_img(t,s,r,c) = 128;//thresholded_label;
//        //                        while(!temp_list.IsEmpty())
//        //                        {
//        //                            bdVoxel v = temp_list.GetLeftEnd();
//        //
//        //                            int rnn, cnn;
//        //                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//        //                            {
//        //                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
//        //                                {
//        //                                    //cout<<" first_RG ";
//        //                                    bdVoxel v_temp;
//        //                                    v_temp(s,rnn,cnn);
//        //                                    temp_list.AddToRightEnd(v_temp);
//        //                                    thresholded_img(t,s,rnn,cnn) = 128;//thresholded_label;
//        //                                }
//        //                            }
//        //
//        //                            temp_list.DeleteLeftEnd();
//        //                        }
//        //                    }
//        //                }
//        //            }
//        //            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
//        
//        //            // Go through the segmented region and grow non-segmented regions adjacent to inner edges. Store these regions to output and thresholded_img,
//        //            // so that these regions are not re-visited in future iterations.
//        //            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
//        //            {
//        //                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
//        //                {
//        //                    unsigned int r = sorted_array[i][j].R();
//        //                    unsigned int c = sorted_array[i][j].C();
//        //
//        //                    if(thresholded_img(t,s,r,c)==edge_label)// if it's an edge pixel...
//        //                    {
//        //                        //... find the background voxel (if any)...
//        //                        int rn, cn;
//        //                        for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
//        //                        {
//        //                            if(thresholded_img(t,s,rn,cn)==0)//if background found...
//        //                            {
//        //                                //cout<<" second_RG ";
//        //                                //... perform region growing on bckground label.
//        //                                bdList< bdVoxel > temp_list;
//        //                                bdVoxel v;
//        //                                v(rn,cn);
//        //                                temp_list.AddToRightEnd(v);
//        //                                thresholded_img(t,s,rn,cn) = thresholded_label;
//        //                                output(t,s,rn,cn) = thresholded_label;
//        //                                while(!temp_list.IsEmpty())
//        //                                {
//        //                                    v = temp_list.GetLeftEnd();
//        //
//        //                                    int rnn, cnn;
//        //                                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
//        //                                    {
//        //                                        if(thresholded_img(t,s,rnn,cnn)==0)//if background pixel
//        //                                        {
//        //                                            bdVoxel v_temp;
//        //                                            v_temp(s,rnn,cnn);
//        //                                            temp_list.AddToRightEnd(v_temp);
//        //                                            thresholded_img(t,s,rnn,cnn) = thresholded_label;
//        //                                            output(t,s,rnn,cnn) = thresholded_label;
//        //                                        }
//        //                                    }
//        //                                    
//        //                                    temp_list.DeleteLeftEnd();
//        //                                }
//        //                            }
//        //                        }
//        //                        
//        //                    }
//        //                }
//        //            }
//        //        }
//        //output.CopyFrom(thresholded_img);
//    }
//    
//    return 1;
//}



int bdMelanoma::SegmentMelanomaInnerStructure5(bdImage &original_img, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    
    unsigned int range_min = 0;
    unsigned int range_max = 100;
    
    bdImage original;
    {
        bdGIP gip;
        gip.RescaleWholeRange(original_img, range_min, range_max, original);
    }
    
    
    output.SetSizeAndPropertiesAs(original);
    output.FillInWith(0);
    
    unsigned short edge_label = 350;//65535;
    unsigned short thresholded_label = 255;
    
    unsigned int t = 0;
    //    for(unsigned int t = 0; t<original.GetNumberOfTimeSeries(); t++)
    {
        bdArray<int> histogram;
        //unsigned int range_min = 0;
        //unsigned int range_max = 65535;
        bdMelanoma::Histogram(original, segmented_melanoma, histogram, range_min, range_max, t, 0, 0, 0);
        
        // create array of sorted voxel values. This will be used to access voxels in order of their values.
        bdArray< bdArray< bdVoxel > > sorted_array;
        sorted_array.Set(histogram.GetNumberOfElements());
        for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
        {
            if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
        }
        
        //cout<<" Creating list ";
        // put the masked voxels from the original image to the sorted_array
        bdArray<int> index_array;
        index_array.Set(sorted_array.GetNumberOfElements());
        index_array.FillInWith(0);
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(t,s,r,c)!=0)
                {
                    sorted_array[original(t,s,r,c)][(index_array[original(t,s,r,c)])](r,c);
                    (index_array[original(t,s,r,c)])++;
                }
            }
        }
        
        //cout<<" Segmenting ";
        
        // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
        // label value. For a pixel that already has a neighbor with a label value give that label value.
        // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
        bdGeometry g;
        g.SetDimensions(1, original.GetNumberOfRows(), original.GetNumberOfColumns());
        
        bdImage thresholded_img;
        thresholded_img.SetSizeAndPropertiesAs(original);
        thresholded_img.FillInWith(0);
        
        int step = 6;
        unsigned int thr=range_min+step; //30+step;
        
        // Perform initial thresholding to fill in parts that do not fall in range of examined threshold values
        for(unsigned int i=0; i<thr-step; i++)// 'i' is the threshold value
        {
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                
                // perform trhesholdeing
                thresholded_img(t,s,r,c) = thresholded_label;
            }
        }
        
        for( ; thr<range_max; thr=thr+step) //for( ; thr<150; thr=thr+step)
        {
            cout<<" thr="<<thr<<" ";
            
            // create the thresholded image
            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
            {
                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
                {
                    
                    unsigned int r = sorted_array[i][j].R();
                    unsigned int c = sorted_array[i][j].C();
                    
                    // perform trhesholdeing
                    thresholded_img(t,s,r,c) = thresholded_label;
                }
            }
            
            // find edges of the thresholded image, add pixels to list for fast access
            bdList< bdVoxel > list_of_edge_pixels;
            for(unsigned int i=thr-step; i<thr; i++)//for(unsigned int i=40; i<thr; i++)// 'i' is the threshold value
            {
                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
                {
                    unsigned int r = sorted_array[i][j].R();
                    unsigned int c = sorted_array[i][j].C();
                    int rn, cn;
                    for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
                    {
                        if(thresholded_img(t,s,rn,cn)==0)//if edge found
                        {
                            list_of_edge_pixels.AddToRightEnd(sorted_array[i][j]);
                            thresholded_img(t,s,r,c) = edge_label;
                            break;
                        }
                    }
                    
                }
            }
            // At this point: edge pixels are in the list_of_edge_pixels and are marked with value edge_label in the thresholded_img.
            
            // Filter the list of edge pixels to maintin only a single pixel per edge.
            bdList< bdVoxel > list_of_filtered_edge_pixels;
            {
                bdListIterator< bdVoxel > it;
                for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
                {
                    unsigned int r = it.GetElementPointer()->R();
                    unsigned int c = it.GetElementPointer()->C();
                    
                    if(thresholded_img(t,s,r,c) == edge_label)
                    {
                        list_of_filtered_edge_pixels.AddToRightEnd(it.GetElement());
                        
                        bdList< bdVoxel > temp_list;
                        temp_list.AddToRightEnd(it.GetElement());
                        thresholded_img(t,s,r,c) = thresholded_label;
                        while(!temp_list.IsEmpty())
                        {
                            bdVoxel v = temp_list.GetLeftEnd();
                            
                            int rnn, cnn;
                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
                            {
                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
                                {
                                    bdVoxel v_temp;
                                    v_temp(s,rnn,cnn);
                                    temp_list.AddToRightEnd(v_temp);
                                    thresholded_img(t,s,rnn,cnn) = thresholded_label;
                                }
                            }
                            
                            temp_list.DeleteLeftEnd();
                        }
                    }
                }
            }
            // At this point: list_of_filtered_edge_pixels contains list of pixels where each pixel represents one contour (taken from the contour).
            // Edge pixels are no longer marked in the thresholded_img.
            
            // Go through the list_of_filtered_edge_pixels and check if the pixel belongs to the outer edge by growing a region
            // and checking if it touches the border of the segmented region (pixel value 0 in segmented image).
            bdListIterator< bdVoxel > it;
            for(it.SetLeftEnd(list_of_filtered_edge_pixels); it.IsValid(); it.MoveRight())
            {
                //unsigned int r = it.GetElementPointer()->R();
                //unsigned int c = it.GetElementPointer()->C();
                
                // OVDE STAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                // list of pixels belonging to the region
                bdList< bdVoxel > region_pixels_list;
                //region_pixels_list.AddToRightEnd(it.GetElement());
                
                bdList< bdVoxel > temp_list;
                temp_list.AddToRightEnd(it.GetElement());
                //thresholded_img(t,s,r,c) = 128;//thresholded_label;
                while(!temp_list.IsEmpty())
                {
                    bdVoxel v = temp_list.GetLeftEnd();
                    
                    int rnn, cnn;
                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
                    {
                        if(thresholded_img(t,s,rnn,cnn)==0)//if edge found
                        {
                            // if the segmented image is 0 it means that we have reached the end of the segmented melanoma, so the edge is outer...
                            if(segmented_melanoma(t,s,rnn,cnn)==0)
                            {
                                //... and since the edge is outer, we need to undo the region growing...
                                while(!region_pixels_list.IsEmpty())
                                {
                                    bdVoxel v2 = region_pixels_list.GetLeftEnd();
                                    thresholded_img(t,s,v2.R(),v2.C()) = 0;
                                    output(t,s,v2.R(),v2.C()) = 0;
                                    region_pixels_list.DeleteLeftEnd();
                                }
                                
                                //region_pixels_list.Reset();
                                temp_list.Reset();
                                break;
                            }
                            // if the segmented image is non-zero we consider this edge an inner egde (untill proven otherwise)...
                            else
                            {
                                //... so perform segmentation.
                                thresholded_img(t,s,rnn,cnn) = 128;
                                output(t,s,rnn,cnn) = 255;
                                
                                bdVoxel v_temp;
                                v_temp(s,rnn,cnn);
                                temp_list.AddToRightEnd(v_temp);
                                region_pixels_list.AddToRightEnd(v_temp);
                            }
                        }
                    }
                    
                    temp_list.DeleteLeftEnd();
                }
                
                
            }
            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
        }
        
        
        
        //            // Go through the list_of_edge_pixels and check if the pixel belongs to the outer edge
        //            // the idea here is to project horizontal and vertical liines from the pixel to the edges of the image.
        //            // In case there is no object pixel between the current pixel and the edge, the pixel belongs to the outer edge,
        //            // and we can extract the outer edges by region growing. The remaining edges are inner edges.
        //            bdListIterator< bdVoxel > it;
        //            for(it.SetLeftEnd(list_of_edge_pixels); it.IsValid(); it.MoveRight())
        //            {
        //                unsigned int r = it.GetElementPointer()->R();
        //                unsigned int c = it.GetElementPointer()->C();
        //
        //                if(thresholded_img(t,s,r,c) == edge_label)
        //                {
        //
        //                    int is_object_found = 0;
        //                    for(int rn = r-1; rn>=0; rn-- )
        //                    {
        //                        if(thresholded_img(t,s,rn,c)!=0)
        //                        {
        //                            is_object_found = 1;
        //                            break;
        //                        }
        //                    }
        //
        //                    if(is_object_found)
        //                    {
        //                        is_object_found = 0;
        //                        for(int rn = r+1; rn<thresholded_img.GetNumberOfRows(); rn++ )
        //                        {
        //                            if(thresholded_img(t,s,rn,c)!=0)
        //                            {
        //                                is_object_found = 1;
        //                                break;
        //                            }
        //                        }
        //                    }
        //
        //                    if(is_object_found)
        //                    {
        //                        is_object_found = 0;
        //                        for(int cn = c-1; cn>=0; cn-- )
        //                        {
        //                            if(thresholded_img(t,s,r,cn)!=0)
        //                            {
        //                                is_object_found = 1;
        //                                break;
        //                            }
        //                        }
        //                    }
        //
        //                    if(is_object_found)
        //                    {
        //                        is_object_found = 0;
        //                        for(int cn = c+1; cn<thresholded_img.GetNumberOfColumns(); cn++ )
        //                        {
        //                            if(thresholded_img(t,s,r,cn)!=0)
        //                            {
        //                                is_object_found = 1;
        //                                break;
        //                            }
        //                        }
        //                    }
        //
        //                    // if no object is found, the current pixel belongs to the outer edge, so we have to track other pixels by region growing.
        //                    if(!is_object_found)
        //                    {
        //                        bdList< bdVoxel > temp_list;
        //                        temp_list.AddToRightEnd(it.GetElement());
        //                        thresholded_img(t,s,r,c) = 128;//thresholded_label;
        //                        while(!temp_list.IsEmpty())
        //                        {
        //                            bdVoxel v = temp_list.GetLeftEnd();
        //
        //                            int rnn, cnn;
        //                            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
        //                            {
        //                                if(thresholded_img(t,s,rnn,cnn)==edge_label)//if edge found
        //                                {
        //                                    //cout<<" first_RG ";
        //                                    bdVoxel v_temp;
        //                                    v_temp(s,rnn,cnn);
        //                                    temp_list.AddToRightEnd(v_temp);
        //                                    thresholded_img(t,s,rnn,cnn) = 128;//thresholded_label;
        //                                }
        //                            }
        //
        //                            temp_list.DeleteLeftEnd();
        //                        }
        //                    }
        //                }
        //            }
        //            // At this point: only the inner edge pixels have edge_label values in the thresholded_img.
        
        //            // Go through the segmented region and grow non-segmented regions adjacent to inner edges. Store these regions to output and thresholded_img,
        //            // so that these regions are not re-visited in future iterations.
        //            for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
        //            {
        //                for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
        //                {
        //                    unsigned int r = sorted_array[i][j].R();
        //                    unsigned int c = sorted_array[i][j].C();
        //
        //                    if(thresholded_img(t,s,r,c)==edge_label)// if it's an edge pixel...
        //                    {
        //                        //... find the background voxel (if any)...
        //                        int rn, cn;
        //                        for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
        //                        {
        //                            if(thresholded_img(t,s,rn,cn)==0)//if background found...
        //                            {
        //                                //cout<<" second_RG ";
        //                                //... perform region growing on bckground label.
        //                                bdList< bdVoxel > temp_list;
        //                                bdVoxel v;
        //                                v(rn,cn);
        //                                temp_list.AddToRightEnd(v);
        //                                thresholded_img(t,s,rn,cn) = thresholded_label;
        //                                output(t,s,rn,cn) = thresholded_label;
        //                                while(!temp_list.IsEmpty())
        //                                {
        //                                    v = temp_list.GetLeftEnd();
        //
        //                                    int rnn, cnn;
        //                                    for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
        //                                    {
        //                                        if(thresholded_img(t,s,rnn,cnn)==0)//if background pixel
        //                                        {
        //                                            bdVoxel v_temp;
        //                                            v_temp(s,rnn,cnn);
        //                                            temp_list.AddToRightEnd(v_temp);
        //                                            thresholded_img(t,s,rnn,cnn) = thresholded_label;
        //                                            output(t,s,rnn,cnn) = thresholded_label;
        //                                        }
        //                                    }
        //
        //                                    temp_list.DeleteLeftEnd();
        //                                }
        //                            }
        //                        }
        //
        //                    }
        //                }
        //            }
        //        }
        //output.CopyFrom(thresholded_img);
    }
    
    return 1;
}




int bdMelanoma::SegmentMelanomaInnerStructure6(bdImage &original, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int s = 0;
    
    unsigned short edge_label = 350;//65535;
    unsigned short thresholded_label = 255;
    
    
    // We continue to work on the slice which is a combination of the RGB cahnnels.
    bdImage slice;
    slice.SetSize(1,1, original.GetNumberOfRows(), original.GetNumberOfColumns());
    slice.SetVisualizationPropertiesToMatchInput(original);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // compute combined gray image as AVERAGE of RGB channels and rescale to 255 values
            slice(r,c) = (original(0,0,r,c) + original(1,0,r,c) + original(2,0,r,c)) / 3;
            
            
            //            // MEDIAN OF ALL CHANNELS
            //            if(input(1,0,r,c)>input(2,0,r,c))
            //            {
            //                if(input(2,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(2,0,r,c);
            //                else extracted_slice(r,c) = input(0,0,r,c);
            //            }
            //            else
            //            {
            //                if(input(1,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
            //                else extracted_slice(r,c) = input(0,0,r,c);
            //            }
            //
            //            // MAX OF ALL CHANNELS
            //            if(input(1,0,r,c)>input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
            //            else extracted_slice(r,c) = input(2,0,r,c);
            //            if(input(0,0,r,c)>extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);
            //
            //            // MIDRANGE OF ALL CHANNELS
            //            unsigned int min;
            //            if(input(1,0,r,c)<input(2,0,r,c)) min = input(1,0,r,c);
            //            else min = input(2,0,r,c);
            //            if(input(0,0,r,c)<min) min = input(0,0,r,c);
            //            extracted_slice(r,c) = (extracted_slice(r,c) + min) / 2;
            //
            //            // MIN OF ALL CHANNELS
            //            if(input(1,0,r,c)<input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
            //            else extracted_slice(r,c) = input(2,0,r,c);
            //            if(input(0,0,r,c)<extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);
        }
    }
    
    bdGIP gip;
    bdImage extracted_slice;
    unsigned int range_min = 0;
    unsigned int range_max = 255;
    gip.RescaleWholeRange(slice, range_min, range_max, extracted_slice);
    
    output.SetSizeAndPropertiesAs(extracted_slice);
    output.FillInWith(0);

    
    // Histogram of extracted slice
    bdArray<int> histogram;
    bdMelanoma::Histogram(extracted_slice, segmented_melanoma, histogram, 0, range_max, 0, 0, 0, 0);
    
    // Create array of sorted voxel values. This will be used to access voxels in order of their values.
    bdArray< bdArray< bdVoxel > > sorted_array;
    sorted_array.Set(histogram.GetNumberOfElements());
    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
    {
        if(histogram[i]!=0) sorted_array[i].Set(histogram[i]);
    }

    //cout<<" Creating list ";
    
    // Put the masked voxels from the original image to the sorted_array
    bdArray<int> index_array;
    index_array.Set(sorted_array.GetNumberOfElements());
    index_array.FillInWith(0);
    for(unsigned int r=0; r< extracted_slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<extracted_slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0)
            {
                sorted_array[extracted_slice(r,c)][(index_array[extracted_slice(r,c)])](r,c);
                (index_array[extracted_slice(r,c)])++;
            }
        }
    }
    
    
    // Go through the sorted array and for a pixel that has no neighbors in the output image give a new
    // label value. For a pixel that already has a neighbor with a label value give that label value.
    // For a pixel that has two or more pixels with DIFFERENT label values assign value 0.
    bdGeometry g;
    g.SetDimensions(1, extracted_slice.GetNumberOfRows(), extracted_slice.GetNumberOfColumns());
    
    bdImage thresholded_img;
    thresholded_img.SetSizeAndPropertiesAs(extracted_slice);
    thresholded_img.FillInWith(0);
    
    int step = 6;
    unsigned int thr = 6+step;
    
    // Perform initial thresholding to fill in parts that do not fall in range of examined threshold values
    for(unsigned int i=0; i<thr-step; i++)// 'i' is the threshold value
    {
        for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
        {
            
            unsigned int r = sorted_array[i][j].R();
            unsigned int c = sorted_array[i][j].C();
            
            if(thresholded_img(r,c)!=0) //DODAO!!!!!!!
            // perform trhesholdeing
            thresholded_img(r,c) = thresholded_label;
        }
    }
    
    
    for( ; thr<range_max; thr=thr+step)//for( ; thr<150; thr=thr+step)
    {
        cout<<" thr="<<thr<<" ";
        
        // create the thresholded image
        for(unsigned int i=thr-step; i<thr; i++)// 'i' is the threshold value
        {
            cout<<" i="<<i<<" ";
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                
                if(thresholded_img(r,c)!=0) //DODAO!!!!!!!
                // perform trhesholdeing
                thresholded_img(r,c) = thresholded_label;
            }
        }
        
        // find edges of the thresholded image, add pixels to list for fast access
        bdList< bdVoxel > list_of_edge_pixels;
        for(unsigned int i=thr-step; i<thr; i++)//for(unsigned int i=40; i<thr; i++)// 'i' is the threshold value
        {
            //cout<<" i="<<i<<" ";
            for(unsigned int j=0; j<sorted_array[i].GetNumberOfElements(); j++)
            {
                unsigned int r = sorted_array[i][j].R();
                unsigned int c = sorted_array[i][j].C();
                int rn, cn;
                for(g.ForCoordinates_4_Neighbors(r,c); g.Get_4_Neighbors(rn,cn); )
                {
                    if(thresholded_img(rn,cn)==0)//if edge found
                    {
                        //cout<<" & ";
                        list_of_edge_pixels.AddToRightEnd(sorted_array[i][j]);
                        thresholded_img(r,c) = edge_label;
                        break;
                    }
                }
                
            }
        }
        // At this point: edge pixels are in the list_of_edge_pixels and are marked with value edge_label in the thresholded_img.

        // Process the list of edge pixels by finding individual edges and looking if they are inner or outer edges.
        {
            //cout<<" 1 ";

            bdListIterator< bdVoxel > it;
            int debug_number_of_iterations = 0;
            for(it.SetLeftEnd(list_of_edge_pixels), debug_number_of_iterations=0; it.IsValid() && debug_number_of_iterations<2; it.MoveRight(), debug_number_of_iterations++)
            {
                //cout<<" 2 ";
                
                
                
                unsigned int r = it.GetElementPointer()->R();
                unsigned int c = it.GetElementPointer()->C();
                
               //output(r,c) = 255;
                
                if(thresholded_img(r,c) == edge_label)
                {
                    bdList< bdVoxel > single_edge_pixels;
                    
                    bdList< bdVoxel > temp_list;
                    temp_list.AddToRightEnd(it.GetElement());
                    thresholded_img(r,c) = thresholded_label;
                    while(!temp_list.IsEmpty())
                    {
                        bdVoxel v = temp_list.GetLeftEnd();
                        
                        single_edge_pixels.AddToRightEnd(v);
                        
                        int rnn, cnn;
                        for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rnn,cnn); )
                        {
                            if(thresholded_img(rnn,cnn)==edge_label)//if edge found
                            {
                                bdVoxel v_temp;
                                v_temp(s,rnn,cnn);
                                temp_list.AddToRightEnd(v_temp);
                                thresholded_img(rnn,cnn) = thresholded_label;
                            }
                        }
                        
                        temp_list.DeleteLeftEnd();
                    }
                    // At this point: list 'single_edge_pixels' is a list of pixels of a single contour.
                    // Its pixels are no longer marked in the thresholded_img (they were marked with 'edge_label' value).
                    
                    // Now we have to check if the edge is an inner edge by finding its edge of background voxels
                    // and looking if edge of background voxels is shorter than the edge (if it is, it is an inner edge).
                    unsigned int single_edge_size = single_edge_pixels.GetNumberOfElements();
                    
//                    while(!single_edge_pixels.IsEmpty())
//                    {
//                        bdVoxel v = single_edge_pixels.GetLeftEnd();
//
//                        //output(v.R(),v.C()) = 128;
//                        
//                        int rnn, cnn;
//                        for(g.ForCoordinates_4_Neighbors(v.R(),v.C()); g.Get_4_Neighbors(rnn,cnn); )
//                        {
//                            if(thresholded_img(rnn,cnn)==0 && segmented_melanoma(rnn,cnn)!=0) //if edge found
//                            {
//                                bdVoxel v_temp;
//                                v_temp(0,rnn,cnn);
//                                temp_list.AddToRightEnd(v_temp);
//                                thresholded_img(rnn,cnn) = single_edge_size;// !!!!!!!!!!!!!!!!!! 400;
//                                //output(v_temp.R(),v_temp.C()) = 255;
//
//                            }
//                        }
//                        
//                        single_edge_pixels.DeleteLeftEnd();
//                    }
//                    // At this point: temp_list contains the background pixels, which are marked in the thresholded_img with value 400.
                    
                    cout<<" single_edge_size="<<single_edge_size<<"  ";
                    cout<<" temp_list_size="<<temp_list.GetNumberOfElements()<<"  ";
                    
                    if(single_edge_size < 100)
                    {
                        while(!single_edge_pixels.IsEmpty())
                        {
                            bdVoxel v = single_edge_pixels.GetLeftEnd();
    
                            output(v.R(),v.C()) = 128;
    
                            single_edge_pixels.DeleteLeftEnd();
                        }
                        
                        

                    }
                    
                    
                    
                    
                    
//                    // If the size of the temp_list is smaller than single_edge_size, the edge is inner, so record the region to output.
//                    if(temp_list.GetNumberOfElements() < single_edge_size && temp_list.GetNumberOfElements()>0)
//                    {
//                        cout<<" * ";
//                        while(!temp_list.IsEmpty())
//                        {
//                            bdVoxel v = temp_list.GetLeftEnd();
//                            output(v.R(),v.C()) = 50*(debug_number_of_iterations+1);
//                            
//                            int rnn, cnn;
//                            for(g.ForCoordinates_4_Neighbors(v.R(),v.C()); g.Get_4_Neighbors(rnn,cnn); )
//                            {
//                                if(thresholded_img(rnn,cnn)==0)//if edge found
//                                {
//                                    bdVoxel v_temp;
//                                    v_temp(0,rnn,cnn);
//                                    temp_list.AddToRightEnd(v_temp);
//                                    thresholded_img(rnn,cnn) = 400;
//                                }
//                            }
//                            
//                            temp_list.DeleteLeftEnd();
//                        }
//                        
//                    }
//                    else
//                    {
//                        // The edge is outer edge, so go thorugh temp_list to reset the values in thresholded_img to 0.
//                        bdListIterator< bdVoxel > it2;
//                        for(it2.SetLeftEnd(temp_list); it2.IsValid(); it2.MoveRight())
//                        {
//                            thresholded_img(it.GetElementPointer()->R(),it.GetElementPointer()->C()) = 0;
//                        }
//                    }
                }
            }
        }
    }
    
    output.CopyFrom(thresholded_img);
    
    return 1;
}


int bdMelanoma::SegmentMelanomaInnerStructure7(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
{
    if(input_img.IsEmpty()) return 0;
    if(segmented_melanoma.IsEmpty()) return 0;
    
    bdGIP gip;
    
    bdImage input;
    gip.RescaleWholeRange(input_img, 0,255, input);
    
    
    bdImage original;// negative of the input;
    gip.Negative(input,original);
    
    cout<<" 1 ";
    
    
    // We continue to work on the slice which is a combination of the RGB cahnnels.
    bdImage extracted_slice;
    extracted_slice.SetSize(1,1, original.GetNumberOfRows(), original.GetNumberOfColumns());
    extracted_slice.SetVisualizationPropertiesToMatchInput(original);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            //// compute combined gray image as AVERAGE of RGB channels.
            //extracted_slice(0,0,r,c) = (original(0,0,r,c) + original(1,0,r,c) + original(2,0,r,c)) / 3;
            
                        // MEDIAN OF ALL CHANNELS
                        if(input(1,0,r,c)>input(2,0,r,c))
                        {
                            if(input(2,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(2,0,r,c);
                            else extracted_slice(r,c) = input(0,0,r,c);
                        }
                        else
                        {
                            if(input(1,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
                            else extracted_slice(r,c) = input(0,0,r,c);
                        }
            
//                        // MAX OF ALL CHANNELS
//                        if(input(1,0,r,c)>input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
//                        else extracted_slice(r,c) = input(2,0,r,c);
//                        if(input(0,0,r,c)>extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);
//            //
            //            // MIDRANGE OF ALL CHANNELS
            //            unsigned int min;
            //            if(input(1,0,r,c)<input(2,0,r,c)) min = input(1,0,r,c);
            //            else min = input(2,0,r,c);
            //            if(input(0,0,r,c)<min) min = input(0,0,r,c);
            //            extracted_slice(r,c) = (extracted_slice(r,c) + min) / 2;
            
            
            //            // MIN OF ALL CHANNELS
            //            if(input(1,0,r,c)<input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
            //            else extracted_slice(r,c) = input(2,0,r,c);
            //            if(input(0,0,r,c)<extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);
            
            
//            // BLUE CHANNEL
//            extracted_slice(r,c) = input(2,0,r,c);
//            
            
        }
    }
    
    cout<<" 2a ";
    
    
    // The algorithms which will be used work by extracting rigdes (local maxima), so we will perform
    // some morphology to make the regions more uniform (using median) and smooth
    /// the transition between the regions (mean).
    //    bdImage extracted_slice_max;
    //    gip.MaximumCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 36,extracted_slice_max);//25,extracted_slice_max);
    
    cout<<" 2b ";
    
    bdImage slice;
    //    gip.MeanCircleFor2DSlice(extracted_slice_max, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);
    
//    gip.MeanCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);
    
    
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVER-RIDDEN HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    slice.CopyFrom(extracted_slice);
    
    cout<<" 3 ";
    
    // Get range of values found in the masked region of the image
    unsigned int gray_value_max = 0;
    unsigned int gray_value_min = 65535;
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0)
            {
                if(slice(r,c)>gray_value_max) gray_value_max = slice(r,c);
                if(slice(r,c)<gray_value_min) gray_value_min = slice(r,c);
            }
        }
    }
    
    
    cout<<"  gray_value_max="<<gray_value_max<<", gray_value_min="<<gray_value_min<<"  ";
    
    
    // Calculate the histograms of the masked area in the morphologically modified image for a single component
    bdArray<int> histogram_of_slice;
    this->Histogram(slice,segmented_melanoma,histogram_of_slice,0,gray_value_max,0,0,0,0);
    
    
    // create array of sorted voxel values. This will be used to access voxels in order of their values.
    bdArray< bdArray< bdVoxel > > sorted_array;
    sorted_array.Set(histogram_of_slice.GetNumberOfElements());
    
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        if(histogram_of_slice[i]!=0) sorted_array[i].Set(histogram_of_slice[i]);
    }
    
    
    // Put the masked voxels from the original image to the sorted_array
    bdArray<int> index_array;
    index_array.Set(sorted_array.GetNumberOfElements());
    index_array.FillInWith(0);
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0 && slice(r,c)!=0)
            {
                sorted_array[slice(0,0,r,c)] [(index_array[slice(0,0,r,c)])] (r,c);
                (index_array[slice(0,0,r,c)])++;
            }
        }
    }
    
    
    cout<<" 4 ";
    
    
    output.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(slice);
    output.FillInWith(0);
    
    bdImage temp;
    temp.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    temp.SetVisualizationPropertiesToMatchInput(slice);
    temp.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    
    
    //   unsigned int min, max;
    //   slice.GetVoxelMinimumAndMaximumValue(&min, &max);
    
    unsigned int t=0, s=0;
    
    unsigned short index_of_new_CCs = 53;//we start AT LEAST with value 3 because value 1 is reserved for a new voxel not yet assigned to a region and 2 is the region boundary.
    //    for(int threshold = max; threshold>min; threshold--)
    for(int threshold = gray_value_max; threshold>gray_value_min; threshold--)
        
    {
        //cout<<"thr="<<threshold<<",noel="<<sorted_array[threshold].GetNumberOfElements()<<"  ";
        bdList<bdVoxel> added_pixels;
        //for(unsigned int r=0; r<temp.GetNumberOfRows(); r++)
        //{
        //for(unsigned int c=0; c<temp.GetNumberOfColumns(); c++)
        for(unsigned int j=0; j<sorted_array[threshold].GetNumberOfElements(); j++)
        {
            unsigned int r = sorted_array[threshold][j].R();
            unsigned int c = sorted_array[threshold][j].C();
            
            if(temp(t,s,r,c)==0)
            {
                if(slice(t,s,r,c)==threshold)//if(slice(t,s,r,c)>=threshold)
                {
                    int rn,cn;
                    int is_previous_region_found = 0;
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(100,rn,cn); )
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(9,rn,cn); )
                    for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(4,rn,cn); )
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(25,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(16,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(36,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(81,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                        
                    {
                        if(temp(t,s,rn,cn)!=0)
                        {
                            temp(t,s,r,c) = temp(t,s,rn,cn);
                            is_previous_region_found = 1;
                            break;
                        }
                    }
                    if(!is_previous_region_found)
                    {
                        temp(t,s,r,c) = index_of_new_CCs;
                        index_of_new_CCs++;
                    }
                    bdVoxel v; v(r,c);
                    added_pixels.AddToRightEnd(v);
                }
            }
        }
        //}
        
        //Go through the list of added pixels and find pixels that border with pixels of other value (that is greater than 1), and make borders for these regions.
        //bdListIterator<bdVoxel> it;
        while(!added_pixels.IsEmpty())
        {
            bdVoxel *pv = added_pixels.GetLeftEndPointerToElement();
            unsigned short region_value = temp(t,s,pv->R(),pv->C());
            
            //output(0,0,pv->R(),pv->C()) = region_value;
            
            //check if it is a border with other region...
            if(region_value>1)
            {
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(pv->R(),pv->C()); g.Get_8_Neighbors(rn,cn); )
                {
                    if(temp(t,s,rn,cn)!=region_value && temp(t,s,rn,cn)>0)
                    {
                        //... if it is a border, grow a bound for the region (convert its edges into bounds)...
                        bdRegion2D reg;//, edge;
                        //reg.CreateRegion_8_FromSeedAndValueHigherThan(temp,threshold,pv->C(),pv->R(),0,0);
                        reg.CreateRegion_8_FromSeedPoint(temp,region_value,t,s,pv->R(),pv->C());
                        bdRegion2DIterator itr;
                        for(itr.SetBegin(&reg); itr.IsValid(); itr.MoveToNext())
                        {
                            output(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = region_value;
                            temp(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = 1;
                        }
                        
                        break;
                    }
                }
            }
            
            added_pixels.DeleteLeftEnd();
        }
    }
    
    
    //output.CopyFrom(extracted_slice);
    
    return 1;
}


//int bdMelanoma::SegmentMelanomaInnerStructure8(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
//{
//    unsigned int SE_radius = 6;
//    unsigned short segmented_label = 65535;
//    
//    if(!this->SegmentStructures(input_img, 0, segmented_melanoma, (SE_radius*SE_radius), 80, output, segmented_label)) return 0;
//    
//    //----- Erosion of the inner structures (Comment out if not needed) -----
//    {
//        cout<<"  Erosion of segmented inner bright structure...";
//        bdBIP bip;
//        bdImage temp;
//        bdGIP gip;
//        gip.Negative(output, temp);
//        bip.ErosionByMapping_Circle(temp, segmented_melanoma, output, 1);
//        cout<<" completed."<<endl;
//    }
//    //----------
//    
//    
//    //----- Label connected components -----
//    {
//        unsigned short label = 1;
//        
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(output(r,c)!=0)
//                {
//                    if(output(r,c)==segmented_label)
//                    {
//                        //perform region growing
//                        bdList< bdVoxel > temp_list;
//                        bdVoxel v;
//                        temp_list.AddToLeftEnd(v(r,c));
//                        output(v.R(),v.C()) = label;
//                        
//                        bdGeometry g;
//                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//                        while(!temp_list.IsEmpty())
//                        {
//                            v = temp_list.GetLeftEnd();
//                            
//                            int rn, cn;
//                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
//                            {
//                                if(output(rn,cn)==segmented_label)
//                                {
//                                    bdVoxel vn;
//                                    temp_list.AddToRightEnd(vn(rn,cn));
//                                    output(vn.R(),vn.C()) = label;
//                                }
//                            }
//                            
//                            temp_list.DeleteLeftEnd();
//                        }
//                        
//                        label = (label+1) % segmented_label;
//                    }
//                }
//            }
//        }
//    }
//    //----------
//    bdImage output2;
//    output2.SetSizeAndPropertiesAs(output);
//    output2.FillInWith(0);
//    
//    {
//        bdGeometry g;
//        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//        
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(segmented_melanoma(r,c)!=0)
//                {
//                    if(output(r,c)==0)
//                    {
//                        int number_of_CCs = 0;
//                        unsigned short CC_values[3];// 3 different CCs are enough to proclaim network structure.
//                        int rn, cn;
//                        for(g.ForCoordinates_Circle(r,c,SE_radius*SE_radius); g.Get_Circle((SE_radius+1)*(SE_radius+1), rn, cn); )
//                        {
//                            if(output(rn,cn)!=0)
//                            {
//                                int is_CC_already_recorded = 0;
//                                for(unsigned int i=0; i<number_of_CCs && i<3; i++)
//                                {
//                                    if(CC_values[i]==output(rn,cn))
//                                    {
//                                        is_CC_already_recorded = 1;
//                                        break;
//                                    }
//                                }
//                                
//                                if(!is_CC_already_recorded)
//                                {
//                                    CC_values[number_of_CCs] = output(rn,cn);
//                                    number_of_CCs++;
//                                }
//                                
//                                if(number_of_CCs>=3)
//                                {
//                                    output2(r,c) = 255;
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    output.CopyFrom(output2);
//    
//    
//    return 1;
//}



//int bdMelanoma::SegmentMelanomaInnerStructure8(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
//{
//    unsigned int SE_radius = 6;
//    unsigned short segmented_label = 65535;
//    unsigned int t_input = 0;
//    unsigned int threshold = 0;
//    
//    if(!this->SegmentStructures(input_img, t_input, segmented_melanoma, (SE_radius*SE_radius), 70, output, segmented_label)) return 0;
//    //if(!this->SegmentStructures(input_img, t_input, segmented_melanoma, (SE_radius*SE_radius), 60, output, segmented_label)) return 0;
//    
//    bdArray<int> histogram;
//    this->Histogram(input_img, segmented_melanoma, histogram, 0, 255, t_input, 0, 0, 0);
//    //this->Histogram(input_img, input_img, histogram, 0, 255, t_input, 0, 0, 0);
//    int max_value = histogram[0];
//    int max_value_index = 0;
//    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
//    {
//        if(max_value<histogram[i])
//        {
//            max_value = histogram[i];
//            max_value_index = i;
//        }
//    }
//    cout<<" max_value_index="<<max_value_index<<"  ";
//    
//    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
//    unsigned int last_non_zero_value = 1;
//    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
//    {
//        if(histogram[i]==0)
//        {
//            histogram[i] = last_non_zero_value;
//        }
//        else
//        {
//            last_non_zero_value = histogram[i];
//        }
//    }
//    
//    //----- Determine threshold value 1 -----
//    {
//        // Now that we have the maximum peak, search for the bending (foot of the curve) between the maxima in the histogram.
//        
//        double line_pos1_src[3], line_pos2_src[3];
//        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[max_value_index]; line_pos1_src[2] = max_value_index;
//        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 0;
//        double max_squared_distance = 0;
//        unsigned int index_of_max_distance = max_value_index-1;
//        for(int i=max_value_index-1; i>1; i--)// for(unsigned int i=index_of_max_value1-1; i>0; i--)
//        {
//            double position_to_project_src[3], projected_position_src[3], n;
//            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
//            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
//            if(n>0 && n<1)
//            {
//                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
//                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
//                
//                if(d>max_squared_distance)
//                {
//                    max_squared_distance = d;
//                    index_of_max_distance = i;
//                }
//            }
//        }
//        
//        // The threshold value is the index in of the crossing with x axis
//        threshold = index_of_max_distance;
//        
//        cout<<" threshold="<<threshold<<"  ";
//    }
//    //------------------------------
//
//    
//    
//    
//    
//
//    
//    //----- Erosion of the inner structures (Comment out if not needed) -----
//    {
//        cout<<"  Erosion of segmented inner bright structure...";
//        bdBIP bip;
//        bdImage temp;
//        bdGIP gip;
//        gip.Negative(output, temp);
//        //bip.ErosionByMapping_Circle(temp, segmented_melanoma, output, 1);
//        bip.ErosionByMapping_Circle2(temp, segmented_melanoma, output, 3, 60);
//        cout<<" completed."<<endl;
//    }
//    //----------
//    
//    
//    //----- Label connected components -----
//    {
//        unsigned short label = 1;
//        
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(output(r,c)!=0)
//                {
//                    if(output(r,c)==segmented_label)
//                    {
//                        //perform region growing
//                        bdList< bdVoxel > temp_list;
//                        bdVoxel v;
//                        temp_list.AddToLeftEnd(v(r,c));
//                        output(v.R(),v.C()) = label;
//                        
//                        bdGeometry g;
//                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//                        while(!temp_list.IsEmpty())
//                        {
//                            v = temp_list.GetLeftEnd();
//                            
//                            int rn, cn;
//                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
//                            {
//                                if(output(rn,cn)==segmented_label)
//                                {
//                                    bdVoxel vn;
//                                    temp_list.AddToRightEnd(vn(rn,cn));
//                                    output(vn.R(),vn.C()) = label;
//                                }
//                            }
//                            
//                            temp_list.DeleteLeftEnd();
//                        }
//                        
//                        label = (label+1) % segmented_label;
//                    }
//                }
//            }
//        }
//    }
//    //----------
//    
//    cout<<"  Extraction of network...";
//    
//    bdImage output2;
//    output2.SetSizeAndPropertiesAs(output);
//    output2.FillInWith(0);
//    
//    {
//        bdGeometry g;
//        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//        
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(segmented_melanoma(r,c)!=0)
//                {
//                    if(output(r,c)==0)
//                    {
//                        bdList<bdVoxel> seeds;
//                        bdList<unsigned short> CC_values;
//                        
//                        
//                        int rn, cn;
//                        for(g.ForCoordinates_Circle(r,c,((SE_radius*SE_radius)/4)); g.Get_Circle((SE_radius/2+1)*(SE_radius/2+1), rn, cn); )
//                        {
//                            if(output(rn,cn)!=0)
//                            {
//                                if(!bdListing::HasElement(CC_values, output(rn,cn)))
//                                {
//                                    CC_values.AddToRightEnd(output(rn,cn));
//                                    bdVoxel v;
//                                    seeds.AddToRightEnd(v(rn,cn));
//                                }                                
//                            }
//                        }
//                        
//                        if(CC_values.GetNumberOfElements()>=3)
//                        {
//                            bdGeometry g;
//                            g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//                            while(!seeds.IsEmpty())
//                            {
//                                //cout<<" * ";
//                                bdVoxel v = seeds.GetLeftEnd();
//                                
//                                bdList<bdVoxel> temp_list;
//                                temp_list.AddToRightEnd(v);
//                                
//                                int n_of_pixels_in_CC = 0;
//                                int sum_in_CC = 0;
//                                int n_of_pixels_in_border = 0;
//                                int sum_in_border = 0;
//                                int max_in_border = 0;
//
//                                bdList<bdVoxel> CC_list;
//                                
//                                while(!temp_list.IsEmpty())
//                                {
//                                    //cout<<" 2 ";
//                                    v = temp_list.GetLeftEnd();
//                                    
//                                    CC_list.AddToRightEnd(v);
//                                
//                                    int rn, cn;
//                                    for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
//                                    {
//                                        if(output(rn,cn)==0)
//                                        {
//                                            n_of_pixels_in_border++;
//                                            sum_in_border += input_img(t_input,0,rn,cn);
//                                            if(max_in_border<input_img(t_input,0,rn,cn)) max_in_border = input_img(t_input,0,rn,cn);
//                                        }
//                                        
//                                        if(output(rn,cn)!=0 && output2(rn,cn)==0)
//                                        {
//                                            bdVoxel vn;
//                                            temp_list.AddToRightEnd(vn(rn,cn));
//                                            output2(vn.R(),vn.C()) = 255;
//                                            
//                                            n_of_pixels_in_CC++;
//                                            sum_in_CC += input_img(t_input,0,rn,cn);
//                                        }
//                                    }
//                                    
//                                    temp_list.DeleteLeftEnd();
//                                }
//                                
//                                if(sum_in_CC/(n_of_pixels_in_CC+1) <= sum_in_border/(n_of_pixels_in_border+1) || (sum_in_border/(n_of_pixels_in_border+1))>threshold)//(max_value_index*8)/10)//100)//160)
//                                //if(!(sum_in_CC/(n_of_pixels_in_CC+1) > max_in_border))
//                                
//                                {
//                                    while(!CC_list.IsEmpty())
//                                    {
//                                        //cout<<" 3 ";
//
//                                        v = CC_list.GetLeftEnd();
//                                        
//                                        output2(v.R(),v.C()) = 0;
//                                        
//                                        CC_list.DeleteLeftEnd();
//                                    }
//                                }
//                                
//                                seeds.DeleteLeftEnd();
//                            }
//                        }
//
//                    }
//                }
//            }
//        }
//    }
//    cout<<" completed."<<endl;
//    
//    output.CopyFrom(output2);
//    
//    
//    return 1;
//}


int bdMelanoma::SegmentMelanomaInnerStructure8(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int SE_radius = 6;
    unsigned short segmented_label = 65535;
    unsigned int t_input = 0;
//    unsigned int threshold = 0;

    if(!this->SegmentStructures(input_img, t_input, segmented_melanoma, (SE_radius*SE_radius), 70, output, segmented_label)) return 0; //if(!this->SegmentStructures(input_img, t_input, segmented_melanoma, (SE_radius*SE_radius), 60, output, segmented_label)) return 0;

    bdArray<int> histogram;
    this->Histogram(input_img, segmented_melanoma, histogram, 0, 255, t_input, 0, 0, 0); //this->Histogram(input_img, input_img, histogram, 0, 255, t_input, 0, 0, 0);
    int max_value = histogram[0];
    int max_value_index = 0;
    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
    {
        if(max_value<histogram[i])
        {
            max_value = histogram[i];
            max_value_index = i;
        }
    }
    cout<<" max_value_index="<<max_value_index<<"  ";



    //----- Erosion of the inner structures (Comment out if not needed) -----
    {
        cout<<"  Erosion of segmented inner bright structure...";
        bdBIP bip;
        bdImage temp;
        bdGIP gip;
        gip.Negative(output, temp);
        //bip.ErosionByMapping_Circle(temp, segmented_melanoma, output, 1);
        bip.ErosionByMapping_Circle2(temp, segmented_melanoma, output, 3, 60);
        cout<<" completed."<<endl;
    }
    //----------


//    //----- Label connected components -----
//    {
//        unsigned short label = 1;
//
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(output(r,c)!=0)
//                {
//                    if(output(r,c)==segmented_label)
//                    {
//                        //perform region growing
//                        bdList< bdVoxel > temp_list;
//                        bdVoxel v;
//                        temp_list.AddToLeftEnd(v(r,c));
//                        output(v.R(),v.C()) = label;
//
//                        bdGeometry g;
//                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//                        while(!temp_list.IsEmpty())
//                        {
//                            v = temp_list.GetLeftEnd();
//
//                            int rn, cn;
//                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
//                            {
//                                if(output(rn,cn)==segmented_label)
//                                {
//                                    bdVoxel vn;
//                                    temp_list.AddToRightEnd(vn(rn,cn));
//                                    output(vn.R(),vn.C()) = label;
//                                }
//                            }
//
//                            temp_list.DeleteLeftEnd();
//                        }
//
//                        label = (label+1) % segmented_label;
//                    }
//                }
//            }
//        }
//    }
//    //----------
    
    
    //----- Label connected components -----
    {
        unsigned short label = 1;
        
        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                if(output(r,c)!=0)
                {
                    if(output(r,c)==segmented_label)
                    {
                        //perform region growing
                        bdList< bdVoxel > temp_list;
                        bdVoxel v;
                        temp_list.AddToLeftEnd(v(r,c));
                        output(v.R(),v.C()) = label;
                        
                        bdList< bdVoxel > CC_list;
                        CC_list.AddToLeftEnd(v);
                        
                        bdVoxel CC_min_r, CC_min_c, CC_max_r, CC_max_c;// borders of the CC (left, right, bottom and top of the rectangular border)
                        CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
                        CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
                        CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
                        CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
                        
                        
                        
                        bdGeometry g;
                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
                        while(!temp_list.IsEmpty())
                        {
                            v = temp_list.GetLeftEnd();
                            
                            
                            // Update rectangular border
                            if(v.R() > CC_max_r.R())
                            {
                                CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
                            }
                            if(v.C() > CC_max_c.C())
                            {
                                CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
                            }
                            if(v.R() < CC_min_r.R())
                            {
                                CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
                            }
                            if(v.C() < CC_min_c.C())
                            {
                                CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
                            }
                            
                            int rn, cn;
                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
                            {
                                if(output(rn,cn)==segmented_label)
                                {
                                    bdVoxel vn;
                                    temp_list.AddToRightEnd(vn(rn,cn));
                                    output(vn.R(),vn.C()) = label;
                                    
                                    CC_list.AddToLeftEnd(vn);
                                }
                            }
                            
                            temp_list.DeleteLeftEnd();
                        }
                        
                        // Now the CC is labeled with value 'label' and is contained in CC_list,
                        // determine the area of rectangular ROI and compare to area of CC.
                        unsigned int rectangular_border_side_r = (CC_max_r.R()-CC_min_r.R());
                        unsigned int rectangular_border_side_c = (CC_max_c.C()-CC_min_c.C());
                        unsigned int rectangular_border_area;
                        if(rectangular_border_side_r>rectangular_border_side_c)
                        {
                            rectangular_border_area = rectangular_border_side_r * rectangular_border_side_r;
                        }
                        else
                        {
                            rectangular_border_area = rectangular_border_side_c * rectangular_border_side_c;
                        }
                        
                        // If the CC is not long circular/rectangular enough, erase it.
                        if(CC_list.GetNumberOfElements() < (rectangular_border_area*60)/100) // 3.14/4 = 0.785
                        {
                            while(!CC_list.IsEmpty())
                            {
                                v = CC_list.GetLeftEnd();
                                
                                output(v.R(),v.C()) = 0;
                                
                                CC_list.DeleteLeftEnd();
                            }
                        }
                        
                        label = (label+1) % segmented_label;
                    }
                }
            }
        }
    }
    //----------

    
    
    //----- Extract network -----
    cout<<"  Extraction of network...";

    bdImage output2;
    output2.SetSizeAndPropertiesAs(output);
    output2.FillInWith(0);

    {
        bdGeometry g;
        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());

        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(r,c)!=0)
                {
                    if(output(r,c)==0)
                    {
                        bdList<bdVoxel> seeds;
                        bdList<unsigned short> CC_values;


                        int rn, cn;
                        for(g.ForCoordinates_Circle(r,c,((SE_radius*SE_radius)/4)); g.Get_Circle((SE_radius/2+1)*(SE_radius/2+1), rn, cn); )
                        {
                            if(output(rn,cn)!=0)
                            {
                                if(!bdListing::HasElement(CC_values, output(rn,cn)))
                                {
                                    CC_values.AddToRightEnd(output(rn,cn));
                                    bdVoxel v;
                                    seeds.AddToRightEnd(v(rn,cn));
                                }
                            }
                        }

                        if(CC_values.GetNumberOfElements()>=3)
                        {
                            bdGeometry g;
                            g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
                            while(!seeds.IsEmpty())
                            {
                                //cout<<" * ";
                                bdVoxel v = seeds.GetLeftEnd();

                                bdList<bdVoxel> temp_list;
                                temp_list.AddToRightEnd(v);

                                int n_of_pixels_in_CC = 0;
                                int sum_in_CC = 0;
                                int n_of_pixels_in_border = 0;
                                int sum_in_border = 0;
                                int max_in_border = 0;

                                bdList<bdVoxel> CC_list;

                                while(!temp_list.IsEmpty())
                                {
                                    //cout<<" 2 ";
                                    v = temp_list.GetLeftEnd();

                                    CC_list.AddToRightEnd(v);

                                    int rn, cn;
                                    for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
                                    {
                                        if(output(rn,cn)==0)
                                        {
                                            n_of_pixels_in_border++;
                                            sum_in_border += input_img(t_input,0,rn,cn);
                                            if(max_in_border<input_img(t_input,0,rn,cn)) max_in_border = input_img(t_input,0,rn,cn);
                                        }

                                        if(output(rn,cn)!=0 && output2(rn,cn)==0)
                                        {
                                            bdVoxel vn;
                                            temp_list.AddToRightEnd(vn(rn,cn));
                                            output2(vn.R(),vn.C()) = 255;

                                            n_of_pixels_in_CC++;
                                            sum_in_CC += input_img(t_input,0,rn,cn);
                                        }
                                    }

                                    temp_list.DeleteLeftEnd();
                                }

////!!!!!!!!!!!!! THIS PARAMETER IS CURRENTLY SET AD-HOC !!!!!!!!!!!!!!!!!!!
//                                if(sum_in_CC/(n_of_pixels_in_CC+1) <= sum_in_border/(n_of_pixels_in_border+1) || (sum_in_border/(n_of_pixels_in_border+1))>(max_value_index*85)/100)//threshold)//100)//160)
//                                //if(!(sum_in_CC/(n_of_pixels_in_CC+1) > max_in_border))
//                                //if(sum_in_CC/(n_of_pixels_in_CC+1) - sum_in_border/(n_of_pixels_in_border+1) <= 4)
//                                {
//                                    while(!CC_list.IsEmpty())
//                                    {
//                                        //cout<<" 3 ";
//
//                                        v = CC_list.GetLeftEnd();
//
//                                        output2(v.R(),v.C()) = 0;
//
//                                        CC_list.DeleteLeftEnd();
//                                    }
//                                }

                                seeds.DeleteLeftEnd();
                            }
                        }

                    }
                }
            }
        }
    }
    cout<<" completed."<<endl;

    output.CopyFrom(output2);
    //----------


    return 1;
}





int bdMelanoma::SegmentMelanomaInnerDarkStructure(bdImage &input, bdImage &segmented_melanoma, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    if(segmented_melanoma.IsEmpty()) return 0;
    
    bdGIP gip;
    bdImage original;// negative of the input;
    gip.Negative(input,original);

    cout<<" 1 ";
    
    
    // Get range of values found in the masked region of the image
    unsigned int range_min = 65535, range_max = 0;
    for(unsigned int t=0; t<original.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int r=0; r< original.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
            {
                if(segmented_melanoma(r,c)!=0)
                {
                    if(original(t,0,r,c)<range_min) range_min = original(t,0,r,c);
                    if(original(t,0,r,c)>range_max) range_max = original(t,0,r,c);
                }
            }
        }
    }
    
    cout<<"  range="<<range_min<<","<<range_max<<"  ";
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVER-RIDDEN HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    range_min = 0;

    // Calculate the histograms of the masked area for each of the components (R,G,B)
    bdArray< bdArray<int> > histograms_per_RGB;
    histograms_per_RGB.Set(original.GetNumberOfTimeSeries());
    //unsigned int range_min = 0, range_max = 65535;
    for(unsigned int t=0; t<original.GetNumberOfTimeSeries(); t++)
    {
        this->Histogram(original,segmented_melanoma,histograms_per_RGB[t],range_min,range_max,t,0,0,0);
    }
    
    
    // For each histogram calculate variance and find the index of the max variance
    double max_variance = -1;//initial non-valid value.
    unsigned int max_variance_index = 0;
    
    for(unsigned int t=0; t<original.GetNumberOfTimeSeries(); t++)
    {
        // calculate mean
        double mean = 0, sum = 0, n=0;
        for(unsigned int i=1; i<histograms_per_RGB[t].GetNumberOfElements(); i++)
        {
            sum += histograms_per_RGB[t][i] * i;
            n += histograms_per_RGB[t][i];
        }
        mean = sum / n;
        
        // calculate variance
        double variance = 0, sum_for_variance = 0;
        for(unsigned int i=1; i<histograms_per_RGB[t].GetNumberOfElements(); i++)
        {
            sum_for_variance += (i - mean) * (i - mean) * histograms_per_RGB[t][i];
        }
        variance = sum_for_variance / n;
        
        cout<<"["<<t<<"],variance="<<variance<<",mean="<<mean<<"  ";
        
        if(variance>max_variance || max_variance<0)
        {
            max_variance_index = t;
            max_variance = variance;
        }
    }
 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVER-RIDDEN HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //max_variance_index = 0;//for debug purposes
    
    cout<<" 2 ";
    
    // We continue to work on the slice that has the highest variance, so extract the slice
    bdImage extracted_slice;
    extracted_slice.SetSize(1,1, original.GetNumberOfRows(), original.GetNumberOfColumns());
    extracted_slice.SetVisualizationPropertiesToMatchInput(original);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            extracted_slice(0,0,r,c) = original(max_variance_index,0,r,c);
        }
    }
    
    cout<<" 2a ";
    
    
    // The algorithms which will be used work by extracting rigdes (local maxima), so we will perform
    // some morphology to make the regions more uniform (using median) and smooth
    /// the transition between the regions (mean).
    bdImage extracted_slice_max;
    gip.MaximumCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 25,extracted_slice_max);
    
    cout<<" 2b ";
    
    bdImage slice;
    gip.MeanCircleFor2DSlice(extracted_slice_max, segmented_melanoma, 0,0,0,0, 36,slice);
    
    
    cout<<" 3 ";
    
    
    // Calculate the histograms of the masked area in the morphologically modified image for a single component
    bdArray<int> histogram_of_slice;
    this->Histogram(slice,segmented_melanoma,histogram_of_slice,range_min,range_max,0,0,0,0);
    
    
    // create array of sorted voxel values. This will be used to access voxels in order of their values.
    bdArray< bdArray< bdVoxel > > sorted_array;
    //sorted_array.Set(histograms_per_RGB[max_variance_index].GetNumberOfElements());
    sorted_array.Set(histogram_of_slice.GetNumberOfElements());

    //for(unsigned int i=0; i<histograms_per_RGB[max_variance_index].GetNumberOfElements(); i++)
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        //if(histograms_per_RGB[max_variance_index][i]!=0) sorted_array[i].Set(histograms_per_RGB[max_variance_index][i]);
        if(histogram_of_slice[i]!=0) sorted_array[i].Set(histogram_of_slice[i]);

    }
    

    // Put the masked voxels from the original image to the sorted_array
    bdArray<int> index_array;
    index_array.Set(sorted_array.GetNumberOfElements());
    index_array.FillInWith(0);
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0 && slice(r,c)!=0)
            {
                sorted_array[slice(0,0,r,c)] [(index_array[slice(0,0,r,c)])] (r,c);
                (index_array[slice(0,0,r,c)])++;
            }
        }
    }
    
    
    cout<<" 4 ";
    
    
    output.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(slice);
    output.FillInWith(0);
    
    bdImage temp;
    temp.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    temp.SetVisualizationPropertiesToMatchInput(slice);
    temp.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    
    
    unsigned int min, max;
    slice.GetVoxelMinimumAndMaximumValue(&min, &max);
    
    unsigned int t=0, s=0;

    unsigned short index_of_new_CCs = 53;//we start AT LEAST with value 3 because value 1 is reserved for a new voxel not yet assigned to a region and 2 is the region boundary.
    for(int threshold = max; threshold>min; threshold--)
    {
        //cout<<"thr="<<threshold<<",noel="<<sorted_array[threshold].GetNumberOfElements()<<"  ";
        bdList<bdVoxel> added_pixels;
        //for(unsigned int r=0; r<temp.GetNumberOfRows(); r++)
        //{
            //for(unsigned int c=0; c<temp.GetNumberOfColumns(); c++)
            for(unsigned int j=0; j<sorted_array[threshold].GetNumberOfElements(); j++)
            {
                unsigned int r = sorted_array[threshold][j].R();
                unsigned int c = sorted_array[threshold][j].C();

                if(temp(t,s,r,c)==0)
                {
                    if(slice(t,s,r,c)==threshold)//if(slice(t,s,r,c)>=threshold)
                    {
                        int rn,cn;
                        int is_previous_region_found = 0;
                        //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(100,rn,cn); )
                        //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(9,rn,cn); )
                        //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(25,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                        //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(16,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                        //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(36,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                        for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(81,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
  
                        {
                            if(temp(t,s,rn,cn)!=0)
                            {
                                temp(t,s,r,c) = temp(t,s,rn,cn);
                                is_previous_region_found = 1;
                                break;
                            }
                        }
                        if(!is_previous_region_found)
                        {
                            temp(t,s,r,c) = index_of_new_CCs;
                            index_of_new_CCs++;
                        }
                        bdVoxel v; v(r,c);
                        added_pixels.AddToRightEnd(v);
                    }
                }
            }
        //}
        
        //Go through the list of added pixels and find pixels that border with pixels of other value (that is greater than 1), and make borders for these regions.
        //bdListIterator<bdVoxel> it;
        while(!added_pixels.IsEmpty())
        {
            bdVoxel *pv = added_pixels.GetLeftEndPointerToElement();
            unsigned short region_value = temp(t,s,pv->R(),pv->C());
            
            //output(0,0,pv->R(),pv->C()) = region_value;
            
            //check if it is a border with other region...
            if(region_value>1)
            {
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(pv->R(),pv->C()); g.Get_8_Neighbors(rn,cn); )
                {
                    if(temp(t,s,rn,cn)!=region_value && temp(t,s,rn,cn)>0)
                    {
                        //... if it is a border, grow a bound for the region (convert its edges into bounds)...
                        bdRegion2D reg;//, edge;
                        //reg.CreateRegion_8_FromSeedAndValueHigherThan(temp,threshold,pv->C(),pv->R(),0,0);
                        reg.CreateRegion_8_FromSeedPoint(temp,region_value,t,s,pv->R(),pv->C());
                        bdRegion2DIterator itr;
                        for(itr.SetBegin(&reg); itr.IsValid(); itr.MoveToNext())
                        {
                            output(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = region_value;
                            temp(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = 1;
                        }
                        
                        break;
                    }
                }
            }
            
            added_pixels.DeleteLeftEnd();
        }
    }
    
    
    //output.CopyFrom(extracted_slice);
    
    return 1;
}



int bdMelanoma::SegmentMelanomaInnerDarkStructure2(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
{
    if(input_img.IsEmpty()) return 0;
    if(segmented_melanoma.IsEmpty()) return 0;
    
    bdGIP gip;
    
    bdImage input;
    gip.RescaleWholeRange(input_img, 0,255, input);
    
    bdImage original;// negative of the input;
    gip.Negative(input,original);
    
    cout<<" 1 ";
    
    
    // We continue to work on the slice which is a combination of the RGB cahnnels.
    bdImage extracted_slice;
    extracted_slice.SetSize(1,1, original.GetNumberOfRows(), original.GetNumberOfColumns());
    extracted_slice.SetVisualizationPropertiesToMatchInput(original);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
//            // compute combined gray image as AVERAGE of RGB channels.
//            extracted_slice(0,0,r,c) = (original(0,0,r,c) + original(1,0,r,c) + original(2,0,r,c)) / 3;
            
            // compute gray image as only RED channel.
            extracted_slice(0,0,r,c) = original(0,0,r,c);


//            // MEDIAN OF ALL CHANNELS
//            if(input(1,0,r,c)>input(2,0,r,c))
//            {
//                if(input(2,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(2,0,r,c);
//                else extracted_slice(r,c) = input(0,0,r,c);
//            }
//            else
//            {
//                if(input(1,0,r,c)>input(0,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
//                else extracted_slice(r,c) = input(0,0,r,c);
//            }
//
//            // MAX OF ALL CHANNELS
//            if(input(1,0,r,c)>input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
//            else extracted_slice(r,c) = input(2,0,r,c);
//            if(input(0,0,r,c)>extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);
//
//            // MIDRANGE OF ALL CHANNELS
//            unsigned int min;
//            if(input(1,0,r,c)<input(2,0,r,c)) min = input(1,0,r,c);
//            else min = input(2,0,r,c);
//            if(input(0,0,r,c)<min) min = input(0,0,r,c);
//            extracted_slice(r,c) = (extracted_slice(r,c) + min) / 2;

            
//            // MIN OF ALL CHANNELS
//            if(input(1,0,r,c)<input(2,0,r,c)) extracted_slice(r,c) = input(1,0,r,c);
//            else extracted_slice(r,c) = input(2,0,r,c);
//            if(input(0,0,r,c)<extracted_slice(r,c)) extracted_slice(r,c) = input(0,0,r,c);


        }
    }
    
    cout<<" 2a ";
    
    
    // The algorithms which will be used work by extracting rigdes (local maxima), so we will perform
    // some morphology to make the regions more uniform (using median) and smooth
    /// the transition between the regions (mean).
//    bdImage extracted_slice_max;
//    gip.MaximumCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 36,extracted_slice_max);//25,extracted_slice_max);
    
    cout<<" 2b ";
    
    bdImage slice;
//    gip.MeanCircleFor2DSlice(extracted_slice_max, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);
    
    //gip.MeanCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);

    
    
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVER-RIDDEN HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    slice.CopyFrom(extracted_slice);
    
    cout<<" 3 ";
    
    // Get range of values found in the masked region of the image
    unsigned int gray_value_max = 0;
    unsigned int gray_value_min = 65535;
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0)
            {
                if(slice(r,c)>gray_value_max) gray_value_max = slice(r,c);
                if(slice(r,c)<gray_value_min) gray_value_min = slice(r,c);
            }
        }
    }

    
    cout<<"  gray_value_max="<<gray_value_max<<", gray_value_min="<<gray_value_min<<"  ";


    // Calculate the histograms of the masked area in the morphologically modified image for a single component
    bdArray<int> histogram_of_slice;
    this->Histogram(slice,segmented_melanoma,histogram_of_slice,0,gray_value_max,0,0,0,0);
    
    
    // create array of sorted voxel values. This will be used to access voxels in order of their values.
    bdArray< bdArray< bdVoxel > > sorted_array;
    sorted_array.Set(histogram_of_slice.GetNumberOfElements());
    
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        if(histogram_of_slice[i]!=0) sorted_array[i].Set(histogram_of_slice[i]);
    }
    
    
    // Put the masked voxels from the original image to the sorted_array
    bdArray<int> index_array;
    index_array.Set(sorted_array.GetNumberOfElements());
    index_array.FillInWith(0);
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0 && slice(r,c)!=0)
            {
                sorted_array[slice(0,0,r,c)] [(index_array[slice(0,0,r,c)])] (r,c);
                (index_array[slice(0,0,r,c)])++;
            }
        }
    }
    
    
    cout<<" 4 ";
    
    
    output.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(slice);
    output.FillInWith(0);
    
    bdImage temp;
    temp.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    temp.SetVisualizationPropertiesToMatchInput(slice);
    temp.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    
    
 //   unsigned int min, max;
 //   slice.GetVoxelMinimumAndMaximumValue(&min, &max);
    
    unsigned int t=0, s=0;
    
    unsigned short index_of_new_CCs = 53;//we start AT LEAST with value 3 because value 1 is reserved for a new voxel not yet assigned to a region and 2 is the region boundary.
//    for(int threshold = max; threshold>min; threshold--)
    for(int threshold = gray_value_max; threshold>gray_value_min; threshold--)

    {
        //cout<<"thr="<<threshold<<",noel="<<sorted_array[threshold].GetNumberOfElements()<<"  ";
        bdList<bdVoxel> added_pixels;
        //for(unsigned int r=0; r<temp.GetNumberOfRows(); r++)
        //{
        //for(unsigned int c=0; c<temp.GetNumberOfColumns(); c++)
        for(unsigned int j=0; j<sorted_array[threshold].GetNumberOfElements(); j++)
        {
            unsigned int r = sorted_array[threshold][j].R();
            unsigned int c = sorted_array[threshold][j].C();
            
            if(temp(t,s,r,c)==0)
            {
                if(slice(t,s,r,c)==threshold)//if(slice(t,s,r,c)>=threshold)
                {
                    int rn,cn;
                    int is_previous_region_found = 0;
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(100,rn,cn); )
                    for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(9,rn,cn); )
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(25,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(16,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(36,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                    //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(81,rn,cn); )// THIS SHOULD BE AN INPUT PARAMETER !!!!!!!!!!!!!!!!!!!!!
                        
                    {
                        if(temp(t,s,rn,cn)!=0)
                        {
                            temp(t,s,r,c) = temp(t,s,rn,cn);
                            is_previous_region_found = 1;
                            break;
                        }
                    }
                    if(!is_previous_region_found)
                    {
                        temp(t,s,r,c) = index_of_new_CCs;
                        index_of_new_CCs++;
                    }
                    bdVoxel v; v(r,c);
                    added_pixels.AddToRightEnd(v);
                }
            }
        }
        //}
        
        //Go through the list of added pixels and find pixels that border with pixels of other value (that is greater than 1), and make borders for these regions.
        //bdListIterator<bdVoxel> it;
        while(!added_pixels.IsEmpty())
        {
            bdVoxel *pv = added_pixels.GetLeftEndPointerToElement();
            unsigned short region_value = temp(t,s,pv->R(),pv->C());
            
            //output(0,0,pv->R(),pv->C()) = region_value;
            
            //check if it is a border with other region...
            if(region_value>1)
            {
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(pv->R(),pv->C()); g.Get_8_Neighbors(rn,cn); )
                {
                    if(temp(t,s,rn,cn)!=region_value && temp(t,s,rn,cn)>0)
                    {
                        //... if it is a border, grow a bound for the region (convert its edges into bounds)...
                        bdRegion2D reg;//, edge;
                        //reg.CreateRegion_8_FromSeedAndValueHigherThan(temp,threshold,pv->C(),pv->R(),0,0);
                        reg.CreateRegion_8_FromSeedPoint(temp,region_value,t,s,pv->R(),pv->C());
                        bdRegion2DIterator itr;
                        for(itr.SetBegin(&reg); itr.IsValid(); itr.MoveToNext())
                        {
                            output(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = region_value;
                            temp(t,s,itr.GetIndexRow(),itr.GetIndexColumn()) = 1;
                        }
                        
                        break;
                    }
                }
            }
            
            added_pixels.DeleteLeftEnd();
        }
    }
    
    
    //output.CopyFrom(extracted_slice);
    
    return 1;
}


int bdMelanoma::SegmentMelanomaInnerDarkStructure3(bdImage &input_img, bdImage &segmented_melanoma, bdImage &output)
{
    unsigned int SE_radius = 6;
    unsigned short segmented_label = 65535;
    unsigned int t_input = 0;
    
    if(!this->SegmentStructures(input_img, t_input, segmented_melanoma, (SE_radius*SE_radius), 40, output, segmented_label)) return 0;
    
    
    
//    // We continue to work on the slice which is a grayscale combination of green and red channel.
//    bdImage slice;
//    slice.SetSize(1,1, input_img.GetNumberOfRows(), input_img.GetNumberOfColumns());
//    slice.SetVisualizationPropertiesToMatchInput(input_img);
//    for(unsigned int r=0; r<input_img.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<input_img.GetNumberOfColumns(); c++)
//        {
//            slice(r,c) = (input_img(0,0,r,c) + input_img(0,0,r,c))/2;
//        }
//    }
//    
//    if(!this->SegmentStructures(slice, 0, segmented_melanoma, (SE_radius*SE_radius), 40, output, segmented_label)) return 0;
    
    
    
    //----- Erosion of the inner dark structures (Comment out if not needed) -----
    {
        cout<<"Erosion of segmented inner dark structure...";
        bdBIP bip;
        bdImage temp;
        temp.CopyFrom(output);
        //bip.ErosionByMapping_Circle(temp, segmented_melanoma, output, 2);
        bip.ErosionByMapping_Circle2(temp, segmented_melanoma, output, 2, 95);
        cout<<" completed."<<endl;
    }
    //----------

    


    
    //----- Label connected components -----
    {
        //unsigned short label = 1;
        unsigned short label_globules = 255;
        unsigned short label_dots = 127;

        
        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                if(output(r,c)!=0)
                {
                    if(output(r,c)==segmented_label)
                    {
                        //perform region growing
                        bdList< bdVoxel > temp_list;
                        bdVoxel v;
                        temp_list.AddToLeftEnd(v(r,c));
                        output(v.R(),v.C()) = label_globules;//label;
                        
                        bdList< bdVoxel > CC_list;
                        CC_list.AddToLeftEnd(v);
                        
                        unsigned short min_value = 65535;//, max_value = 0;
                        
                        bdGeometry g;
                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
                        while(!temp_list.IsEmpty())
                        {
                            v = temp_list.GetLeftEnd();
                            
                            //if(input_img(t_input,0,v.R(),v.C())>max_value) max_value = input_img(t_input,0,v.R(),v.C());
                            unsigned int gray_value = ( input_img(0,0,v.R(),v.C())+input_img(1,0,v.R(),v.C())+input_img(2,0,v.R(),v.C()) )/3;
                            if(gray_value<min_value) min_value = gray_value;
//                            if(input_img(0,0,v.R(),v.C())<min_value) min_value = input_img(0,0,v.R(),v.C());
//                            if(input_img(1,0,v.R(),v.C())<min_value) min_value = input_img(1,0,v.R(),v.C());
//                            if(input_img(2,0,v.R(),v.C())<min_value) min_value = input_img(2,0,v.R(),v.C());

                            
                            int rn, cn;
                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
                            {
                                if(output(rn,cn)==segmented_label)
                                {
                                    bdVoxel vn;
                                    temp_list.AddToRightEnd(vn(rn,cn));
                                    output(vn.R(),vn.C()) = label_globules;//label;
                                    
                                    CC_list.AddToLeftEnd(vn);
                                }
                            }
                            
                            temp_list.DeleteLeftEnd();
                        }
                        
                        // Now the CC is labeled with value 'label' and is contained in CC_list.
                        
                        
                        
                        
                        bdList< bdVoxel > CC_list2;//to record the filtered region.
                        
                        {
                            while(!CC_list.IsEmpty())
                            {
                                v = CC_list.GetLeftEnd();
                                
                                //unsigned int gray_value = ( input_img(0,0,v.R(),v.C())+input_img(1,0,v.R(),v.C())+input_img(2,0,v.R(),v.C()) )/3;
                                unsigned int gray_value = ( input_img(0,0,v.R(),v.C())+input_img(1,0,v.R(),v.C()) ) /2;
                                
//!!!!!!!!!!!!! PARAMETERS SET AD-HOC !!!!!!!!!!!!!!!!!!!!!!!!!!
                                if(gray_value>100) output(v.R(),v.C()) = 0;
                                else
                                {
                                    CC_list2.AddToRightEnd(v);
                                    //if(CC_n_of_pixels<13) output(v.R(),v.C()) = label_dots;
                                }
                                
                                CC_list.DeleteLeftEnd();
                            }
                        }
                        
                        
                        // filter by the area to determine if the structures are dots or globules
                        double CC_n_of_pixels = CC_list2.GetNumberOfElements();
                        while(!CC_list2.IsEmpty())
                        {
                            v = CC_list2.GetLeftEnd();
                            
                            if(CC_n_of_pixels<13) output(v.R(),v.C()) = label_dots;
                            
                            CC_list2.DeleteLeftEnd();
                        }
                        
                        //label = (label+1) % segmented_label;
                    }
                }
            }
        }
    }
    //----------
    
    
    
    
    
    
//    //----- Label connected components -----
//    {
//        unsigned short label = 1;
//        
//        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//            {
//                if(output(r,c)!=0)
//                {
//                    if(output(r,c)==segmented_label)
//                    {
//                        //perform region growing
//                        bdList< bdVoxel > temp_list;
//                        bdVoxel v;
//                        temp_list.AddToLeftEnd(v(r,c));
//                        output(v.R(),v.C()) = label;
//                        
//                        bdList< bdVoxel > CC_list;
//                        CC_list.AddToLeftEnd(v);
//                        
//                        bdVoxel CC_min_r, CC_min_c, CC_max_r, CC_max_c;// borders of the CC (left, right, bottom and top of the rectangular border)
//                        CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
//                        CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
//                        CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
//                        CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
//                        
//                        
//                        
//                        bdGeometry g;
//                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
//                        while(!temp_list.IsEmpty())
//                        {
//                            v = temp_list.GetLeftEnd();
//                            
//                            
//                            // Update rectangular border
//                            if(v.R() > CC_max_r.R())
//                            {
//                                CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
//                            }
//                            if(v.C() > CC_max_c.C())
//                            {
//                                CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
//                            }
//                            if(v.R() < CC_min_r.R())
//                            {
//                                CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
//                            }
//                            if(v.C() < CC_min_c.C())
//                            {
//                                CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
//                            }
//                            
//                            int rn, cn;
//                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
//                            {
//                                if(output(rn,cn)==segmented_label)
//                                {
//                                    bdVoxel vn;
//                                    temp_list.AddToRightEnd(vn(rn,cn));
//                                    output(vn.R(),vn.C()) = label;
//                                    
//                                    CC_list.AddToLeftEnd(vn);
//                                }
//                            }
//                            
//                            temp_list.DeleteLeftEnd();
//                        }
//                        
//                        // Now the CC is labeled with value 'label' and is contained in CC_list,
//                        // determine the area of rectangular ROI and compare to area of CC.
//                        unsigned int rectangular_border_side_r = (CC_max_r.R()-CC_min_r.R());
//                        unsigned int rectangular_border_side_c = (CC_max_c.C()-CC_min_c.C());
//                        unsigned int rectangular_border_area;
//                        if(rectangular_border_side_r>rectangular_border_side_c)
//                        {
//                            rectangular_border_area = rectangular_border_side_r * rectangular_border_side_r;
//                        }
//                        else
//                        {
//                            rectangular_border_area = rectangular_border_side_c * rectangular_border_side_c;
//                        }
//                        
//                        // If the CC is not long circular/rectangular enough, erase it.
//                        if(CC_list.GetNumberOfElements() < (rectangular_border_area*60)/100) // 3.14/4 = 0.785
//                        {
//                            while(!CC_list.IsEmpty())
//                            {
//                                v = CC_list.GetLeftEnd();
//                                
//                                output(v.R(),v.C()) = 0;
//                                
//                                CC_list.DeleteLeftEnd();
//                            }
//                        }
//                        
//                        label = (label+1) % segmented_label;
//                    }
//                }
//            }
//        }
//    }
//    //----------
    

    
    
    
    return 1;
}




int bdMelanoma::SegmentStructures(bdImage &original, unsigned int t_original, bdImage &mask, unsigned int SE_squared_radius, unsigned int percentage, bdImage &output, unsigned short segmented_label)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    bdGIP gip;
    
    bdImage input;
    gip.RescaleWholeRange(original, 0,255, input);
    
    bdImage input_negative;// negative of the input;
    gip.Negative(input,input_negative);
    
    cout<<" 1 ";
    
    
    // We continue to work on the slice with index 't'
    bdImage extracted_slice;
    extracted_slice.SetSize(1,1, original.GetNumberOfRows(), original.GetNumberOfColumns());
    extracted_slice.SetVisualizationPropertiesToMatchInput(original);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            extracted_slice(0,0,r,c) = input_negative(t_original,0,r,c);
        }
    }
    
    cout<<" 2a ";
    
    
    // The algorithms which will be used work by extracting rigdes (local maxima), so we will perform
    // some morphology to make the regions more uniform (using median) and smooth
    /// the transition between the regions (mean).
    //    bdImage extracted_slice_max;
    //    gip.MaximumCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 36,extracted_slice_max);//25,extracted_slice_max);
    
    //cout<<" 2b ";
    
    bdImage slice;
    //    gip.MeanCircleFor2DSlice(extracted_slice_max, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);
    
    //gip.MeanCircleFor2DSlice(extracted_slice, segmented_melanoma, 0,0,0,0, 49,slice);//36,slice);
    
    
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVER-RIDDEN HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    slice.CopyFrom(extracted_slice);
    
    cout<<" 3 ";
    
    // Get range of values found in the masked region of the image
    unsigned int gray_value_max = 0;
    unsigned int gray_value_min = 65535;
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                if(slice(r,c)>gray_value_max) gray_value_max = slice(r,c);
                if(slice(r,c)<gray_value_min) gray_value_min = slice(r,c);
            }
        }
    }
    
    
    cout<<"  gray_value_max="<<gray_value_max<<", gray_value_min="<<gray_value_min<<"  ";
    
    
    // Calculate the histograms of the masked area in the morphologically modified image for a single component
    bdArray<int> histogram_of_slice;
    this->Histogram(slice,mask,histogram_of_slice,0,gray_value_max,0,0,0,0);

    
    
    // create array of sorted voxel values. This will be used to access voxels in order of their values.
    bdArray< bdArray< bdVoxel > > sorted_array;
    sorted_array.Set(histogram_of_slice.GetNumberOfElements());
    
    for(unsigned int i=0; i<histogram_of_slice.GetNumberOfElements(); i++)
    {
        if(histogram_of_slice[i]!=0) sorted_array[i].Set(histogram_of_slice[i]);
    }
    
    
    // Put the masked voxels from the original image to the sorted_array
    bdArray<int> index_array;
    index_array.Set(sorted_array.GetNumberOfElements());
    index_array.FillInWith(0);
    for(unsigned int r=0; r< slice.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<slice.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0 && slice(r,c)!=0)
            {
                sorted_array[slice(0,0,r,c)] [(index_array[slice(0,0,r,c)])] (r,c);
                (index_array[slice(0,0,r,c)])++;
            }
        }
    }
    
    
    cout<<" 4 ";
    
    
    output.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(slice);
    output.FillInWith(0);
    
    bdImage temp;
    temp.SetSize(slice.GetNumberOfTimeSeries(),slice.GetNumberOfSlices(),slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    temp.SetVisualizationPropertiesToMatchInput(slice);
    temp.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(slice.GetNumberOfRows(),slice.GetNumberOfColumns());
    
    
    unsigned int t=0, s=0;
    
    unsigned short index_of_new_CCs = 53;//we start AT LEAST with value 3 because value 1 is reserved for a new voxel not yet assigned to a region and 2 is the region boundary.
    
    //unsigned int SE_squared_radius = 16;
    
    unsigned int n_of_pixels_in_SE = 0;
    int rn,cn;
    for(g.ForCoordinates_Circle(slice.GetNumberOfRows()/2,slice.GetNumberOfColumns()/2,0); g.Get_Circle(SE_squared_radius,rn,cn); )
    {
        n_of_pixels_in_SE++;
    }
    
    int step = 1;
    for(int threshold_main = gray_value_max; threshold_main>gray_value_min; threshold_main = threshold_main-step)
    {
        bdList<bdVoxel> added_pixels;
        
        for(int threshold = threshold_main; threshold>gray_value_min && threshold>threshold_main-step; threshold--)
        {
            for(unsigned int j=0; j<sorted_array[threshold].GetNumberOfElements(); j++)
            {
                unsigned int r = sorted_array[threshold][j].R();
                unsigned int c = sorted_array[threshold][j].C();
                
                if(temp(t,s,r,c)==0)
                {
                    if(slice(t,s,r,c)==threshold)
                    {
                        temp(t,s,r,c) = index_of_new_CCs;
                        bdVoxel v; v(r,c);
                        added_pixels.AddToRightEnd(v);
                    }
                }
            }
        }
        
        
        
        //Go through the list of added pixels and find pixels that border with pixels of other value (that is greater than 1), and make borders for these regions.
        while(!added_pixels.IsEmpty())
        {
            bdVoxel *pv = added_pixels.GetLeftEndPointerToElement();
            
            //output(0,0,pv->R(),pv->C()) = region_value;
            
            //check if it is a border with other region...
            
            int n_of_non_zero_pixels = 0;
            //int rn,cn;
            for(g.ForCoordinates_Circle(pv->R(),pv->C(),0); g.Get_Circle(SE_squared_radius,rn,cn); )
            {
                if(temp(t,s,rn,cn)>0)
                {
                    n_of_non_zero_pixels++;
                }
            }
            
            if( n_of_non_zero_pixels < ((n_of_pixels_in_SE*percentage)/100) ) output(t,s,pv->R(),pv->C()) = segmented_label;
            
            added_pixels.DeleteLeftEnd();
        }
        
    }
    
    
    //output.CopyFrom(extracted_slice);
    
    return 1;
}



int bdMelanoma::SegmentStreaks(bdImage &original, bdImage &segmented_melanoma, unsigned int percentage, bdImage &output)
{
    unsigned int SE_radius = 6;
    unsigned short segmented_label = 65535;
    unsigned int t_input = 0;
    unsigned int d_threshold = 170;
    
    if(!this->SegmentStructures(original, t_input, segmented_melanoma, (SE_radius*SE_radius), percentage, output, segmented_label)) return 0;


    //----- Label connected components -----
    {
        unsigned short label = 1;

        for(unsigned int r=0; r< output.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
            {
                if(output(r,c)!=0)
                {
                    if(output(r,c)==segmented_label)
                    {
                        
                        
                        //perform region growing
                        bdList< bdVoxel > temp_list;
                        bdVoxel v;
                        temp_list.AddToLeftEnd(v(r,c));
                        output(v.R(),v.C()) = label;

                        bdList< bdVoxel > CC_list;
                        CC_list.AddToLeftEnd(v);
                        
                        bdVoxel CC_min_r, CC_min_c, CC_max_r, CC_max_c;// borders of the CC (left, right, bottom and top of the rectangular border)
                        CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
                        CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
                        CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
                        CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
                        

                        
                        bdGeometry g;
                        g.SetDimensions(1, output.GetNumberOfRows(), output.GetNumberOfColumns());
                        while(!temp_list.IsEmpty())
                        {
                            v = temp_list.GetLeftEnd();

                            
                            // Update rectangular border
                            if(v.R() > CC_max_r.R())
                            {
                                CC_max_r.R() = v.R(); CC_max_r.C() = v.C();
                            }
                            if(v.C() > CC_max_c.C())
                            {
                                CC_max_c.R() = v.R(); CC_max_c.C() = v.C();
                            }
                            if(v.R() < CC_min_r.R())
                            {
                                CC_min_r.R() = v.R(); CC_min_r.C() = v.C();
                            }
                            if(v.C() < CC_min_c.C())
                            {
                                CC_min_c.R() = v.R(); CC_min_c.C() = v.C();
                            }

                            int rn, cn;
                            for(g.ForCoordinates_4_Neighbors(v.R(), v.C()); g.Get_4_Neighbors(rn, cn); )
                            {
                                if(output(rn,cn)==segmented_label)
                                {
                                    bdVoxel vn;
                                    temp_list.AddToRightEnd(vn(rn,cn));
                                    output(vn.R(),vn.C()) = label;

                                    CC_list.AddToLeftEnd(vn);
                                }
                            }

                            temp_list.DeleteLeftEnd();
                        }

                        // Now the CC is labeled with value 'label' and is contained in CC_list.
                        // determine the max distance of the recorded voxels on the rectangular border
                        unsigned int d1 = CC_max_r.DistanceEuclideanSquared(CC_max_c);
                        unsigned int d2 = CC_max_r.DistanceEuclideanSquared(CC_min_c);
                        unsigned int d3 = CC_max_r.DistanceEuclideanSquared(CC_min_r);
                        unsigned int d4 = CC_max_c.DistanceEuclideanSquared(CC_min_c);
                        unsigned int d5 = CC_max_c.DistanceEuclideanSquared(CC_min_r);
                        unsigned int d6 = CC_min_c.DistanceEuclideanSquared(CC_min_r);
                        
                        unsigned int d_max = d1;
                        if(d2>d_max) d_max = d2;
                        if(d3>d_max) d_max = d3;
                        if(d4>d_max) d_max = d4;
                        if(d5>d_max) d_max = d5;
                        if(d6>d_max) d_max = d6;
                        
                        
                        unsigned int d_maxR_maxC = d1;
                        if(d2 > d_maxR_maxC) d_maxR_maxC = d2;
                        
                        unsigned int d_maxR_minC = d2;
                        if(d5 > d_maxR_minC) d_maxR_minC = d5;
                        
                        unsigned int d_maxR_minR = d3;
                        if(d4 > d_maxR_minR) d_maxR_minR = d4;
                        
                        unsigned int d_min = d_maxR_maxC;
                        if(d_maxR_minC<d_min) d_min = d_maxR_minC;
                        if(d_maxR_minR<d_min) d_min = d_maxR_minR;


                        
//                        unsigned int d_min = d1;
//                        if(d2<d_min) d_min = d2;
//                        if(d3<d_min) d_min = d3;
//                        if(d4<d_min) d_min = d4;
//                        if(d5<d_min) d_min = d5;
//                        if(d6<d_min) d_min = d6;
                        
                        unsigned int d = d_max - d_min;

                        // If the CC is not long enough, erase it.
                        //if(d_max < d_threshold)
                        if(d < d_threshold)
                        {
                            while(!CC_list.IsEmpty())
                            {
                                v = CC_list.GetLeftEnd();

                                output(v.R(),v.C()) = 0;

                                CC_list.DeleteLeftEnd();
                            }
                        }
                        
                        label = (label+1) % segmented_label;
                    }
                }
            }
        }
    }
    //----------
    
    return 1;
    
}




int bdMelanoma::SegmentStructurelessArea(bdImage &original, bdImage &segmented_melanoma, unsigned int percentage, bdImage &output)
{
    unsigned int t_input = 0;
    unsigned int circle_squared_radius = 2;
    
    output.SetSize(1,1,original.GetNumberOfRows(),original.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(original);
    output.FillInWith(0);

    
    bdGIP gip;
    bdImage blurred_image;
    gip.MeanCircleFor2DSlice(original, segmented_melanoma, t_input, 0, 0, 0, circle_squared_radius, blurred_image);
    
    for(unsigned int r=0; r<segmented_melanoma.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_melanoma.GetNumberOfColumns(); c++)
        {
            if(segmented_melanoma(r,c)!=0)
            {
                int diff = original(t_input,0,r,c) - blurred_image(0,0,r,c);
                if(diff<0) diff = -diff;
                output(r,c) = diff;
            }
        }
    }
    
    return 1;
}



int bdMelanoma::SegmentStructurelessArea2(bdImage &mask, bdImage &border_image, bdImage &dark_structures_image, bdImage &bright_structures_image, bdImage &streaks_image, bdImage &output)
{
    output.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(mask);
    output.FillInWith(0);
    
    // Calculate area of the lesion
    double area_in_pixels = 0;
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                area_in_pixels++;
            }
        }
    }

    
    bdGeometry g;
    g.SetDimensions(1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
    
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                if(border_image(0,0,r,c)==0 && dark_structures_image(r,c)==0 && bright_structures_image(r,c)==0 && streaks_image(r,c)==0)
                {
                    int is_radius_found = 0;
                    unsigned int found_squared_radius = 0;
                    for(unsigned int squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        int rn, cn;
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius+1,rn,cn); )
                        {
                            if(border_image(rn,cn)!=0)
                            {
                                is_radius_found = 1;
                                found_squared_radius = 0;
                                break;
                            }
                            if(dark_structures_image(rn,cn)!=0 || bright_structures_image(rn,cn)!=0 || streaks_image(rn,cn)!=0)
                            {
                                is_radius_found = 1;
                                found_squared_radius = squared_radius;
                            }
                        }
                    }
//!!!!!!! PARAMETER SET AD_HOC !!!!!!!!!!
                    if(area_in_pixels/10.0 <= found_squared_radius*3.14 * 3) output(r,c) = found_squared_radius;
                    else output(r,c) = 0;
                    
                }
            }
        }
    }
    
    bdImage temp;
    temp.CopyFrom(output);
    
    //Go through the temp image and grow a sphere for every non-zero pixel
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0 && temp(r,c)!=0)
            {
                //int is_radius_found = 0;
                //for(unsigned int squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                //{
                    int rn, cn;
                    for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(temp(r,c),rn,cn); )
                    {
                        output(rn,cn) = 255;
                    }
                //}
            }
        }
    }
    
    return 1;
}





int bdMelanoma::AsymmetryOfMelanoma(bdImage &original, bdImage &mask, unsigned int number_of_symmetry_axes, bdArray<double> &output_min_symmetry_coeff, bdArray<double> &output_max_symmetry_coeff, bdArray<double> &output_mean_symmetry_coeff)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    if(number_of_symmetry_axes<2) return 0;
    
    
    // Calculate the center of mass for segmented melanoma (masked region in mask image)
    double center_of_mass_r = 0, center_of_mass_c = 0, n_of_pixels = 0;
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                center_of_mass_c += c;
                center_of_mass_r += r;
                n_of_pixels += 1;
            }
        }
    }
    center_of_mass_c = center_of_mass_c / n_of_pixels;
    center_of_mass_r = center_of_mass_r / n_of_pixels;
    
    //cout<<"  center_of_mass=("<<center_of_mass_r<<","<<center_of_mass_c<<")"<<endl;
    
    
    // Set the size of output arrays to the number of time series of original image and initialize the values to incorrect (impossible) values
    output_min_symmetry_coeff.Set(original.GetNumberOfTimeSeries());
    output_max_symmetry_coeff.Set(original.GetNumberOfTimeSeries());
    output_mean_symmetry_coeff.Set(original.GetNumberOfTimeSeries());
    output_min_symmetry_coeff.FillInWith(65535);
    output_max_symmetry_coeff.FillInWith(-1);
    output_mean_symmetry_coeff.FillInWith(0);
    
    
    // For each image time instance (modality)...
    for(unsigned int t=0; t<original.GetNumberOfTimeSeries(); t++)
    {
        // initialize a sum of all coefficients for this time instance to 0.
        double sum_of_asymmetry_coeff_for_this_t = 0;
        
        //... and for each angle check symmetry
        for(unsigned int a=0; a<number_of_symmetry_axes; a++)
        {
            // Calculate the parameters for the line equation: y = k*x + n or in our case: c = k*r + n.
            double angle = (180.0 / ((double)number_of_symmetry_axes)) * a;
            double k = (sin(angle)) / (cos(angle));// delta_y / delta_x.
            double n = center_of_mass_c - (k*center_of_mass_r);// because c = k*r + n --> n = c - k*r.
            
            //cout<<"  angle="<<angle<<", k="<<k<<", n="<<n<<endl;

            
            // Create indicator image to see if the pixel symmetry has already been tested.
            bdImage indicator;
            indicator.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
            indicator.FillInWith(0);
            
            
            double n_of_symmetric_pixels = 0;//number of pixels for which we were able to find a symmetric pair.
            double asymmetry_coeff_sum = 0;
            for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0 && indicator(r,c)==0)
                    {
                        // set the indicator that we have processed this pixel.
                        indicator(r,c) = 1;
                        
                        // calculate the symmetric indexes
                        double d = (r + (c - n)*k)/(1 + k*k); // d = (x + (y - c)*k)/(1 + k^2)
                        double r_symmetric_real = 2*d - r; // x' = 2*d - x
                        double c_symmetric_real = 2*d*k - c + 2*n; // y' = 2*d*k - y + 2c

                        int r_symmetric_index = (int) r_symmetric_real;
                        int c_symmetric_index = (int) c_symmetric_real;
                        
                        // Check if the symmetric index falls within the image...
                        if(r_symmetric_index>=0 && r_symmetric_index<mask.GetNumberOfRows() && c_symmetric_index>=0 && c_symmetric_index<mask.GetNumberOfColumns())
                        {
                            //... and if so, calculate the symmetry for the given pixel and its symmetric pixel.
                            if(mask(r_symmetric_index,c_symmetric_index)==0)
                            {
                                asymmetry_coeff_sum += 1;
                                original(t,0,r,c) = 255;
                            }
                            n_of_symmetric_pixels += 1;
                        }
                    }
                }
            }
            
            //cout<<"  asymmetry_coeff_sum="<<asymmetry_coeff_sum<<"  ";
            
            // calculate the asymmetry coeff for this angle and set the output array values (if need be)
            double asymmetry_coeff = asymmetry_coeff_sum / n_of_symmetric_pixels;
            sum_of_asymmetry_coeff_for_this_t += asymmetry_coeff;
            
            //cout<<"("<<t<<","<<a<<")="<<asymmetry_coeff<<"   ";
            
            if(asymmetry_coeff>output_max_symmetry_coeff[t]) { output_max_symmetry_coeff[t] = asymmetry_coeff; }
            if(asymmetry_coeff<output_min_symmetry_coeff[t]) { output_min_symmetry_coeff[t] = asymmetry_coeff; }
        }
        
        output_mean_symmetry_coeff[t] = sum_of_asymmetry_coeff_for_this_t / ((double) number_of_symmetry_axes);
        
        
        //cout<<" ["<<t<<"]_min="<<output_min_symmetry_coeff[t]<<",max="<<output_max_symmetry_coeff[t]<<",avrg="<<output_mean_symmetry_coeff[t]<<" ";
        
    }
    
    
    return 1;
}



int bdMelanoma::AsymmetryOfColorsInRegionLists_Helper(bdList<bdColorCode> &list_of_colors_in_region1, bdList<bdColorCode> &list_of_colors_in_region2, double n_of_pixels_in_region1, double n_of_pixels_in_region2,  double &max_diff, const char *text_for_output)
{
    double multiplication_coeff = 10; // equals 10 because a color/structure has to be at least 10% of the whole region to be considered important.
    
    
    bdListIterator<bdColorCode> itcc;
    bdList<bdColorCode> list_of_accessed_colors;//list of accessed colors, to keep track which ones were processed.
    
    for(itcc.SetLeftEnd(list_of_colors_in_region1); itcc.IsValid(); itcc.MoveRight())
    {
        if(!bdListing::GetNodeWithElement(list_of_accessed_colors, itcc.GetElement()))
        {
            double n_of_pixels_of_current_structure = 1;
            bdListNode<bdColorCode> *lncc = bdListing::GetNodeWithElement(list_of_colors_in_region2, itcc.GetElement());
            if(lncc)
            {
                n_of_pixels_of_current_structure = lncc->GetElementPointer()->m_intensity;
            }
            
            double ratio_coeff;
            if(n_of_pixels_of_current_structure < itcc.GetElementPointer()->m_intensity)
            {
                ratio_coeff = (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / n_of_pixels_of_current_structure;
            }
            else
            {
                ratio_coeff = (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / itcc.GetElementPointer()->m_intensity;
            }
            if(ratio_coeff<-1 || ratio_coeff>1) ratio_coeff = 1;
            
            double diff_current = ratio_coeff * multiplication_coeff * (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / ((double)((n_of_pixels_in_region1+n_of_pixels_in_region2)/2));
            if(diff_current<0) diff_current = -diff_current;
            
            if(diff_current>max_diff) max_diff = diff_current;
            
            cout<<endl<<"  diff of "<<text_for_output<<" ("<<itcc.GetElementPointer()->m_red<<","<<itcc.GetElementPointer()->m_green<<","<<itcc.GetElementPointer()->m_blue<<") is "<<diff_current;
            
            list_of_accessed_colors.AddToLeftEnd(itcc.GetElement());
        }
    }
    for(itcc.SetLeftEnd(list_of_colors_in_region2); itcc.IsValid(); itcc.MoveRight())
    {
        if(!bdListing::GetNodeWithElement(list_of_accessed_colors, itcc.GetElement()))
        {
            double n_of_pixels_of_current_structure = 1;
            bdListNode<bdColorCode> *lncc = bdListing::GetNodeWithElement(list_of_colors_in_region1, itcc.GetElement());
            if(lncc)
            {
                n_of_pixels_of_current_structure = lncc->GetElementPointer()->m_intensity;
            }
            
            double ratio_coeff;
            if(n_of_pixels_of_current_structure < itcc.GetElementPointer()->m_intensity)
            {
                ratio_coeff = (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / n_of_pixels_of_current_structure;
            }
            else
            {
                ratio_coeff = (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / itcc.GetElementPointer()->m_intensity;
            }
            if(ratio_coeff<-1 || ratio_coeff>1) ratio_coeff = 1;
            
            double diff_current = ratio_coeff * multiplication_coeff * (itcc.GetElementPointer()->m_intensity - n_of_pixels_of_current_structure) / ((double)((n_of_pixels_in_region1+n_of_pixels_in_region2)/2));
            if(diff_current<0) diff_current = -diff_current;
            
            if(diff_current>max_diff) max_diff = diff_current;
            
            cout<<endl<<"  diff of "<<text_for_output<<" ("<<itcc.GetElementPointer()->m_red<<","<<itcc.GetElementPointer()->m_green<<","<<itcc.GetElementPointer()->m_blue<<") is "<<diff_current;
            
            list_of_accessed_colors.AddToLeftEnd(itcc.GetElement());
        }
    }
    
    return 1;
}



int bdMelanoma::AsymmetryOfMelanoma2(bdImage &original, bdImage &mask, bdImage &mapped_chart_colors_image, bdImage &combined_structure_image, unsigned int number_of_symmetry_axes, bdImage &output)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    // Image containing the masked area devided into 4 subareas by calculated symmetry axes.
    bdImage quartiles_image;
    //quartiles_image.SetSize(3,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
    //quartiles_image.SetVisualizationPropertiesToMatchInput(mask);
    quartiles_image.SetSizeAndPropertiesAs(original);
    quartiles_image.FillInWith(0);
    
    // Calculate the center of mass for segmented melanoma (masked region in mask image)
    double center_of_mass_r = 0, center_of_mass_c = 0;
    unsigned int lesion_area_in_pixels = 0;
    for(unsigned int r=0; r< mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                lesion_area_in_pixels++;
                center_of_mass_c += c;
                center_of_mass_r += r;
            }
        }
    }
    center_of_mass_c = center_of_mass_c / ((double)lesion_area_in_pixels);
    center_of_mass_r = center_of_mass_r / ((double)lesion_area_in_pixels);
    
    // Put the masked voxels into an array for faster access.
    bdArray<bdVoxel> voxel_array;
    voxel_array.Set(lesion_area_in_pixels);
    unsigned int i=0;
    for(unsigned int r=0; r< mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                quartiles_image(0,0,r,c) = mask(r,c);
                quartiles_image(1,0,r,c) = mask(r,c);
                quartiles_image(2,0,r,c) = mask(r,c);
                voxel_array[i](r,c);
                i++;
            }
        }
    }
    
    
    //----- Find the most symmetric axis in terms of shape -----
    double min_asymmetry_coeff = -1;
    unsigned int index_of_min_asymmetry_coeff = 0;
    
    //... and for each angle check symmetry
    for(unsigned int a=0; a<number_of_symmetry_axes; a++)
    {
        // Calculate the parameters for the line equation: y = k*x + n or in our case: c = k*r + n.
        double angle = (180.0 / ((double)number_of_symmetry_axes)) * a;
        double k = (sin(angle)) / (cos(angle));// delta_y / delta_x.
        double n = center_of_mass_c - (k*center_of_mass_r);// because c = k*r + n --> n = c - k*r.
        
        //cout<<"  angle="<<angle<<", k="<<k<<", n="<<n<<endl;
        
        // Create indicator image to see if the pixel symmetry has already been tested.
        bdImage indicator;
        indicator.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
        indicator.FillInWith(0);
        
        double n_of_mapped_pixels = 0;//number of pixels for which we were able to find a symmetric pair.
        double asymmetry_coeff_sum = 0;
        for(unsigned int i=0; i<voxel_array.GetNumberOfElements(); i++)
        {
            unsigned int r = voxel_array[i].R();
            unsigned int c = voxel_array[i].C();
            if(indicator(r,c)==0)
            {
                // set the indicator that we have processed this pixel.
                indicator(r,c) = 1;

                // calculate the symmetric indexes
                double d = (r + (c - n)*k)/(1 + k*k); // d = (x + (y - c)*k)/(1 + k^2)
                double r_symmetric_real = 2*d - r; // x' = 2*d - x
                double c_symmetric_real = 2*d*k - c + 2*n; // y' = 2*d*k - y + 2c

                int r_symmetric_index = (int) r_symmetric_real;
                int c_symmetric_index = (int) c_symmetric_real;

                // Check if the symmetric index falls within the image...
                if(r_symmetric_index>=0 && r_symmetric_index<mask.GetNumberOfRows() && c_symmetric_index>=0 && c_symmetric_index<mask.GetNumberOfColumns())
                {
                    //... and if so, calculate the symmetry for the given pixel and its symmetric pixel.
                    if(mask(r_symmetric_index,c_symmetric_index)==0)
                    {
                        asymmetry_coeff_sum += 1;
                        //original(t,0,r,c) = 255;
                    }
                    n_of_mapped_pixels += 1;
                }
            }
        }
        
        //cout<<"  asymmetry_coeff_sum="<<asymmetry_coeff_sum<<"  ";
        
        // calculate the asymmetry coeff for this angle
        double asymmetry_coeff = asymmetry_coeff_sum / n_of_mapped_pixels;
        

        cout<<"asymmetry_coeff["<<a<<"]="<<asymmetry_coeff<<"   ";

        
        if(asymmetry_coeff<min_asymmetry_coeff || min_asymmetry_coeff<0)
        {
            index_of_min_asymmetry_coeff = a;
            min_asymmetry_coeff = asymmetry_coeff;
        }
        
    }
    //----------
    
    
    
    cout<<"min_asymmetry_coeff["<<index_of_min_asymmetry_coeff<<"]="<<min_asymmetry_coeff<<"   ";
    
    //----- Determine asymetry axes -----
    int normal_R[2], normal_C[2]; // normals will store directions of asymmetry axes.
    int is_R_normal_determined = 0;
    int is_C_normal_determined = 0;
    {
        unsigned int a = index_of_min_asymmetry_coeff;
        
        // Calculate the parameters for the line equation: y = k*x + n or in our case: c = k*r + n.
        double angle = (180.0 / ((double)number_of_symmetry_axes)) * a;
        double k = (sin(angle)) / (cos(angle));// delta_y / delta_x.
        double n = center_of_mass_c - (k*center_of_mass_r);// because c = k*r + n --> n = c - k*r.
        
        //output.CopyFrom(mask);
        
        // Iterate over rows
        for(unsigned int r=0; r< mask.GetNumberOfRows(); r++)
        {
            int c = (int) (k*((double)r) + n);
            if(c>=0 && c<mask.GetNumberOfColumns())
            {
//                output(0,0,r,c) = 0;
//                output(1,0,r,c) = 0;
//                output(2,0,r,c) = 0;
                quartiles_image(0,0,r,c) = 0;
                quartiles_image(1,0,r,c) = 0;
                quartiles_image(2,0,r,c) = 0;
                
                if(!is_R_normal_determined)
                {
                    normal_R[0] = ((int)center_of_mass_r) - r;
                    normal_R[1] = ((int)center_of_mass_c) - c;
                }
            }
            
            int dr = r - ((int)center_of_mass_r);
            int dc = c - ((int)center_of_mass_c);
            
            int rn = ((int)center_of_mass_r) - dc;
            int cn = ((int)center_of_mass_c) + dr;
            if(rn>=0 && rn<mask.GetNumberOfRows() && cn>=0 && cn<mask.GetNumberOfColumns())
            {
                //output(0,0,rn,cn) = 0;
                //output(1,0,rn,cn) = 0;
                //output(2,0,rn,cn) = 0;
                quartiles_image(0,0,rn,cn) = 0;
                quartiles_image(1,0,rn,cn) = 0;
                quartiles_image(2,0,rn,cn) = 0;
                
                if(!is_C_normal_determined)
                {
                    normal_C[0] = ((int)center_of_mass_r) - rn;
                    normal_C[1] = ((int)center_of_mass_c) - cn;
                }

            }
            
        }
        
        
        // Iterate over columns
        for(unsigned int c=0; c< mask.GetNumberOfColumns(); c++)
        {
            // c = k*r + n  =>  r = (c - n) / k;
            if(k!=0)
            {
                int r = (int) ((((double)c) - n) / k);
                
                if(r>=0 && r<mask.GetNumberOfRows())
                {
                    //                output(0,0,r,c) = 0;
                    //                output(1,0,r,c) = 0;
                    //                output(2,0,r,c) = 0;
                    quartiles_image(0,0,r,c) = 0;
                    quartiles_image(1,0,r,c) = 0;
                    quartiles_image(2,0,r,c) = 0;
                    
                    if(!is_R_normal_determined)
                    {
                        normal_R[0] = ((int)center_of_mass_r) - r;
                        normal_R[1] = ((int)center_of_mass_c) - c;
                    }
                }
                
                int dr = r - ((int)center_of_mass_r);
                int dc = c - ((int)center_of_mass_c);
                
                int rn = ((int)center_of_mass_r) - dc;
                int cn = ((int)center_of_mass_c) + dr;
                if(rn>=0 && rn<mask.GetNumberOfRows() && cn>=0 && cn<mask.GetNumberOfColumns())
                {
                    //output(0,0,rn,cn) = 0;
                    //output(1,0,rn,cn) = 0;
                    //output(2,0,rn,cn) = 0;
                    quartiles_image(0,0,rn,cn) = 0;
                    quartiles_image(1,0,rn,cn) = 0;
                    quartiles_image(2,0,rn,cn) = 0;
                    
                    if(!is_C_normal_determined)
                    {
                        normal_C[0] = ((int)center_of_mass_r) - rn;
                        normal_C[1] = ((int)center_of_mass_c) - cn;
                    }
                    
                }
            }
            
        }
        
        //cout<<"  angle="<<angle<<", k="<<k<<", n="<<n<<endl;
    }
    //-----------
    
    output.FillInWith(0);
    
    //----- Mark quartiles in the quartiles_image by projecting each pixel to axes -----
    for(unsigned int i=0; i<voxel_array.GetNumberOfElements(); i++)
    {
        int r = voxel_array[i].R();
        int c = voxel_array[i].C();
        if(mask(r,c)!=0)
        {
            output(0,0,r,c) = original(0,0,r,c); // 255;
            output(1,0,r,c) = original(1,0,r,c); //255;
            output(2,0,r,c) = original(2,0,r,c); //255;
        }
        
        // project coordinates to normals
        //if(output(r,c)!=0)
        if(quartiles_image(0,0,r,c)!=0 || quartiles_image(1,0,r,c)!=0 || quartiles_image(2,0,r,c)!=0)
        {
            double p1 = (r-((int)center_of_mass_r))*normal_R[0] + (c-((int)center_of_mass_c))*normal_R[1];
            
            if(p1<0)
            {
                //output(0,0,r,c) = 255;
                //output(1,0,r,c) = 0;
                quartiles_image(0,0,r,c) = 255;
                quartiles_image(1,0,r,c) = 0;
            }
            else
            {
                //output(0,0,r,c) = 127;
                //output(1,0,r,c) = 0;
                quartiles_image(0,0,r,c) = 127;
                quartiles_image(1,0,r,c) = 0;
            }
            
            double p2 = (r-((int)center_of_mass_r))*normal_C[0] + (c-((int)center_of_mass_c))*normal_C[1];
            
            if(p2<0)
            {
                //output(2,0,r,c) = 255;
                //output(1,0,r,c) = 0;
                quartiles_image(2,0,r,c) = 255;
                quartiles_image(1,0,r,c) = 0;
            }
            else
            {
                //output(2,0,r,c) = 127;
                //output(1,0,r,c) = 0;
                quartiles_image(2,0,r,c) = 127;
                quartiles_image(1,0,r,c) = 0;
            }
            
        }
        else
        {
            // mark the axes in the output image.
            output(0,0,r,c) = 0;
            output(1,0,r,c) = 0;
            output(2,0,r,c) = 0;
        }
    }
    //------------
    
    //----- Put quartiles from the quartiles_image to arrays for faster access -----
    bdArray<bdVoxel> q00_voxel_array, q01_voxel_array, q10_voxel_array, q11_voxel_array;
    int n_of_pixels_q00 = 0, n_of_pixels_q01 = 0, n_of_pixels_q10 = 0, n_of_pixels_q11 = 0;
    
    // count the number of pixels in each quartile
    for(unsigned int i=0; i<voxel_array.GetNumberOfElements(); i++)
    {
        int r = voxel_array[i].R();
        int c = voxel_array[i].C();
        
        // If q00
        if(quartiles_image(0,0,r,c)==127 && quartiles_image(2,0,r,c)==127)
        {
            n_of_pixels_q00++;
        }
        // If q01
        if(quartiles_image(0,0,r,c)==127 && quartiles_image(2,0,r,c)==255)
        {
            n_of_pixels_q01++;
        }
        // If q10
        if(quartiles_image(0,0,r,c)==255 && quartiles_image(2,0,r,c)==127)
        {
            n_of_pixels_q10++;
        }
        // If q11
        if(quartiles_image(0,0,r,c)==255 && quartiles_image(2,0,r,c)==255)
        {
            n_of_pixels_q11++;
        }
    }
    // Put the pixels from each quartile to corresponding array.
    {
        q00_voxel_array.Set(n_of_pixels_q00);
        q01_voxel_array.Set(n_of_pixels_q01);
        q10_voxel_array.Set(n_of_pixels_q10);
        q11_voxel_array.Set(n_of_pixels_q11);
        int q00 = 0, q01 = 0, q10 = 0, q11 = 0;
        for(unsigned int i=0; i<voxel_array.GetNumberOfElements(); i++)
        {
            int r = voxel_array[i].R();
            int c = voxel_array[i].C();
            
            // If q00
            if(quartiles_image(0,0,r,c)==127 && quartiles_image(2,0,r,c)==127)
            {
                q00_voxel_array[q00].R() = r;
                q00_voxel_array[q00].C() = c;
                q00++;
            }
            // If q01
            if(quartiles_image(0,0,r,c)==127 && quartiles_image(2,0,r,c)==255)
            {
                q01_voxel_array[q01].R() = r;
                q01_voxel_array[q01].C() = c;
                q01++;
            }
            // If q10
            if(quartiles_image(0,0,r,c)==255 && quartiles_image(2,0,r,c)==127)
            {
                q10_voxel_array[q10].R() = r;
                q10_voxel_array[q10].C() = c;
                q10++;
            }
            // If q11
            if(quartiles_image(0,0,r,c)==255 && quartiles_image(2,0,r,c)==255)
            {
                q11_voxel_array[q11].R() = r;
                q11_voxel_array[q11].C() = c;
                q11++;
            }
        }
    }
    //----------
    
    //----- For each quartile perform analysis of colors and structures -----
    //--- populate Q00 structures list ---
    bdList<bdColorCode> q00_structures_list;
    bdList<bdColorCode> q00_colors_list;
    {
        for(unsigned int q=0; q<q00_voxel_array.GetNumberOfElements(); q++)
        {
            // Structures
            {
                bdColorCode cc;
                cc.Set(combined_structure_image(0,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()),
                       combined_structure_image(1,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()),
                       combined_structure_image(2,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()));
                
                //output(0,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()) = 160;
                //output(1,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()) = 160;
                //output(2,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()) = 160;
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q00_structures_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q00_structures_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
            
            // Colors
            {
                bdColorCode cc;
                cc.Set(mapped_chart_colors_image(0,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()),
                       mapped_chart_colors_image(1,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()),
                       mapped_chart_colors_image(2,0,q00_voxel_array[q].R(),q00_voxel_array[q].C()));
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q00_colors_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q00_colors_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
        }
        // Print out the Q_00 structures.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_00 structures: ";
            for(it.SetLeftEnd(q00_structures_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
        // Print out the Q_00 colors.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_00 colors: ";
            for(it.SetLeftEnd(q00_colors_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
    }
    //---
    //--- populate Q01 structures list---
    bdList<bdColorCode> q01_structures_list;
    bdList<bdColorCode> q01_colors_list;
    {
        for(unsigned int q=0; q<q01_voxel_array.GetNumberOfElements(); q++)
        {
            // Structures
            {
                bdColorCode cc;
                cc.Set(combined_structure_image(0,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()),
                       combined_structure_image(1,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()),
                       combined_structure_image(2,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()));
                
                //output(0,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()) = 160;
                //output(1,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()) = 160;
                //output(2,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()) = 160;
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q01_structures_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q01_structures_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
            
            // Colors
            {
                bdColorCode cc;
                cc.Set(mapped_chart_colors_image(0,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()),
                       mapped_chart_colors_image(1,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()),
                       mapped_chart_colors_image(2,0,q01_voxel_array[q].R(),q01_voxel_array[q].C()));
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q01_colors_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q01_colors_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
        }
        // Print out the Q_01 structures.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_01 structures: ";
            for(it.SetLeftEnd(q01_structures_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
        // Print out the Q_01 colors.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_01 colors: ";
            for(it.SetLeftEnd(q01_colors_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
    }
    //---
    //--- populate Q11 structures list---
    bdList<bdColorCode> q11_structures_list;
    bdList<bdColorCode> q11_colors_list;
    {
        for(unsigned int q=0; q<q11_voxel_array.GetNumberOfElements(); q++)
        {
            // Structures
            {
                bdColorCode cc;
                cc.Set(combined_structure_image(0,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()),
                       combined_structure_image(1,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()),
                       combined_structure_image(2,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()));
                
                //output(0,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()) = 160;
                //output(1,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()) = 160;
                //output(2,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()) = 160;
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q11_structures_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q11_structures_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
            
            // Colors
            {
                bdColorCode cc;
                cc.Set(mapped_chart_colors_image(0,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()),
                       mapped_chart_colors_image(1,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()),
                       mapped_chart_colors_image(2,0,q11_voxel_array[q].R(),q11_voxel_array[q].C()));
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q11_colors_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q11_colors_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
        }
        // Print out the Q_11 structures.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_11 structures: ";
            for(it.SetLeftEnd(q11_structures_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
        // Print out the Q_11 colors.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_11 colors: ";
            for(it.SetLeftEnd(q11_colors_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
    }
    //---
    //--- populate Q10 structures list ---
    bdList<bdColorCode> q10_structures_list;
    bdList<bdColorCode> q10_colors_list;
    {
        for(unsigned int q=0; q<q10_voxel_array.GetNumberOfElements(); q++)
        {
            // Structures
            {
                bdColorCode cc;
                cc.Set(combined_structure_image(0,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()),
                       combined_structure_image(1,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()),
                       combined_structure_image(2,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()));
                
                //output(0,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()) = 160;
                //output(1,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()) = 160;
                //output(2,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()) = 160;
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q10_structures_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q10_structures_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
            
            // Colors
            {
                bdColorCode cc;
                cc.Set(mapped_chart_colors_image(0,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()),
                       mapped_chart_colors_image(1,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()),
                       mapped_chart_colors_image(2,0,q10_voxel_array[q].R(),q10_voxel_array[q].C()));
                
                bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(q10_colors_list, cc);
                if(ln)
                {
                    ln->GetElementPointer()->m_intensity++;
                }
                else
                {
                    q10_colors_list.AddNewToRightEnd()->Set(cc.m_red, cc.m_green, cc.m_blue, 1);
                }
            }
        }
        // Print out the Q_10 structures.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_10 structures: ";
            for(it.SetLeftEnd(q10_structures_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
        // Print out the Q_10 colors.
        {
            bdListIterator<bdColorCode> it;
            cout<<endl<<" Q_10 colors: ";
            for(it.SetLeftEnd(q10_colors_list); it.IsValid(); it.MoveRight())
            {
                cout<<" ("<<it.GetElementPointer()->m_red<<","<<it.GetElementPointer()->m_green<<","<<it.GetElementPointer()->m_blue<<")["<<it.GetElementPointer()->m_intensity<<"]  ";
            }
            cout<<endl;
        }
    }
    //---
    //-----
    
    
    //----- Analyze asymmetry per axis -----
    //
    // 00 01
    // 10 11
    //
    //--------------------------------------
    
    
    //--- Axis 1 analysis: q00,q01 against q10,q11 ---
    {
        
        double max_diff_structure = 0;// max value of normalized difference in structures for the examined qurtiles.
        
        double max_diff_color = 0;// max value of normalized difference in colors for the examined qurtiles.
        
        //--- Compare structures -----
        {
            // Compare q00 to q10
            cout<<endl<<" Comparing structures of q00 nd q10: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q00_structures_list, q10_structures_list, n_of_pixels_q00, n_of_pixels_q10,  max_diff_structure, "structure");
            
            //Compare q11 to q01
            cout<<endl<<" Comparing structures of q11 nd q01: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q11_structures_list, q01_structures_list, n_of_pixels_q11, n_of_pixels_q01,  max_diff_structure, "structure");
        }
        //-----
        
        //--- Compare colors -----
        {
            // Compare q00 to q10
            cout<<endl<<" Comparing colors of q00 nd q10: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q00_colors_list, q10_colors_list, n_of_pixels_q00, n_of_pixels_q10,  max_diff_color, "color");
            
            //Compare q11 to q01
            cout<<endl<<" Comparing structures of q11 nd q01: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q11_colors_list, q01_colors_list, n_of_pixels_q11, n_of_pixels_q01,  max_diff_color, "color");
        }
        //-----
        
        
        cout<<endl<<" Max normalized diff of sructure is "<<max_diff_structure<<endl;
        cout<<endl<<" Max normalized diff of color is "<<max_diff_color<<endl;
        
        if(max_diff_structure>max_diff_color) this->m_asymmetry_soft_score_axis_1 = max_diff_structure;
        else this->m_asymmetry_soft_score_axis_1 = max_diff_color;
        
        if(this->m_asymmetry_soft_score_axis_1 >= 1.0)
        {
            this->m_asymmetry_soft_score_axis_1 = 1;
            this->m_asymmetry_hard_score_axis_1 = 1;
        }
        else
        {
            this->m_asymmetry_hard_score_axis_1 = 0;
        }
        
    }
    //---
    //--- Axis 2 analysis: q00,q10 against q01,q11 ---
    {
        
        double max_diff_structure = 0;// max value of normalized difference in structures for the examined qurtiles.
        
        double max_diff_color = 0;// max value of normalized difference in colors for the examined qurtiles.
        
        //--- Compare structures -----
        {
            // Compare q00 to q01
            cout<<endl<<" Comparing structures of q00 nd q01: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q00_structures_list, q01_structures_list, n_of_pixels_q00, n_of_pixels_q01,  max_diff_structure, "structure");
            
            //Compare q11 to q10
            cout<<endl<<" Comparing structures of q11 nd q10: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q11_structures_list, q10_structures_list, n_of_pixels_q11, n_of_pixels_q10,  max_diff_structure, "structure");
        }
        //-----
        
        //--- Compare colors -----
        {
            // Compare q00 to q01
            cout<<endl<<" Comparing colors of q00 nd q01: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q00_colors_list, q01_colors_list, n_of_pixels_q00, n_of_pixels_q01,  max_diff_color, "color");
            
            //Compare q11 to q10
            cout<<endl<<" Comparing structures of q11 nd q10: "<<endl;
            this->AsymmetryOfColorsInRegionLists_Helper(q11_colors_list, q10_colors_list, n_of_pixels_q11, n_of_pixels_q10,  max_diff_color, "color");
        }
        //-----
        
        
        cout<<endl<<" Max normalized diff of sructure is "<<max_diff_structure<<endl;
        cout<<endl<<" Max normalized diff of color is "<<max_diff_color<<endl;
        
        if(max_diff_structure>max_diff_color) this->m_asymmetry_soft_score_axis_2 = max_diff_structure;
        else this->m_asymmetry_soft_score_axis_2 = max_diff_color;
        
        if(this->m_asymmetry_soft_score_axis_2 >= 1.0)
        {
            this->m_asymmetry_soft_score_axis_2 = 1;
            this->m_asymmetry_hard_score_axis_2 = 1;
        }
        else
        {
            this->m_asymmetry_hard_score_axis_2 = 0;
        }
        
    }
    //---
    
    m_asymmetry_hard_score = m_asymmetry_hard_score_axis_1 + m_asymmetry_hard_score_axis_2;
    
    // If both soft scores are less than 1, take the maximum of the scores to be the final soft score. This way we solve the problem
    // when the cumulative soft score is larger than 1, while the hard score is 0. (e.g. soft = 0.6 + 0.8 = 1.4, but hard = 0).
    if(m_asymmetry_soft_score_axis_1<1 && m_asymmetry_soft_score_axis_2<1)
    {
        if(m_asymmetry_soft_score_axis_1>m_asymmetry_soft_score_axis_2)
        {
            m_asymmetry_soft_score = m_asymmetry_soft_score_axis_1;
        }
        else
        {
            m_asymmetry_soft_score = m_asymmetry_soft_score_axis_2;
        }
    }
    else
    {
        m_asymmetry_soft_score = m_asymmetry_soft_score_axis_1 + m_asymmetry_soft_score_axis_2;
    }
    
    //----------
    

    
    return 1;
}



int bdMelanoma::BorderOfMelanoma(bdImage &original, bdImage &mask, int squared_radius, bdArray<int> &output_contrast_along_border)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    if(squared_radius<2) return 0;
    
    bdGeometry g;
    g.SetDimensions(1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
    
    // Find all borders of the melanoma and store them to an image
    bdImage border_image;
    border_image.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
    border_image.SetVisualizationPropertiesToMatchInput(mask);
    border_image.FillInWith(0);
    
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            // if pixels belongs to the mask...
            if(mask(r,c)!=0 )
            {
                //... check if it is a border
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
                {
                    if(mask(rn,cn)==0)
                    {
                        border_image(r,c) = 255;
                        break;
                    }
                }
            }
        }
    }
    
    // Extract the largest border image.
    bdBIP bip;
    bdImage largest_border_image;
    bip.ExtractLargest_26_ConnectedComponent(border_image, 1, largest_border_image);
    
    // Find the first non-zero voxel in the largest_border_image
    int is_first_pixel_found = 0;
    bdVoxel first_voxel;
    for(unsigned int r=0; r<largest_border_image.GetNumberOfRows() && !is_first_pixel_found; r++)
    {
        for(unsigned int c=0; c<largest_border_image.GetNumberOfColumns(); c++)
        {
            if(largest_border_image(r,c)!=0)
            {
                is_first_pixel_found = 1;
                first_voxel.R() = r;
                first_voxel.C() = c;
                break;
            }
        }
    }
    
    // List of border voxels with calculated contrast
    bdList< bdVoxel > border_voxels_with_contrast;
    
    // Follow the border from the first found border pixel and calculate contrast
    bdList< bdVoxel > temp_list;
    bdImage temp_image;
    temp_image.SetSizeAndPropertiesAs(mask);
    temp_image.FillInWith(0);
    temp_list.AddToLeftEnd(first_voxel);
    temp_image(first_voxel.R(),first_voxel.C()) = 255;
    while(!temp_list.IsEmpty())
    {
        bdVoxel v = temp_list.GetLeftEnd();
        
        //----- Calculate contrast -----
        {
            int rn,cn;
            int sum_object = 0;
            int n_of_pixels_object = 0;
            int sum_background = 0;
            int n_of_pixels_background = 0;

            for(g.ForCoordinates_Circle(v.R(),v.C(),0); g.Get_Circle(squared_radius, rn, cn); )
            {
                if(mask(rn,cn)!=0)
                {
                    sum_object += original(rn,cn);
                    n_of_pixels_object++;
                }
                else
                {
                    sum_background += original(rn,cn);
                    n_of_pixels_background++;
                }
            }
            
            int contrast = (sum_object/n_of_pixels_object) - (sum_background/n_of_pixels_background);
            if(contrast<0) contrast = -contrast;
            
            bdVoxel voxel;
            voxel(v.R(),v.C());
            voxel.SetValue(contrast);
            border_voxels_with_contrast.AddToRightEnd(voxel);
            
        }
        //----------
        
        // Find the next border pixel and enter it to the list
        {
        int rn,cn;
        for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rn,cn); )
        {
            if(largest_border_image(rn,cn)!=0 && temp_image(rn,cn)==0)
            {
                temp_image(rn,cn)=255;
                bdVoxel v1;
                temp_list.AddToRightEnd(v1(rn,cn));
                break;
            }
        }
        }
        
        temp_list.DeleteLeftEnd();
    }
    
    // Finally, copy the list of border voxels with conntrast tooutput array.
    output_contrast_along_border.Set(border_voxels_with_contrast.GetNumberOfElements());
    bdListIterator< bdVoxel > it;
    int i = 0;
    for(it.SetLeftEnd(border_voxels_with_contrast), i = 0; it.IsValid() && i<output_contrast_along_border.GetNumberOfElements(); it.MoveRight(), i++)
    {
        output_contrast_along_border[i] = it.GetElement().V();
    }
    
    mask.CopyFrom(largest_border_image);
    
    return 1;
    
}


//int bdMelanoma::BorderOfMelanoma2(bdImage &mask, bdArray<int> &output_contrast_along_border, bdImage &output)
//{
//    if(mask.IsEmpty()) return 0;
//    
//    bdGeometry g;
//    g.SetDimensions(1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
//    
//    // Find all borders of the melanoma and store them to an image
//    bdImage border_image;
//    border_image.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
//    border_image.SetVisualizationPropertiesToMatchInput(mask);
//    border_image.FillInWith(0);
//    
//    output.SetSize(1,1,mask.GetNumberOfRows(),mask.GetNumberOfColumns());
//    output.SetVisualizationPropertiesToMatchInput(mask);
//    output.FillInWith(0);
//    
//    
//    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
//        {
//            // if pixels belongs to the mask...
//            if(mask(0,0,r,c)!=0 || mask(1,0,r,c)!=0 || mask(2,0,r,c)!=0)
//            {
//                //... check if it is a border
//                int rn,cn;
//                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
//                {
//                    if(mask(0,0,rn,cn)==0 && mask(1,0,rn,cn)==0 && mask(2,0,rn,cn)==0)
//                    {
//                        border_image(r,c) = 255;
//                        break;
//                    }
//                }
//            }
//        }
//    }
//    
//    // Extract the largest border image.
//    bdBIP bip;
//    bdImage largest_border_image;
//    bip.ExtractLargest_26_ConnectedComponent(border_image, 1, largest_border_image);
//    
//    // Find the first non-zero voxel in the largest_border_image
//    int is_first_pixel_found = 0;
//    bdVoxel first_voxel;
//    for(unsigned int r=0; r<largest_border_image.GetNumberOfRows() && !is_first_pixel_found; r++)
//    {
//        for(unsigned int c=0; c<largest_border_image.GetNumberOfColumns(); c++)
//        {
//            if(largest_border_image(r,c)!=0)
//            {
//                is_first_pixel_found = 1;
//                first_voxel.R() = r;
//                first_voxel.C() = c;
//                break;
//            }
//        }
//    }
//    
//    // List of border voxels with calculated contrast
//    bdList< bdVoxel > border_voxels_with_contrast;
//    
//    // Follow the border from the first found border pixel and calculate contrast
//    bdList< bdVoxel > temp_list;
//    bdImage temp_image;
//    temp_image.SetSizeAndPropertiesAs(mask);
//    temp_image.FillInWith(0);
//    temp_list.AddToLeftEnd(first_voxel);
//    temp_image(first_voxel.R(),first_voxel.C()) = 255;
//    while(!temp_list.IsEmpty())
//    {
//        bdVoxel v = temp_list.GetLeftEnd();
//        
//        //----- Calculate contrast -----
//        {
//            int rn,cn;
//            
//            unsigned int squared_radius = 0;
//            int is_radius_found = 0;
//            //int contrast = 0;
//            for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
//            {
//                for(g.ForCoordinates_Circle(v.R(),v.C(),squared_radius); g.Get_Circle(squared_radius, rn, cn); )
//                {
//                    if(mask(1,0,rn,cn)!=0 && mask(2,0,rn,cn)!=0)
//                    {
//                        is_radius_found = 1;
//                        //contrast = squared_radius;
//                        break;
//                    }
//                }
//            }
//            
//            output(0,0,v.R(),v.C()) = squared_radius;
//
//            
//            bdVoxel voxel;
//            voxel(v.R(),v.C());
//            voxel.SetValue(squared_radius);
//            border_voxels_with_contrast.AddToRightEnd(voxel);
//            
//        }
//        //----------
//        
//        // Find the next border pixel and enter it to the list
//        {
//            int rn,cn;
//            for(g.ForCoordinates_8_Neighbors(v.R(),v.C()); g.Get_8_Neighbors(rn,cn); )
//            {
//                if(largest_border_image(rn,cn)!=0 && temp_image(rn,cn)==0)
//                {
//                    temp_image(rn,cn)=255;
//                    bdVoxel v1;
//                    temp_list.AddToRightEnd(v1(rn,cn));
//                    break;
//                }
//            }
//        }
//        
//        temp_list.DeleteLeftEnd();
//    }
//    
//    //output.CopyFrom(largest_border_image);
//
//    
//    // Finally, copy the list of border voxels with conntrast tooutput array.
//    output_contrast_along_border.Set(border_voxels_with_contrast.GetNumberOfElements());
//    bdListIterator< bdVoxel > it;
//    int i = 0;
//    for(it.SetLeftEnd(border_voxels_with_contrast), i = 0; it.IsValid() && i<output_contrast_along_border.GetNumberOfElements(); it.MoveRight(), i++)
//    {
//        output_contrast_along_border[i] = it.GetElement().V();
//        //output(0,0,it.GetElement().R(),it.GetElement().C()) = it.GetElement().V();
//    }
//    
//    //mask.CopyFrom(largest_border_image);
//    
//    return 1;
//    
//
//}



int bdMelanoma::BorderOfMelanoma2(bdImage &segmented_melanoma, bdArray<int> &output_contrast_along_border, bdImage &output)
{
    if(segmented_melanoma.IsEmpty()) return 0;
    
    bdGeometry g;
    g.SetDimensions(1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    
    // Find all borders of the melanoma and store them to an image
    bdImage border_image;
    border_image.SetSize(1,1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    border_image.SetVisualizationPropertiesToMatchInput(segmented_melanoma);
    border_image.FillInWith(0);
    
    output.SetSize(3,1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(segmented_melanoma);
    output.FillInWith(0);
    
    
    for(unsigned int r=0; r<segmented_melanoma.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_melanoma.GetNumberOfColumns(); c++)
        {
            // if pixels belongs to the segmentation (any channel)...
            if(segmented_melanoma(0,0,r,c)!=0 || segmented_melanoma(1,0,r,c)!=0 || segmented_melanoma(2,0,r,c)!=0)
            {
                //... check if it is a border
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
                {
                    if(segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)==0)
                    {
                        border_image(r,c) = 255;
                        break;
                    }
                }
            }
        }
    }
    
    // Extract the largest border image.
    bdBIP bip;
    bdImage largest_border_image;
    bip.ExtractLargest_26_ConnectedComponent(border_image, 1, largest_border_image);
    
    
    
    bdList< bdVoxel > list_of_largest_border_voxels;
    for(unsigned int r=0; r<segmented_melanoma.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_melanoma.GetNumberOfColumns(); c++)
        {
            if(largest_border_image(r,c)!=0)
            {
                unsigned int contrast1 = 0;
                
                //----- Calculate contrast 1-----
                {
                    int rn,cn;
                    
                    unsigned int squared_radius = 0;
                    int is_radius_found = 0;
                    for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius, rn, cn); )
                        {
                            if( (segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)==0) ||
                                (segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)!=0) ||
                                (segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)!=0) )
                            {
                                is_radius_found = 1;
                                contrast1 = squared_radius;
                                break;
                            }
                        }
                    }
                    if(!is_radius_found) contrast1 = g.GetNumberOfCircleElements()-1;
                    
                    if(contrast1>255) output(0,0,r,c) = 255;
                    else output(0,0,r,c) = contrast1;
                    
                    
                }
                //----------
                
                unsigned int contrast2 = 0;
                
                //----- Calculate contrast 2 -----
                {
                    int rn,cn;
                    
                    unsigned int squared_radius = 0;
                    int is_radius_found = 0;
                    //int contrast = 0;
                    for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius, rn, cn); )
                        {
                            if(segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)!=0)
                            {
                                is_radius_found = 1;
                                contrast2 = squared_radius;
                                break;
                            }
                        }
                    }
                    if(!is_radius_found) contrast2 = g.GetNumberOfCircleElements()-1;
                    
                    if(contrast2>255) output(1,0,r,c) = 255;
                    else output(1,0,r,c) = contrast2;
                }
                //----------
                
                
                bdVoxel voxel;
                voxel(r,c);
                voxel.SetValue(contrast1+contrast2); //output(0,0,r,c) + output(2,0,r,c));
                list_of_largest_border_voxels.AddToRightEnd(voxel);
                
                
                if(contrast1+contrast2>255) output(2,0,r,c) = 255;
                else output(2,0,r,c) = contrast1 + contrast2;
                
            }
        }
    }
    
    
    // Finally, copy the list of border voxels with conntrast tooutput array.
    output_contrast_along_border.Set(list_of_largest_border_voxels.GetNumberOfElements());
    bdListIterator< bdVoxel > it;
    int i = 0;
    for(it.SetLeftEnd(list_of_largest_border_voxels), i = 0; it.IsValid() && i<output_contrast_along_border.GetNumberOfElements(); it.MoveRight(), i++)
    {
        output_contrast_along_border[i] = it.GetElement().V();
    }

    
    return 1;
    
}



int bdMelanoma::BorderOfMelanoma3(bdImage &segmented_melanoma, bdImage &mask, bdArray<int> &output_contrast_along_border, unsigned int threshold_squared_distance, bdImage &output)
{
    if(segmented_melanoma.IsEmpty()) return 0;
    
    bdGeometry g;
    g.SetDimensions(1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    
    // Find all borders of the melanoma and store them to an image
    bdImage border_image;
    border_image.SetSize(1,1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    border_image.SetVisualizationPropertiesToMatchInput(segmented_melanoma);
    border_image.FillInWith(0);
    
    output.SetSize(3,1,segmented_melanoma.GetNumberOfRows(),segmented_melanoma.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(segmented_melanoma);
    output.FillInWith(0);
    
    
    for(unsigned int r=0; r<segmented_melanoma.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_melanoma.GetNumberOfColumns(); c++)
        {
            // if pixels belongs to the segmentation (any channel)...
            //if(segmented_melanoma(0,0,r,c)!=0 || segmented_melanoma(1,0,r,c)!=0 || segmented_melanoma(2,0,r,c)!=0)
            if(mask(r,c)!=0)
            {
                //... check if it is a border
                int rn,cn;
                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
                {
                    if(mask(rn,cn)==0)
                    {
                        border_image(r,c) = 255;
                        break;
                    }
                }
            }
        }
    }
    
    // Extract the largest border image.
    bdBIP bip;
    bdImage largest_border_image;
    bip.ExtractLargest_26_ConnectedComponent(border_image, 1, largest_border_image);
    
    
    double n_of_border_pixels = 0;
    double n_of_sharp_border_pixels = 0;
    
    bdList< bdVoxel > list_of_largest_border_voxels;
    for(unsigned int r=0; r<segmented_melanoma.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<segmented_melanoma.GetNumberOfColumns(); c++)
        {
            if(largest_border_image(r,c)!=0)
            {
                unsigned int contrast1 = 0;
                
                //----- Calculate contrast 1-----
                {
                    int rn,cn;
                    
                    unsigned int squared_radius = 0;
                    int is_radius_found = 0;
                    for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius, rn, cn); )
                        {
                            if( (segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)==0) ||
                               (segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)!=0) ||
                               (segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)!=0) )
                            {
                                is_radius_found = 1;
                                contrast1 = squared_radius;
                                break;
                            }
                        }
                    }
                    if(!is_radius_found) contrast1 = g.GetNumberOfCircleElements()-1;
                    
                    //if(contrast1>255) output(0,0,r,c) = 255;
                    //else output(0,0,r,c) = contrast1;
                    
                    
                }
                //----------
                
                unsigned int contrast2 = 0;
                
                //----- Calculate contrast 2 -----
                {
                    int rn,cn;
                    
                    unsigned int squared_radius = 0;
                    int is_radius_found = 0;
                    //int contrast = 0;
                    for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius, rn, cn); )
                        {
                            if(segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)!=0)
                            {
                                is_radius_found = 1;
                                contrast2 = squared_radius;
                                break;
                            }
                        }
                    }
                    if(!is_radius_found) contrast2 = g.GetNumberOfCircleElements()-1;
                    
                    //if(contrast2>255) output(1,0,r,c) = 255;
                    //else output(1,0,r,c) = contrast2;
                }
                //----------
                
                unsigned int contrast3 = 0;
                
                //----- Calculate contrast 3-----
                {
                    int rn,cn;
                    
                    unsigned int squared_radius = 0;
                    int is_radius_found = 0;
                    for(squared_radius = 0; squared_radius<g.GetNumberOfCircleElements()-1 && !is_radius_found; squared_radius++)
                    {
                        for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius, rn, cn); )
                        {
                            if( (segmented_melanoma(0,0,rn,cn)!=0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)==0) ||
                               (segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)!=0 && segmented_melanoma(2,0,rn,cn)==0) ||
                               (segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)!=0) ||
                               (segmented_melanoma(0,0,rn,cn)==0 && segmented_melanoma(1,0,rn,cn)==0 && segmented_melanoma(2,0,rn,cn)==0))
                            {
                                is_radius_found = 1;
                                contrast3 = squared_radius;
                                break;
                            }
                        }
                    }
                    if(!is_radius_found) contrast3 = g.GetNumberOfCircleElements()-1;
                    
                    //if(contrast3>255) output(0,0,r,c) = 255;
                    //else output(2,0,r,c) = contrast3;
                    
                    
                }
                //----------
                
                
                bdVoxel voxel;
                voxel(r,c);
                double sqrt_contrast1 = sqrt ((double)contrast1);
                double sqrt_contrast2 = sqrt ((double)contrast2);
                double sqrt_contrast3 = sqrt ((double)contrast3);
                //double contrast = contrast1 + contrast2 + sqrt((double)(4.0*contrast1*contrast2));
                double contrast = (sqrt_contrast1 + sqrt_contrast2 + sqrt_contrast3) * (sqrt_contrast1 + sqrt_contrast2 + sqrt_contrast3);
                voxel.SetValue(contrast);//contrast1+contrast2); //output(0,0,r,c) + output(2,0,r,c));
                list_of_largest_border_voxels.AddToRightEnd(voxel);
                
                
                //if(contrast1+contrast2>255) output(2,0,r,c) = 255;
                //else output(2,0,r,c) = contrast1 + contrast2;
                
                if(contrast>255) output(0,0,r,c) = 255;
                else output(0,0,r,c) = contrast;
                
                if(contrast<threshold_squared_distance)
                {
                    n_of_sharp_border_pixels += 1;
                    output(2,0,r,c) = 255;
                    output(1,0,r,c) = 255;
                }
                else
                {
                    output(2,0,r,c) = 0;
                    output(1,0,r,c) = 0;
                }
                
                n_of_border_pixels += 1;
                
            }
        }
    }
    
    
    // Finally, copy the list of border voxels with conntrast to output array.
    output_contrast_along_border.Set(list_of_largest_border_voxels.GetNumberOfElements());
    bdListIterator< bdVoxel > it;
    int i = 0;
    for(it.SetLeftEnd(list_of_largest_border_voxels), i = 0; it.IsValid() && i<output_contrast_along_border.GetNumberOfElements(); it.MoveRight(), i++)
    {
        output_contrast_along_border[i] = it.GetElement().V();
    }
    
    
    // Compute the border score
    this->m_border_soft_score = (n_of_sharp_border_pixels / n_of_border_pixels) / 0.125; // 0.125 because it is 1/8, since the original method devides lesion to 8 parts.
    this->m_border_hard_score = (int) (round(this->m_border_soft_score));
    
    
    
    
    return 1;
    
    
}




int bdMelanoma::ColorsOfMelanoma(bdImage &original, bdImage &mask, bdImage &output_mapped_colors_image)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_colors_image.FillInWith(0);
    
    int number_of_colors = 8;

    bdArray< bdColorCode > brown_colors;
    brown_colors.Set(number_of_colors);
    brown_colors[0].m_red = 29; brown_colors[0].m_green = 27; brown_colors[0].m_blue = 26;
    brown_colors[1].m_red = 56; brown_colors[1].m_green = 48; brown_colors[1].m_blue = 44;
    brown_colors[2].m_red = 84; brown_colors[2].m_green = 66; brown_colors[2].m_blue = 57;
    brown_colors[3].m_red = 112; brown_colors[3].m_green = 79; brown_colors[3].m_blue = 63;
    brown_colors[4].m_red = 140; brown_colors[4].m_green = 89; brown_colors[4].m_blue = 64;
    brown_colors[5].m_red = 168; brown_colors[5].m_green = 92; brown_colors[5].m_blue = 55;
    brown_colors[6].m_red = 196; brown_colors[6].m_green = 95; brown_colors[6].m_blue = 45;
    brown_colors[7].m_red = 224; brown_colors[7].m_green = 94; brown_colors[7].m_blue = 29;

    bdArray< bdColorCode > blue_colors;
    blue_colors.Set(number_of_colors);
    blue_colors[0].m_red = 25; blue_colors[0].m_green = 25; blue_colors[0].m_blue = 28;
    blue_colors[1].m_red = 44; blue_colors[1].m_green = 44; blue_colors[1].m_blue = 56;
    blue_colors[2].m_red = 56; blue_colors[2].m_green = 56; blue_colors[2].m_blue = 84;
    blue_colors[3].m_red = 62; blue_colors[3].m_green = 62; blue_colors[3].m_blue = 112;
    blue_colors[4].m_red = 64; blue_colors[4].m_green = 64; blue_colors[4].m_blue = 140;
    blue_colors[5].m_red = 58; blue_colors[5].m_green = 58; blue_colors[5].m_blue = 168;
    blue_colors[6].m_red = 44; blue_colors[6].m_green = 44; blue_colors[6].m_blue = 196;
    blue_colors[7].m_red = 25; blue_colors[7].m_green = 25; blue_colors[7].m_blue = 224;


    bdArray< bdColorCode > gray_colors;
    gray_colors.Set(number_of_colors);
    gray_colors[0].m_red = 28; gray_colors[0].m_green = 28; gray_colors[0].m_blue = 28;
    gray_colors[1].m_red = 56; gray_colors[1].m_green = 56; gray_colors[1].m_blue = 56;
    gray_colors[2].m_red = 84; gray_colors[2].m_green = 84; gray_colors[2].m_blue = 84;
    gray_colors[3].m_red = 112; gray_colors[3].m_green = 112; gray_colors[3].m_blue = 112;
    gray_colors[4].m_red = 140; gray_colors[4].m_green = 140; gray_colors[4].m_blue = 140;
    gray_colors[5].m_red = 168; gray_colors[5].m_green = 168; gray_colors[5].m_blue = 168;
    gray_colors[6].m_red = 196; gray_colors[6].m_green = 196; gray_colors[6].m_blue = 196;
    gray_colors[7].m_red = 224; gray_colors[7].m_green = 224; gray_colors[7].m_blue = 224;
    
    
    /// Iterate over all pixels and find the closest predefined color.
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(0,0,r,c)!=0)
            {

                bdColorCode closest_color;
                double min_distance = 2000000;//initial value is set to be larger than any possible calculated value.

                // iterate over brown colors
                for(unsigned int i=0; i<brown_colors.GetNumberOfElements(); i++)
                {
                    //double d = brown_colors[i].L2(pixel_r[0],pixel_g[0],pixel_b[0]);
                    
                    double d = brown_colors[i].L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                    if(d<min_distance)
                    {
                        closest_color.Set(brown_colors[i].m_red, brown_colors[i].m_green, brown_colors[i].m_blue);
                        min_distance = d;
                    }
                }
                // iterate over blue colors
                for(unsigned int i=0; i<blue_colors.GetNumberOfElements(); i++)
                {
                    //double d = blue_colors[i].L2(pixel_r[0],pixel_g[0],pixel_b[0]);

                    double d = blue_colors[i].L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                    if(d<min_distance)
                    {
                        closest_color.Set(blue_colors[i].m_red, blue_colors[i].m_green, blue_colors[i].m_blue);
                        min_distance = d;
                    }
                }
                // iterate over gray colors
                for(unsigned int i=0; i<gray_colors.GetNumberOfElements(); i++)
                {
                    //double d = gray_colors[i].L2(pixel_r[0],pixel_g[0],pixel_b[0]);
                    
                    double d = gray_colors[i].L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                    if(d<min_distance)
                    {
                        closest_color.Set(gray_colors[i].m_red, gray_colors[i].m_green, gray_colors[i].m_blue);
                        min_distance = d;
                    }
                }

                // assign to the output image pixel value the value of the closest found color
                output_mapped_colors_image(0,0,r,c) = closest_color.m_red;
                output_mapped_colors_image(1,0,r,c) = closest_color.m_green;
                output_mapped_colors_image(2,0,r,c) = closest_color.m_blue;
            }
            // if the original pixel was 0, just assign output value 0.
            else
            {
                output_mapped_colors_image(0,0,r,c) = 0;
                output_mapped_colors_image(1,0,r,c) = 0;
                output_mapped_colors_image(2,0,r,c) = 0;
            }
        }
    }
    
    
    return 1;
}



int bdMelanoma::ColorsOfMelanoma2(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &output_mapped_colors_image)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_colors_image.FillInWith(0);
    
//    bdColorCode color_light_brown; color_light_brown.m_red = 244; color_light_brown.m_green = 94; color_light_brown.m_blue = 29;
//    bdColorCode color_dark_brown; color_dark_brown.m_red = 112; color_dark_brown.m_green = 79; color_dark_brown.m_blue = 63;
//    bdColorCode color_blue_gray; color_blue_gray.m_red = 25; color_blue_gray.m_green = 25; color_blue_gray.m_blue = 224;
//    bdColorCode color_black; color_black.m_red = 40; color_black.m_green = 40; color_black.m_blue = 40;
//    bdColorCode color_white; color_white.m_red = 200; color_white.m_green = 200; color_white.m_blue = 200;

    
    /// Iterate over all pixels and find the closest predefined color.
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(0,0,r,c)!=0)
            {
//                if(c==795 && (r==original.GetNumberOfRows()-1052))
//                {
//                    cout<<"orig=("<<original(0,0,r,c)<<","<<original(1,0,r,c)<<","<<original(2,0,r,c)<<")  ";
//                }
                
                bdColorCode closest_color;
                double min_distance = -1;//initial value
                
                // iterate over brown colors
                for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
                {
                    for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
                    {
                        bdColorCode color;
                        color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
                        
                        //double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                        double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));

                        if(d<min_distance || min_distance<0)
                        {
                            closest_color.Set(color.m_red, color.m_green, color.m_blue);
                            min_distance = d;
                        }
                    }
                }
                
                // iterate over red colors
                for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
                {
                    for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
                    {
                        bdColorCode color;
                        color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
                        
                        //double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                        double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));

                        if(d<min_distance)
                        {
                            closest_color.Set(color.m_red, color.m_green, color.m_blue);
                            min_distance = d;
                        }
                    }
                }
                
                // iterate over blue colors
                for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
                {
                    for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
                    {
                        bdColorCode color;
                        color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
                        
                        //double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                        double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));

                        if(d<min_distance)
                        {
                            closest_color.Set(color.m_red, color.m_green, color.m_blue);
                            min_distance = d;
                        }
                    }
                }
                
                // assign to the output image pixel value the value of the closest found color
                output_mapped_colors_image(0,0,r,c) = closest_color.m_red;
                output_mapped_colors_image(1,0,r,c) = closest_color.m_green;
                output_mapped_colors_image(2,0,r,c) = closest_color.m_blue;
                
//                if(c==795 && (r==original.GetNumberOfRows()-1052))
//                {
//                    cout<<" output=("<<output_mapped_colors_image(0,0,r,c)<<","<<output_mapped_colors_image(1,0,r,c)<<","<<output_mapped_colors_image(2,0,r,c)<<")  ";
//                }

            }
            // if the original pixel was 0, just assign output value 0.
            else
            {
                output_mapped_colors_image(0,0,r,c) = 0;
                output_mapped_colors_image(1,0,r,c) = 0;
                output_mapped_colors_image(2,0,r,c) = 0;
            }
        }
    }
    
    output_mapped_colors_image(0,0,0,0) = 0; output_mapped_colors_image(1,0,0,0) = 0; output_mapped_colors_image(2,0,0,0) = 0;
    output_mapped_colors_image(0,0,0,1) = 255; output_mapped_colors_image(1,0,0,1) = 255; output_mapped_colors_image(2,0,0,1) = 255;

    
    return 1;
}



int bdMelanoma::ColorsOfMelanoma3(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    output_mapped_chart_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_chart_colors_image.FillInWith(0);
    
    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_colors_image.FillInWith(0);
    
    // Build histogram to find the most common GRAYSCALE value.
    bdArray<int> histogram;
    histogram.Set(256);
    histogram.FillInWith(0);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
            histogram[v]++;
        }
    }
    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
    unsigned int last_non_zero_value = 1;
    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
    {
        if(histogram[i]==0)
        {
            histogram[i] = last_non_zero_value;
        }
        else
        {
            last_non_zero_value = histogram[i];
        }
    }
    int start_search_index = 25;
    int most_common_gray_value = start_search_index;
    for(unsigned int r=start_search_index+1; r<histogram.GetNumberOfElements()-1; r++)
    {
        if(histogram[r]>histogram[most_common_gray_value]) most_common_gray_value = r;
    }
    
    unsigned int white_threshold;
    
    //----- Determine WHITE threshold value-----
    {
        // Now that we have the maximum peak, search for the bending (foot of the curve) between the peak and the end of the histogram.
        
        double line_pos1_src[3], line_pos2_src[3];
        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[most_common_gray_value]; line_pos1_src[2] = most_common_gray_value;
        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 255;
        double max_squared_distance = 0;
        unsigned int index_of_max_distance = most_common_gray_value+1;//-1;
        //for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)
        for(unsigned int i=254; i>most_common_gray_value; i--)
        {
            double position_to_project_src[3], projected_position_src[3], n;
            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
            if(n>0 && n<1)
            {
                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
                
                if(d>max_squared_distance)
                {
                    max_squared_distance = d;
                    index_of_max_distance = i;
                }
            }
        }
        
        // The threshold value is the index in of the crossing with x axis
        white_threshold = index_of_max_distance;
        
        cout<<" white_threshold="<<white_threshold<<" ";
    }
    //------------------------------
    
    
    
    /// Iterate over all pixels and find the closest predefined color.
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(0,0,r,c)!=0)
            {
                // Check if the pixel is "WHITE"
                int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
                //if(v>most_common_gray_value)
                if(v>white_threshold)
                {
                    output_mapped_chart_colors_image(0,0,r,c) = 255;
                    output_mapped_chart_colors_image(1,0,r,c) = 255;
                    output_mapped_chart_colors_image(2,0,r,c) = 255;
                    
                    output_mapped_colors_image(0,0,r,c) = 255;
                    output_mapped_colors_image(1,0,r,c) = 255;
                    output_mapped_colors_image(2,0,r,c) = 255;
                }
                else
                {
                    bdColorCode closest_color;
                    bdColorCode closest_color_map;
                    double min_distance = -1;//initial value
                    
                    // iterate over brown colors
                    for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            
                            if(d<min_distance || min_distance<0)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(brown_map(0,0,rn,cn),brown_map(1,0,rn,cn),brown_map(2,0,rn,cn));
                                min_distance = d;
                            }
                        }
                    }
                    
                    // iterate over red colors
                    for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));

                            
                            if(d<min_distance)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(red_map(0,0,rn,cn),red_map(1,0,rn,cn),red_map(2,0,rn,cn));

                                min_distance = d;
                            }
                        }
                    }
                    
                    // iterate over blue colors
                    for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));

                            
                            if(d<min_distance)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(blue_map(0,0,rn,cn),blue_map(1,0,rn,cn),blue_map(2,0,rn,cn));
                                min_distance = d;
                            }
                        }
                    }
                    
                    // assign to the output image pixel value the value of the closest found color
                    output_mapped_chart_colors_image(0,0,r,c) = closest_color.m_red;
                    output_mapped_chart_colors_image(1,0,r,c) = closest_color.m_green;
                    output_mapped_chart_colors_image(2,0,r,c) = closest_color.m_blue;
                    
                    output_mapped_colors_image(0,0,r,c) = closest_color_map.m_red;
                    output_mapped_colors_image(1,0,r,c) = closest_color_map.m_green;
                    output_mapped_colors_image(2,0,r,c) = closest_color_map.m_blue;
                }
            }
            // if the original pixel was 0, just assign output value 0.
            else
            {
                output_mapped_colors_image(0,0,r,c) = 0;
                output_mapped_colors_image(1,0,r,c) = 0;
                output_mapped_colors_image(2,0,r,c) = 0;
            }
        }
    }
    
    output_mapped_colors_image(0,0,0,0) = 0; output_mapped_colors_image(1,0,0,0) = 0; output_mapped_colors_image(2,0,0,0) = 0;
    output_mapped_colors_image(0,0,0,1) = 255; output_mapped_colors_image(1,0,0,1) = 255; output_mapped_colors_image(2,0,0,1) = 255;
    
    
    return 1;
    
}



int bdMelanoma::ColorsOfMelanoma4(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    output_mapped_chart_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_chart_colors_image.FillInWith(0);
    
    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_colors_image.FillInWith(0);
    
    
    // Build histogram to find the most common GRAYSCALE value.
    bdArray<int> histogram;
    histogram.Set(256);
    histogram.FillInWith(0);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
            histogram[v]++;
        }
    }
    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
    unsigned int last_non_zero_value = 1;
    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
    {
        if(histogram[i]==0)
        {
            histogram[i] = last_non_zero_value;
        }
        else
        {
            last_non_zero_value = histogram[i];
        }
    }
    int most_common_gray_value = 0;
    for(unsigned int r=1; r<histogram.GetNumberOfElements()-1; r++)
    {
        if(histogram[r]>histogram[most_common_gray_value]) most_common_gray_value = r;
    }
    
    unsigned int white_threshold;
    
    //----- Determine WHITE threshold value-----
    {
        // Now that we have the maximum peak, search for the bending (foot of the curve) between the peak and the end of the histogram.
        
        double line_pos1_src[3], line_pos2_src[3];
        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[most_common_gray_value]; line_pos1_src[2] = most_common_gray_value;
        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 255;
        double max_squared_distance = 0;
        unsigned int index_of_max_distance = most_common_gray_value+1;//-1;
        //for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)
        for(unsigned int i=254; i>most_common_gray_value; i--)
        {
            double position_to_project_src[3], projected_position_src[3], n;
            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
            if(n>0 && n<1)
            {
                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
                
                if(d>max_squared_distance)
                {
                    max_squared_distance = d;
                    index_of_max_distance = i;
                }
            }
        }
        
        // The threshold value is the index in of the crossing with x axis
        white_threshold = index_of_max_distance;
        
        cout<<" white_threshold="<<white_threshold<<" ";
    }
    //------------------------------

    
    
    
    /// Iterate over all pixels and find the closest predefined color.
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(0,0,r,c)!=0)
            {
                
                // Check if the pixel is "WHITE"
                int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
                //if(v>most_common_gray_value)
                if(v>white_threshold)
                {
                    output_mapped_chart_colors_image(0,0,r,c) = 255;
                    output_mapped_chart_colors_image(1,0,r,c) = 255;
                    output_mapped_chart_colors_image(2,0,r,c) = 255;
                    
                    output_mapped_colors_image(0,0,r,c) = 255;
                    output_mapped_colors_image(1,0,r,c) = 255;
                    output_mapped_colors_image(2,0,r,c) = 255;
                }
                else
                {
                    
                    
                    //                if(c==795 && (r==original.GetNumberOfRows()-1052))
                    //                {
                    //                    cout<<"orig=("<<original(0,0,r,c)<<","<<original(1,0,r,c)<<","<<original(2,0,r,c)<<")  ";
                    //                }
                    
                    bdColorCode closest_color;
                    bdColorCode closest_color_map;
                    double min_distance = 2000000;//initial value
                    
                    // if RED channel is the most prominent...
                    //if(original(0,0,r,c)>original(1,0,r,c) && original(0,0,r,c)>original(2,0,r,c))
                    if(original(0,0,r,c)>original(2,0,r,c))
                    {
                        //... if GREEN is higher than blue, search from brown chart.
                        if(original(1,0,r,c)>original(2,0,r,c))
                        {
                            // iterate over brown colors
                            for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
                            {
                                for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
                                {
                                    bdColorCode color;
                                    color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
                                    
                                    double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    
                                    if(d<min_distance || min_distance<0)
                                    {
                                        closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                        closest_color_map.Set(brown_map(0,0,rn,cn),brown_map(1,0,rn,cn),brown_map(2,0,rn,cn));
                                        min_distance = d;
                                    }
                                }
                            }
                        }
                        //... else BLUE is higher than green, so we search from the red chart.
                        else
                        {
                            // iterate over red colors
                            for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
                            {
                                for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
                                {
                                    bdColorCode color;
                                    color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
                                    
                                    double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                    
                                    
                                    if(d<min_distance)
                                    {
                                        closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                        closest_color_map.Set(red_map(0,0,rn,cn),red_map(1,0,rn,cn),red_map(2,0,rn,cn));
                                        
                                        min_distance = d;
                                    }
                                }
                            }
                        }
                    }
                    
                    // if BLUE channel is the most prominent (or no other color has been assigned) search from the blue chart.
                    //if((original(2,0,r,c)>original(0,0,r,c) && original(2,0,r,c)>original(1,0,r,c)) || min_distance<0)
                    else
                    {
                    
                        // iterate over blue colors
                        for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
                        {
                            for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
                            {
                                bdColorCode color;
                                color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
                                
                                double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                
                                
                                if(d<min_distance)
                                {
                                    closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                    closest_color_map.Set(blue_map(0,0,rn,cn),blue_map(1,0,rn,cn),blue_map(2,0,rn,cn));
                                    min_distance = d;
                                }
                            }
                        }
                    }
                    
                    // assign to the output image pixel value the value of the closest found color
                    output_mapped_chart_colors_image(0,0,r,c) = closest_color.m_red;
                    output_mapped_chart_colors_image(1,0,r,c) = closest_color.m_green;
                    output_mapped_chart_colors_image(2,0,r,c) = closest_color.m_blue;
                    
                    output_mapped_colors_image(0,0,r,c) = closest_color_map.m_red;
                    output_mapped_colors_image(1,0,r,c) = closest_color_map.m_green;
                    output_mapped_colors_image(2,0,r,c) = closest_color_map.m_blue;
                    
                    
                    //                if(c==795 && (r==original.GetNumberOfRows()-1052))
                    //                {
                    //                    cout<<" output=("<<output_mapped_colors_image(0,0,r,c)<<","<<output_mapped_colors_image(1,0,r,c)<<","<<output_mapped_colors_image(2,0,r,c)<<")  ";
                    //                }
                }
                
            }
            // if the original pixel was 0, just assign output value 0.
            else
            {
                output_mapped_colors_image(0,0,r,c) = 0;
                output_mapped_colors_image(1,0,r,c) = 0;
                output_mapped_colors_image(2,0,r,c) = 0;
            }
        }
    }
    
    output_mapped_colors_image(0,0,0,0) = 0; output_mapped_colors_image(1,0,0,0) = 0; output_mapped_colors_image(2,0,0,0) = 0;
    output_mapped_colors_image(0,0,0,1) = 255; output_mapped_colors_image(1,0,0,1) = 255; output_mapped_colors_image(2,0,0,1) = 255;
    
    
    return 1;
    
}



//int bdMelanoma::ColorsOfMelanoma5(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image)
//{
//    if(original.IsEmpty()) return 0;
//    if(mask.IsEmpty()) return 0;
//    
//    output_mapped_chart_colors_image.SetSizeAndPropertiesAs(original);
//    output_mapped_chart_colors_image.FillInWith(0);
//    
//    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
//    output_mapped_colors_image.FillInWith(0);
//    
//    // Build histogram to find the most common GRAYSCALE value.
//    bdArray<int> histogram;
//    histogram.Set(256);
//    histogram.FillInWith(0);
//    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//        {
//            int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
//            histogram[v]++;
//        }
//    }
//    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
//    unsigned int last_non_zero_value = 1;
//    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
//    {
//        if(histogram[i]==0)
//        {
//            histogram[i] = last_non_zero_value;
//        }
//        else
//        {
//            last_non_zero_value = histogram[i];
//        }
//    }
//    int start_search_index = 25;
//    int most_common_gray_value = start_search_index;
//    for(unsigned int r=start_search_index+1; r<histogram.GetNumberOfElements()-1; r++)
//    {
//        if(histogram[r]>histogram[most_common_gray_value]) most_common_gray_value = r;
//    }
//    
//    unsigned int white_threshold;
//    
//    //----- Determine WHITE threshold value-----
//    {
//        // Now that we have the maximum peak, search for the bending (foot of the curve) between the peak and the end of the histogram.
//        
//        double line_pos1_src[3], line_pos2_src[3];
//        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[most_common_gray_value]; line_pos1_src[2] = most_common_gray_value;
//        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 255;
//        double max_squared_distance = 0;
//        unsigned int index_of_max_distance = most_common_gray_value+1;//-1;
//        //for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)
//        for(unsigned int i=254; i>most_common_gray_value; i--)
//        {
//            double position_to_project_src[3], projected_position_src[3], n;
//            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
//            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
//            if(n>0 && n<1)
//            {
//                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
//                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
//                
//                if(d>max_squared_distance)
//                {
//                    max_squared_distance = d;
//                    index_of_max_distance = i;
//                }
//            }
//        }
//        
//        // The threshold value is the index in of the crossing with x axis
//        white_threshold = index_of_max_distance;
//        
//        cout<<" white_threshold="<<white_threshold<<" ";
//    }
//    //------------------------------
//    
//    
//    
//    /// Iterate over all pixels and find the closest predefined color.
//    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//        {
//            // Work only on pixels within the mask.
//            if(mask(0,0,r,c)!=0)
//            {
//                for(int execute_once=1; execute_once; execute_once=0)
//                {
//                    //----- Check if the pixel is "WHITE" -----
//                    {
//                        int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
//                        if(v>white_threshold)
//                        {
//                            output_mapped_chart_colors_image(0,0,r,c) = 255;
//                            output_mapped_chart_colors_image(1,0,r,c) = 255;
//                            output_mapped_chart_colors_image(2,0,r,c) = 255;
//                            
//                            output_mapped_colors_image(0,0,r,c) = 255;
//                            output_mapped_colors_image(1,0,r,c) = 255;
//                            output_mapped_colors_image(2,0,r,c) = 255;
//                            
//                            break;
//                        }
//                    }
//                    //--------------------------------------------
//                    
//                    //----- Check if the pixel is "BLACK" -----
//                    {
//                        if(original(0,0,r,c)<33 && original(1,0,r,c)<33 && original(2,0,r,c)<33)
//                        {
//                            output_mapped_chart_colors_image(0,0,r,c) = 128;
//                            output_mapped_chart_colors_image(1,0,r,c) = 128;
//                            output_mapped_chart_colors_image(2,0,r,c) = 128;
//                            
//                            output_mapped_colors_image(0,0,r,c) = 128;
//                            output_mapped_colors_image(1,0,r,c) = 128;
//                            output_mapped_colors_image(2,0,r,c) = 128;
//                            
//                            break;
//                        }
//                    }
//                    //--------------------------------------------
//                    
//                    //----- Check if the pixel is "BLUE_GRAY" -----
//                    if((original(2,0,r,c)>original(0,0,r,c) && original(2,0,r,c)>original(1,0,r,c) && original(2,0,r,c)<163) ||
//                       (original(2,0,r,c)-original(0,0,r,c)<10 && original(2,0,r,c)-original(0,0,r,c)>-10 && original(0,0,r,c)<128 && original(1,0,r,c)<128 && original(2,0,r,c)<128) )
//                    {
//                        output_mapped_chart_colors_image(0,0,r,c) = 0;
//                        output_mapped_chart_colors_image(1,0,r,c) = 0;
//                        output_mapped_chart_colors_image(2,0,r,c) = 255;
//                        
//                        output_mapped_colors_image(0,0,r,c) = 0;
//                        output_mapped_colors_image(1,0,r,c) = 0;
//                        output_mapped_colors_image(2,0,r,c) = 255;
//                        
//                        break;
//                    }
//                    //----------------------------------------------
//                    
//                    
//                    bdColorCode closest_color;
//                    bdColorCode closest_color_map;
//                    double min_distance = -1;//initial value
//                    
//                    // iterate over brown colors
//                    for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
//                    {
//                        for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
//                        {
//                            bdColorCode color;
//                            color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
//                            
//                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            
//                            if(d<min_distance || min_distance<0)
//                            {
//                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
//                                closest_color_map.Set(brown_map(0,0,rn,cn),brown_map(1,0,rn,cn),brown_map(2,0,rn,cn));
//                                min_distance = d;
//                            }
//                        }
//                    }
//                    
//                    // iterate over red colors
//                    for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
//                    {
//                        for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
//                        {
//                            bdColorCode color;
//                            color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
//                            
//                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            
//                            
//                            if(d<min_distance)
//                            {
//                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
//                                closest_color_map.Set(red_map(0,0,rn,cn),red_map(1,0,rn,cn),red_map(2,0,rn,cn));
//                                
//                                min_distance = d;
//                            }
//                        }
//                    }
//                    
//                    // iterate over blue colors
//                    for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
//                    {
//                        for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
//                        {
//                            bdColorCode color;
//                            color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
//
//                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
//
//
//                            if(d<min_distance)
//                            {
//                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
//                                closest_color_map.Set(blue_map(0,0,rn,cn),blue_map(1,0,rn,cn),blue_map(2,0,rn,cn));
//                                min_distance = d;
//                            }
//                        }
//                    }
//                    
//                    // assign to the output image pixel value the value of the closest found color
//                    output_mapped_chart_colors_image(0,0,r,c) = closest_color.m_red;
//                    output_mapped_chart_colors_image(1,0,r,c) = closest_color.m_green;
//                    output_mapped_chart_colors_image(2,0,r,c) = closest_color.m_blue;
//                    
//                    output_mapped_colors_image(0,0,r,c) = closest_color_map.m_red;
//                    output_mapped_colors_image(1,0,r,c) = closest_color_map.m_green;
//                    output_mapped_colors_image(2,0,r,c) = closest_color_map.m_blue;
//
//                    
//                    
//                }
//                
//                
//            }
//            // if the original pixel was 0, just assign output value 0.
//            else
//            {
//                output_mapped_colors_image(0,0,r,c) = 0;
//                output_mapped_colors_image(1,0,r,c) = 0;
//                output_mapped_colors_image(2,0,r,c) = 0;
//            }
//        }
//    }
//    
//    output_mapped_colors_image(0,0,0,0) = 0; output_mapped_colors_image(1,0,0,0) = 0; output_mapped_colors_image(2,0,0,0) = 0;
//    output_mapped_colors_image(0,0,0,1) = 255; output_mapped_colors_image(1,0,0,1) = 255; output_mapped_colors_image(2,0,0,1) = 255;
//    
//    
//    return 1;
//    
//}



int bdMelanoma::ColorsOfMelanoma5(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image)
{
    if(original.IsEmpty()) return 0;
    if(mask.IsEmpty()) return 0;
    
    output_mapped_chart_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_chart_colors_image.FillInWith(0);
    
    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
    output_mapped_colors_image.FillInWith(0);
    
    // Build histogram to find the most common GRAYSCALE value.
    bdArray<int> histogram;
    histogram.Set(256);
    histogram.FillInWith(0);
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
            histogram[v]++;
        }
    }
    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
    unsigned int last_non_zero_value = 1;
    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
    {
        if(histogram[i]==0)
        {
            histogram[i] = last_non_zero_value;
        }
        else
        {
            last_non_zero_value = histogram[i];
        }
    }
    int start_search_index = 25;
    int most_common_gray_value = start_search_index;
    for(unsigned int r=start_search_index+1; r<histogram.GetNumberOfElements()-1; r++)
    {
        if(histogram[r]>histogram[most_common_gray_value]) most_common_gray_value = r;
    }
    
    unsigned int white_threshold;
    
    //----- Determine WHITE threshold value-----
    {
        // Now that we have the maximum peak, search for the bending (foot of the curve) between the peak and the end of the histogram.
        
        double line_pos1_src[3], line_pos2_src[3];
        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[most_common_gray_value]; line_pos1_src[2] = most_common_gray_value;
        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 255;
        double max_squared_distance = 0;
        unsigned int index_of_max_distance = most_common_gray_value+1;//-1;
        //for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)
        for(unsigned int i=254; i>most_common_gray_value; i--)
        {
            double position_to_project_src[3], projected_position_src[3], n;
            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
            if(n>0 && n<1)
            {
                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
                
                if(d>max_squared_distance)
                {
                    max_squared_distance = d;
                    index_of_max_distance = i;
                }
            }
        }
        
        // The threshold value is the index in of the crossing with x axis
        white_threshold = index_of_max_distance;
        
        cout<<" white_threshold="<<white_threshold<<" ";
    }
    //------------------------------
    
    
    
    bdColorCode dark_brown_color_map;
    dark_brown_color_map.Set(136,0,21);
    bdColorCode light_brown_color_map;
    light_brown_color_map.Set(185,122,87);
    bdColorCode blue_gray_color_map;
    blue_gray_color_map.Set(0,162,232);
    bdColorCode red_color_map;
    red_color_map.Set(237,28,36);
    bdColorCode black_color_map;
    black_color_map.Set(127,127,127);
    bdColorCode white_color_map;
    white_color_map.Set(255,255,255);
    

    bdList<bdColorCode> brown_colors_list;
    bdList<bdColorCode> red_colors_list;
    bdList<bdColorCode> blue_gray_colors_list;
    bdList<bdColorCode> black_colors_list;
    bdList<bdColorCode> white_colors_list;


    
    /// Iterate over all pixels and find the closest predefined color.
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(r,c)!=0)
            {
                for(int execute_once=1; execute_once; execute_once=0)
                {
                    //----- Check if the pixel is "WHITE" -----
                    {
                        int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
                        if(v>white_threshold)
                        {
                            output_mapped_chart_colors_image(0,0,r,c) = 255;
                            output_mapped_chart_colors_image(1,0,r,c) = 255;
                            output_mapped_chart_colors_image(2,0,r,c) = 255;
                            
                            output_mapped_colors_image(0,0,r,c) = white_color_map.m_red;
                            output_mapped_colors_image(1,0,r,c) = white_color_map.m_green;
                            output_mapped_colors_image(2,0,r,c) = white_color_map.m_blue;
                            
                            //m_color_mapped_white_volume++;
                            
                            break;
                        }
                    }
                    //--------------------------------------------
                    
                    //----- Check if the pixel is "BLACK" -----
                    {
                        if(original(0,0,r,c)<33 && original(1,0,r,c)<33 && original(2,0,r,c)<33)
                        {
                            output_mapped_chart_colors_image(0,0,r,c) = 33;
                            output_mapped_chart_colors_image(1,0,r,c) = 33;
                            output_mapped_chart_colors_image(2,0,r,c) = 33;
                            
                            output_mapped_colors_image(0,0,r,c) = black_color_map.m_red;
                            output_mapped_colors_image(1,0,r,c) = black_color_map.m_green;
                            output_mapped_colors_image(2,0,r,c) = black_color_map.m_blue;
                            
                            //m_color_mapped_black_volume++;
                            
                            break;
                        }
                    }
                    //--------------------------------------------
                    
                    //----- Check if the pixel is "BLUE_GRAY" -----
                    if((original(2,0,r,c)>original(0,0,r,c) && original(2,0,r,c)>original(1,0,r,c) && original(2,0,r,c)<163) ||
                       (original(2,0,r,c)-original(0,0,r,c)<10 && original(2,0,r,c)-original(0,0,r,c)>-10 && original(0,0,r,c)<128 && original(1,0,r,c)<128 && original(2,0,r,c)<128) )
                    {
                        output_mapped_colors_image(0,0,r,c) = blue_gray_color_map.m_red;
                        output_mapped_colors_image(1,0,r,c) = blue_gray_color_map.m_green;
                        output_mapped_colors_image(2,0,r,c) = blue_gray_color_map.m_blue;
                        
                        // iterate over blue colors just to complete the output_mapped_chart_colors_image
                        double min_distance = -1;//initial value
                        for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
                        {
                            for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
                            {
                                bdColorCode color;
                                color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
                                
                                double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                                
                                if(d<min_distance || min_distance<0)
                                {
                                    output_mapped_chart_colors_image(0,0,r,c) = color.m_red;
                                    output_mapped_chart_colors_image(1,0,r,c) = color.m_green;
                                    output_mapped_chart_colors_image(2,0,r,c) = color.m_blue;
                                    min_distance = d;
                                }
                            }
                        }
                        
                        bdColorCode closest_color;
                        closest_color.Set(output_mapped_chart_colors_image(0,0,r,c), output_mapped_chart_colors_image(1,0,r,c), output_mapped_chart_colors_image(2,0,r,c));
                        bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(brown_colors_list, closest_color);
                        if(ln!=NULL)
                        {
                            ln->GetElementPointer()->m_intensity += 1;
                        }
                        else
                        {
                            bdColorCode *pcc = brown_colors_list.AddNewToRightEnd();
                            pcc->Set(closest_color.m_red, closest_color.m_green, closest_color.m_blue, 1);
                        }
                        
                        //m_color_mapped_blue_gray_volume++;
                        
                        break;
                    }
                    //----------------------------------------------
                    
                    
                    bdColorCode closest_color;
                    bdColorCode closest_color_map;
                    double min_distance = -1;//initial value
                    bdList<bdColorCode> *pl = NULL; //pointer to the list from which the color was mapped.
                    //int indicator_of_mapped_color = 0;// use 1 for brown, 2 for red, 3 for blue-gray.
                    
                    // iterate over brown colors
                    for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            
                            if(d<min_distance || min_distance<0)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(brown_map(0,0,rn,cn),brown_map(1,0,rn,cn),brown_map(2,0,rn,cn));
                                min_distance = d;
                                pl = &brown_colors_list;
                                //indicator_of_mapped_color = 1;
                            }
                        }
                    }
                    
                    // iterate over red colors
                    for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            
                            
                            if(d<min_distance)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(red_map(0,0,rn,cn),red_map(1,0,rn,cn),red_map(2,0,rn,cn));
                                min_distance = d;
                                pl = &red_colors_list;
                                //indicator_of_mapped_color = 2;
                            }
                        }
                    }
                    
                    // iterate over blue colors
                    for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
                    {
                        for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
                        {
                            bdColorCode color;
                            color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
                            
                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
                            
                            
                            if(d<min_distance)
                            {
                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
                                closest_color_map.Set(blue_map(0,0,rn,cn),blue_map(1,0,rn,cn),blue_map(2,0,rn,cn));
                                min_distance = d;
                                pl = &blue_gray_colors_list;
                                //indicator_of_mapped_color = 3;
                            }
                        }
                    }
                    
                    // assign to the output image pixel value the value of the closest found color
                    output_mapped_chart_colors_image(0,0,r,c) = closest_color.m_red;
                    output_mapped_chart_colors_image(1,0,r,c) = closest_color.m_green;
                    output_mapped_chart_colors_image(2,0,r,c) = closest_color.m_blue;
                    
                    output_mapped_colors_image(0,0,r,c) = closest_color_map.m_red;
                    output_mapped_colors_image(1,0,r,c) = closest_color_map.m_green;
                    output_mapped_colors_image(2,0,r,c) = closest_color_map.m_blue;
                    
                    if(pl)
                    {
                        bdListNode<bdColorCode> *ln = bdListing::GetNodeWithElement(*pl, closest_color);
                        if(ln!=NULL)
                        {
                            (ln->GetElementPointer()->m_intensity)++;
                        }
                        else
                        {
                            bdColorCode *pcc = pl->AddNewToRightEnd();
                            pcc->Set(closest_color.m_red, closest_color.m_green, closest_color.m_blue, 1);
                        }
                    }
                    
                }
            }
            // if the original pixel was 0, just assign output value 0.
            else
            {
                output_mapped_colors_image(0,0,r,c) = 0;
                output_mapped_colors_image(1,0,r,c) = 0;
                output_mapped_colors_image(2,0,r,c) = 0;
            }
        }
    }

    
    //----- Calculate the color scores -----
    this->m_color_blue_gray_area = blue_gray_colors_list.GetNumberOfElements();
    this->m_color_brown_area = brown_colors_list.GetNumberOfElements();
    this->m_color_red_area = red_colors_list.GetNumberOfElements();
    bdListIterator<bdColorCode> itcc;
    
    this->m_color_blue_gray_volume = 0;
    for(itcc.SetLeftEnd(blue_gray_colors_list); itcc.IsValid(); itcc.MoveRight())
    {
        this->m_color_blue_gray_volume += itcc.GetElementPointer()->m_intensity;
    }
    
    this->m_color_brown_volume = 0;
    for(itcc.SetLeftEnd(brown_colors_list); itcc.IsValid(); itcc.MoveRight())
    {
        this->m_color_brown_volume += itcc.GetElementPointer()->m_intensity;
    }

    this->m_color_red_volume = 0;
    for(itcc.SetLeftEnd(red_colors_list); itcc.IsValid(); itcc.MoveRight())
    {
        this->m_color_red_volume += itcc.GetElementPointer()->m_intensity;
    }

    // Calculte mapped colors
    m_color_mapped_black_volume = 0;
    m_color_mapped_white_volume = 0;
    m_color_mapped_light_brown_volume = 0;
    m_color_mapped_dark_brown_volume = 0;
    m_color_mapped_red_volume = 0;
    m_color_mapped_blue_gray_volume = 0;
    
    /// Iterate over all pixels and find the closest predefined color.
    double total_area = 0.1;
    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
        {
            // Work only on pixels within the mask.
            if(mask(r,c)!=0)
            {
                total_area++;
                
                if(output_mapped_colors_image(0,0,r,c)==dark_brown_color_map.m_red && output_mapped_colors_image(1,0,r,c)==dark_brown_color_map.m_green && output_mapped_colors_image(2,0,r,c)==dark_brown_color_map.m_blue)
                {
                    m_color_mapped_dark_brown_volume++;
                }
                if(output_mapped_colors_image(0,0,r,c)==light_brown_color_map.m_red && output_mapped_colors_image(1,0,r,c)==light_brown_color_map.m_green && output_mapped_colors_image(2,0,r,c)==light_brown_color_map.m_blue)
                {
                    m_color_mapped_light_brown_volume++;
                }
                if(output_mapped_colors_image(0,0,r,c)==blue_gray_color_map.m_red && output_mapped_colors_image(1,0,r,c)==blue_gray_color_map.m_green && output_mapped_colors_image(2,0,r,c)==blue_gray_color_map.m_blue)
                {
                    m_color_mapped_blue_gray_volume++;
                }
                if(output_mapped_colors_image(0,0,r,c)==red_color_map.m_red && output_mapped_colors_image(1,0,r,c)==red_color_map.m_green && output_mapped_colors_image(2,0,r,c)==red_color_map.m_blue)
                {
                    m_color_mapped_red_volume++;
                }
                if(output_mapped_colors_image(0,0,r,c)==black_color_map.m_red && output_mapped_colors_image(1,0,r,c)==black_color_map.m_green && output_mapped_colors_image(2,0,r,c)==black_color_map.m_blue)
                {
                    m_color_mapped_black_volume++;
                }
                if(output_mapped_colors_image(0,0,r,c)==white_color_map.m_red && output_mapped_colors_image(1,0,r,c)==white_color_map.m_green && output_mapped_colors_image(2,0,r,c)==white_color_map.m_blue)
                {
                    m_color_mapped_white_volume++;
                }
            }
        }
    }
    
    m_color_mapped_dark_brown_volume = m_color_mapped_dark_brown_volume / total_area;
    m_color_mapped_light_brown_volume = m_color_mapped_light_brown_volume / total_area;
    m_color_mapped_blue_gray_volume = m_color_mapped_blue_gray_volume / total_area;
    m_color_mapped_red_volume = m_color_mapped_red_volume / total_area;
    m_color_mapped_black_volume = m_color_mapped_black_volume / total_area;
    m_color_mapped_white_volume = m_color_mapped_white_volume / total_area;
    
    
    m_color_hard_score = 0;
    if(m_color_mapped_dark_brown_volume>=0.1) m_color_hard_score += 1;
    if(m_color_mapped_light_brown_volume>=0.1) m_color_hard_score += 1;
    if(m_color_mapped_blue_gray_volume>=0.1) m_color_hard_score += 1;
    if(m_color_mapped_red_volume>=0.1) m_color_hard_score += 1;
    if(m_color_mapped_black_volume>=0.1) m_color_hard_score += 1;
    if(m_color_mapped_white_volume>=0.1) m_color_hard_score += 1;
    
    m_color_soft_score = 0;
    if(m_color_mapped_dark_brown_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_dark_brown_volume * 10.0;
    if(m_color_mapped_light_brown_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_light_brown_volume * 10.0;
    if(m_color_mapped_blue_gray_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_blue_gray_volume * 10.0;
    if(m_color_mapped_red_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_red_volume * 10.0;
    if(m_color_mapped_black_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_black_volume * 10.0;
    if(m_color_mapped_white_volume>=0.1) m_color_soft_score += 1;
    else m_color_soft_score += m_color_mapped_white_volume * 10.0;
    //----------
    
   
    
    
    return 1;
    
}




//int bdMelanoma::ColorsOfMelanoma6(bdImage &original, bdImage &mask, bdImage &red_chart, bdImage &brown_chart, bdImage &blue_chart, bdImage &red_map, bdImage &brown_map, bdImage &blue_map, bdImage &output_mapped_chart_colors_image, bdImage &output_mapped_colors_image)
//{
//    if(original.IsEmpty()) return 0;
//    if(mask.IsEmpty()) return 0;
//
//    output_mapped_chart_colors_image.SetSizeAndPropertiesAs(original);
//    output_mapped_chart_colors_image.FillInWith(0);
//
//    output_mapped_colors_image.SetSizeAndPropertiesAs(original);
//    output_mapped_colors_image.FillInWith(0);
//
//    // Build histogram to find the most common GRAYSCALE value.
//    bdArray<int> histogram;
//    histogram.Set(256);
//    histogram.FillInWith(0);
//    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//        {
//            int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
//            histogram[v]++;
//        }
//    }
//    // Modify the histogram to NOT contain 0 values after the first non-zero value is found.
//    unsigned int last_non_zero_value = 1;
//    for(unsigned int i=0; i<histogram.GetNumberOfElements(); i++)
//    {
//        if(histogram[i]==0)
//        {
//            histogram[i] = last_non_zero_value;
//        }
//        else
//        {
//            last_non_zero_value = histogram[i];
//        }
//    }
//    int start_search_index = 25;
//    int most_common_gray_value = start_search_index;
//    for(unsigned int r=start_search_index+1; r<histogram.GetNumberOfElements()-1; r++)
//    {
//        if(histogram[r]>histogram[most_common_gray_value]) most_common_gray_value = r;
//    }
//
//    unsigned int white_threshold;
//
//    //----- Determine WHITE threshold value-----
//    {
//        // Now that we have the maximum peak, search for the bending (foot of the curve) between the peak and the end of the histogram.
//
//        double line_pos1_src[3], line_pos2_src[3];
//        line_pos1_src[0] = 0; line_pos1_src[1] = histogram[most_common_gray_value]; line_pos1_src[2] = most_common_gray_value;
//        line_pos2_src[0] = 0; line_pos2_src[1] = 0; line_pos2_src[2] = 255;
//        double max_squared_distance = 0;
//        unsigned int index_of_max_distance = most_common_gray_value+1;//-1;
//        //for(unsigned int i=index_of_max_value1-1; i>index_of_max_value2; i--)
//        for(unsigned int i=254; i>most_common_gray_value; i--)
//        {
//            double position_to_project_src[3], projected_position_src[3], n;
//            position_to_project_src[0] = 0; position_to_project_src[1] = histogram[i]; position_to_project_src[2] = i;
//            this->ProjectPositionToLine(line_pos1_src, line_pos2_src, position_to_project_src, projected_position_src, n);
//            if(n>0 && n<1)
//            {
//                double d = (projected_position_src[1]-position_to_project_src[1])*(projected_position_src[1]-position_to_project_src[1]) +
//                (projected_position_src[2]-position_to_project_src[2])*(projected_position_src[2]-position_to_project_src[2]);
//
//                if(d>max_squared_distance)
//                {
//                    max_squared_distance = d;
//                    index_of_max_distance = i;
//                }
//            }
//        }
//
//        // The threshold value is the index in of the crossing with x axis
//        white_threshold = index_of_max_distance;
//
//        cout<<" white_threshold="<<white_threshold<<" ";
//    }
//    //------------------------------
//
//
//
//    /// Iterate over all pixels and find the closest predefined color.
//    for(unsigned int r=0; r<original.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<original.GetNumberOfColumns(); c++)
//        {
//            // Work only on pixels within the mask.
//            if(mask(0,0,r,c)!=0)
//            {
//                for(int execute_once=1; execute_once; execute_once=0)
//                {
//                    //----- Check if the pixel is "WHITE" -----
//                    {
//                        int v = (original(0,0,r,c)+original(1,0,r,c)+original(2,0,r,c)) / 3;
//                        if(v>white_threshold)
//                        {
//                            output_mapped_chart_colors_image(0,0,r,c) = 255;
//                            output_mapped_chart_colors_image(1,0,r,c) = 255;
//                            output_mapped_chart_colors_image(2,0,r,c) = 255;
//
//                            output_mapped_colors_image(0,0,r,c) = 255;
//                            output_mapped_colors_image(1,0,r,c) = 255;
//                            output_mapped_colors_image(2,0,r,c) = 255;
//
//                            break;
//                        }
//                    }
//                    //--------------------------------------------
//
//                    //----- Check if the pixel is "BLACK" -----
//                    {
//                        if(original(0,0,r,c)<33 && original(1,0,r,c)<33 && original(2,0,r,c)<33)
//                        {
//                            output_mapped_chart_colors_image(0,0,r,c) = 128;
//                            output_mapped_chart_colors_image(1,0,r,c) = 128;
//                            output_mapped_chart_colors_image(2,0,r,c) = 128;
//
//                            output_mapped_colors_image(0,0,r,c) = 128;
//                            output_mapped_colors_image(1,0,r,c) = 128;
//                            output_mapped_colors_image(2,0,r,c) = 128;
//
//                            break;
//                        }
//                    }
//                    //--------------------------------------------
//
//                    //----- Check if the pixel is "BLUE_GRAY" -----
//                    if((original(2,0,r,c)>original(0,0,r,c) && original(2,0,r,c)>original(1,0,r,c) && original(2,0,r,c)<163) ||
//                       (original(2,0,r,c)-original(0,0,r,c)<10 && original(2,0,r,c)-original(0,0,r,c)>-10 && original(0,0,r,c)<128 && original(1,0,r,c)<128 && original(2,0,r,c)<128) )
//                    {
//                        output_mapped_chart_colors_image(0,0,r,c) = 0;
//                        output_mapped_chart_colors_image(1,0,r,c) = 0;
//                        output_mapped_chart_colors_image(2,0,r,c) = 255;
//
//                        output_mapped_colors_image(0,0,r,c) = 0;
//                        output_mapped_colors_image(1,0,r,c) = 0;
//                        output_mapped_colors_image(2,0,r,c) = 255;
//
//                        break;
//                    }
//                    //----------------------------------------------
//                    
//                    
//                    //----- Check if the pixel is "UNCLASSIFIED" (green) -----
//                    if(original(0,0,r,c)>original(1,0,r,c) && original(2,0,r,c)>original(1,0,r,c))
//                    {
//                        output_mapped_chart_colors_image(0,0,r,c) = 0;
//                        output_mapped_chart_colors_image(1,0,r,c) = 255;
//                        output_mapped_chart_colors_image(2,0,r,c) = 0;
//                        
//                        output_mapped_colors_image(0,0,r,c) = 0;
//                        output_mapped_colors_image(1,0,r,c) = 255;
//                        output_mapped_colors_image(2,0,r,c) = 0;
//                        
//                        break;
//                    }
//                    //----------------------------------------------
//                    
//                    
//                    //----- Check if the pixel is "BROWN" -----
//                    if(original(0,0,r,c)>original(1,0,r,c) && original(1,0,r,c)>original(2,0,r,c))
//                    {
//                        int diff_R_B = original(0,0,r,c) - original(2,0,r,c);
//                        int diff_G_B = original(1,0,r,c) - original(2,0,r,c);
//                        int projected_green_value = (diff_R_B+1) / 2;
//                        int allowed_error = (projected_green_value+1) / 2;
//                        
//                        
//                        
//                        //if( (diff_R_B - (2*diff_G_B) > -5) && (diff_R_B - (2*diff_G_B) < 5) )
//                        if( (diff_G_B>projected_green_value-allowed_error) && (diff_G_B<projected_green_value+allowed_error) )
//                        {
//                            output_mapped_chart_colors_image(0,0,r,c) = 160;
//                            output_mapped_chart_colors_image(1,0,r,c) = 80;
//                            output_mapped_chart_colors_image(2,0,r,c) = 0;
//                            
//                            output_mapped_colors_image(0,0,r,c) = 160;
//                            output_mapped_colors_image(1,0,r,c) = 80;
//                            output_mapped_colors_image(2,0,r,c) = 0;
//                            
//                            break;
//                        }
//                    }
//                    else
//                    {
//                        output_mapped_chart_colors_image(0,0,r,c) = 250;
//                        output_mapped_chart_colors_image(1,0,r,c) = 00;
//                        output_mapped_chart_colors_image(2,0,r,c) = 0;
//                        
//                        output_mapped_colors_image(0,0,r,c) = 250;
//                        output_mapped_colors_image(1,0,r,c) = 0;
//                        output_mapped_colors_image(2,0,r,c) = 0;
//                    }
//                    //----------------------------------------------
//
//
//
//
////                    bdColorCode closest_color;
////                    bdColorCode closest_color_map;
////                    double min_distance = -1;//initial value
////
////                    // iterate over brown colors
////                    for(unsigned int rn = 0; rn < brown_chart.GetNumberOfRows(); rn++)
////                    {
////                        for(unsigned int cn = 0; cn < brown_chart.GetNumberOfColumns(); cn++)
////                        {
////                            bdColorCode color;
////                            color.m_red = brown_chart(0,0,rn,cn); color.m_green = brown_chart(1,0,rn,cn); color.m_blue = brown_chart(2,0,rn,cn);
////
////                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////
////                            if(d<min_distance || min_distance<0)
////                            {
////                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
////                                closest_color_map.Set(brown_map(0,0,rn,cn),brown_map(1,0,rn,cn),brown_map(2,0,rn,cn));
////                                min_distance = d;
////                            }
////                        }
////                    }
////
////                    // iterate over red colors
////                    for(unsigned int rn = 0; rn < red_chart.GetNumberOfRows(); rn++)
////                    {
////                        for(unsigned int cn = 0; cn < red_chart.GetNumberOfColumns(); cn++)
////                        {
////                            bdColorCode color;
////                            color.m_red = red_chart(0,0,rn,cn); color.m_green = red_chart(1,0,rn,cn); color.m_blue = red_chart(2,0,rn,cn);
////
////                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////
////
////                            if(d<min_distance)
////                            {
////                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
////                                closest_color_map.Set(red_map(0,0,rn,cn),red_map(1,0,rn,cn),red_map(2,0,rn,cn));
////
////                                min_distance = d;
////                            }
////                        }
////                    }
////
////                    // iterate over blue colors
////                    for(unsigned int rn = 0; rn < blue_chart.GetNumberOfRows(); rn++)
////                    {
////                        for(unsigned int cn = 0; cn < blue_chart.GetNumberOfColumns(); cn++)
////                        {
////                            bdColorCode color;
////                            color.m_red = blue_chart(0,0,rn,cn); color.m_green = blue_chart(1,0,rn,cn); color.m_blue = blue_chart(2,0,rn,cn);
////
////                            double d = color.L2(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.Lmax(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////                            //double d = color.L2order(original(0,0,r,c),original(1,0,r,c),original(2,0,r,c));
////
////
////                            if(d<min_distance)
////                            {
////                                closest_color.Set(color.m_red, color.m_green, color.m_blue);
////                                closest_color_map.Set(blue_map(0,0,rn,cn),blue_map(1,0,rn,cn),blue_map(2,0,rn,cn));
////                                min_distance = d;
////                            }
////                        }
////                    }
////
////                    // assign to the output image pixel value the value of the closest found color
////                    output_mapped_chart_colors_image(0,0,r,c) = closest_color.m_red;
////                    output_mapped_chart_colors_image(1,0,r,c) = closest_color.m_green;
////                    output_mapped_chart_colors_image(2,0,r,c) = closest_color.m_blue;
////
////                    output_mapped_colors_image(0,0,r,c) = closest_color_map.m_red;
////                    output_mapped_colors_image(1,0,r,c) = closest_color_map.m_green;
////                    output_mapped_colors_image(2,0,r,c) = closest_color_map.m_blue;
//
//
//
//                }
//
//
//            }
//            // if the original pixel was 0, just assign output value 0.
//            else
//            {
//                output_mapped_colors_image(0,0,r,c) = 0;
//                output_mapped_colors_image(1,0,r,c) = 0;
//                output_mapped_colors_image(2,0,r,c) = 0;
//            }
//        }
//    }
//
//    output_mapped_colors_image(0,0,0,0) = 0; output_mapped_colors_image(1,0,0,0) = 0; output_mapped_colors_image(2,0,0,0) = 0;
//    output_mapped_colors_image(0,0,0,1) = 255; output_mapped_colors_image(1,0,0,1) = 255; output_mapped_colors_image(2,0,0,1) = 255;
//
//
//    return 1;
//
//}





int bdMelanoma::Quantify_D_Feature(bdImage &mask, bdImage &segmented_melanoma, bdImage &dark_structures_image, bdImage &bright_structures_image, bdImage &streaks_image, bdImage &structureless_area_image)
{
    bdImage dark_structures_image2;
    bdImage bright_structures_image2;
    bdImage streaks_image2;
    bdImage structureless_area_image2;
    
    bdBIP bip;
    bip.Closing_Sphere(dark_structures_image, dark_structures_image2, 25);
    bip.Closing_Sphere(bright_structures_image, bright_structures_image2, 25);
    bip.Closing_Sphere(streaks_image, streaks_image2, 25);
    bip.Closing_Sphere(structureless_area_image, structureless_area_image2, 25);
    
    double area_total = 0;
    double area_globules = 0;
    double area_dots = 0;
    double area_network = 0;
    double area_structureless = 0;
    double area_streaks = 0;
    
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                //We examine the structures only at the part of the lesion that is segmented in all 3 channels.
                if(segmented_melanoma(0,0,r,c)!=0 && segmented_melanoma(1,0,r,c)!=0 && segmented_melanoma(2,0,r,c)!=0)
                {
                    area_total++;
                    if(dark_structures_image2(r,c)==255) area_globules++;
                    if(dark_structures_image2(r,c)==127) area_dots++;
                    if(bright_structures_image2(r,c)!=0) area_network++;
                    if(streaks_image2(r,c)!=0) area_streaks++;
                    if(structureless_area_image2(r,c)!=0) area_structureless++;
                }
            }
        }
    }
    
    m_DS_dots = area_dots / area_total;
    m_DS_globules = area_globules / area_total;
    m_DS_streaks = area_streaks / area_total;
    m_DS_structureless_area = area_structureless / area_total;
    m_DS_network = area_network / area_total;
    
    
    m_DS_soft_score = 0;
    if(m_DS_dots>0.1) m_DS_soft_score += 1;
    else m_DS_soft_score += m_DS_dots;
    if(m_DS_globules>0.1) m_DS_soft_score += 1;
    else m_DS_soft_score += m_DS_globules;
    if(m_DS_streaks>0.1) m_DS_soft_score += 1;
    else m_DS_soft_score += m_DS_streaks;
    if(m_DS_structureless_area>0.1) m_DS_soft_score += 1;
    else m_DS_soft_score += m_DS_structureless_area;
    if(m_DS_network>0.1) m_DS_soft_score += 1;
    else m_DS_soft_score += m_DS_network;

    
    m_DS_hard_score = 0;
    if(m_DS_dots>0.1) m_DS_hard_score += 1;
    if(m_DS_globules>0.1) m_DS_hard_score += 1;
    if(m_DS_streaks>0.1) m_DS_hard_score += 1;
    if(m_DS_structureless_area>0.1) m_DS_hard_score += 1;
    if(m_DS_network>0.1) m_DS_hard_score += 1;
    

    
    return 1;
}





int bdMelanoma::DiameterOfMelanoma(bdImage &border_image, double *output_min_diameter_in_pixels, double *output_max_diameter_in_pixels, double *output_mean_diameter_in_pixels)
{
    if(border_image.IsEmpty()) return 0;

    
    // Calculate the center of mass for border of melanoma
    double center_of_mass_r = 0, center_of_mass_c = 0, n_of_pixels = 0;
    for(unsigned int r=0; r<border_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<border_image.GetNumberOfColumns(); c++)
        {
            if(border_image(r,c)!=0)
            {
                center_of_mass_c += c;
                center_of_mass_r += r;
                n_of_pixels += 1;
            }
        }
    }
    center_of_mass_c = center_of_mass_c / n_of_pixels;
    center_of_mass_r = center_of_mass_r / n_of_pixels;
    
    // Calculate the meinimum, maximum and average distance of each pixel on the border from the center of mass
    *output_min_diameter_in_pixels = 1000000;
    *output_max_diameter_in_pixels = 0;
    *output_mean_diameter_in_pixels = 0;
    for(unsigned int r=0; r<border_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<border_image.GetNumberOfColumns(); c++)
        {
            if(border_image(r,c)!=0)
            {
                double d = sqrt( (r-center_of_mass_r)*(r-center_of_mass_r) + (c-center_of_mass_c)*(c-center_of_mass_c) );
                if(*output_min_diameter_in_pixels > d) *output_min_diameter_in_pixels = d;
                if(*output_max_diameter_in_pixels < d) *output_max_diameter_in_pixels = d;
                *output_mean_diameter_in_pixels = *output_mean_diameter_in_pixels + d;
            }
        }
    }

    *output_mean_diameter_in_pixels = *output_mean_diameter_in_pixels / n_of_pixels;
    
    
    *output_min_diameter_in_pixels = 2* (*output_min_diameter_in_pixels);
    *output_max_diameter_in_pixels = 2* (*output_max_diameter_in_pixels);
    *output_mean_diameter_in_pixels = 2* (*output_mean_diameter_in_pixels);
    
    return 1;

}


int bdMelanoma::SubtractImages(bdImage &image1, bdImage &image2, bdImage &mask, bdImage &output_image)
{
    if(!image1.IsEqualSizeAs_2D(image2)) return 0;
    
    output_image.SetSizeAndPropertiesAs(image1);

    for(unsigned int t=0; t<image1.GetNumberOfTimeSeries() && t<image2.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<image1.GetNumberOfSlices() && s<image2.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<image1.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<image1.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0)
                    {
                        int diff = image1(t,s,r,c) - image2(t,s,r,c);
                        if(diff<0) diff = -diff;
                        output_image(t,s,r,c) = diff;
                    }
                }
            }
        }
    }
    return 1;
}



int bdMelanoma::MakeCombinedSegmentedImage(bdImage &mask, bdImage &melanoma_bright_structure, bdImage &melanoma_dark_structure, bdImage &melanoma_dark_streaks, bdImage &output_combined_image)
{
    if(!mask.IsEqualSizeAs_2D(melanoma_bright_structure)) return 0;
    if(!mask.IsEqualSizeAs_2D(melanoma_dark_structure)) return 0;
    if(!mask.IsEqualSizeAs_2D(melanoma_dark_streaks)) return 0;
    
    output_combined_image.SetSize(3, 1, mask.GetNumberOfRows(), mask.GetNumberOfColumns());
    output_combined_image.SetVisualizationPropertiesToMatchInput(mask);
    output_combined_image.FillInWith(0);
    
    for(unsigned int r=0; r<mask.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<mask.GetNumberOfColumns(); c++)
        {
            if(mask(r,c)!=0)
            {
                
                if(melanoma_bright_structure(r,c)!=0) output_combined_image(1,0,r,c) = 255;

                if(melanoma_dark_structure(r,c)!=0) output_combined_image(2,0,r,c) = melanoma_dark_structure(r,c);
                
                if(melanoma_bright_structure(r,c)==0 && melanoma_dark_structure(r,c)==0) output_combined_image(0,0,r,c) = 255;

                if(melanoma_dark_streaks(r,c)!=0)
                {
                    output_combined_image(0,0,r,c) = 255; output_combined_image(1,0,r,c) = 255; output_combined_image(2,0,r,c) = 255;
                }
//                if(melanoma_bright_structure(r,c)!=0) output_combined_image(1,0,r,c) = 255;
//                else
//                {
//                    if(melanoma_dark_structure(r,c)!=0) output_combined_image(2,0,r,c) = 255;
//                    else output_combined_image(0,0,r,c) = 255;
//                }
            }
        }
    }
    
    return 1;

}


int bdMelanoma::Extract_RedImage(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            output_image(r,c) = white_image(0,0,r,c);
        }
    }
    
    return 1;
}


int bdMelanoma::Extract_GreenImage(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            output_image(r,c) = white_image(1,0,r,c);
        }
    }
    
    return 1;
}


int bdMelanoma::Extract_BlueImage(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            output_image(r,c) = white_image(2,0,r,c);
        }
    }
    
    return 1;
}


int bdMelanoma::Extract_MaxOfGreenAndBlueImage(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            if(white_image(1,0,r,c) > white_image(2,0,r,c)) output_image(r,c) = white_image(1,0,r,c);
            else output_image(r,c) = white_image(2,0,r,c);
        }
    }
    
    return 1;
}


int bdMelanoma::Extract_MaxOfRGBImage(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            // MAX
            if(white_image(1,0,r,c) > white_image(2,0,r,c)) output_image(r,c) = white_image(1,0,r,c);
            else output_image(r,c) = white_image(2,0,r,c);
            if(white_image(0,0,r,c) > output_image(r,c)) output_image(r,c) = white_image(0,0,r,c);
            
            
            // MEDIAN
//            if(white_image(0,0,r,c) > white_image(1,0,r,c))
//            {
//                if(white_image(0,0,r,c) > white_image(2,0,r,c)) output_image(r,c) = white_image(2,0,r,c);
//                else output_image(r,c) = white_image(0,0,r,c);
//            }
//            else
//            {
//                if(white_image(1,0,r,c) > white_image(2,0,r,c)) output_image(r,c) = white_image(2,0,r,c);
//                else output_image(r,c) = white_image(1,0,r,c);
//            }
            
        }
    }
    
    return 1;
}


// THIS VERSION DOES NOT GIVE SATISFACTORY RESULT
//int bdMelanoma::Extract_MaxOfRGBImage(bdImage &white_image, bdImage &output_image)
//{
//    if(white_image.IsEmpty()) return 0;
//    
//    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
//    output_image.SetVisualizationPropertiesToMatchInput(white_image);
//    
//    bdGeometry g;
//    g.SetDimensions(1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
//    
//    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r=r+2)//r++)
//    {
//        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c=c+2)//++)
//        {
//            int max_diff = -1;
//            int index_of_max_diff = 0;
////            bdArray<int> sum, n_of_pixels;
////            sum.Set(white_image.GetNumberOfTimeSeries());
////            n_of_pixels.Set(white_image.GetNumberOfTimeSeries());
////            sum.FillInWith(0);
////            n_of_pixels.FillInWith(0);
//            for(unsigned int t=0; t<white_image.GetNumberOfTimeSeries(); t++)
//            {
//                int sum = 0;
//                int n_of_pixels = 0;
//                
//                int rn,cn;
//                for(g.ForCoordinates_8_Neighbors(r,c); g.Get_8_Neighbors(rn,cn); )
//                {
////                    sum[t] += white_image(t,0,rn,cn);
////                    (n_of_pixels[t])++;
//                    sum += white_image(t,0,rn,cn);
//                    n_of_pixels++;
//                }
//                int diff = white_image(t,0,r,c) - (sum/n_of_pixels);
//                if(diff<0) diff = -diff;
//                if(diff>max_diff)
//                {
//                    max_diff = diff;
//                    index_of_max_diff = t;
//                }
//            }
//            
//            //output_image(r,c) = white_image(index_of_max_diff,0,r,c);
//            
//            {
//                int rn,cn;
//                for(g.ForCoordinates_9_Neighbors(r,c); g.Get_9_Neighbors(rn,cn); )
//                {
//                    output_image(rn,cn) = white_image(index_of_max_diff,0,rn,cn);
//                }
//            }
//
//        }
//    }
//    
//    return 1;
//}




int bdMelanoma::Extract_DifferenceOfRGB_Image(bdImage &white_image, bdImage &output_image)
{
    if(white_image.IsEmpty()) return 0;
    
    output_image.SetSize(1, 1, white_image.GetNumberOfRows(), white_image.GetNumberOfColumns());
    output_image.SetVisualizationPropertiesToMatchInput(white_image);
    
    for(unsigned int r=0; r<white_image.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<white_image.GetNumberOfColumns(); c++)
        {
            // new value is |red-green| + |green-blue| + |blue-red|.
            int d_rg = white_image(0,0,r,c) - white_image(1,0,r,c);
            if(d_rg<0) d_rg = -d_rg;
            int d_gb = white_image(1,0,r,c) - white_image(2,0,r,c);
            if(d_gb<0) d_gb = -d_gb;
            int d_br = white_image(2,0,r,c) - white_image(0,0,r,c);
            if(d_br<0) d_br = -d_br;

            output_image(r,c) = d_rg;
            if(output_image(r,c) < d_gb) output_image(r,c) = d_gb;
            if(output_image(r,c) < d_br) output_image(r,c) = d_br;
            
            //output_image(r,c) = d_rg + d_gb + d_br;
        }
    }
    
    return 1;
}



//int bdMelanoma::Extract2DSlice(bdImage &input, unsigned int t, unsigned int s, bdImage &output)
//{
//    if(input.IsEmpty()) return 0;
//    if(!(t<input.GetNumberOfTimeSeries())) return 0;
//    if(!(s<input.GetNumberOfSlices())) return 0;
//    
//    output.SetSize(1,1,input.GetNumberOfRows(),input.GetNumberOfColumns());
//    output.SetVisualizationPropertiesToMatchInput(input);
//    
//    for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
//    {
//        for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
//        {
//            output(r,c) = input(t,s,r,c);
//        }
//    }
//    output.Tag(0,0) = input.Tag(t,s);
//    
//    return 1;
//}



int bdMelanoma::PrintResults(const char *file_name)
{
    if(!file_name) return 0;
    
    bdString bds, bds_dir, bds_file_name, bds_ext;
    bds.Assign(file_name);
    bds.ExtractFileNamePathParts(bds_dir, bds_file_name, bds_ext);
    if(!(bds_ext=="txt" || bds_ext=="TXT"))
    {
        bds.Append(".txt");
    }
    
    ofstream file;
    file.open(bds.C_String(),ios::binary);

    
    cout<<endl;
    cout<<"----- ASYMMETRY -----"<<endl<<endl;
    cout<<"  axis 1 soft score = "<<m_asymmetry_soft_score_axis_1<<endl;
    cout<<"  axis 2 soft score = "<<m_asymmetry_soft_score_axis_2<<endl;
    cout<<"  axis 1 hard score = "<<m_asymmetry_hard_score_axis_1<<endl;
    cout<<"  axis 2 hard score = "<<m_asymmetry_hard_score_axis_2<<endl;
    cout<<endl;
    cout<<"  Asymmetry soft score = "<<m_asymmetry_soft_score<<endl;
    cout<<"  Asymmetry hard score = "<<m_asymmetry_hard_score<<endl;
    cout<<"---------------------"<<endl<<endl;
    
    
    file<<endl;
    file<<"----- ASYMMETRY -----"<<endl<<endl;
    file<<"  axis 1 soft score = "<<m_asymmetry_soft_score_axis_1<<endl;
    file<<"  axis 2 soft score = "<<m_asymmetry_soft_score_axis_2<<endl;
    file<<"  axis 1 hard score = "<<m_asymmetry_hard_score_axis_1<<endl;
    file<<"  axis 2 hard score = "<<m_asymmetry_hard_score_axis_2<<endl;
    file<<endl;
    file<<"  Asymmetry soft score = "<<m_asymmetry_soft_score<<endl;
    file<<"  Asymmetry hard score = "<<m_asymmetry_hard_score<<endl;
    file<<"---------------------"<<endl<<endl;
    
    
    cout<<endl;
    cout<<"----- BORDER -----"<<endl<<endl;
    cout<<"  Border soft score = "<<m_border_soft_score<<endl;
    cout<<"  Border hard score = "<<m_border_hard_score<<endl;
    cout<<"------------------"<<endl<<endl;
    
    
    file<<endl;
    file<<"----- BORDER -----"<<endl<<endl;
    file<<"  Border soft score = "<<m_border_soft_score<<endl;
    file<<"  Border hard score = "<<m_border_hard_score<<endl;
    file<<"------------------"<<endl<<endl;

    
    cout<<endl;
    cout<<"----- COLOR -----"<<endl<<endl;
    cout<<"  Brown color area = "<<m_color_brown_area<<endl;
    cout<<"  Red color area = "<<m_color_red_area<<endl;
    cout<<"  Blue-gray color area = "<<m_color_blue_gray_area<<endl;
    cout<<endl;
    cout<<"  Brown color volume = "<<m_color_brown_volume<<endl;
    cout<<"  Red color volume = "<<m_color_red_volume<<endl;
    cout<<"  Blue-gray color volume = "<<m_color_blue_gray_volume<<endl;
    cout<<endl;
    cout<<"  Mapped Dark Brown color volume = "<<m_color_mapped_dark_brown_volume<<endl;
    cout<<"  Mapped Light Brown color volume = "<<m_color_mapped_light_brown_volume<<endl;
    cout<<"  Mapped Red color volume = "<<m_color_mapped_red_volume<<endl;
    cout<<"  Mapped Blue-Gray color volume = "<<m_color_mapped_blue_gray_volume<<endl;
    cout<<"  Mapped Black color volume = "<<m_color_mapped_black_volume<<endl;
    cout<<"  Mapped White color volume = "<<m_color_mapped_white_volume<<endl;
    cout<<endl;
    cout<<"  Color soft score = "<<m_color_soft_score<<endl;
    cout<<"  Color hard score = "<<m_color_hard_score<<endl;
    cout<<"-----------------"<<endl<<endl;

    
    file<<endl;
    file<<"----- COLOR -----"<<endl<<endl;
    file<<"  Brown color area = "<<m_color_brown_area<<endl;
    file<<"  Red color area = "<<m_color_red_area<<endl;
    file<<"  Blue-gray color area = "<<m_color_blue_gray_area<<endl;
    file<<endl;
    file<<"  Brown color volume = "<<m_color_brown_volume<<endl;
    file<<"  Red color volume = "<<m_color_red_volume<<endl;
    file<<"  Blue-gray color volume = "<<m_color_blue_gray_volume<<endl;
    file<<endl;
    file<<"  Mapped Dark Brown color volume = "<<m_color_mapped_dark_brown_volume<<endl;
    file<<"  Mapped Light Brown color volume = "<<m_color_mapped_light_brown_volume<<endl;
    file<<"  Mapped Red color volume = "<<m_color_mapped_red_volume<<endl;
    file<<"  Mapped Blue-Gray color volume = "<<m_color_mapped_blue_gray_volume<<endl;
    file<<"  Mapped Black color volume = "<<m_color_mapped_black_volume<<endl;
    file<<"  Mapped White color volume = "<<m_color_mapped_white_volume<<endl;
    file<<endl;
    file<<"  Color soft score = "<<m_color_soft_score<<endl;
    file<<"  Color hard score = "<<m_color_hard_score<<endl;
    file<<"-----------------"<<endl<<endl;
  
    
    
    cout<<endl;
    cout<<"----- DIFFERENTIAL STRUCTURES -----"<<endl<<endl;
    cout<<"  Normalized area of dots = "<<m_DS_dots<<endl;
    cout<<"  Normalized area of globules = "<<m_DS_globules<<endl;
    cout<<"  Normalized area of streaks = "<<m_DS_streaks<<endl;
    cout<<"  Normalized area of structureless area = "<<m_DS_structureless_area<<endl;
    cout<<"  Normalized area of network = "<<m_DS_network<<endl;
    cout<<endl;
    cout<<"  D soft score = "<<m_DS_soft_score<<endl;
    cout<<"  D hard score = "<<m_DS_hard_score<<endl;
    cout<<"-----------------------------------"<<endl<<endl;

    
    
    file<<endl;
    file<<"----- DIFFERENTIAL STRUCTURES -----"<<endl<<endl;
    file<<"  Normalized area of dots = "<<m_DS_dots<<endl;
    file<<"  Normalized area of globules = "<<m_DS_globules<<endl;
    file<<"  Normalized area of streaks = "<<m_DS_streaks<<endl;
    file<<"  Normalized area of structureless area = "<<m_DS_structureless_area<<endl;
    file<<"  Normalized area of network = "<<m_DS_network<<endl;
    file<<endl;
    file<<"  D soft score = "<<m_DS_soft_score<<endl;
    file<<"  D hard score = "<<m_DS_hard_score<<endl;
    file<<"-----------------------------------"<<endl<<endl;

    


    
    
    
    file.close();
    
    return 1;
}





