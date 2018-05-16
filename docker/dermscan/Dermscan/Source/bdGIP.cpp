/**************************************************************************
 Template class for processing gray images.

 Author: Danilo Babin
 File name: "bdGIP.h"
***************************************************************************/


#include "bdGIP.h"




int bdGIP::MaximumCircleFor2DSlice(bdImage &input, bdImage &mask, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask, unsigned int circle_squared_radius, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(),input.GetNumberOfRows(),input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    output.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(1,input.GetNumberOfRows(),input.GetNumberOfColumns());
    
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            if(mask(t_mask,s_mask,r,c)!=0)
            {
                unsigned short max = 0;
                int rn,cn;
                for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(circle_squared_radius,rn,cn); )
                {
                    if(max<input(t_input,s_input,rn,cn)) max = input(t_input,s_input,rn,cn);
                }
                
                output(t_input,s_input,r,c) = max;
            }
        }
    }

    return 1;
}


int bdGIP::MaskPer_2D_Slice(bdImage &input, bdImage &mask, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    if(input.GetNumberOfRows()!=mask.GetNumberOfRows() || input.GetNumberOfColumns()!=mask.GetNumberOfColumns()) return 0;
    
    output.SetSize(input.GetNumberOfTimeSeries(), input.GetNumberOfSlices(), input.GetNumberOfRows(), input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    
    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<input.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0) output(t,s,r,c) = input(t,s,r,c);
                    else output(t,s,r,c) = 0;
                }
            }
        }
    }
    
    return 1;
}





int bdGIP::MeanCircleFor2DSlice(bdImage &input, bdImage &mask, unsigned int t_input, unsigned int s_input, unsigned int t_mask, unsigned int s_mask, unsigned int circle_squared_radius, bdImage &output)
{
    if(input.IsEmpty()) return 0;
    output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(),input.GetNumberOfRows(),input.GetNumberOfColumns());
    output.SetVisualizationPropertiesToMatchInput(input);
    output.FillInWith(0);
    
    bdGeometry g;
    g.SetDimensions(1,input.GetNumberOfRows(),input.GetNumberOfColumns());
    
    for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
    {
        for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
        {
            if(mask(t_mask,s_mask,r,c)!=0)
            {
                unsigned int sum = 0;
                unsigned int n_of_voxels = 0;
                int rn,cn;
                for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(circle_squared_radius,rn,cn); )
                {
                    sum += input(t_input,s_input,rn,cn);
                    n_of_voxels++;
                }
                
                if(n_of_voxels) output(t_input,s_input,r,c) = sum / n_of_voxels;
                else output(t_input,s_input,r,c) = 0;
            }
        }
    }
    
    return 1;
}




int bdGIP::Negative(bdImage &input, bdImage &output)
{
	if(input.IsEmpty()) return 0;
	output.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices(),input.GetNumberOfRows(),input.GetNumberOfColumns());
	output.SetVisualizationPropertiesToMatchInput(input);

	//Find the min and max value
	unsigned int min, max;
	input.GetVoxelMinimumAndMaximumValue(&min,&max);
	int range = max-min;

	for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<input.GetNumberOfSlices(); s++)
		{
			for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
			{
				if(this->IsAbortRequested()) return 0;
				//int progress = ((t*s*input.GetNumberOfRows() + r)*100) / (input.GetNumberOfTimeSeries()*input.GetNumberOfSlices()*input.GetNumberOfRows());
				//this->SetProgressCounterRelativeValue(progress);

				for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
				{
					double coeff = 1.0 - ((double)(input(t,s,r,c)-min)) / ((double)range);
					output(t,s,r,c) = (unsigned short)(coeff * ((double)max));
				}
			}
		}
	}
	return 1;
}






int bdGIP::RescaleWholeRange(bdImage &input_image, unsigned int lower_value, unsigned int upper_value, bdImage &output_image)
{
	if(input_image.IsEmpty()) return 0;
	unsigned int full_range_min, full_range_max;
	
	input_image.GetVoxelValueFullRange(&full_range_min, &full_range_max);
	if(upper_value>full_range_max) upper_value = full_range_max;
	if(lower_value<full_range_min) lower_value = full_range_min;

	output_image.SetSize(input_image.GetNumberOfTimeSeries(),input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
	output_image.SetVisualizationPropertiesToMatchInput(input_image);

	//Find min and max value
	unsigned int min = 0, max = 1;
	input_image.GetVoxelMinimumAndMaximumValue(&min, &max);
	unsigned int range = max-min;

	for(unsigned int t=0; t<input_image.GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
		{
			for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
			{
				for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
				{
					output_image(t,s,r,c) = input_image(t,s,r,c) - min;
					output_image(t,s,r,c) = (unsigned short)( lower_value + (((double)(output_image(t,s,r,c))) * (upper_value-lower_value)) / ((double)range) );				
				}
			}
		}
	}

	return 1;
}




