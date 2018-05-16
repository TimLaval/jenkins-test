/**************************************************************************
 Class of operations defined on binary images

 Author: Danilo Babin
 File name: "bdBIP.cpp"
***************************************************************************/


#include "bdBIP.h"








int bdBIP::CenterOfObjectVolume(bdImage &image, unsigned int &output_center_s, unsigned int &output_center_r, unsigned int &output_center_c, unsigned int t)
{
	if(image.IsEmpty()) return 0;
    if(t>=image.GetNumberOfTimeSeries()) return 0;

	unsigned int n=0;
	for(unsigned int s=0; s<image.GetNumberOfSlices(); s++)
	{
		for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
			{
				if(image(t,s,r,c)!=0)
				{
					output_center_s += s;
                    output_center_r += r;
                    output_center_c += c;
                    n++;
				}
			}
		}
	}

	output_center_s = output_center_s/n;
	output_center_r = output_center_r/n;
	output_center_c = output_center_c/n;

	return 1;
}





int bdBIP::AreaOfObject_3D(bdImage &image, unsigned int t)
{
	if(t>=image.GetNumberOfTimeSeries()) return -1;
	int n = 0;
	for(unsigned int s=0; s<image.GetNumberOfSlices(); s++)
	{
		for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
			{
				if(image(t,s,r,c)!=0) n++;
			}
		}
	}
	return n;
}


int bdBIP::AreaOfObject_2D(bdImage &image, unsigned int sclice_index, unsigned int t)
{
	if(t>=image.GetNumberOfTimeSeries()) return -1;
	int n = 0;
	for(unsigned int r=0; r<image.GetNumberOfRows(); r++)
	{
		for(unsigned int c=0; c<image.GetNumberOfColumns(); c++)
		{
			if(image(t,sclice_index,r,c)!=0) n++;
		}
	}
	return n;
}


int bdBIP::Threshold(bdImage &input_image, int threshold_value, bdImage &output_image)
{
	if(input_image.IsEmpty()) return 0;
	output_image.SetSize(input_image.GetNumberOfTimeSeries(),input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
	output_image.SetVisualizationPropertiesToMatchInput(input_image);

	for(unsigned int t=0; t<input_image.GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
		{
			for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
			{
				for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
				{
					if(input_image(t,s,r,c)>threshold_value) output_image(t,s,r,c) = input_image(t,s,r,c);
					else output_image(t,s,r,c) = 0;
				}
			}
		}
	}
	return 1;
}



int bdBIP::NumberOfNonZeroVoxelsIn26Neighborhood(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	//if(image.IsEmpty()) return -1;
	if(!(t<image.GetNumberOfTimeSeries() && s<image.GetNumberOfSlices() && r<image.GetNumberOfRows() && c<image.GetNumberOfColumns())) return -1;

	bdGeometry g;
	g.SetDimensions(image.GetNumberOfSlices(), image.GetNumberOfRows(), image.GetNumberOfColumns());

	int number_of_neighbors = 0;
	int sn,rn,cn;
	for(g.ForCoordinates_26_Neighbors(s,r,c); g.Get_26_Neighbors(sn,rn,cn); )
	{
		//if(g.IsValidPosition(sn,rn,cn)) { if(image(sn,rn,cn)!=0) number_of_neighbors++; }
		if(image(t,sn,rn,cn)!=0) number_of_neighbors++; 
	}

	return number_of_neighbors;
}


int bdBIP::NumberOfNonZeroVoxelsIn26Neighborhood(bdImage &image, unsigned int s, unsigned int r, unsigned int c)
{
	return this->NumberOfNonZeroVoxelsIn26Neighborhood(image,0,s,r,c);
}



int bdBIP::NumberOfNonZeroVoxelsIn27Neighborhood(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
	//if(image.IsEmpty()) return -1;
	if(!(t<image.GetNumberOfTimeSeries() && s<image.GetNumberOfSlices() && r<image.GetNumberOfRows() && c<image.GetNumberOfColumns())) return -1;

	bdGeometry g;
	g.SetDimensions(image.GetNumberOfSlices(), image.GetNumberOfRows(), image.GetNumberOfColumns());

	int number_of_neighbors = 0;
	int sn,rn,cn;
	for(g.ForCoordinates_27_Neighbors(s,r,c); g.Get_27_Neighbors(sn,rn,cn); )
	{ 
		//if(g.IsValidPosition(sn,rn,cn)) { if(image(sn,rn,cn)!=0) number_of_neighbors++; }
		if(image(t,sn,rn,cn)!=0) number_of_neighbors++;
	}
	return number_of_neighbors;
}


int bdBIP::NumberOfNonZeroVoxelsIn27Neighborhood(bdImage &image, unsigned int s, unsigned int r, unsigned int c)
{
	return this->NumberOfNonZeroVoxelsIn27Neighborhood(image,0,s,r,c);
}


int bdBIP::GetClosestNonZeroVoxelToSeedIndexes(bdImage &image, unsigned int t, unsigned int s, unsigned int r, unsigned int c, unsigned int &out_s, unsigned int &out_r, unsigned int &out_c)
{
    if(image.IsEmpty()) return 0;
    
    if(image(t,s,r,c)!=0)
    {
        out_s = s; out_r = r; out_c = c;
        return 1;
    }
    
    bdGeometry g;
    g.SetDimensions(image.GetNumberOfSlices(), image.GetNumberOfRows(), image.GetNumberOfColumns());
    
    for(int squared_radius = 0; squared_radius<g.GetNumberOfSphereElements(); squared_radius++)
    {
        int sn, rn, cn;
        for(g.ForCoordinates_Sphere(s,r,c, squared_radius); g.Get_Sphere(squared_radius+1, sn, rn, cn); )
        {
            if(image(t,sn,rn,cn)!=0)
            {
                out_s = sn; out_r = rn; out_c = cn;
                return 1;
            }
        }
    }
    // if it came to here, no non-zero voxel was found.
    return 0;
}

int bdBIP::GetClosestNonZeroVoxelToSeedIndexes(bdImage &image, unsigned int s, unsigned int r, unsigned int c, unsigned int &out_s, unsigned int &out_r, unsigned int &out_c)
{
    return this->GetClosestNonZeroVoxelToSeedIndexes(image, 0, s, r, c, out_s, out_r, out_c);
}




int bdBIP::ExtractLargest_26_ConnectedComponent(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t)
{
	if(input_image.IsEmpty()) return 0;

	//Initialize output image with the current binary image
	output_image.CopyFrom(input_image);
	
	bdVoxel p;
	int sn, rn, cn;
	bdGeometry g;
	g.SetDimensions(input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
	bdList<bdVoxel> l, l1, l2, *p_l_previous, *p_l_new, *p_l_temp;
	p_l_previous = &l1;
	p_l_new = &l2;
	for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
	{
		for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
		{
			for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
			{
				//Find the whole 3D region
				if(output_image(t,s,r,c)>threshold)
				{
					l.AddToRightEnd(p(s,r,c));
					p_l_new->AddToRightEnd(p);
					output_image(t,s,r,c) = 0;//We delete the pixel from the image after placing it in the list

					while(!l.IsEmpty())
					{
						p = l.GetLeftEnd();
						for(g.ForCoordinates_26_Neighbors(p.S(),p.R(),p.C()); g.Get_26_Neighbors(sn,rn,cn); )
						{
							if(output_image(t,sn,rn,cn)>threshold)
							{
								l.AddToRightEnd(p(sn,rn,cn));
								p_l_new->AddToRightEnd(p);
								output_image(t,sn,rn,cn) = 0;//Delete the pixel from output_image
							}
						}
						l.DeleteLeftEnd();
					}
				}
				else output_image(t,s,r,c) = 0;
				//After finding the whole 3D region, check if it is larger than the previous one...
				if(p_l_new->GetNumberOfElements()>p_l_previous->GetNumberOfElements())
				{
					//...if so, replace the previous region with the new one
					p_l_previous->Reset();//DeleteAll();
					p_l_temp = p_l_previous;
					p_l_previous = p_l_new;
					p_l_new = p_l_temp;
				}
				else
				{
					//...if not, just delete the new region
					p_l_new->Reset();//DeleteAll();
				}
			}
		}
	}
	//After this part the 'p_l_previous' points to the list of the largest region and output_image is empty,
	// so we need to copy the list to the image

	while(!p_l_previous->IsEmpty())
	{
		p = p_l_previous->GetLeftEnd();
		output_image(t,p.S(),p.R(),p.C()) = input_image(t,p.S(),p.R(),p.C());//(*this)[p.s][p.r][p.c];
		p_l_previous->DeleteLeftEnd();
	}
	return 1;
}


//int bdBIP::Extract_CC26_LargestAndClosestToSeed(bdImage &input_image, int threshold, bdImage &output_image, unsigned int seed_s, unsigned int seed_r, unsigned int seed_c, double max_distance, unsigned int t)
//{
///// NOT FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!
//    
//    if(input_image.IsEmpty()) return 0;
//    
//    //Initialize output image with the current binary image
//    output_image.CopyFrom(input_image);
//    
//    double max_score = 0;//maximum score of distance/max_distance*area over all CCs up to now.
//    bdVoxel p;
//    int sn, rn, cn;
//    bdGeometry g;
//    g.SetDimensions(input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
//    bdList<bdVoxel> l, l1, l2, *p_l_previous, *p_l_new, *p_l_temp;
//    p_l_previous = &l1;
//    p_l_new = &l2;
//    for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
//    {
//        for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
//        {
//            for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
//            {
//                //Find the whole 3D region
//                if(output_image(t,s,r,c)>threshold)
//                {
//                    l.AddToRightEnd(p(s,r,c));
//                    p_l_new->AddToRightEnd(p);
//                    output_image(t,s,r,c) = 0;//We delete the pixel from the image after placing it in the list
//                    
//                    while(!l.IsEmpty())
//                    {
//                        p = l.GetLeftEnd();
//                        for(g.ForCoordinates_26_Neighbors(p.S(),p.R(),p.C()); g.Get_26_Neighbors(sn,rn,cn); )
//                        {
//                            if(output_image(t,sn,rn,cn)>threshold)
//                            {
//                                l.AddToRightEnd(p(sn,rn,cn));
//                                p_l_new->AddToRightEnd(p);
//                                output_image(t,sn,rn,cn) = 0;//Delete the pixel from output_image
//                            }
//                        }
//                        l.DeleteLeftEnd();
//                    }
//                }
//                else output_image(t,s,r,c) = 0;
//                //After finding the whole 3D region, check if it is larger than the previous one...
//                if(p_l_new->GetNumberOfElements()>p_l_previous->GetNumberOfElements())
//                {
//                    //...if so, replace the previous region with the new one
//                    p_l_previous->Reset();//DeleteAll();
//                    p_l_temp = p_l_previous;
//                    p_l_previous = p_l_new;
//                    p_l_new = p_l_temp;
//                }
//                else
//                {
//                    //...if not, just delete the new region
//                    p_l_new->Reset();//DeleteAll();
//                }
//            }
//        }
//    }
//    //After this part the 'p_l_previous' points to the list of the largest region and output_image is empty,
//    // so we need to copy the list to the image
//    
//    while(!p_l_previous->IsEmpty())
//    {
//        p = p_l_previous->GetLeftEnd();
//        output_image(t,p.S(),p.R(),p.C()) = input_image(t,p.S(),p.R(),p.C());//(*this)[p.s][p.r][p.c];
//        p_l_previous->DeleteLeftEnd();
//    }
//    return 1;
//
//}


int bdBIP::ConnectedComponents_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t, int *output_number_of_components)
{
    if(input_image.IsEmpty()) return 0;
    
    bdImage thresholded_image;
    this->Threshold(input_image, threshold, thresholded_image);

    
    //Initialize output image.
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
    
    *output_number_of_components = n_of_regions-1;
    
    return 1;
    
//	if(input_image.IsEmpty()) return 0;
//
//	//Initialize output image.
//	output_image.SetSize(input_image.GetNumberOfTimeSeries(),input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());
//	output_image.FillInWith(0);
//
//	unsigned short n_of_regions = 1;
//	for(unsigned int s=0; s<input_image.GetNumberOfSlices(); s++)
//	{
//		for(unsigned int r=0; r<input_image.GetNumberOfRows(); r++)
//		{
//			for(unsigned int c=0; c<input_image.GetNumberOfColumns(); c++)
//			{
//				//Find the 3D region
//				if(input_image(t,s,r,c)>threshold && output_image(t,s,r,c)==0)
//				{
//					bdRegion3D reg;
//					reg.CreateRegion_26_FromSeed(input_image,c,r,s,t);
//
//					bdRegion3DIterator itr;
//					for(itr.SetBegin(&reg); itr.IsValid(); itr.MoveToNext())
//					{
//						output_image(itr.GetIndexSlice(),itr.GetIndexRow(),itr.GetIndexColumn()) = n_of_regions;
//					}
//					n_of_regions++;
//                    cout<<" n:"<<n_of_regions<<" ";
//				}
//			}
//		}
//	}
//    
//    *output_number_of_components = n_of_regions-1;
//    
//	return 1;
}


int bdBIP::Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, bdList<bdDiscreteCoordinates3D> &seeds, unsigned int t)
{
    if(input_image.IsEmpty()) return 0;
    
    //Initialize output image with the current binary image
    output_image.CopyFrom(input_image);
    output_image.FillInWith(0);
    
    
    bdList<bdDiscreteCoordinates3D> l;
    
    bdListIterator<bdDiscreteCoordinates3D> it;
    for(it.SetLeftEnd(seeds); it.IsValid(); it.MoveRight())
    {
        if(input_image(t,it.GetElementPointer()->S(),it.GetElementPointer()->R(),it.GetElementPointer()->C()) > threshold )
        {
            cout<<" seed: "<<it.GetElementPointer()->S()<<","<<it.GetElementPointer()->R()<<","<<it.GetElementPointer()->C()<<" ";
            output_image(t,it.GetElementPointer()->S(),it.GetElementPointer()->R(),it.GetElementPointer()->C()) = input_image(t,it.GetElementPointer()->S(),it.GetElementPointer()->R(),it.GetElementPointer()->C());
            l.AddToLeftEnd(it.GetElement());
        }
    }
    
    bdGeometry g;
    g.SetDimensions(input_image.GetNumberOfSlices(),input_image.GetNumberOfRows(),input_image.GetNumberOfColumns());

    while(!l.IsEmpty())
    {
        bdDiscreteCoordinates3D p = l.GetLeftEnd();
        int sn, rn, cn;
        for(g.ForCoordinates_26_Neighbors(p.S(),p.R(),p.C()); g.Get_26_Neighbors(sn,rn,cn); )
        {
            if(input_image(t,sn,rn,cn)>threshold && output_image(t,sn,rn,cn) == 0)
            {
                l.AddToRightEnd(p(sn,rn,cn));
                output_image(t,sn,rn,cn) = input_image(t,sn,rn,cn);
            }
        }
        l.DeleteLeftEnd();
    }
    
    return 1;
}


int bdBIP::Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int t, unsigned int s, unsigned int r, unsigned int c)
{
    bdList<bdDiscreteCoordinates3D> seed;
    bdDiscreteCoordinates3D *dc = seed.AddNewToRightEnd();
    dc->S() = s; dc->R() = r; dc->C() = c;
    return this->Extract_26_ConnectedComponent_3D(input_image, threshold, output_image, seed, t);
}


int bdBIP::Extract_26_ConnectedComponent_3D(bdImage &input_image, int threshold, bdImage &output_image, unsigned int s, unsigned int r, unsigned int c)
{
    return this->Extract_26_ConnectedComponent_3D(input_image, threshold, output_image, 0, s, r, c);
}





int bdBIP::RadiusSquared_Sphere(bdImage &input_skeleton_image, bdImage &input_segmented_image, bdImage &output_image)
{
	if(input_skeleton_image.IsEmpty()) return 0;
	if(!input_segmented_image.IsEqualSizeAs(input_skeleton_image)) return 0;

	output_image.SetSize(input_skeleton_image.GetNumberOfTimeSeries(),input_skeleton_image.GetNumberOfSlices(),input_skeleton_image.GetNumberOfRows(),input_skeleton_image.GetNumberOfColumns());
	output_image.SetVisualizationPropertiesToMatchInput(input_segmented_image);
	output_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(input_skeleton_image.GetNumberOfSlices(),input_skeleton_image.GetNumberOfRows(),input_skeleton_image.GetNumberOfColumns());

	for(unsigned int t=0; t<output_image.GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<output_image.GetNumberOfSlices(); s++)
		{
			for(unsigned int r=0; r<output_image.GetNumberOfRows(); r++)
			{
				for(unsigned int c=0; c<output_image.GetNumberOfColumns(); c++)
				{
					if(input_skeleton_image(t,s,r,c)!=0)
					{
						int is_radius_found = 0;
						unsigned int squared_radius = 0;
						while(!is_radius_found)
						{
							squared_radius++;
							int sn,rn,cn;
							for(g.ForCoordinates_Sphere(s,r,c,squared_radius); g.Get_Sphere(squared_radius,sn,rn,cn); )
							{
								if(input_segmented_image(t,sn,rn,cn)==0)
								{
									is_radius_found = 1;
									break;
								}
							}
						
							if(squared_radius>=g.GetNumberOfSphereElements()) is_radius_found = 1;
						}
						output_image(t,s,r,c) = squared_radius;
					}
				}
			}
		}
	}
	return 1;
}


int bdBIP::RadiusSquared_CirclePer2DSlice(bdImage &input_mask_image, bdImage &input_segmented_image, bdImage &output_image)
{
	if(input_mask_image.IsEmpty()) return 0;
	if(!input_segmented_image.IsEqualSizeAs(input_mask_image)) return 0;

	output_image.SetSize(input_mask_image.GetNumberOfTimeSeries(),input_mask_image.GetNumberOfSlices(),input_mask_image.GetNumberOfRows(),input_mask_image.GetNumberOfColumns());
	output_image.SetVisualizationPropertiesToMatchInput(input_segmented_image);
	output_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(1,input_mask_image.GetNumberOfRows(),input_mask_image.GetNumberOfColumns());

	for(unsigned int t=0; t<output_image.GetNumberOfTimeSeries(); t++)
	{
		for(unsigned int s=0; s<output_image.GetNumberOfSlices(); s++)
		{
			for(unsigned int r=0; r<output_image.GetNumberOfRows(); r++)
			{
				for(unsigned int c=0; c<output_image.GetNumberOfColumns(); c++)
				{
					if(input_mask_image(t,s,r,c)!=0)
					{
						int is_radius_found = 0;
						unsigned int squared_radius = 0;
						while(!is_radius_found)
						{
							squared_radius++;
							int rn,cn;
							for(g.ForCoordinates_Circle(r,c,squared_radius); g.Get_Circle(squared_radius,rn,cn); )
							{
								if(input_segmented_image(t,s,rn,cn)==0)
								{
									is_radius_found = 1;
									break;
								}
							}
						
							if(squared_radius>=g.GetNumberOfSphereElements()) is_radius_found = 1;
						}
						output_image(t,s,r,c) = squared_radius;
					}
				}
			}
		}
	}
	return 1;
}




int bdBIP::ProfileVolume_Sphere(bdImage &input_mask_image, bdImage &input_segmented_image, bdImage &output_image, unsigned int offset_percent)
{
	if(offset_percent<1 || offset_percent>100) offset_percent = 100;
	if(input_mask_image.IsEmpty()) return 0;
	if(!input_segmented_image.IsEqualSizeAs(input_mask_image)) return 0;

	output_image.SetSizeAndPropertiesAs(input_mask_image);
	output_image.FillInWith(0);

	bdGeometry g;
	g.SetDimensions(input_mask_image.GetNumberOfSlices(),input_mask_image.GetNumberOfRows(),input_mask_image.GetNumberOfColumns());

    for(unsigned int t=0; t<output_image.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output_image.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output_image.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output_image.GetNumberOfColumns(); c++)
                {
                    if(input_mask_image(t,s,r,c)!=0)
                    {
                        int is_radius_found = 0;
                        int number_of_non_zero_voxels = 0;

                        int sn,rn,cn;
                        unsigned int previous_index = 0;
                        for(g.ForCoordinates_Sphere(s,r,c,0); g.Get_Sphere(g.GetNumberOfSphereElements(),sn,rn,cn); )
                        {
                            //Check if the condition for exiting is fulfulled.
                            if(is_radius_found)
                            {
                                if(g.GetCurrentFORLoopIndex() != previous_index) { break; }
                            }
                            
                            if(input_segmented_image(t,sn,rn,cn)==0) { is_radius_found = 1; }
                            else number_of_non_zero_voxels++;

                            //Update the previous index value.
                            previous_index = g.GetCurrentFORLoopIndex();
                        }
                        output_image(t,s,r,c) = number_of_non_zero_voxels;
                    }
                }
            }
        }
    }
    
	return 1;
}









int bdBIP::Dilation_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	if(input.IsEmpty()) return 0;

	output.SetSizeAndPropertiesAs(input);
	output.FillInWith(0);
	
	bdGeometry g;
	g.SetDimensions(output.GetNumberOfSlices(), output.GetNumberOfRows(), output.GetNumberOfColumns());

    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    if(input(t,s,r,c)!=0)
                    {
                        int sn,rn,cn;
                        for(g.ForCoordinates_Sphere(s,r,c,0); g.Get_Sphere(sphere_squared_radius,sn,rn,cn); )
                        {
                            output(t,sn,rn,cn) = input(t,s,r,c);
                        }
                    }
                }
            }
        }
    }
	return 1;
}



int bdBIP::DilationWithImageResize_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	if(input.IsEmpty()) return 0;
	int sphere_radius = (int)( sqrt((double)sphere_squared_radius));

	//Resize by adding 2*(sphere_radius+2) to each dimension
	bdImage resized_image;
	resized_image.SetSize(input.GetNumberOfTimeSeries(),input.GetNumberOfSlices()+2*(sphere_radius+2),input.GetNumberOfRows()+2*(sphere_radius+2),input.GetNumberOfColumns()+2*(sphere_radius+2));
	resized_image.FillInWith(0);

    for(unsigned int t=0; t<input.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<input.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<input.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<input.GetNumberOfColumns(); c++)
                {
                    resized_image(t,sphere_radius+2+s,sphere_radius+2+r,sphere_radius+2+c) = input(t,s,r,c);
                }
            }
        }
    }

	//Now just dilate the resized image
	return this->Dilation_Sphere(resized_image,output,sphere_squared_radius);
}



int bdBIP::Erosion_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	if(input.IsEmpty()) return 0;

	output.CopyFrom(input);

	bdGeometry g;
	g.SetDimensions(output.GetNumberOfSlices(), output.GetNumberOfRows(), output.GetNumberOfColumns());

    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    if(input(t,s,r,c)!=0)
                    {
                        int sn,rn,cn;
                        for(g.ForCoordinates_Sphere(s,r,c,0); g.Get_Sphere(sphere_squared_radius,sn,rn,cn); )
                        {
                            if(input(t,sn,rn,cn)==0)
                            {
                                output(t,s,r,c) = 0;
                                break;
                            }
                        }
                    }
                    //else output(s,r,c) = 0;
                }
            }
        }
    }
	return 1;
}



int bdBIP::Erosion_Circle(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius)
{
    if(input.IsEmpty()) return 0;
    
    output.CopyFrom(input);
    
    bdGeometry g;
    g.SetDimensions(output.GetNumberOfSlices(), output.GetNumberOfRows(), output.GetNumberOfColumns());
    
    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0)
                    {
                        if(input(t,s,r,c)!=0)
                        {
                            int rn,cn;
                            for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn); )
                            {
                                if(input(t,s,rn,cn)==0)
                                {
                                    output(t,s,r,c) = 0;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}



int bdBIP::ErosionByMapping_Circle(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius)
{
    if(input.IsEmpty()) return 0;
    
    output.SetSizeAndPropertiesAs(input);
    output.FillInWith(0);
    
    double radius = sqrt ((double)squared_radius);
    
    bdGeometry g;
    g.SetDimensions(output.GetNumberOfSlices(), output.GetNumberOfRows(), output.GetNumberOfColumns());
    
    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0)
                    {
                        if(input(t,s,r,c)!=0)
                        {
                            int is_zero_voxel_found = 0;
                            int rn,cn;
                            for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn); )
                            {
                                if(input(t,s,rn,cn)==0)
                                {
                                    is_zero_voxel_found = 1;
                                    break;
                                }
                            }
                            if(!is_zero_voxel_found)
                            {
                                unsigned int squared_radius_larger =  (unsigned int)((radius+1)*(radius+1));
                                //for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn); )
                                for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius_larger,rn,cn); )
                                {
                                    output(t,s,rn,cn) = input(t,s,rn,cn);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}



int bdBIP::ErosionByMapping_Circle2(bdImage &input, bdImage &mask, bdImage &output, unsigned int squared_radius, unsigned int percentage)
{
    if(input.IsEmpty()) return 0;
    
    output.SetSizeAndPropertiesAs(input);
    output.FillInWith(0);
    
    //double radius = sqrt ((double)squared_radius);
    
    bdGeometry g;
    g.SetDimensions(output.GetNumberOfSlices(), output.GetNumberOfRows(), output.GetNumberOfColumns());
    
    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    if(mask(r,c)!=0)
                    {
                        if(input(t,s,r,c)!=0)
                        {
                            int n_of_non_zero_pixels = 0;
                            int n_of_pixels = 0;
                            
                            int rn,cn;
                            for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn); )
                            {
                                if(input(t,s,rn,cn)!=0)
                                {
                                    n_of_non_zero_pixels++;
                                }
                                n_of_pixels++;
                            }
                            
                            if(n_of_non_zero_pixels >= ((n_of_pixels*percentage)/100))
                            {
                                int rn,cn;
                                for(g.ForCoordinates_Circle(r,c,0); g.Get_Circle(squared_radius,rn,cn); )
                                {
                                    if(input(t,s,rn,cn)!=0)
                                    {
                                        output(t,s,rn,cn) = input(t,s,rn,cn);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}



int bdBIP::Opening_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	bdImage temp;
	if(!this->Erosion_Sphere(input,temp,sphere_squared_radius)) return 0;
	return this->Dilation_Sphere(temp,output,sphere_squared_radius);
}



int bdBIP::Closing_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	bdImage temp;
	if(!this->Dilation_Sphere(input,temp,sphere_squared_radius)) return 0;
	return this->Erosion_Sphere(temp,output,sphere_squared_radius);
}



int bdBIP::ClosingAccurate_Sphere(bdImage &input, bdImage &output, unsigned int sphere_squared_radius)
{
	bdImage temp1,temp2;
	if(!this->DilationWithImageResize_Sphere(input,temp1,sphere_squared_radius)) return 0;
	if(!this->Erosion_Sphere(temp1,temp2,sphere_squared_radius)) return 0;

	//Now just resize the image to its original value
	int sphere_radius = (int)( sqrt((double)sphere_squared_radius));
	output.SetSizeAndPropertiesAs(input);
    for(unsigned int t=0; t<output.GetNumberOfTimeSeries(); t++)
    {
        for(unsigned int s=0; s<output.GetNumberOfSlices(); s++)
        {
            for(unsigned int r=0; r<output.GetNumberOfRows(); r++)
            {
                for(unsigned int c=0; c<output.GetNumberOfColumns(); c++)
                {
                    output(t,s,r,c) = temp2(t,sphere_radius+2+s,sphere_radius+2+r,sphere_radius+2+c);
                }
            }
        }
    }
	return 1;
}



