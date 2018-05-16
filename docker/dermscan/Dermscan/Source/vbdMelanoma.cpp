/**************************************************************************
 Methods loading/saving and visual manipulation of images of melanoma.
 
 Author: Danilo Babin
 File name: "vbdMelanoma.h"
 ***************************************************************************/



#include "vbdMelanoma.h"




int vbdMelanoma::LoadGrayScaleImageFromPNGFile(const char *file_name, bdImage &output_loaded_image)
{
    // Read the image
    vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
    reader->SetFileName(file_name);
    reader->Update();
    
    // Extract R, G, B components
    vtkSmartPointer<vtkImageExtractComponents> extract_red_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract_red_filter->SetInputConnection(reader->GetOutputPort());
    extract_red_filter->SetComponents(0);
    extract_red_filter->Update();
    
//    vtkSmartPointer<vtkImageExtractComponents> extract_green_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
//    extract_green_filter->SetInputConnection(reader->GetOutputPort());
//    extract_green_filter->SetComponents(1);
//    extract_green_filter->Update();
//    
//    vtkSmartPointer<vtkImageExtractComponents> extract_blue_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
//    extract_blue_filter->SetInputConnection(reader->GetOutputPort());
//    extract_blue_filter->SetComponents(2);
//    extract_blue_filter->Update();
    
    vtkSmartPointer<vtkImageData> red_image = vtkSmartPointer<vtkImageData>::New();
    red_image->DeepCopy(extract_red_filter->GetOutput());
    
//    vtkSmartPointer<vtkImageData> green_image = vtkSmartPointer<vtkImageData>::New();
//    green_image->DeepCopy(extract_green_filter->GetOutput());
//    
//    vtkSmartPointer<vtkImageData> blue_image = vtkSmartPointer<vtkImageData>::New();
//    blue_image->DeepCopy(extract_blue_filter->GetOutput());
    
    // Set image dimensions.
    output_loaded_image.SetSize(1, 1, red_image->GetDimensions()[1], red_image->GetDimensions()[0]);
    
    // Depending on the scalar size we will have to address the pixels in a different way (unsigned char or unsigned short)
    // vtkImageData::GetScalarSize returns the number of Bytes of the scalar.
    
    // If image has unsigned char elements...
    if(red_image->GetScalarSize()==1)
    {
        /// ...iterate over all pixels and copy to output image.
        for(unsigned int r=0; r<output_loaded_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output_loaded_image.GetNumberOfColumns(); c++)
            {
                unsigned char *pixel_r = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
//                unsigned char *pixel_g = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
//                unsigned char *pixel_b = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));
//                
                output_loaded_image(0,0,r,c) = pixel_r[0];
//                output_loaded_image(1,0,r,c) = pixel_g[0];
//                output_loaded_image(2,0,r,c) = pixel_b[0];
            }
        }
    }
    
    // If image has unsigned short elements...
    if(red_image->GetScalarSize()==2)
    {
        /// ...iterate over all pixels and copy to output image.
        for(unsigned int r=0; r<output_loaded_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output_loaded_image.GetNumberOfColumns(); c++)
            {
                unsigned short *pixel_r = static_cast<unsigned short*>(red_image->GetScalarPointer(c,r,0));
//                unsigned short *pixel_g = static_cast<unsigned short*>(green_image->GetScalarPointer(c,r,0));
//                unsigned short *pixel_b = static_cast<unsigned short*>(blue_image->GetScalarPointer(c,r,0));
//                
                output_loaded_image(0,0,r,c) = pixel_r[0];
//                output_loaded_image(1,0,r,c) = pixel_g[0];
//                output_loaded_image(2,0,r,c) = pixel_b[0];
            }
        }
    }
    
    return 1;


}


int vbdMelanoma::LoadWhiteLightRGBImageFromPNGFile(const char *file_name, bdImage &output_loaded_image)
{
    // Read the image
    vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
    reader->SetFileName(file_name);
    reader->Update();
    
    // Extract R, G, B components
    vtkSmartPointer<vtkImageExtractComponents> extract_red_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract_red_filter->SetInputConnection(reader->GetOutputPort());
    extract_red_filter->SetComponents(0);
    extract_red_filter->Update();
    
    vtkSmartPointer<vtkImageExtractComponents> extract_green_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract_green_filter->SetInputConnection(reader->GetOutputPort());
    extract_green_filter->SetComponents(1);
    extract_green_filter->Update();
    
    vtkSmartPointer<vtkImageExtractComponents> extract_blue_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extract_blue_filter->SetInputConnection(reader->GetOutputPort());
    extract_blue_filter->SetComponents(2);
    extract_blue_filter->Update();
    
    vtkSmartPointer<vtkImageData> red_image = vtkSmartPointer<vtkImageData>::New();
    red_image->DeepCopy(extract_red_filter->GetOutput());
    
    vtkSmartPointer<vtkImageData> green_image = vtkSmartPointer<vtkImageData>::New();
    green_image->DeepCopy(extract_green_filter->GetOutput());
    
    vtkSmartPointer<vtkImageData> blue_image = vtkSmartPointer<vtkImageData>::New();
    blue_image->DeepCopy(extract_blue_filter->GetOutput());
    
    // Set image dimensions.
    output_loaded_image.SetSize(3, 1, red_image->GetDimensions()[1], red_image->GetDimensions()[0]);
    
    // Depending on the scalar size we will have to address the pixels in a different way (unsigned char or unsigned short)
    // vtkImageData::GetScalarSize returns the number of Bytes of the scalar.
    
    // If image has unsigned char elements...
    if(red_image->GetScalarSize()==1)
    {
        /// ...iterate over all pixels and copy to output image.
        for(unsigned int r=0; r<output_loaded_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output_loaded_image.GetNumberOfColumns(); c++)
            {
                unsigned char *pixel_r = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
                unsigned char *pixel_g = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
                unsigned char *pixel_b = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));

                output_loaded_image(0,0,r,c) = pixel_r[0];
                output_loaded_image(1,0,r,c) = pixel_g[0];
                output_loaded_image(2,0,r,c) = pixel_b[0];
            }
        }
    }
    
    // If image has unsigned short elements...
    if(red_image->GetScalarSize()==2)
    {
        /// ...iterate over all pixels and copy to output image.
        for(unsigned int r=0; r<output_loaded_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<output_loaded_image.GetNumberOfColumns(); c++)
            {
                unsigned short *pixel_r = static_cast<unsigned short*>(red_image->GetScalarPointer(c,r,0));
                unsigned short *pixel_g = static_cast<unsigned short*>(green_image->GetScalarPointer(c,r,0));
                unsigned short *pixel_b = static_cast<unsigned short*>(blue_image->GetScalarPointer(c,r,0));
                
                output_loaded_image(0,0,r,c) = pixel_r[0];
                output_loaded_image(1,0,r,c) = pixel_g[0];
                output_loaded_image(2,0,r,c) = pixel_b[0];
            }
        }
    }

    return 1;
}



int vbdMelanoma::SaveWhiteLightRGBImageToPNGFile_8bit(bdImage &decomposed_RGB_image, const char *output_file_name)
{
    bdGIP gip;
    bdImage decomposed_RGB_image_rescaled;
    gip.RescaleWholeRange(decomposed_RGB_image, 0, 255, decomposed_RGB_image_rescaled);
    
    vtkSmartPointer<vtkImageData> red_image = vtkSmartPointer<vtkImageData>::New();
    red_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    red_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//red_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);

    vtkSmartPointer<vtkImageData> green_image = vtkSmartPointer<vtkImageData>::New();
    green_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    green_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//green_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);

    vtkSmartPointer<vtkImageData> blue_image = vtkSmartPointer<vtkImageData>::New();
    blue_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    blue_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//blue_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);

    // If input image is really a color image...
    if(decomposed_RGB_image.GetNumberOfTimeSeries()>2)
    {
        for(unsigned int r=0; r<decomposed_RGB_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<decomposed_RGB_image.GetNumberOfColumns(); c++)
            {
//                unsigned short *pixel;
//                pixel = static_cast<unsigned short*>(red_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
//
//                pixel = static_cast<unsigned short*>(green_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(1,0,r,c);
//
//                pixel = static_cast<unsigned short*>(blue_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(2,0,r,c);
                
                unsigned char *pixel;
                pixel = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(1,0,r,c);
                
                pixel = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(2,0,r,c);

            }
        }
    }
    // If input image is really a gray-scale image...
    else
    {
        for(unsigned int r=0; r<decomposed_RGB_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<decomposed_RGB_image.GetNumberOfColumns(); c++)
            {
//                unsigned short *pixel;
//                pixel = static_cast<unsigned short*>(red_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
//                
//                pixel = static_cast<unsigned short*>(green_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
//                
//                pixel = static_cast<unsigned short*>(blue_image->GetScalarPointer(c,r,0));
//                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
                
                unsigned char *pixel;
                pixel = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image_rescaled(0,0,r,c);

            }
        }
    }
    
    
    vtkSmartPointer<vtkImageAppendComponents> append_filter = vtkSmartPointer<vtkImageAppendComponents>::New();
    append_filter->SetInputDataObject(0, red_image);
    append_filter->AddInputDataObject(0, green_image);
    append_filter->AddInputDataObject(0, blue_image);
    append_filter->Update();
    
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(output_file_name);
    writer->SetInputDataObject(append_filter->GetOutput());
    writer->Write();

    return 1;
}



int vbdMelanoma::SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(bdImage &decomposed_RGB_image, const char *output_file_name)
{
    //bdGIP gip;
    //bdImage decomposed_RGB_image_rescaled;
    //gip.RescaleWholeRange(decomposed_RGB_image, 0, 255, decomposed_RGB_image_rescaled);
    
    vtkSmartPointer<vtkImageData> red_image = vtkSmartPointer<vtkImageData>::New();
    red_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    red_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//red_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);
    
    vtkSmartPointer<vtkImageData> green_image = vtkSmartPointer<vtkImageData>::New();
    green_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    green_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//green_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);
    
    vtkSmartPointer<vtkImageData> blue_image = vtkSmartPointer<vtkImageData>::New();
    blue_image->SetDimensions(decomposed_RGB_image.GetNumberOfColumns(),decomposed_RGB_image.GetNumberOfRows(),1);
    blue_image->AllocateScalars(VTK_UNSIGNED_CHAR,1);//blue_image->AllocateScalars(VTK_UNSIGNED_SHORT,1);
    
    // If input image is really a color image...
    if(decomposed_RGB_image.GetNumberOfTimeSeries()>2)
    {
        for(unsigned int r=0; r<decomposed_RGB_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<decomposed_RGB_image.GetNumberOfColumns(); c++)
            {
                unsigned char *pixel;
                pixel = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(0,0,r,c); //decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(1,0,r,c); //decomposed_RGB_image_rescaled(1,0,r,c);
                
                pixel = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(2,0,r,c);//decomposed_RGB_image_rescaled(2,0,r,c);
                
            }
        }
    }
    // If input image is really a gray-scale image...
    else
    {
        for(unsigned int r=0; r<decomposed_RGB_image.GetNumberOfRows(); r++)
        {
            for(unsigned int c=0; c<decomposed_RGB_image.GetNumberOfColumns(); c++)
            {
                unsigned char *pixel;
                pixel = static_cast<unsigned char*>(red_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(0,0,r,c); //decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(green_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(0,0,r,c); //decomposed_RGB_image_rescaled(0,0,r,c);
                
                pixel = static_cast<unsigned char*>(blue_image->GetScalarPointer(c,r,0));
                pixel[0] = decomposed_RGB_image(0,0,r,c); //decomposed_RGB_image_rescaled(0,0,r,c);
                
            }
        }
    }
    
    
    vtkSmartPointer<vtkImageAppendComponents> append_filter = vtkSmartPointer<vtkImageAppendComponents>::New();
    append_filter->SetInputDataObject(0, red_image);
    append_filter->AddInputDataObject(0, green_image);
    append_filter->AddInputDataObject(0, blue_image);
    append_filter->Update();
    
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(output_file_name);
    writer->SetInputDataObject(append_filter->GetOutput());
    writer->Write();
    
    return 1;
}


