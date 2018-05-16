/**************************************************************************
 Example of using VTK over libraries.

 Author: Danilo Babin
 File name: "DS_FusedRGB.cpp"
***************************************************************************/



#include <vtkSmartPointer.h>
#include <vtkAxes.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPropPicker.h>
#include <vtkPicker.h>
#include <vtkCellPicker.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPNGReader.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageAppendComponents.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkImageWeightedSum.h>
#include <vtkImageCast.h>


//#include <conio.h>
//#include <math.h>
#include <iostream>
#include <sstream>

#include "bdArray.h"
#include "bdMatrix.h"
#include "bdString.h"



int main(int argc, char** argv)
{
    
	//===== THIS IS GENERAL VTK PART - NEEDED FOR ALL EXAMPLES =====
	vtkSmartPointer<vtkRenderer> renderer = vtkRenderer::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);	
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);
	//renWin->SetSize(900, 900);
	//iren->Initialize();//Initialize the interactor here in order for the SplineWighet to work properly
	//==============================================================

    
    
//    // Read the image
//    vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
//    reader->SetFileName("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000001/ISIC_0000001.png");
//    reader->Update();
//    
//    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New(); //reader->GetOutput();
//    image->DeepCopy(reader->GetOutput());
//    
//    
//    vtkSmartPointer<vtkImageExtractComponents> extractRedFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
//    //extractRedFilter->SetInputConnection(reader->GetOutputPort());
//    extractRedFilter->SetInputDataObject(image);
//    extractRedFilter->SetComponents(0);
//    extractRedFilter->Update();
//    
//    vtkSmartPointer<vtkImageExtractComponents> extractGreenFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
//    //extractGreenFilter->SetInputConnection(reader->GetOutputPort());
//    extractGreenFilter->SetInputDataObject(image);
//    extractGreenFilter->SetComponents(1);
//    extractGreenFilter->Update();
//    
//    vtkSmartPointer<vtkImageExtractComponents> extractBlueFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
//    //extractBlueFilter->SetInputConnection(reader->GetOutputPort());
//    extractBlueFilter->SetInputDataObject(image);
//    extractBlueFilter->SetComponents(2);
//    extractBlueFilter->Update();
//    
//    vtkSmartPointer<vtkImageData> red_image = vtkSmartPointer<vtkImageData>::New();
//    red_image->DeepCopy(extractRedFilter->GetOutput());
//    
//    vtkSmartPointer<vtkImageData> green_image = vtkSmartPointer<vtkImageData>::New();
//    green_image->DeepCopy(extractGreenFilter->GetOutput());
//    
//    vtkSmartPointer<vtkImageData> blue_image = vtkSmartPointer<vtkImageData>::New();
//    blue_image->DeepCopy(extractBlueFilter->GetOutput());
//    
//    //vtkImageData *red_image = extractRedFilter->GetOutput();
//    //vtkImageData *green_image = extractGreenFilter->GetOutput();
//    //vtkImageData *blue_image = extractBlueFilter->GetOutput();
//    
//    vtkSmartPointer<vtkImageAppendComponents> appendFilter = vtkSmartPointer<vtkImageAppendComponents>::New();
//    appendFilter->SetInputDataObject(0, red_image);
//    appendFilter->AddInputDataObject(0, green_image);
//    appendFilter->AddInputDataObject(0, blue_image);
//    appendFilter->Update();
//    
//    //vtkImageData *combined_image = appendFilter->GetOutput();
//    
//    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
//    writer->SetFileName("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000001/cobined.png");
//    writer->SetInputDataObject(appendFilter->GetOutput());//appendFilter->GetOutputPort());
//    writer->Write();
    
    int number_of_files_to_use = 7;
    int number_of_files = 7;
    
    bdArray<bdString> file_names;
    file_names.Set(number_of_files);
    
    file_names[0].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_006_inner_royalblue.png");
    file_names[1].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_007_inner_uv.png");
    file_names[2].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_007_outer_green.png");
    file_names[3].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_007_outer_royalblue.png");
    file_names[4].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_008_outer_deepred.png");
    file_names[5].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_008_outer_farred.png");
    file_names[6].Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/selected/image_009_outer_uv.png");
    
    // Read the images
    bdArray< vtkSmartPointer<vtkImageData> > loaded_images;
    loaded_images.Set(number_of_files_to_use);
    for(unsigned int i=0; i<loaded_images.GetNumberOfElements(); i++)
    {
        vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
        reader->SetFileName(file_names[i].C_String());
        reader->Update();
        loaded_images[i] = vtkSmartPointer<vtkImageData>::New();
        loaded_images[i]->DeepCopy(reader->GetOutput());
    }
    
    bdArray< vtkSmartPointer<vtkImageData> > red_images, green_images, blue_images;
    red_images.Set(number_of_files_to_use);
    green_images.Set(number_of_files_to_use);
    blue_images.Set(number_of_files_to_use);
    for(unsigned int i=0; i<red_images.GetNumberOfElements(); i++)
    {
        {
            vtkSmartPointer<vtkImageExtractComponents> extract_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
            extract_filter->SetInputDataObject(loaded_images[i]);
            extract_filter->SetComponents(0);
            extract_filter->Update();
            red_images[i] = vtkSmartPointer<vtkImageData>::New();
            red_images[i]->DeepCopy(extract_filter->GetOutput());
        }
        {
            vtkSmartPointer<vtkImageExtractComponents> extract_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
            extract_filter->SetInputDataObject(loaded_images[i]);
            extract_filter->SetComponents(1);
            extract_filter->Update();
            green_images[i] = vtkSmartPointer<vtkImageData>::New();
            green_images[i]->DeepCopy(extract_filter->GetOutput());
        }
        {
            vtkSmartPointer<vtkImageExtractComponents> extract_filter = vtkSmartPointer<vtkImageExtractComponents>::New();
            extract_filter->SetInputDataObject(loaded_images[i]);
            extract_filter->SetComponents(2);
            extract_filter->Update();
            blue_images[i] = vtkSmartPointer<vtkImageData>::New();
            blue_images[i]->DeepCopy(extract_filter->GetOutput());
        }
    }
    
    int dimension_rows = loaded_images[0]->GetDimensions()[1];
    int dimension_columns = loaded_images[0]->GetDimensions()[0];
    
    bdMatrix<double> combined_red, combined_green, combined_blue;
    combined_red.SetSize(dimension_rows, dimension_columns);
    combined_green.SetSize(dimension_rows, dimension_columns);
    combined_blue.SetSize(dimension_rows, dimension_columns);
    combined_red.FillInWith(0);
    combined_green.FillInWith(0);
    combined_blue.FillInWith(0);
    
    
    for(unsigned int r=0; r<dimension_rows; r++)
    {
        for(unsigned int c=0; c<dimension_columns; c++)
        {
            for(unsigned int i=0; i<red_images.GetNumberOfElements(); i++)
            {
                {
                    unsigned char *pixel = static_cast<unsigned char*>(red_images[i]->GetScalarPointer(c,r,0));
                    combined_red(r,c) += pixel[0];
                }
                {
                    unsigned char *pixel = static_cast<unsigned char*>(green_images[i]->GetScalarPointer(c,r,0));
                    combined_green(r,c) += pixel[0];
                }
                {
                    unsigned char *pixel = static_cast<unsigned char*>(blue_images[i]->GetScalarPointer(c,r,0));
                    combined_blue(r,c) += pixel[0];
                }
            }
            
            combined_red(r,c) = combined_red(r,c) / ((double)loaded_images.GetNumberOfElements());
            combined_green(r,c) = combined_green(r,c) / ((double)loaded_images.GetNumberOfElements());
            combined_blue(r,c) = combined_blue(r,c) / ((double)loaded_images.GetNumberOfElements());
            
            {
                unsigned char *pixel = static_cast<unsigned char*>(red_images[0]->GetScalarPointer(c,r,0));
                pixel[0] = (unsigned char) combined_red(r,c);
            }
            {
                unsigned char *pixel = static_cast<unsigned char*>(green_images[0]->GetScalarPointer(c,r,0));
                pixel[0] = (unsigned char) combined_green(r,c);
            }
            {
                unsigned char *pixel = static_cast<unsigned char*>(blue_images[0]->GetScalarPointer(c,r,0));
                pixel[0] = (unsigned char) combined_blue(r,c);
            }
        }
    }

    
    vtkSmartPointer<vtkImageAppendComponents> appendFilter = vtkSmartPointer<vtkImageAppendComponents>::New();
    appendFilter->SetInputDataObject(0, red_images[0]);
    appendFilter->AddInputDataObject(0, green_images[0]);
    appendFilter->AddInputDataObject(0, blue_images[0]);
    appendFilter->Update();


    
    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2016_12_09/combined.png");
    writer->SetInputDataObject(appendFilter->GetOutput());//appendFilter->GetOutput());
    writer->Write();
    
    
    
    

	//===== START THE RENDERING =====
	renWin->Render();
	iren->Start();
	//=====
    
    
    return 1;
}
