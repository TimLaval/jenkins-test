/**************************************************************************
 DermScan Focal point analysis.

 Author: Danilo Babin
 File name: "DS_FocusAnalysis.cpp"
 ***************************************************************************/




#include <iostream>
#include <sstream>

#include "bdString.h"
#include "bdMelanoma.h"
#include "vbdMelanoma.h"
#include "tinydir.h"


int GetFileNamesForExtension(bdList<bdString> &file_names, const char *directory_path, const char *extension)
{
    
    tinydir_dir dir;
    if(tinydir_open_sorted(&dir,directory_path) == -1)
    {
        cout<<"GetFileNamesForExtension(): Error opening directory:'"<<directory_path<<"'.";
        tinydir_close(&dir);
        return 0;
    }
    
    file_names.Reset();
    
    for(int i=0; i<dir.n_files; i++) //for(int i=dir.n_files-1; i>=0; i--)
    {
        tinydir_file file;
        if(tinydir_readfile_n(&dir, &file, i) == -1)
        {
            cout<<endl<<"GetFileNamesForExtension(): Error getting file.";
            tinydir_close(&dir);
            return 0;
        }
        else
        {
            if(!file.is_dir)
            {
                bdString name_of_file;
                name_of_file.Assign(file.path);
                
                bdString file_dir, file_name_root, file_extension;
                name_of_file.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
                if(file_extension==extension)
                {
                    file_names.AddToRightEnd(name_of_file);
                }
            }
            //cout<<file.name<<endl;
        }
    }
    
    
    return 1;
}



int GetFileNamesAndZeroFocusReferenceImageName(bdList<bdString> &file_names, bdString &reference_image_file_name, const char *directory_path, const char *extension)
{
    
    tinydir_dir dir;
    if(tinydir_open_sorted(&dir,directory_path) == -1)
    {
        cout<<"GetFileNamesForExtension(): Error opening directory:'"<<directory_path<<"'.";
        tinydir_close(&dir);
        return 0;
    }
    
    file_names.Reset();
    
    for(int i=0; i<dir.n_files; i++) //for(int i=dir.n_files-1; i>=0; i--)
    {
        tinydir_file file;
        if(tinydir_readfile_n(&dir, &file, i) == -1)
        {
            cout<<endl<<"GetFileNamesForExtension(): Error getting file.";
            tinydir_close(&dir);
            return 0;
        }
        else
        {
            if(!file.is_dir)
            {
                bdString name_of_file;
                name_of_file.Assign(file.path);
                
                bdString file_dir, file_name_root, file_extension;
                name_of_file.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
                if(file_extension==extension)
                {
                    
                    file_names.AddToRightEnd(name_of_file);
                    
                    // check if this file will be the reference image
                    bdString ref_image_name;
                    ref_image_name.Assign(file.name);
                    bdList<double> list_of_extracted_doubles;
                    ref_image_name.ExtractAllDoubleNumbersToList(list_of_extracted_doubles);
                    if(list_of_extracted_doubles.GetLeftEnd()>-0.005 && list_of_extracted_doubles.GetLeftEnd()<0.005)
                    {
                        reference_image_file_name.Assign(name_of_file);
                        //cout<<endl<<"ref_image = '"<<ref_image_name<<"'"<<endl;
                    }
//                    // this part should be active only if you want to exclude the reference image from the
//                    // list of images (in that case comment the above line adding the file name to the list).
//                    else
//                    {
//                        file_names.AddToRightEnd(name_of_file);
//                    }
                }
            }
            //cout<<file.name<<endl;
        }
    }
    
    
    return 1;
}



int GetFileNamesAndZeroFocusReferenceImageName(bdList<bdString> &file_names, bdString &reference_image_file_name, const char *directory_path, const char *extension, const char *prefix)
{
    tinydir_dir dir;
    if(tinydir_open_sorted(&dir,directory_path) == -1)
    {
        cout<<"GetFileNamesForExtension(): Error opening directory:'"<<directory_path<<"'.";
        tinydir_close(&dir);
        return 0;
    }
    
    file_names.Reset();
    
    for(int i=0; i<dir.n_files; i++) //for(int i=dir.n_files-1; i>=0; i--)
    {
        tinydir_file file;
        if(tinydir_readfile_n(&dir, &file, i) == -1)
        {
            cout<<endl<<"GetFileNamesForExtension(): Error getting file.";
            tinydir_close(&dir);
            return 0;
        }
        else
        {
            if(!file.is_dir)
            {
                bdString name_of_file_without_path;
                name_of_file_without_path.Assign(file.name);
                if(name_of_file_without_path.HasPrefix(prefix))
                {
                
                    bdString name_of_file;
                    name_of_file.Assign(file.path);
                    
                    bdString file_dir, file_name_root, file_extension;
                    name_of_file.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
                    if(file_extension==extension)
                    {
                        
                        file_names.AddToRightEnd(name_of_file);
                        
                        // check if this file will be the reference image
                        bdString ref_image_name;
                        ref_image_name.Assign(file.name);
                        bdList<double> list_of_extracted_doubles;
                        ref_image_name.ExtractAllDoubleNumbersToList(list_of_extracted_doubles);
                        if(list_of_extracted_doubles.GetLeftEnd()>-0.005 && list_of_extracted_doubles.GetLeftEnd()<0.005)
                        {
                            reference_image_file_name.Assign(name_of_file);
                            //cout<<endl<<"ref_image = '"<<ref_image_name<<"'"<<endl;
                        }
                        //                    // this part should be active only if you want to exclude the reference image from the
                        //                    // list of images (in that case comment the above line adding the file name to the list).
                        //                    else
                        //                    {
                        //                        file_names.AddToRightEnd(name_of_file);
                        //                    }
                    }
                }
            }
            //cout<<file.name<<endl;
        }
    }
    
    
    return 1;
}






int main(int argc, char** argv)
{
    bdString directory_path;
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_test002_20170607092053");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_test003_20170607094448");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_test004_20170607121159");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_Test005_20170608095136");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_Test006_20170608113501");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_06_22/capture_test007_20170608122203");
    //directory_path.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_07_04/capture_g125346_20170429084909/png");
    
    if ( argc != 2 )
    {
        cout<<"usage: "<< argv[0] <<" <directory_path>\n";
        return 0;
    }
    else
    {
        directory_path.Assign(argv[1]);
    }
    
    
    bdArray<bdString> file_name_prefixes;
    file_name_prefixes.Set(5);
    {
        int c = 0;
        file_name_prefixes[c].Assign("output_color_blue_polarized"); c++;
        file_name_prefixes[c].Assign("output_color_blue_unpolarized"); c++;
        file_name_prefixes[c].Assign("output_color_deep_red_polarized"); c++;
        file_name_prefixes[c].Assign("output_color_far_red_polarized"); c++;
        file_name_prefixes[c].Assign("output_color_green_polarized"); c++;
    }
    
    
    // Segmented reference image is defined here as it will be the same for all bands.
    bdImage segmented_reference_image;
    bdImage ROI_image; // contains region of interest for focus analysis
    
    
    for(unsigned int file_name_prefix_counter=0; file_name_prefix_counter<file_name_prefixes.GetNumberOfElements(); file_name_prefix_counter++)
    {
    
    
        bdString reference_image_file_name;
        bdList<bdString> file_names;
        GetFileNamesAndZeroFocusReferenceImageName(file_names, reference_image_file_name, directory_path.C_String(), "png", file_name_prefixes[file_name_prefix_counter].C_String());
        
        //----- List the files being loaded -----
        bdListIterator<bdString> it;
        cout<<"Loading files:"<<endl;
        for(it.SetLeftEnd(file_names); it.IsValid(); it.MoveRight())
        {
            //it.GetElement().
            cout<<it.GetElement()<<endl;
        }
        cout<<"Reference image: "<<reference_image_file_name<<endl;
        //----------
        
        
        //----- Data used for various steps (global data) -----
        bdImage reference_image;
        bdArray<bdImage> focus_images;
        bdArray<bdImage> subtracted_images;
        bdArray<bdImage> normalized_subtracted_images;
        bdArray<bdImage> blurred_images;
        //----------
        
        
        
        //----- File names for writing images based on input file name -----
        bdString file_dir, file_name_root, file_extension;
        reference_image_file_name.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
        bdString dir_name_processed_images;
        dir_name_processed_images.Assign(file_dir);
        dir_name_processed_images.Append("processed");
        dir_name_processed_images.AppendDirSeparator();
        
        bdString file_name_segmented_reference_image;
        file_name_segmented_reference_image.Assign(dir_name_processed_images);
        file_name_segmented_reference_image.Append(file_name_prefixes[file_name_prefix_counter]);
        file_name_segmented_reference_image.Append("_segmented.png");
        
        bdString file_name_ROI_reference_image;
        file_name_ROI_reference_image.Assign(dir_name_processed_images);
        file_name_ROI_reference_image.Append(file_name_prefixes[file_name_prefix_counter]);
        file_name_ROI_reference_image.Append("_ROI.png");
        
        bdArray<bdString> file_names_of_subtracted_images;
        file_names_of_subtracted_images.Set(file_names.GetNumberOfElements());
        for(unsigned int i=0; i<file_names_of_subtracted_images.GetNumberOfElements(); i++)
        {
            file_names_of_subtracted_images[i].Assign(file_dir);
            file_names_of_subtracted_images[i].Append("processed/");
            bdString temp;
            temp.Assign(file_names[i]);
            temp.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
            file_names_of_subtracted_images[i].Append(file_name_root);
            file_names_of_subtracted_images[i].Append("_subtracted.png");
        }
        
        bdArray<bdString> file_names_of_normalized_subtracted_images;
        file_names_of_normalized_subtracted_images.Set(file_names.GetNumberOfElements());
        for(unsigned int i=0; i<file_names_of_normalized_subtracted_images.GetNumberOfElements(); i++)
        {
            file_names_of_normalized_subtracted_images[i].Assign(file_dir);
            file_names_of_normalized_subtracted_images[i].Append("processed/");
            bdString temp;
            temp.Assign(file_names[i]);
            temp.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
            file_names_of_normalized_subtracted_images[i].Append(file_name_root);
            file_names_of_normalized_subtracted_images[i].Append("_subtracted_normalized.png");
        }
        //----------
        
        
        //----- Load the reference image as a white light image -----
        {
            cout<<"Loading image reference image...";
            vbdMelanoma v_melanoma;
            v_melanoma.LoadWhiteLightRGBImageFromPNGFile(reference_image_file_name.C_String(), reference_image);
            cout<<" completed."<<endl;
        }
        //-----------
        
        
        //----- Segment melanoma from reference image and make ROI image ----- (we perform this only on BLUE POLARIZED band)
        if(file_name_prefix_counter==0)
        {
            {
                cout<<"Segmenting melanoma from reference image...";
                bdMelanoma melanoma;
                melanoma.SegmentMelanoma_Thresholding(reference_image, segmented_reference_image);
                cout<<" completed."<<endl;
            }
            //-----------
            //----- Save the images -----
            {
                cout<<"Saving segmented image...";
                bdGIP gip;
                bdImage segmented_image_rescaled;//vbdBasicImage segmented_image_rescaled;
                gip.RescaleWholeRange(segmented_reference_image, 0, 255, segmented_image_rescaled);
                vbdMelanoma v_melanoma;
                cout<<" "<<file_name_segmented_reference_image.C_String()<<" ";
                v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(segmented_image_rescaled, file_name_segmented_reference_image.C_String());
                cout<<" completed."<<endl;
            }
            //-----------
        
        
            //----- Dilate the segmented image to obtain larger ROI -----
            {
                cout<<"ROI extraction... ";
                bdMelanoma melanoma;
                melanoma.DilateBinaryImageWithCircleSE(segmented_reference_image, ROI_image, 400);
                cout<<" completed."<<endl;
            }
            //-----------
            //----- Save the images -----
            {
                cout<<"Saving ROI image...";
                bdGIP gip;
                bdImage ROI_image_rescaled;
                gip.RescaleWholeRange(ROI_image, 0, 255, ROI_image_rescaled);
                vbdMelanoma v_melanoma;
                v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(ROI_image_rescaled, file_name_ROI_reference_image.C_String());
                cout<<" completed."<<endl;
            }
            //-----------
        }
        //-----------

        
        //----- Load the images of different focus values -----
        {
            cout<<"Loading images...";
            vbdMelanoma v_melanoma;
            focus_images.Set(file_names.GetNumberOfElements());
            for(unsigned int i=0; i<focus_images.GetNumberOfElements(); i++)
            {
                v_melanoma.LoadWhiteLightRGBImageFromPNGFile(file_names[i].C_String(), focus_images[i]);
            }
            cout<<" completed."<<endl;
        }
        //-----------
        
        
        
        //----- Compute blurred images from reference image -----
        {
            cout<<"Computing blurred images from reference image...";
            bdMelanoma melanoma;
            unsigned int number_of_blur_steps = 11;
            blurred_images.Set(number_of_blur_steps);
            for(unsigned int blur_step = 0; blur_step<number_of_blur_steps; blur_step++)
            {
                cout<<"["<<blur_step<<"/"<<number_of_blur_steps-1<<"]";
                melanoma.MeanOperatorPerRGBComponent(reference_image,ROI_image,(blur_step*blur_step), blurred_images[blur_step]);
            }
            cout<<" completed."<<endl;
        }
        //-----------
        
        
        //----- Subtract images -----
        {
            cout<<"Subtracting images...";
            bdMelanoma melanoma;
            subtracted_images.Set(focus_images.GetNumberOfElements());
            for(unsigned int i=0; i<focus_images.GetNumberOfElements(); i++)
            {
                // Each of the focus images subtract with all blurred images and save the result that yields the lowest score.
                double minimum_subtract_score = -1;
                int index_of_minimum_subtract_score = -1;
                
                for(unsigned int b = 0; b<blurred_images.GetNumberOfElements(); b++)
                {
                    bdImage subtracted_image_temp;
                    melanoma.SubtractImages(blurred_images[b], focus_images[i], ROI_image, subtracted_image_temp);
                    double subtract_score = melanoma.SumOfVoxelValuesForAllRGBComponents(subtracted_image_temp,ROI_image);
                    if(minimum_subtract_score<0 || minimum_subtract_score>subtract_score)
                    {
                        subtracted_images[i].CopyFrom(subtracted_image_temp);
                        minimum_subtract_score = subtract_score;
                        index_of_minimum_subtract_score = b;
                    }
                    cout<<" subtract_score("<<i<<","<<b<<")="<<subtract_score;
                }
                cout<<endl<<" min_subtract_score("<<i<<") for blur SE radius "<<index_of_minimum_subtract_score<<" is "<<minimum_subtract_score<<endl;
            }
            cout<<" completed."<<endl;
        }
        //-----------
        //----- Save the images -----
        {
            cout<<"Saving subtracted images...";
            bdGIP gip;
            for(unsigned int i=0; i<subtracted_images.GetNumberOfElements(); i++)
            {
                bdImage subtracted_image_rescaled;
                gip.RescaleWholeRange(subtracted_images[i], 0, 255, subtracted_image_rescaled);
                vbdMelanoma v_melanoma;
                v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(subtracted_image_rescaled, file_names_of_subtracted_images[i].C_String());
            }
            cout<<" completed."<<endl;
        }
        //-----------


        //----- Normalize subtracted images for better visualization -----
        {
            cout<<"Normalizing/leveling subtracted images...";
            bdMelanoma melanoma;
            normalized_subtracted_images.Set(focus_images.GetNumberOfElements());
            for(unsigned int i=0; i<focus_images.GetNumberOfElements(); i++)
            {
                melanoma.LevelComponentsOfRGBImage(subtracted_images[i], 25, normalized_subtracted_images[i]);
                //melanoma.NormalizeWholeRGBRangeImage(subtracted_images[i], 0, 25, normalized_subtracted_images[i]);
            }
            cout<<" completed."<<endl;
        }
        //-----------
        //----- Save the images -----
        {
            cout<<"Saving normalized/leveled subtracted images...";
            bdGIP gip;
            for(unsigned int i=0; i<subtracted_images.GetNumberOfElements(); i++)
            {
                bdImage subtracted_normalized_image_rescaled;
                gip.RescaleWholeRange(normalized_subtracted_images[i], 0, 255, subtracted_normalized_image_rescaled);
                vbdMelanoma v_melanoma;
                v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(subtracted_normalized_image_rescaled, file_names_of_normalized_subtracted_images[i].C_String());
            }
            cout<<" completed."<<endl;
        }
        //-----------

    }
    

    cout<<endl<<"END."<<endl;

    

    
    return 1;
}
