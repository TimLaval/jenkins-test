/**************************************************************************
 DermScan ABCD(E) feature quantification.

 Author: Danilo Babin
 File name: "DS_ABCD.cpp"
 ***************************************************************************/




#include <iostream>
#include <sstream>

#include "bdString.h"
#include "bdMelanoma.h"
#include "vbdMelanoma.h"



//#define PATH_TO_BROWN_CHART "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/brown_red9x9.png"
//#define PATH_TO_RED_CHART "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/red9x9.png"
//#define PATH_TO_BLUE_CHART "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/blue9x9.png"

//#define PATH_TO_BROWN_MAP "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/brown9x9map.png"
//#define PATH_TO_BLUE_MAP "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/blue9x9map.png"
//#define PATH_TO_RED_MAP "/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/red9x9map.png"



int main(int argc, char** argv)
{
    bdString file_name_white_light_image;
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000002/ISIC_0000002.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000001/ISIC_0000001.png");// Streaks.
    
    
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_g125346_20170429084909_white_unpolarized_focus_0_with_csc.png");//LARGE
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test002_20170607092053_white_unpolarized_focus_0_with_csc.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test003_20170607094448_white_unpolarized_focus_0_with_csc.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test004_20170607121159_white_unpolarized_focus_0_with_csc.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_Test005_20170608095136_white_unpolarized_focus_0_with_csc.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_Test006_20170608113501_white_unpolarized_focus_0_with_csc.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test007_20170608122203_white_unpolarized_focus_0_with_csc.png");//HAIRS
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test007_20170608122706_white_unpolarized_focus_0_with_csc.png");//HAIRS
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_test007_20170608123158_white_unpolarized_focus_0_with_csc.png");//HAIRS
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/white_light_images/unpolarized/capture_Test008_20170608131129_white_unpolarized_focus_0_with_csc.png");
    
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/capture_Test012_20170629101442/output_color_white_unpolarized_focus_0_exposure_66700_intensity_0.07_with_csc.png");//Blue Marker
    
    
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/capture_63087607_20170728084233/output_color_white_unpolarized_focus_0_exposure_66700_intensity_0.07_with_csc.png");// Net structure with hairs
    
    //---------- D feature test images ----------
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_61277491_20170802115723.png");// Streaks?
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_62587217_20170705064239.png");// Dark melanoma with hardly visible network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_63268684_20170906143249.png");// Network
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_63268684_20170906143322.png");// Hardly visibile network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_63268684_20170906143352.png");// NO visibile network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_63268684_20170906143523.png");// NO visibile network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_67416909_20170712125217.png");// Very very little network
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_70638960_20170828090924.png");// Network
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_73893430_20170906071225.png");// Dots with a little bit of network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_78804945_20170829134343.png");// Globules - Ok, but a bit larger.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_81022170_20170830082903.png");// Globules.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/capture_86361540_20170828075901.png");// Dots.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature/ISIC_0000000.png");// Dots.

        //-------------------------------------------
    
    //---------- C feature test images ----------
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_60924040_20170712143129.png");// Red lesion
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_62587217_20170705064506.png");// Melanoma D/L brown
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_73893430_20170906071409.png");// Melanoma D/L brown
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_78379450_20170705095646.png");// Red and brown D/L
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_80518590_20170919114308.png");// Black spot (some brown)
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_82079542_20170706125837.png");// Red lesion
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_86250594_20170802134149.png");// Red and brown lesion (hardly visible)
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_test010_20170621103229.png");// Melanoma D/L brown
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_test016_20170705101416.png");// Red lesion with dark (red/brown) spots
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/capture_test016_20170705101627.png");// L brown lesion with some red part and D brown spots
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000127.png");// White veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000140.png");//  Black, brown, blue veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000142.png");//  Black, brown, gray veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000181.png");//  Brown melanoma with some blue/gray veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000192.png");//  Blue-gray veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000235.png");// Blue-gray veil
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_C_Feature/ISIC_0000477.png");// Blue-gray veil
    //-------------------------------------------
    
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature_G/capture_Test006_20170608113805.png");// Globules.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature_G/capture_70638960_20170828090924.png");// Network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature_G/capture_61277491_20170802115723.png");// Nothing special, a bit of Network.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_D_Feature_G/capture_73893430_20170906071441.png"); //Some dots and a streak
    
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_B_Feature_G/capture_test125_20171116122806.png");// Melanoma partly sharp edges.
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Barco_Images_2017_10_17_16bit/-Test_B_Feature_G/capture_unidentified_20170906095446.png"); //Melanoma sharp edges

    

    
    
    if ( argc != 2 )
    {
        cout<<"usage: "<< argv[0] <<" <filename>\n";
        return 0;
    }
    else
    {
        file_name_white_light_image.Assign(argv[1]);
    }
    
    
    //----- Data used for various steps (global data) -----
    bdImage white_light_image;
    bdImage extracted_channel_red_image;
    bdImage extracted_channel_green_image;
    bdImage extracted_channel_blue_image;
    bdImage extracted_channel_max_RGB_image;
    bdImage segmented_image;
    bdImage mask_image;
    bdImage asymmetry_image;
    bdImage border_image;
    bdImage color_mapped_image;
    bdImage color_chart_mapped_image;
    bdImage inner_bright_structure_image;
    bdImage inner_dark_structure_image;
    bdImage inner_dark_streaks_image;
    bdImage inner_structureless_image;
    bdImage inner_combined_structure_image;
    //----------
    
    
    
    //----- File names for writing images based on input file name -----
    bdString file_dir, file_name_root, file_extension;
    file_name_white_light_image.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
    
    bdString file_name_extracted_channel_red_image;
    file_name_extracted_channel_red_image.Assign(file_dir);
    file_name_extracted_channel_red_image.Append(file_name_root);
    file_name_extracted_channel_red_image.Append("_extracted_channel_red.");
    file_name_extracted_channel_red_image.Append(file_extension);
    bdString file_name_extracted_channel_green_image;
    file_name_extracted_channel_green_image.Assign(file_dir);
    file_name_extracted_channel_green_image.Append(file_name_root);
    file_name_extracted_channel_green_image.Append("_extracted_channel_green.");
    file_name_extracted_channel_green_image.Append(file_extension);
    bdString file_name_extracted_channel_blue_image;
    file_name_extracted_channel_blue_image.Assign(file_dir);
    file_name_extracted_channel_blue_image.Append(file_name_root);
    file_name_extracted_channel_blue_image.Append("_extracted_channel_blue.");
    file_name_extracted_channel_blue_image.Append(file_extension);
    bdString file_name_extracted_channel_max_RGB_image;
    file_name_extracted_channel_max_RGB_image.Assign(file_dir);
    file_name_extracted_channel_max_RGB_image.Append(file_name_root);
    file_name_extracted_channel_max_RGB_image.Append("_extracted_channel_RGB_max.");
    file_name_extracted_channel_max_RGB_image.Append(file_extension);
    
    bdString file_name_modified_white_light_image;
    file_name_modified_white_light_image.Assign(file_dir);
    file_name_modified_white_light_image.Append(file_name_root);
    file_name_modified_white_light_image.Append("_modified_white_light.");
    file_name_modified_white_light_image.Append(file_extension);

    
    
    bdString file_name_extracted_channel_red_histogram;
    file_name_extracted_channel_red_histogram.Assign(file_dir);
    file_name_extracted_channel_red_histogram.Append(file_name_root);
    file_name_extracted_channel_red_histogram.Append("_extracted_channel_red_histogram.");
    file_name_extracted_channel_red_histogram.Append(".m");
    bdString file_name_extracted_channel_green_histogram;
    file_name_extracted_channel_green_histogram.Assign(file_dir);
    file_name_extracted_channel_green_histogram.Append(file_name_root);
    file_name_extracted_channel_green_histogram.Append("_extracted_channel_green_histogram.");
    file_name_extracted_channel_green_histogram.Append(".m");
    bdString file_name_extracted_channel_blue_histogram;
    file_name_extracted_channel_blue_histogram.Assign(file_dir);
    file_name_extracted_channel_blue_histogram.Append(file_name_root);
    file_name_extracted_channel_blue_histogram.Append("_extracted_channel_blue_histogram.");
    file_name_extracted_channel_blue_histogram.Append(".m");
    bdString file_name_extracted_channel_max_RGB_histogram;
    file_name_extracted_channel_max_RGB_histogram.Assign(file_dir);
    file_name_extracted_channel_max_RGB_histogram.Append(file_name_root);
    file_name_extracted_channel_max_RGB_histogram.Append("_extracted_channel_max_RGB_histogram.");
    file_name_extracted_channel_max_RGB_histogram.Append(".m");
    
    bdString file_name_extracted_channel_red_histogram2;
    file_name_extracted_channel_red_histogram2.Assign(file_dir);
    file_name_extracted_channel_red_histogram2.Append(file_name_root);
    file_name_extracted_channel_red_histogram2.Append("_extracted_channel_red_histogram2.");
    file_name_extracted_channel_red_histogram2.Append(".m");
    bdString file_name_extracted_channel_green_histogram2;
    file_name_extracted_channel_green_histogram2.Assign(file_dir);
    file_name_extracted_channel_green_histogram2.Append(file_name_root);
    file_name_extracted_channel_green_histogram2.Append("_extracted_channel_green_histogram2.");
    file_name_extracted_channel_green_histogram2.Append(".m");
    bdString file_name_extracted_channel_blue_histogram2;
    file_name_extracted_channel_blue_histogram2.Assign(file_dir);
    file_name_extracted_channel_blue_histogram2.Append(file_name_root);
    file_name_extracted_channel_blue_histogram2.Append("_extracted_channel_blue_histogram2.");
    file_name_extracted_channel_blue_histogram2.Append(".m");
    

    
    bdString file_name_extracted_channel_red_masked_histogram;
    file_name_extracted_channel_red_masked_histogram.Assign(file_dir);
    file_name_extracted_channel_red_masked_histogram.Append(file_name_root);
    file_name_extracted_channel_red_masked_histogram.Append("_extracted_channel_red_masked_histogram.");
    file_name_extracted_channel_red_masked_histogram.Append(".m");
    bdString file_name_extracted_channel_green_masked_histogram;
    file_name_extracted_channel_green_masked_histogram.Assign(file_dir);
    file_name_extracted_channel_green_masked_histogram.Append(file_name_root);
    file_name_extracted_channel_green_masked_histogram.Append("_extracted_channel_green_masked_histogram.");
    file_name_extracted_channel_green_masked_histogram.Append(".m");
    bdString file_name_extracted_channel_blue_masked_histogram;
    file_name_extracted_channel_blue_masked_histogram.Assign(file_dir);
    file_name_extracted_channel_blue_masked_histogram.Append(file_name_root);
    file_name_extracted_channel_blue_masked_histogram.Append("_extracted_channel_blue_masked_histogram.");
    file_name_extracted_channel_blue_masked_histogram.Append(".m");
    bdString file_name_extracted_channel_max_RGB_masked_histogram;
    file_name_extracted_channel_max_RGB_masked_histogram.Assign(file_dir);
    file_name_extracted_channel_max_RGB_masked_histogram.Append(file_name_root);
    file_name_extracted_channel_max_RGB_masked_histogram.Append("_extracted_channel_max_RGB_masked_histogram.");
    file_name_extracted_channel_max_RGB_masked_histogram.Append(".m");


    
    bdString file_name_segmented_image;
    file_name_segmented_image.Assign(file_dir);
    file_name_segmented_image.Append(file_name_root);
    file_name_segmented_image.Append("_segmented.");
    file_name_segmented_image.Append(file_extension);
    bdString file_name_mask_image;
    file_name_mask_image.Assign(file_dir);
    file_name_mask_image.Append(file_name_root);
    file_name_mask_image.Append("_mask.");
    file_name_mask_image.Append(file_extension);
    bdString file_name_asymmetry_image;
    file_name_asymmetry_image.Assign(file_dir);
    file_name_asymmetry_image.Append(file_name_root);
    file_name_asymmetry_image.Append("_asymmetry.");
    file_name_asymmetry_image.Append(file_extension);
    bdString file_name_border_image;
    file_name_border_image.Assign(file_dir);
    file_name_border_image.Append(file_name_root);
    file_name_border_image.Append("_border.");
    file_name_border_image.Append(file_extension);
    bdString file_name_color_mapped_image;
    file_name_color_mapped_image.Assign(file_dir);
    file_name_color_mapped_image.Append(file_name_root);
    file_name_color_mapped_image.Append("_color_mapped.");
    file_name_color_mapped_image.Append(file_extension);
    bdString file_name_color_chart_mapped_image;
    file_name_color_chart_mapped_image.Assign(file_dir);
    file_name_color_chart_mapped_image.Append(file_name_root);
    file_name_color_chart_mapped_image.Append("_color_chart_mapped.");
    file_name_color_chart_mapped_image.Append(file_extension);
    bdString file_name_inner_bright_structure_image;
    file_name_inner_bright_structure_image.Assign(file_dir);
    file_name_inner_bright_structure_image.Append(file_name_root);
    file_name_inner_bright_structure_image.Append("_inner_bright_structure.");
    file_name_inner_bright_structure_image.Append(file_extension);
    bdString file_name_inner_dark_structure_image;
    file_name_inner_dark_structure_image.Assign(file_dir);
    file_name_inner_dark_structure_image.Append(file_name_root);
    file_name_inner_dark_structure_image.Append("_inner_dark_structure.");
    file_name_inner_dark_structure_image.Append(file_extension);
    bdString file_name_inner_dark_streaks_image;
    file_name_inner_dark_streaks_image.Assign(file_dir);
    file_name_inner_dark_streaks_image.Append(file_name_root);
    file_name_inner_dark_streaks_image.Append("_inner_dark_streaks_image.");
    file_name_inner_dark_streaks_image.Append(file_extension);
    bdString file_name_inner_structureless_image;
    file_name_inner_structureless_image.Assign(file_dir);
    file_name_inner_structureless_image.Append(file_name_root);
    file_name_inner_structureless_image.Append("_inner_structureless_image.");
    file_name_inner_structureless_image.Append(file_extension);
    bdString file_name_inner_combined_structure_image;
    file_name_inner_combined_structure_image.Assign(file_dir);
    file_name_inner_combined_structure_image.Append(file_name_root);
    file_name_inner_combined_structure_image.Append("_inner_combined_structure.");
    file_name_inner_combined_structure_image.Append(file_extension);
    //----------

    
    
    
    bdMelanoma melanoma;

    
    
    //----- Load the white light image -----
    {
        cout<<"Loading image...";
        vbdMelanoma v_melanoma;
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(file_name_white_light_image.C_String(), white_light_image);
        cout<<" completed."<<endl;
    }
    //-----------
    
    
    //----- Perform downsampling (Comment this part out if you don't want downsampling) -----
    {
        bdImage temp;
        temp.CopyFrom(white_light_image);
        //bdMelanoma melanoma;
        //melanoma.Downsample_Mean(temp, white_light_image, 5);
        melanoma.Downsample_Max(temp, white_light_image, 4);
    }
    //-----------
    
    
    //----- Perform value rescaling (Comment this part out if you don't want rescaling) -----
    {
        bdImage temp;
        temp.CopyFrom(white_light_image);
        bdGIP gip;
        gip.RescaleWholeRange(temp,0, 255, white_light_image);
    }
    //-----------
    
    
//    //----- Extract images and save histograms for Debug -----
//    {
//        cout<<"Extracting images for debug...";
//        bdMelanoma melanoma;
//        melanoma.Extract_RedImage(white_light_image, extracted_channel_red_image);
//        melanoma.Extract_GreenImage(white_light_image, extracted_channel_green_image);
//        melanoma.Extract_BlueImage(white_light_image, extracted_channel_blue_image);
//        melanoma.Extract_MaxOfRGBImage(white_light_image, extracted_channel_max_RGB_image);  //Extract_MaxOfGreenAndBlueImage(white_light_image, extracted_channel_max_RGB_image);
//        
//        bdCurveXY histogram_R, histogram_G, histogram_B, histogram_max_GB;
//        unsigned int min, max;
//        
//        extracted_channel_red_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_red_image, histogram_R, min, max);
//        histogram_R.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_red_histogram2.C_String());
//        
//        extracted_channel_green_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_green_image, histogram_G, min, max);
//        histogram_G.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_green_histogram2.C_String());
//        
//        extracted_channel_blue_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_blue_image, histogram_B, min, max);
//        histogram_B.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_blue_histogram2.C_String());
//          
//        cout<<" completed."<<endl;
//    }
//    //-----------


    
    //----- Perform equalization (Comment this part out if you don't want equalization) -----
    {
        bdImage temp;
        temp.CopyFrom(white_light_image);
        //bdMelanoma melanoma;
        melanoma.EqualizeBrightnessOfChannels(temp, white_light_image);
    }
    //-----------
//    //----- Perform value rescaling (Comment this part out if you don't want rescaling) -----
//    {
//        bdImage temp;
//        temp.CopyFrom(white_light_image);
//        bdGIP gip;
//        gip.RescaleWholeRange(temp,0, 255, white_light_image);
//    }
//    //-----------


    
    //----- Save the modified white light image -----
    {
        cout<<"Saving modified white light image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(white_light_image, file_name_modified_white_light_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

    
    
    
    //----- Segment melanoma -----
    {
        cout<<"Segmenting melanoma...";
        //bdMelanoma melanoma;
        //melanoma.SegmentMelanoma_Thresholding(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding2(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding3(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding4(white_light_image, segmented_image, segmented_image2);
        //melanoma.SegmentMelanoma_Thresholding5(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding7(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding8(white_light_image, segmented_image);
        melanoma.SegmentMelanoma_Thresholding9(white_light_image, segmented_image);

        cout<<" completed."<<endl;
    }
    //-----------
    //----- Save the images -----
    {
        cout<<"Saving segmented image(s)...";
        bdGIP gip;
        bdImage segmented_image_rescaled;
        gip.RescaleWholeRange(segmented_image, 0, 255, segmented_image_rescaled);
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(segmented_image_rescaled, file_name_segmented_image.C_String());
        
        //gip.RescaleWholeRange(segmented_image2, 0, 255, segmented_image_rescaled);
        //v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(segmented_image_rescaled, file_name_segmented_image2.C_String());
        //cout<<" completed."<<endl;
    }
    //-----------
    
    
    
    
    //----- Create Mask image from segmented melanoma image -----
    {
        cout<<"Creating Mask image...";
        //bdMelanoma melanoma;
        melanoma.CreateMaskImage(segmented_image, mask_image);
        cout<<" completed."<<endl;
    }
    //-----------
    //----- Save the images -----
    {
        cout<<"Saving mask image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(mask_image, file_name_mask_image.C_String());
    }
    //-----------

    
    
    
    //----- Extract images and save histograms for Debug -----
    {
        cout<<"Extracting images for debug...";
        //bdMelanoma melanoma;
        melanoma.Extract_RedImage(white_light_image, extracted_channel_red_image);
        melanoma.Extract_GreenImage(white_light_image, extracted_channel_green_image);
        melanoma.Extract_BlueImage(white_light_image, extracted_channel_blue_image);
        melanoma.Extract_MaxOfRGBImage(white_light_image, extracted_channel_max_RGB_image);
        //melanoma.Extract_DifferenceOfRGB_Image(white_light_image, extracted_channel_max_GB_image);
        
        
//        bdCurveXY histogram_R, histogram_G, histogram_B, histogram_max_GB;
//        bdCurveXY histogram_masked_R, histogram_masked_G, histogram_masked_B, histogram_masked_max_GB;
//        unsigned int min, max;
        
//        extracted_channel_red_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_red_image, histogram_R, min, max);
//        histogram_R.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_red_histogram.C_String());
//        melanoma.Histogram(extracted_channel_red_image, segmented_image, histogram_masked_R, min, max);
//        histogram_masked_R.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_red_masked_histogram.C_String());

//        extracted_channel_green_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_green_image, histogram_G, min, max);
//        histogram_G.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_green_histogram.C_String());
//        melanoma.Histogram(extracted_channel_green_image, segmented_image, histogram_masked_G, min, max);
//        histogram_masked_G.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_green_masked_histogram.C_String());

//        extracted_channel_blue_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_blue_image, histogram_B, min, max);
//        histogram_B.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_blue_histogram.C_String());
//        melanoma.Histogram(extracted_channel_blue_image, segmented_image, histogram_masked_B, min, max);
//        histogram_masked_B.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_blue_masked_histogram.C_String());
        
//        extracted_channel_max_RGB_image.GetVoxelMinimumAndMaximumValue(&min, &max);
//        melanoma.Histogram(extracted_channel_max_RGB_image, histogram_max_GB, min, max);
//        histogram_max_GB.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_max_RGB_histogram.C_String());
//        melanoma.Histogram(extracted_channel_max_RGB_image, segmented_image, histogram_masked_max_GB, min, max);
//        histogram_masked_max_GB.SaveToSingleArrayMatlab_M_File(file_name_extracted_channel_max_RGB_masked_histogram.C_String());
        
        cout<<" completed."<<endl;
    }
    //-----------
    //----- Save the extracted images (Comment out if you do not want any images to be saved) -----
    {
        cout<<"Saving extracted images for debug....";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(extracted_channel_red_image, file_name_extracted_channel_red_image.C_String());
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(extracted_channel_green_image, file_name_extracted_channel_green_image.C_String());
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(extracted_channel_blue_image, file_name_extracted_channel_blue_image.C_String());
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(extracted_channel_max_RGB_image, file_name_extracted_channel_max_RGB_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

    
    
    
    
    
    
    //----- Feature B (Border) -----
    {
        cout<<"Border calculation...";
        //int squared_radius = 16;
        bdArray<int> contrast_along_border;
        //bdMelanoma melanoma;
        //border_image.CopyFrom(segmented_image);
        //melanoma.BorderOfMelanoma(white_light_image, border_image, squared_radius, contrast_along_border);
        //melanoma.BorderOfMelanoma2(segmented_image, contrast_along_border, border_image);
        melanoma.BorderOfMelanoma3(segmented_image, mask_image, contrast_along_border, 24, border_image);// 24 as 5^2-1.
        cout<<" completed."<<endl;
        cout<<" Contrast values along border: ";
        for(unsigned int i=0; i<contrast_along_border.GetNumberOfElements(); i++)
        {
            cout<<" "<<contrast_along_border[i];
        }
        cout<<endl;
    }
    //----------
    //----- Save border image -----
    {
        cout<<"Saving border image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(border_image, file_name_border_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

//
//    
//    
    //----- Feature C (Color) -----
    {
        cout<<"Color image computation...";
        //bdMelanoma melanoma;
        //melanoma.ColorsOfMelanoma(white_light_image, white_light_image, color_mapped_image); //melanoma.ColorsOfMelanoma(white_light_image, segmented_image, color_mapped_image);
        
//        //--- v.2 ---
//        cout<<" Loading charts...";
//        vbdMelanoma v_melanoma;
//        bdImage brown_chart, blue_chart, red_chart;
//        v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/brown9x9.png", brown_chart);
//        v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/blue9x9.png", blue_chart);
//        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/violet_red9x9.png", red_chart);
//        v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/red9x9.png", red_chart);
//        cout<<" done...";
//        melanoma.ColorsOfMelanoma2(white_light_image, white_light_image, red_chart, brown_chart, blue_chart, color_mapped_image);
//        //-----------
        
        //--- v.3,4,5 ---
        cout<<" Loading charts...";
        vbdMelanoma v_melanoma;
        bdImage brown_chart, blue_chart, red_chart;
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/brown_red9x9.png", brown_chart);
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/blue9x9.png", blue_chart);
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/red9x9.png", red_chart);

        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_BROWN_CHART, brown_chart);
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_BLUE_CHART, blue_chart);
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_RED_CHART, red_chart);
        
        bdImage brown_map, blue_map, red_map;
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/brown9x9map.png", brown_map);
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/blue9x9map.png", blue_map);
        //v_melanoma.LoadWhiteLightRGBImageFromPNGFile("/Users/danilobabin/Code/DermScan_CMAKE/ColorCharts/red9x9map.png", red_map);
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_BROWN_MAP, brown_map);
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_BLUE_MAP, blue_map);
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(PATH_TO_RED_MAP, red_map);

        
        cout<<" done...";
        melanoma.ColorsOfMelanoma5(white_light_image, mask_image, red_chart, brown_chart, blue_chart, red_map, brown_map, blue_map, color_chart_mapped_image, color_mapped_image);
        //melanoma.ColorsOfMelanoma5(white_light_image, white_light_image, red_chart, brown_chart, blue_chart, red_map, brown_map, blue_map, color_chart_mapped_image, color_mapped_image);
        //melanoma.ColorsOfMelanoma3(white_light_image, white_light_image, red_chart, brown_chart, blue_chart, red_map, brown_map, blue_map, color_chart_mapped_image, color_mapped_image);
        //melanoma.ColorsOfMelanoma4(white_light_image, white_light_image, red_chart, brown_chart, blue_chart, red_map, brown_map, blue_map, color_chart_mapped_image, color_mapped_image);
        //-----------
        
        cout<<" completed."<<endl;
    }
    //----------
    //----- Save the images -----
    {
        cout<<"Saving color mapped image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(color_mapped_image, file_name_color_mapped_image.C_String());
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(color_chart_mapped_image, file_name_color_chart_mapped_image.C_String());

        //v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(color_mapped_image, file_name_color_mapped_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

//
//
//    //----- Feature D1 (Diameter only) -----
//    {
//        cout<<"Diameter of melanoma computation...";
//        bdMelanoma melanoma;
//        double min_diameter_in_pixels, max_diameter_in_pixels, mean_diameter_in_pixels;
//        melanoma.DiameterOfMelanoma(border_image, &min_diameter_in_pixels, &max_diameter_in_pixels, &mean_diameter_in_pixels);
//        cout<<" completed."<<endl;
//        cout<<" Diameters (in pixels): min="<<min_diameter_in_pixels<<", max="<<max_diameter_in_pixels<<", mean="<<mean_diameter_in_pixels<<endl;
//    }
//    //----------


    
    //----- Feature D2a Inner structure segmentation -----
    {
        cout<<"Segmenting inner bright structure of melanoma...";
        //bdMelanoma melanoma;
        //melanoma.SegmentMelanomaInnerStructure3(white_light_image, segmented_image, inner_bright_structure_image);
        //melanoma.SegmentMelanomaInnerStructure5(white_light_image, segmented_image, inner_bright_structure_image);
        //melanoma.SegmentMelanomaInnerStructure6(white_light_image, segmented_image, inner_bright_structure_image);//Not working yet
        //melanoma.SegmentMelanomaInnerStructure7(white_light_image, segmented_image, inner_bright_structure_image);
        //melanoma.SegmentStreaks(white_light_image, segmented_image, 80, inner_bright_structure_image);
        melanoma.SegmentMelanomaInnerStructure8(white_light_image, mask_image, inner_bright_structure_image);
        cout<<" completed."<<endl;
    }
    //----------
//    //----- Erosion of the inner structures (Comment out if not needed) -----
//    {
//        cout<<"Erosion of segmented inner bright structure...";
//        bdBIP bip;
//        bdImage temp;
//        bdGIP gip;
//        gip.Negative(inner_bright_structure_image, temp);
//        //temp.CopyFrom(inner_bright_structure_image);
//        bip.ErosionByMapping_Circle(temp, segmented_image, inner_bright_structure_image, 1);
//        cout<<" completed."<<endl;
//    }
//    //----------
    //----- Save the images -----
    {
        cout<<"Saving melanoma inner bright structure image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(inner_bright_structure_image, file_name_inner_bright_structure_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

    
    
    //----- Feature D2b Inner dark structure segmentation -----
    {
        cout<<"Segmenting inner dark structure of melanoma...";
        //bdMelanoma melanoma;
        //melanoma.SegmentMelanomaInnerDarkStructure(white_light_image, segmented_image, inner_dark_structure_image);
        //melanoma.SegmentMelanomaInnerDarkStructure2(white_light_image, segmented_image, inner_dark_structure_image);
        //melanoma.SegmentStreaks(white_light_image, segmented_image, 40, inner_dark_structure_image);
        melanoma.SegmentMelanomaInnerDarkStructure3(white_light_image, mask_image, inner_dark_structure_image);
        cout<<" completed."<<endl;
    }
    //----------
//    //----- Erosion of the inner dark structures (Comment out if not needed) -----
//    {
//        cout<<"Erosion of segmented inner dark structure...";
//        bdBIP bip;
//        bdImage temp;
//        temp.CopyFrom(inner_dark_structure_image);
//        bip.ErosionByMapping_Circle(temp, segmented_image, inner_dark_structure_image, 1);
//        cout<<" completed."<<endl;
//    }
//    //----------
    //----- Save the images -----
    {
        cout<<"Saving melanoma inner dark structure image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(inner_dark_structure_image, file_name_inner_dark_structure_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

    
    
    
    
    //----- Feature D2c Inner dark streaks segmentation -----
    {
//        bdImage temp;
//        bdGIP gip;
//        gip.MeanCircleFor2DSlice(white_light_image, segmented_image, 0,0,0,0, 3, temp);
//        white_light_image.CopyFrom(temp);
        
        cout<<"Segmenting inner dark streaks of melanoma...";
        //bdMelanoma melanoma;
        melanoma.SegmentStreaks(white_light_image, mask_image, 30, inner_dark_streaks_image);
        cout<<" completed."<<endl;
    }
    //----------
//    //----- Erosion of the inner dark streaks (Comment out if not needed) -----
//    {
//        cout<<"Erosion of segmented inner dark streaks...";
//        bdImage temp;
//        bdGIP gip;
//        gip.Negative(inner_dark_streaks_image, temp);
//        //temp.CopyFrom(inner_dark_streaks_image);
//        bdBIP bip;
//        bip.ErosionByMapping_Circle(temp, segmented_image, inner_dark_streaks_image, 1);
//        cout<<" completed."<<endl;
//    }
//    //----------
    //----- Save the images -----
    {
        cout<<"Saving melanoma inner dark streaks image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(inner_dark_streaks_image, file_name_inner_dark_streaks_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------
    
    
    //----- Feature D2d Inner structureless area segmentation -----
    {
        cout<<"Segmenting inner structureless area of melanoma...";
        //bdMelanoma melanoma;
        //melanoma.SegmentStructurelessArea(white_light_image, mask_image, 30, inner_structureless_image);
        melanoma.SegmentStructurelessArea2(mask_image,border_image, inner_dark_structure_image, inner_bright_structure_image, inner_dark_streaks_image, inner_structureless_image);
        cout<<" completed."<<endl;
    }
    //----- Save the images -----
    {
        cout<<"Saving melanoma inner structureless area image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(inner_structureless_image, file_name_inner_structureless_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------

    
    
    
    
    //----- Make combined inner structure image -----
    {
        cout<<"Making combined inner structure image...";
        //bdMelanoma melanoma;
        melanoma.MakeCombinedSegmentedImage(mask_image, inner_bright_structure_image, inner_dark_structure_image, inner_dark_streaks_image, inner_combined_structure_image);
        cout<<" completed."<<endl;
    }
    //----------
    //----- Save the images -----
    {
        cout<<"Saving combined inner structure image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(inner_combined_structure_image, file_name_inner_combined_structure_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------
    
    
    
    //----- Feature A (Asymmetry) -----
    {
        cout<<"Asymmetry calculation...";
        unsigned int number_of_symmetry_axes = 8;
        //bdArray<double> min_symmetry_coeff;
        //bdArray<double> max_symmetry_coeff;
        //bdArray<double> mean_symmetry_coeff;
        //bdMelanoma melanoma;
        asymmetry_image.CopyFrom(white_light_image);
        //melanoma.AsymmetryOfMelanoma(asymmetry_image, segmented_image, number_of_symmetry_axes, min_symmetry_coeff, max_symmetry_coeff, mean_symmetry_coeff);
        melanoma.AsymmetryOfMelanoma2(white_light_image, mask_image, color_chart_mapped_image, inner_combined_structure_image, number_of_symmetry_axes, asymmetry_image);

        cout<<" completed."<<endl;
    }
    //----------
    //----- Save asymmetry image -----
    {
        cout<<"Saving asymmetry image...";
        vbdMelanoma v_melanoma;
        v_melanoma.SaveWhiteLightRGBImageToPNGFile_8bit(asymmetry_image, file_name_asymmetry_image.C_String());
        cout<<" completed."<<endl;
    }
    //-----------
    
    
    //----- Quantify D -----
    {
        melanoma.Quantify_D_Feature(mask_image, segmented_image, inner_dark_structure_image, inner_bright_structure_image, inner_dark_streaks_image, inner_structureless_image);
    }


    //----- PRINT RESULTS -----
    {
        melanoma.PrintResults(file_name_white_light_image.C_String());
    }
    
    


    cout<<endl<<"END."<<endl;

    

    
    return 1;
}
