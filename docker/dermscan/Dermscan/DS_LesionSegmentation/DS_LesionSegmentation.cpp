/**************************************************************************
 DermScan lesion segmentation feature quantification.

 Author: Danilo Babin
 File name: "DS_LesionSegmentation.cpp"
 ***************************************************************************/




#include <iostream>
#include <sstream>

#include "bdString.h"
#include "bdMelanoma.h"
#include "vbdMelanoma.h"




int main(int argc, char** argv)
{
    bdString file_name_white_light_image;
    
////----- Uncomment this to pre-set the file name (comment out the part for command line)-----
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000002/ISIC_0000002.png");
    //file_name_white_light_image.Assign("/Users/danilobabin/-DIP_IMAGES/-MELANOMA/Images_test_RGB/ISIC_0000001/ISIC_0000001.png");
    
    
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
////---------------------------------------------------
    
    
    
    
//----- Uncomment this to use from command line -----
    if ( argc != 2 )
    {
        cout<<"usage: "<< argv[0] <<" <filename>\n";
        return 0;
    }
    else
    {
        file_name_white_light_image.Assign(argv[1]);
    }
//----------------------------------------------------
    
    
    //----- Data used for various steps (global data) -----
    bdImage white_light_image;
    bdImage segmented_image;
    //bdImage segmented_image2;
    //----------
    
    
    
    //----- File names for writing images based on input file name -----
    bdString file_dir, file_name_root, file_extension;
    file_name_white_light_image.ExtractFileNamePathParts(file_dir,file_name_root,file_extension);
    bdString file_name_segmented_image;
    file_name_segmented_image.Assign(file_dir);
    file_name_segmented_image.Append(file_name_root);
    file_name_segmented_image.Append("_segmented.");
    file_name_segmented_image.Append(file_extension);
    //bdString file_name_segmented_image2;
    //file_name_segmented_image2.Assign(file_dir);
    //file_name_segmented_image2.Append(file_name_root);
    //file_name_segmented_image2.Append("_segmented_0002.");
    //file_name_segmented_image2.Append(file_extension);


    
    
    //----- Load the white light image -----
    {
        cout<<"Loading image...";
        vbdMelanoma v_melanoma;
        v_melanoma.LoadWhiteLightRGBImageFromPNGFile(file_name_white_light_image.C_String(), white_light_image);
        cout<<" completed."<<endl;
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

    
    
    //    //----- Perform equalization (Comment this part out if you don't want equalization) -----
    //    {
    //        bdImage temp;
    //        temp.CopyFrom(white_light_image);
    //        bdMelanoma melanoma;
    //        melanoma.EqualizeBrightnessOfChannels(temp, white_light_image);
    //    }
    //    //-----------
    
    
    
    //----- Segment melanoma -----
    {
        cout<<"Segmenting melanoma...";
        bdMelanoma melanoma;
        //melanoma.SegmentMelanoma_Thresholding(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding2(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding3(white_light_image, segmented_image);
        //melanoma.SegmentMelanoma_Thresholding4(white_light_image, segmented_image, segmented_image2);
        //melanoma.SegmentMelanoma_Thresholding5(white_light_image, segmented_image);
        ////melanoma.SegmentMelanoma_Thresholding6(white_light_image, segmented_image);
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
    
    


    cout<<endl<<"END."<<endl;

    

    
    return 1;
}
