/**************************************************************************
 Methods loading/saving and visual manipulation of images of melanoma.

 Author: Danilo Babin
 File name: "vbdMelanoma.h"
 ***************************************************************************/



#ifndef VBD_MELANOMA_DEF
    #define VBD_MELANOMA_DEF



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


#include "bdImage.h"
#include "bdObject.h"
#include "bdBIP.h"
#include "bdGIP.h"





///  Algorithms for loading, saving and visual manipulation of images of melanoma.

class vbdMelanoma
{
public:
    
    /// Load grayscale image from a PNG file.
    int LoadGrayScaleImageFromPNGFile(const char *file_name, bdImage &output_loaded_image);
    
    /// Load white light RGB image from a PNG file. The image will be decomposed into R, G, B bands, each of which will be stored at a different time (t) index.
    int LoadWhiteLightRGBImageFromPNGFile(const char *file_name, bdImage &output_loaded_image);
    
    
    /// Save white light RGB image to a PNG file. The input is a 2d+Time image where each time instance is also image of different channel (channels from a RGB image stored at a different time (t) index).
    /// Note: this method cas sace 2D grayscale images, too. Performs rescaling to 8-bit range before saving!
    int SaveWhiteLightRGBImageToPNGFile_8bit(bdImage &decomposed_RGB_image, const char *output_file_name);
    
    /// Save white light RGB image to a PNG file. The input is a 2d+Time image where each time instance is also image of different channel (channels from a RGB image stored at a different time (t) index).
    /// Note: this method cas sace 2D grayscale images too. rescaling is not performed!
    int SaveWhiteLightRGBImageToPNGFile_8bit_NO_RESCALE(bdImage &decomposed_RGB_image, const char *output_file_name);

    
    
//    /// Save gray scale image to a PNG file. The input is a 2D image.
//    int SaveGrayScaleImageFromPNGFile(bdImage &gray_image, const char *output_file_name);


    
};

#endif