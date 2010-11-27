/*=========================================================================
  
  Copyright 2002, 2006-2010, Dr. Sandra Black
  Linda C. Campbell Cognitive Neurology Unit
  Sunnybrook Health Sciences Center
  
  This file is part of the Sunnybrook Image Software Processing (SIPS) package

  SIPS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

//Last Modified: 10/12/2008

#include "itkImage.h"

#include <iostream>
#include <string>

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"
#include "../Common/itk_imgmorph_util_templates.h"


int usage( char* argv[] )
{
	std::cout << std::endl;
	std::cout << "Usage: " << argv[0] << std::endl;
	std::cout << "       <seg.img>" << std::endl;
	std::cout << "       <t1csf.img>" << std::endl;
	std::cout << "       <outfn_prefix>" << std::endl;
  std::cout << "       <t1csf_thresh>" << std::endl;
	std::cout << std::endl;
	return 1;
}

int main( int argc, char* argv[] )
{

	if ( argc != 5 )
		return usage( argv );

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          // typedef all required image types //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	//    ImageType and ImageType2DType are used for all internal 3D and 2D calculations
	//    ShortImageType and UcharImageType are used just for saving the output
	typedef double PixelType;
  const int dims = 3;

  typedef itk::Image< PixelType, 3 > ImageType;
	typedef itk::Image< PixelType, 2 > ImageType2D;
	typedef itk::Image< short, 3 > ShortImageType;
	typedef itk::Image< unsigned char, 3 > UcharImageType;

	ImageType::Pointer t1csf = ImageType::New();
	ImageType::Pointer seg = ImageType::New();


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	           // Load the 3 input images //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	seg       = io_utils::ReadItkVol< PixelType, dims >( argv[1] );
	t1csf     = io_utils::ReadItkVol< PixelType, dims >( argv[2] );

  const PixelType T1CSF_THRESH = atof( argv[4] );


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	            // Extract vcsf from seg //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ImageType::Pointer vcsf = ImageType::New();

  vcsf = itk_iops::ThresholdImageToBinary< PixelType, PixelType, 3>
                   ( seg, 1, 0, 7, 7 );
	            

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	          // Dilate vcsf and fill holes //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // 2D dilation with Ball SE (radius 1) 
  vcsf = itk_imorph::DilateBinaryVol2D< ImageType, ImageType >
         ( vcsf, itk_imorph::CreateCrossSE< PixelType, 2 > ( 1 ) );

 
  // 2D Morphological closing operation with Ball SE (radius 2)
  // ** closing must be 2D ** //
 
  vcsf = itk_imorph::DilateBinaryVol2D< ImageType, ImageType >
         ( vcsf, itk_imorph::CreateBallSE< PixelType, 2 >( 2 ) ); 
 // vcsf = itk_imorph::ErodeBinaryVol2D< ImageType, ImageType >
 //        ( vcsf, itk_imorph::CreateBallSE< PixelType, 2 >( 2 ) );
  
	// Fill in any holes (2D)
  
  vcsf = itk_iops::FillHolesInBinaryVol2D< PixelType, PixelType >
         ( vcsf );
  

  

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	             // Find perimeter of brain //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ImageType::Pointer perim = ImageType::New();
  perim = itk_imorph::GetObjectPerimeter< PixelType, 3> 
          ( &*seg );
  

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	                 // Extract scsf //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  ImageType::Pointer scsf = ImageType::New();

  // Threshold csf fuzzy cluster membership map
  scsf = itk_iops::ThresholdImageToBinary< PixelType, PixelType, 3 >
         ( t1csf, 1, 0, T1CSF_THRESH, itk::NumericTraits< PixelType >::max() ); 
  
  // Add perimeter
  scsf = itk_imath::ImageAddImage< PixelType, PixelType, PixelType, 3>
         ( perim, scsf );
  
  scsf = itk_iops::ThresholdImageToBinary< PixelType, PixelType, 3 >
         ( scsf, 1, 0, 1, 2 );
  
  // Remove vcsf ( ?? might want to dilate vcsf even more just for this step... ?? )
  scsf = itk_imath::ImageAddImage< PixelType, PixelType, PixelType, 3 >
         ( scsf, scsf ); // scsf voxels = 2
  scsf = itk_imath::ImageAddImage< PixelType, PixelType, PixelType, dims >
         ( scsf, vcsf ); // vcsf voxels = 1, add, good scsf voxels will still equal 2
  scsf = itk_iops::ThresholdImageToBinary< PixelType, PixelType, dims >
         ( scsf, 1, 0, 2, 2 );
 
           
  // Remove anything not connected to perimeter (using 3D connectivity)
    // relabelConnectedComp labels the largest objects starting from 1,
    // therefore, thresholding the output at 1 should return all csf
    // connected to perimeter, provided that it makes up the largest 
    //  csf object in the labelled image
  scsf = itk_iops::RelabelConnectedComponents< PixelType, PixelType, 3 >
         ( scsf, 1 );
  scsf = itk_iops::ThresholdImageToBinary< PixelType, PixelType, 3 >
         ( scsf, 1, 0, 2, 2 );
 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    // Estimate scsfgm //
 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ImageType::Pointer scsfgm = ImageType::New();

  // 3D dilation with Ball SE ( radius = 1 ) from scsf
  scsfgm = itk_imorph::DilateBinaryImage< PixelType, PixelType, dims > 
          ( scsf, itk_imorph::CreateBallSE< PixelType, dims >( 1 ) );
 
  // apply 2D median filter on dilation result
  scsfgm = itk_iops::MedianFilterVol2D< PixelType, PixelType  >
          ( scsfgm, 1);

  // 2D dilation with Cross SE ( radius = 1 ) 
  scsfgm = itk_imorph::DilateBinaryVol2D< ImageType, ImageType >
          ( scsfgm, itk_imorph::CreateCrossSE< PixelType, 2>( 1 ) );


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create image with scsfgm and vcsf  //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //vcsf = itk_imorph::ErodeBinaryVol2D< ImageType, ImageType >
  //       ( vcsf, itk_imorph::CreateCrossSE< PixelType, 2 >( 2 ) );

  // Create output image with lesions overlayed onto the estimated scsfg and vcsf compartments.
  ImageType::Pointer csfgm = ImageType::New();
  // Brain voxels set to 2. 
  csfgm = itk_imath::ImageFindAndSet<PixelType, PixelType, dims>
                     ( seg, seg, 2 );

  // SCSFGM exclude voxels set to 1.
  csfgm = itk_imath::ImageFindAndSet<PixelType, PixelType, dims>
                     ( csfgm, scsfgm, 1 );
  // Dilated vcsf voxels set to 7.
  csfgm = itk_imath::ImageFindAndSet<PixelType, PixelType, dims>
                     ( csfgm, vcsf, 7 );
  csfgm = itk_imath::ImageFindValueAndSet<PixelType, PixelType, dims>
                     ( csfgm, seg, 0, 0 );


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    // Write output //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  std::string outfn;	outfn = argv[3]; outfn += "_csfgm.img";
  io_utils::WriteItkVol< PixelType, char, dims >
            ( outfn.c_str(), csfgm );

	return 0;
}

