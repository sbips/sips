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

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"

#include <string>
#include <vector>

#include "FCM_FLAIR_histo.h"

int usage( char* argv[] )
{
  std::cout << std::endl;
	std::cout << "Usage: " << argv[0] << std::endl;
	std::cout << "       <fl_fc_m1.img>" << std::endl;
	std::cout << "       <fl.img>" << std::endl;
  std::cout << "       <outfn_prefix>" << std::endl;
  std::cout << "       <thresh_low>" << std::endl;
  std::cout << "       <thresh_high>" << std::endl;
  std::cout << std::endl;
  return 1;
}

int main( int argc, char* argv[] )
{
  if ( argc != 6 )
    return usage( argv );
  
  // Set input variables and declare typedefs.
  std::string outfn;  outfn = argv[3]; outfn += "_flhyper.img";
  const double threshLow = atof( argv[4] );  // 70
  const double threshHigh = atof( argv[5] ); // 500 / 1000

  const int dims = 3;

  typedef double PixelType;
  typedef itk::Image< PixelType, dims > ImageType;


  // Load FLAIR image and FLAIR fuzzy cluster results
  //                     (background membership class)
  ImageType::Pointer fcfl = ImageType::New();
  ImageType::Pointer fl = ImageType::New();
  fcfl = io_utils::ReadItkVol< PixelType, dims >( argv[1] );
  fl   = io_utils::ReadItkVol< PixelType, dims >( argv[2] );


  // Threshold fuzzy cluster membership map.
  fcfl = itk_iops::ThresholdImageToBinary< PixelType, PixelType, dims >
         ( fcfl, 1, 0, threshLow, threshHigh );


  
  // Mask out FLAIR image with thresholded fuzzy cluster results.
  ImageType::Pointer flThresh = ImageType::New();
  flThresh = itk_imath::ImageFindValueAndSet< PixelType, PixelType, dims >
       ( fl, fcfl, 0, 0 );

  // Create histograms for each slice and find cutoffs
  std::vector<int> flLesionCutoffs = fl_hist::histo( fl, flThresh );
 

  // Threshold FLAIR slices using cutoffs from histogram
  fl = fl_hist::GetFlairLesions( fl, flLesionCutoffs );



  // Convert to binary.
  fl = itk_iops::ThresholdImageToBinary< PixelType, PixelType, dims >
       ( fl, 1, 0, 1, itk::NumericTraits< PixelType >::max() );



  // Save output binary FLAIR lesion image.
  io_utils::WriteItkVol< PixelType, char, dims > 
            ( outfn.c_str(), fl );



  return 0;

}
