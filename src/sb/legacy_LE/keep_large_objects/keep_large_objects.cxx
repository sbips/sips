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

#include "itkImage.h"

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"


int usage( char* argv[] )
{
	std::cout << std::endl << std::endl;
	std::cout << "Usage: " << argv[0] << std::endl;
	std::cout << "       <input.img>" << std::endl;
        std::cout << "       <min_object_size>" << std::endl;
	std::cout << "       <output.img>" << std::endl;
	std::cout << std::endl;
	std::cout << "Notes: " << std::endl;
	std::cout << "  - the .img file suffix must be provided" << std::endl;
	std::cout << "  - any object >= to the min_object_size is kept" << std::endl;
	std::cout << "  - any object < the min_object_size is removed" << std::endl;
	std::cout << "  - objects are defined in 3D" << std::endl;	
	std::cout << std::endl;
  return 1;
}

int main( int argc, char* argv[] )
{

  // Check number of arguments.
  if ( argc != 4 )
    return usage( argv );
  

  // Convenient typedefs.
  typedef double PixelType;
  const int dims = 3;

  typedef itk::Image< PixelType, dims > ImageType;

  ImageType::Pointer lesionmap = ImageType::New();


  // Read input.
  lesionmap = io_utils::ReadItkVol< PixelType, dims >( argv[1] );
  

  // Define constants.  
  const int min_lesion_size = atoi( argv[2] );


  // Remove small lesions.
  //ImageType::Pointer lesionmap_lab = ImageType::New();
  lesionmap = itk_iops::RelabelConnectedComponents< double, double, dims >
                  ( lesionmap, min_lesion_size );

  

  // ( Threshold labelled result and assign to lesionmap, 
  //     in case small lesions were removed )
  lesionmap = itk_iops::ThresholdImageToBinary< PixelType, PixelType, dims >
              ( lesionmap, 1, 0, 2, itk::NumericTraits< PixelType >::max() );



  // Save output.
  std::string outfn;
  outfn = argv[3]; 
  io_utils::WriteItkVol< PixelType, char, dims >
            ( outfn.c_str(), lesionmap );

  
}