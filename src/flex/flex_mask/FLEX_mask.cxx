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

#include <iostream>


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageIOBase.h"

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"

int usage( char* argv[] )
{
	std::cout << std::endl << std::endl;
	std::cout << "Usage: " << argv[0] << std::endl;
	std::cout << "       <input.img>" << std::endl;
    std::cout << "       <mask.img>" << std::endl; 
	std::cout << "       <output.img>" << std::endl;
	std::cout << std::endl;
	std::cout << "Notes: " << std::endl;
	std::cout << "  - the .img file suffix must be provided" << std::endl;
	std::cout << std::endl;
 
	return 1;
}

int main( int argc, char* argv[] )
{

  // Check arguments
	if ( argc != 4 )
		return usage(argv);

  // typedefs
  typedef double PixelType;
  const int dims = 3;

  typedef itk::Image< PixelType, dims > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;


  itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
 
    ReaderType::Pointer reader = ReaderType::New();

 
  // Read mask image
  ImageType::Pointer mask = ImageType::New();
  reader->SetImageIO( io );
  reader->SetFileName( argv[2] );
  reader->Update();
  mask = reader->GetOutput();
  mask->DisconnectPipeline();

  // Read input image ( read input after mask so input type determines the output type) 

  ImageType::Pointer input = ImageType::New();
  reader->SetImageIO( io );
  reader->SetFileName( argv[1] ); 
  reader->Update();
  input = reader->GetOutput();
  input->DisconnectPipeline();

  // Do masking
  input = itk_imath::ImageMask< PixelType, PixelType, dims >( input, mask );
  
  // Save output.
  std::string outfn;
  outfn = argv[3]; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned char ) )
  {
    io_utils::WriteItkVol< double, unsigned char, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( char ) )
  {   
    io_utils::WriteItkVol< double, char, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned short ) )
  {   
    io_utils::WriteItkVol< double, unsigned short, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid (  short ) )
  {   
    io_utils::WriteItkVol< double, short, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned long ) )
  {
    io_utils::WriteItkVol< double, unsigned long, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( long ) )
  {   
    io_utils::WriteItkVol< double, long, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( float ) )
  {   
    io_utils::WriteItkVol< double, float, 3 > ( argv[3], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( double ) )
  {   
   io_utils::WriteItkVol< double, double, 3 > ( argv[3], input );
  }
  else 

	return 0;
}
