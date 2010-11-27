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
#include "itkImageFileReader.h"
#include "../Common/itk_imgio_util_templates.h"

void Usage( char* argv[] )
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " imgIN" << std::endl;
  std::cout << std::endl;
  std::cout << "Datatypes:" << std::endl;
  std::cout << "  1:  unsigned  char" << std::endl;
  std::cout << "  2:  signed    char" << std::endl;
  std::cout << "  3:  unsigned  short" << std::endl;
  std::cout << "  4:  signed    short" << std::endl;
  std::cout << "  5:  unsigned  int" << std::endl;
  std::cout << "  6:  signed    int" << std::endl;
  std::cout << "  7:  unsigned  long" << std::endl;
  std::cout << "  8:  signed    long" << std::endl;
  std::cout << "  9:  float" << std::endl;
  std::cout << "  10: double" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << " -file extension required (.img/.nii/nii.gz)" << std::endl;
  std::cout << " -returns the datatype of the input image" << std::endl;
  std::cout << std::endl;
  exit( 1 );
}

int main( int argc, char* argv[] )
{
  // Error checking
  if ( argc < 2 )
    Usage( argv );

  const int Dims = 3; 
  typedef double PixelType;
  typedef itk::Image< PixelType, Dims >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  // Read image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
 
  if ( !egitk::ReadImage< PixelType, Dims >( reader ) )
    return 1;

  // Report datatype
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned char ) )
    std::cout << "1"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( signed char ) )
	  std::cout << "2"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned short ) )
	  std::cout << "3"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( signed short ) )
	  std::cout << "4"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( signed int ) )
	  std::cout << "5"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned int ) )
	  std::cout << "6"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned long ) )
	  std::cout << "7"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( signed long ) )
	  std::cout << "8"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( float ) )
	  std::cout << "9"; 
  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid( double ) )
	  std::cout << "10"; 

  std::cout << std::endl;
  return 0;
}