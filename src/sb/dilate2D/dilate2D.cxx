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
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgmorph_util_templates.h"

void Usage( char* argv[] )
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " imgIN imgOUT voxel_value dilate_radius structuring_element" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << " imgIN:               input image" << std::endl;
  std::cout << " imgOUT:              output image" << std::endl;
  std::cout << " voxel_value:         floating point value in the input image for dilation" << std::endl;
  std::cout << " dilate_radius:       integer dilation radius" << std::endl;
  std::cout << " structuring_element: \"1\" for cross \"2\" for ball" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << " -performs 2D dilation on voxels in the input image that have a value" << std::endl;
  std::cout << "  equal to the the voxel_value using a cross (1) or ball (2)" << std::endl;
  std::cout << "  structuring element with radius equal to dilate_radius" << std::endl;
  std::cout << " -values in the input image unaffected by the dilation operation remain" << std::endl;
  std::cout << "  unchanged" << std::endl;
  std::cout << " -file extension required (.img/.nii/nii.gz)" << std::endl;
  std::cout << std::endl;
  exit( 1 );
}

int main( int argc, char* argv[] )
{
  // Error checking
  if ( argc < 6 )
    Usage( argv );

  const int Dims = 3; 
  typedef double PixelType;
  typedef itk::Image< PixelType, Dims >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  float dilate_value = atof( argv[3] );
  float dilate_radius = atof( argv[4] );

  // Read image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );


  if ( !egitk::ReadImage< PixelType, Dims >( reader ) )
    return 1;

  ImageType::Pointer dilate = ImageType::New();
  dilate = egitk::CreateEmptyImage< PixelType, Dims >( reader->GetOutput() );

  dilate = egitk::ImageFindValueAndSet< PixelType, Dims >( dilate, reader->GetOutput(), dilate_value, 1 );


  
  // THIS CODE NEEDS UPDATING (IT IS BUGGY!!)
  if ( atoi( argv[5] ) == 2 )
  {
    typedef itk::BinaryBallStructuringElement<PixelType, 2 > SeType;
    SeType se = egitk::CreateBallSE< PixelType, 2 >( dilate_radius );
    se.SetRadius( dilate_radius ); // necessary because CreateCrossSE returns 3x3 se for all radii
    dilate = egitk::DilateBinaryVol2D< ImageType, ImageType, SeType >( dilate, se );
  }
  else
  {
    typedef itk::BinaryCrossStructuringElement<PixelType, 2 > SeType;
    SeType se = egitk::CreateCrossSE< PixelType, 2 >( dilate_radius );
    se.SetRadius( dilate_radius ); // necessary because CreateCrossSE returns 3x3 se for all radii
    dilate = egitk::DilateBinaryVol2D< ImageType, ImageType, SeType >( dilate, se );
  }
  
 

  typedef itk::ImageRegionIterator< ImageType > ItType;
  ItType dilate_it( dilate, dilate->GetLargestPossibleRegion() );
  ItType input_it( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() ); 


  
  for ( input_it.GoToBegin(), dilate_it.GoToBegin(); !input_it.IsAtEnd(); ++dilate_it, ++input_it )
  {
    if ( dilate_it.Get() != 0 )
      input_it.Set( dilate_value );
  }

  egitk::WriteImageSamePixelType< PixelType, Dims >( argv[2], reader, reader->GetOutput() );
  
  return 0;
}