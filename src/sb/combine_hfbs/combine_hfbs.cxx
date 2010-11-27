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

#include <string>
#include <iostream>

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"

void Usage( char* argv[] )
{
	std::cout << "Usage:" << std::endl;
	std::cout << "  " << argv[0] << " imgHFBauto imgHFBtemplate imgHFBcombined" << std::endl;
	std::cout << std::endl;
	std::cout << "Notes:" << std::endl;
	std::cout << "  -file extension required (.img/.nii/.nii.gz)" << std::endl;
 	std::cout << std::endl;
	exit(1);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] )
{
	
  if ( argc < 4 )
    Usage( argv );

	const int Dim = 3;
	typedef double PixelType;
	typedef itk::Image< PixelType, Dim > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

  ReaderType::Pointer ahfb = ReaderType::New();
  ahfb->SetFileName( argv[1] );

	ReaderType::Pointer thfb = ReaderType::New();
	thfb->SetFileName( argv[2] );

  if ( !egitk::ReadImage< PixelType, Dim >( ahfb ) || !egitk::ReadImage< PixelType, Dim >( thfb ) )
    exit(1);

  ImageType::SizeType ahfbSize;
  ahfbSize = ahfb->GetOutput()->GetLargestPossibleRegion().GetSize();

  ImageType::SizeType thfbSize;
  thfbSize = thfb->GetOutput()->GetLargestPossibleRegion().GetSize();

  for ( int i=0; i < Dim; ++ i )
  {
    if ( ahfbSize[i] != thfbSize[i] )
    {
      std::cout << "Error: Unequal voxel sizes (Dimension " << i << ": " << ahfbSize[i] <<" and " << thfbSize[i] << ")" << std::endl;
      exit(1);
    }
  }

  
  typedef itk::ImageRegionIterator< ImageType > ItType;

  ItType ahfb_it( ahfb->GetOutput(), ahfb->GetOutput()->GetLargestPossibleRegion() );
  ItType thfb_it( thfb->GetOutput(), thfb->GetOutput()->GetLargestPossibleRegion() );

  for ( ahfb_it.GoToBegin(), thfb_it.GoToBegin(); !ahfb_it.IsAtEnd(); ++ahfb_it, ++thfb_it )
  {
    if (  thfb_it.Get() > 0 )
      ahfb_it.Set( 8 );
  }

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  
  egitk::WriteImageSamePixelType< PixelType, Dim >( argv[3], ahfb, ahfb->GetOutput() );

	return 0;

}

