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
#include "itkImageRegionIterator.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkScalarConnectedComponentImageFilter.h"

#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgio_util_templates.h"

void Usage( char* argv[] )
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " -lseg imgLESIONSEG -t1seg imgT1SEG -out imgOUT [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options" << std::endl;
  std::cout << " -lseg <imgLESIONSEG>: LESION segmentation, required " << std::endl;
  std::cout << " -t1seg <imgT1SEG>:    T1 segmentation, required" << std::endl;
  std::cout << " -out <imgOUT>:        Output image filename, required" << std::endl;
  std::cout << " -bh <float>:          Black hole value to be assigned, default = 2" << std::endl;
  std::cout << " -csf <float>:         CSF value on T1 segmentation, default = 5" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << " -file extension required (.img/.nii/nii.gz)" << std::endl;
  std::cout << " -CSF and VCSF must be uniquely labeled on T1 segmentation" << std::endl;
  std::cout << " -any non-zero value on the -lseg image that corresponds to csf voxel on the" << std::endl;
  std::cout << "  -t1seg image will be relabelled as a black hole (i.e. assigned the -bh value)" << std::endl;
  std::cout << std::endl;
  exit( 1 );
}

int main( int argc, char* argv[] )
{

  const int Dims = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dims >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

 
  // Error checking
  if ( argc < 6 )
    Usage( argv );

  // Set default values
  double bh_value = 2;
  double csf_value = 5;

  // Parse command line
  std::string t1segfn;
  std::string lsegfn;
  std::string outfn;

  for ( int i=1; i < argc; ++i )
  {
    std::string cl( argv[i] );
    if ( i+1 < argc )
    {
      if ( cl == "-bh" )
        bh_value = atof( argv[ i+1 ] );
      if ( cl == "-csf" )
        csf_value = atof( argv[ i+1 ] );
      if ( cl == "-t1seg" )
        t1segfn = argv[ i+1 ];
      if ( cl == "-lseg" )
        lsegfn = argv[ i+1 ];
      if ( cl == "-out" )
        outfn = argv[ i+1 ];
    }     
  }

  std::cout << "Black hole valule:  " << bh_value << std::endl;
  std::cout << "CSF value:  " << csf_value << std::endl;

 
  // Read T1 seg and FLAIR seg
  ReaderType::Pointer wmhSegReader = ReaderType::New();
  ReaderType::Pointer t1SegReader = ReaderType::New();

  wmhSegReader->SetFileName( lsegfn );
  t1SegReader->SetFileName( t1segfn );
 
  if ( !egitk::ReadImage< PixelType, Dims >( wmhSegReader ) ) return 1;;
  if ( !egitk::ReadImage< PixelType, Dims >( t1SegReader ) ) return 1;
  
  typedef itk::ImageRegionIterator< ImageType > ItType;
  ItType wmh_it( wmhSegReader->GetOutput(), wmhSegReader->GetOutput()->GetLargestPossibleRegion() );
  ItType t1_it( t1SegReader->GetOutput(), t1SegReader->GetOutput()->GetLargestPossibleRegion() );
 
  //  if any lesion voxel segments as csf on the T1 segmentation, relabel as black hole
  for ( wmh_it.GoToBegin(), t1_it.GoToBegin(); !wmh_it.IsAtEnd(); ++wmh_it, ++t1_it )
  {
     if ( t1_it.Get() == csf_value && wmh_it.Get() != 0 )  
       wmh_it.Set( bh_value );
  }

  // Write Output
  egitk::WriteImageSamePixelType< PixelType, 3 >( outfn.c_str(), wmhSegReader, wmhSegReader->GetOutput() );
  
  return 0;
}