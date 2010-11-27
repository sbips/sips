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

#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgio_util_templates.h"

void Usage( char* argv[] )
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " -t1seg imgT1SEG -lseg imgLESIONSEG -out imgOUT [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options" << std::endl;
  std::cout << " -t1seg <img>:    T1 segmentation image, required" << std::endl;
  std::cout << " -lseg <img>:     lesion segmentation image, required" << std::endl;
  std::cout << " -out <img>:      output image, required" << std::endl;
  std::cout << " -vcsf   <float>: value of VCSF on T1 segmentation, default = 7" << std::endl;
  std::cout << " -offset <float>: offset value for relabeled WMH, default = 2" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << " -file extension required (.img/.nii/nii.gz)" << std::endl;
  std::cout << " -any lesion on the -lseg image connected in 3D to VCSF on the" << std::endl;
  std::cout << "  on the -t1seg image is relabelled as periventricular by adding the" << std::endl;
  std::cout << "  offset value to each connected voxel, for example, using an offset of 2 with"  << std::endl;
  std::cout << "  a -fl image containing 1s for WMH and 2s for black holes would produce:"  << std::endl << std::endl;
  std::cout << "    LESION                     BEFORE   AFTER" << std::endl;
  std::cout << "    Deep WMH                     1        1    " << std::endl;
  std::cout << "    Deep Black Hole              2        2    " << std::endl;
  std::cout << "    Periventricular WMH          1        3    " << std::endl;
  std::cout << "    Periventriculalr Black Hole  2        4    " << std::endl;
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
  if ( argc < 7 )
    Usage( argv );

  double vcsf_value = 7;
  double offset_value  = 2;
  std::string flfn;
  std::string t1fn;
  std::string outfn;

  for ( int i=1; i<argc; ++i )
	{
		std::string cl( argv[i] );

    if ( i < argc - 1 )
    {
		  if ( cl == "-lseg" )
        flfn = argv[i+1];
      if ( cl == "-t1seg" )
        t1fn = argv[i+1];
      if ( cl == "-out" )
        outfn = argv[i+1];
      if ( cl == "-vcsf" )
        vcsf_value = atof( argv[i+1] );
      if ( cl == "-offset_value" )
        offset_value = atof( argv[i+1] );
    }
  }
  std::cout << "Offset value: " << offset_value << std::endl;
  std::cout << "VCSF value:   " << vcsf_value << std::endl;

  if ( flfn.empty() )  { std::cout << "Error: must supply -fl option" << std::endl; return 1; }
  if ( t1fn.empty() )  { std::cout << "Error: must supply -t1 option" << std::endl; return 1; }
  if ( outfn.empty() ) { std::cout << "Error must supply -out option" << std::endl; return 1; }
  
  // Read T1 seg and FLAIR seg
  ReaderType::Pointer flSegReader = ReaderType::New();
  ReaderType::Pointer t1SegReader = ReaderType::New();

  t1SegReader->SetFileName( t1fn );
  flSegReader->SetFileName( flfn );
  
  if ( !egitk::ReadImage< PixelType, Dims >( flSegReader ) || !egitk::ReadImage< PixelType, Dims >( t1SegReader ) )
    return 1;



  // Isolate CSF on T1 seg
  ImageType::Pointer vcsf = ImageType::New();
  vcsf = egitk::CreateEmptyImage< PixelType, Dims >( t1SegReader->GetOutput() );
  vcsf = egitk::ImageFindValueAndSet< PixelType, Dims >( vcsf, t1SegReader->GetOutput(), vcsf_value, 1  ); 

  // Convert FLAIR seg to binary
  ImageType::Pointer fl = ImageType::New();
  fl = flSegReader->GetOutput();
  fl = egitk::ImageFindAndSet< PixelType, Dims >( fl, fl, 1 );
  
  // Keep WMH that is connected to VCSF
  ImageType::Pointer vwmh = ImageType::New();
  vwmh = egitk::KeepOnAConnectedToB< PixelType, Dims >(fl, vcsf, 1, 0, 1, itk::NumericTraits< double >::max() );

  typedef itk::ImageRegionIterator< ImageType > ItType;
  ItType vwmh_it( vwmh, vwmh->GetLargestPossibleRegion() );
  ItType fl_it( flSegReader->GetOutput(), flSegReader->GetOutput()->GetLargestPossibleRegion() ); 

  for ( fl_it.GoToBegin(), vwmh_it.GoToBegin(); !fl_it.IsAtEnd(); ++fl_it, ++vwmh_it )
  {
    if ( vwmh_it.Get() != 0 )
      fl_it.Set( fl_it.Get() + offset_value );
  }

  // Write Output
  egitk::WriteImageSamePixelType< PixelType, 3 >( outfn.c_str(), flSegReader, flSegReader->GetOutput() );
  
  return 0;
}