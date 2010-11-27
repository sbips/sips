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
  std::cout << argv[0] << " -fl imgFLSEG -t1 imgT1SEG -out imgOUT [options]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options" << std::endl;
  std::cout << " -fl <imgFLSEG>:  FLAIR segmentation, required " << std::endl;
  std::cout << " -t1 <imgT1SEG>:  T1 segmentation, required" << std::endl;
  std::cout << " -out <imgOUT>:   Output image filename, required" << std::endl;
  std::cout << " -gm <float>:     Gray Matter value on T1 segmentation, default = 4" << std::endl;
  std::cout << " -csf <float>:    CSF value on T1 segmentation, default = 5" << std::endl;
  std::cout << " -min <float>:    Minimum CSF object size, default = 1" << std::endl;
  std::cout << " -max <float>:    Maximum CSF object size, default = 10000" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << " -file extension required (.img/.nii/nii.gz)" << std::endl;
  std::cout << " -CSF and VCSF must be uniquely labeled on T1 segmentation" << std::endl;
  std::cout << " -FLAIR segmentation should have zeros for non-lesion and non-zeros for lesions" << std::endl;
  std::cout << " -false positives should be removed from FLAIR segmentation" << std::endl;
  std::cout << " -any CSF voxels connected in 3D to a non-zero voxel on the FLAIR" << std::endl;
  std::cout << "  is reclaimed as a black hole (2s), and dilated (1 voxel), and if any" << std::endl;
  std::cout << "  GM is found in the dilated area, it is reclaimed as WMH (1s)" << std::endl; 
  std::cout << std::endl;
  exit( 1 );
}

int main( int argc, char* argv[] )
{

  const int Dims = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dims >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  // Set default values
  double accept_value = 1;
  double reject_value = 0;
  double min_object_size = 1;
  double max_object_size = 10000;
  double gm_value = 4;
  double csf_value = 5;

  // Parse command line
  std::string t1segfn;
  std::string flsegfn;
  std::string outfn;

  for ( int i=1; i < argc; ++i )
  {
    std::string cl( argv[i] );
    if ( i+1 < argc )
    {
      if ( cl == "-gm" )
        gm_value = atof( argv[ i+1 ] );
      if ( cl == "-csf" )
        csf_value = atof( argv[ i+1 ] );
      if ( cl == "-min" )
        min_object_size = atof( argv[ i+1 ] );
      if ( cl == "-max" )
        max_object_size = atof( argv[ i+1 ] );
      if ( cl == "-t1" )
        t1segfn = argv[ i+1 ];
      if ( cl == "-fl" )
        flsegfn = argv[ i+1 ];
      if ( cl == "-out" )
        outfn = argv[ i+1 ];
    }     
  }

  std::cout << "GM valule:  " << gm_value << std::endl;
  std::cout << "CSF value:  " << csf_value << std::endl;
  std::cout << "Maximum object size:  " << max_object_size << std::endl;
  std::cout << "Minimum object size:  " << min_object_size << std::endl;

  // Error checking
  if ( argc < 7 )
    Usage( argv );

  if ( outfn.empty() )    { std::cout << "Error:  Must supply -out option"; return 1; }
  if ( flsegfn.empty() )  { std::cout << "Error:  Must supply -fl option";  return 1; }
  if ( t1segfn.empty() )  { std::cout << "Error:  Must supply -t1 option";  return 1; }

  // Read T1 seg and FLAIR seg
  ReaderType::Pointer flSegReader = ReaderType::New();
  ReaderType::Pointer t1SegReader = ReaderType::New();

  flSegReader->SetFileName( flsegfn.c_str() );
  t1SegReader->SetFileName( t1segfn.c_str() );
 
  egitk::ReadImage< PixelType, Dims >( flSegReader );
  egitk::ReadImage< PixelType, Dims >( t1SegReader );
  
  // Isolate CSF on T1 seg
  ImageType::Pointer csf = ImageType::New();
  csf = egitk::CreateEmptyImage< PixelType, Dims >( t1SegReader->GetOutput() );
  csf = egitk::ImageFindValueAndSet< PixelType, Dims >( csf, t1SegReader->GetOutput(), csf_value, 1  ); 

  // Convert FLAIR seg to binary
  ImageType::Pointer fl = ImageType::New();
  fl = flSegReader->GetOutput();
  fl = egitk::ImageFindAndSet< PixelType, Dims >( fl, fl, 1 );
  
  // Keep CSF that is connected to FLAIR seg
  csf = egitk::KeepOnAConnectedToB< PixelType, Dims >(csf, fl, accept_value, reject_value, min_object_size, max_object_size );

  // Dilate remaining CSF 
  int sevals[] = {0,0,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0};
	typedef itk::Neighborhood< PixelType, Dims > NeighborhoodType;
	NeighborhoodType se;
	se.SetRadius( 1 );
	int i=0;
	for ( itk::Neighborhood< PixelType, Dims>::Iterator it = se.Begin(); it != se.End(); ++it, ++i )
		*it = sevals[i];

  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, NeighborhoodType > DilateFilterType;
  DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
	dilateFilter->SetInput( csf );
	dilateFilter->SetKernel( se );
	dilateFilter->SetDilateValue( 1 );
	dilateFilter->Update();

  // Add back any remaining CSF/GM voxels found in the dilated CSF 
  typedef itk::ImageRegionIterator< ImageType > ItType;
  ItType csf_it( dilateFilter->GetOutput(), dilateFilter->GetOutput()->GetLargestPossibleRegion() );
  ItType seg_it( t1SegReader->GetOutput(), t1SegReader->GetOutput()->GetLargestPossibleRegion() );
  ItType fl_it(  fl, fl->GetLargestPossibleRegion() );

  // Iterate over FL seg (and dilated CSF and T1 seg)
  for ( fl_it.GoToBegin(), csf_it.GoToBegin(), seg_it.GoToBegin(); !csf_it.IsAtEnd(); ++csf_it, ++seg_it, ++fl_it )
  {
     if ( csf_it.Get() != 0 )  // if found CSF voxel
     {
       if ( seg_it.Get() == gm_value )  // if GM on T1 seg
         fl_it.Set( 1 ); // set as FLAIR WMH
       if ( seg_it.Get() == csf_value ) // if CSF on T1 seg 
         fl_it.Set( 2 ); // set as FLAIR black hole
     }
     else if ( fl_it.Get() != 0 && seg_it.Get() == csf_value )  // if WMH on FLAIR seg and CSF on T1 seg
       fl_it.Set(2);  // set as FLAIR black hole

  }

  // Write Output
  egitk::WriteImage< PixelType, unsigned char, 3 >( outfn.c_str(), fl );
  
  return 0;
}