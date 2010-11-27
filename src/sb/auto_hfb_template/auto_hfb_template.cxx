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
#include <numeric>
#include <math.h>

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
//==================================================================
// COMMAND LINE FUNCTIONS
//==================================================================
//-----------------------------------------------------------------------------
class Parameters
{
public:
  Parameters():
      n_pd_threshold(-1),
      n_t2_threshold(-1),
      pd_starting_right_cutoff(-1),
      t2_starting_right_cutoff(-1),
      cutoff_value(50)
      {};

  bool IsInitialized()
  {
    bool initialized = true;
    if ( n_pd_threshold == -1 )
    {
      std::cout << "Error: must supply -npdt option" << std::endl;
      initialized = false;
    }
    if ( n_t2_threshold == -1 )
    {
      std::cout << "Error: must supply -nt2t option" << std::endl;
      initialized = false;
    }
    if ( t2_filename.length() == 0 )
    {
      std::cout << "Error: must supply -t2 option" << std::endl;
      initialized = false;
    }
    if ( pd_filename.length() == 0  )
    {
      std::cout << "Error: must supply -pd option" << std::endl;
      initialized = false;
    }
    if ( avg_filename.length() == 0  )
    {
      std::cout << "Error: must supply -avg option" << std::endl;
      initialized = false;
    }
    if ( out_filename.length() == 0  )
    {
      std::cout << "Error: must supply -out option" << std::endl;
      initialized = false;
    }
    return initialized;
  }

  double      n_pd_threshold;
  double      n_t2_threshold;
  int         pd_starting_right_cutoff;
  int         t2_starting_right_cutoff;
  int         cutoff_value;
  std::string pd_filename;
  std::string t2_filename;
  std::string avg_filename;
  std::string out_filename;
};
//-----------------------------------------------------------------------------
void GetCommandLineParameters( int c, char* argv[], Parameters& parameters )
{
	for ( int i=1; i < c; ++i )
	{
		std::string cl( argv[i] );
    
		if ( cl == "-npdt" )
    { 
      ++i; if ( i >= c ) break; 
      parameters.n_pd_threshold = atof( argv[i] );
      std::cout << "normalized pd threshold = " << parameters.n_pd_threshold << std::endl;
    }
		if ( cl == "-nt2t" )
    { 
      ++i; if ( i >= c ) break;
      parameters.n_t2_threshold = atof( argv[i] );
      std::cout << "normalized t2 threshold = " << parameters.n_t2_threshold << std::endl;
    }
		if ( cl == "-t2" )
    { 
      ++i; if ( i >= c ) break;
      parameters.t2_filename = argv[i];
      std::cout << "T2 filename = " << parameters.t2_filename << std::endl;
    }
		if ( cl == "-pd" )
    { 
      ++i; if ( i >= c ) break;
      parameters.pd_filename = argv[i];
      std::cout << "PD filename = " << parameters.pd_filename << std::endl;
    }
    if ( cl == "-avg" )
    { 
      ++i; if ( i >= c ) break;
      parameters.avg_filename = argv[i];
      std::cout << "Average template filename = " << parameters.avg_filename << std::endl;
    }
		if ( cl == "-out" )
    { 
      ++i; if ( i >= c ) break;
      parameters.out_filename = argv[i];
      std::cout << "Output filename = " << parameters.out_filename << std::endl;
    }
 		if ( cl == "-cut" )
    { 
      ++i; if ( i >= c ) break;
      parameters.cutoff_value = atoi( argv[i] );
      std::cout << "Cutoff value = " << parameters.cutoff_value << std::endl;
    }   
    
  }// end for loop
} // end GetCommandLineParameters(...)
//-----------------------------------------------------------------------------
void Usage( char* argv[] )
{
	std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " -t2 <T2> -pd <PD> -avg <TEMPLATE> -nt2t <t2_threshold> -npdt <pd_threshold> [options]" << std::endl;
	std::cout << std::endl;
	std::cout << " -t2 <vol>:    (T2 image, usually in T1 acquisition space, required)" << std::endl;
	std::cout << " -pd <vol>:    (PD image, usually in T1 acqisition space, required)" << std::endl;
	std::cout << " -avg <vol>:    (coregistered HfB template, 1s for brain, 0s for bkgd, required)" << std::endl;
	std::cout << " -nt2t <float>: (normalized T2 threshold, between 0 and 1, required)" << std::endl;
	std::cout << " -npdt <float>: (normalized PD threshold, between 0 and 1, required)" << std::endl;
	std::cout << " -cut <int>:    (cutoff value, optional, default = 50)" << std::endl;
	std::cout << std::endl;
	std::cout << "Notes:" << std::endl;
	std::cout << " -good cut value would be the number of voxels expected in the tail" << std::endl << "  of the PD and T2 histograms" << std::endl;
	std::cout << " -file extensions are required (.img/.nii/.nii.gz)" << std::endl;
	std::cout << std::endl;
	exit(1);
}

//==================================================================
// FUNCTIONS
//==================================================================

template < typename TPixelType, int TDim >
int GetCutoff( typename itk::Image< TPixelType, TDim >::Pointer img, int cutoff_count )
{
  // This function will find a good value with which to normalize the PD/T2 images

  // Calculate image maximum
  double max_value = egitk::ImageMax< TPixelType, TDim >( img );
  std::cout << "Image Maximum = " << max_value << std::endl;
   
  // Calculate histogram
  std::vector<double> histogram( (int)max_value + 1, 0 ); 
  
  typename itk::ImageRegionIterator< itk::Image< TPixelType, TDim > > img_it( img, img->GetLargestPossibleRegion() );
  for ( img_it.GoToBegin(); !img_it.IsAtEnd(); ++img_it )
    ++histogram[ img_it.Get() ];

  // Smooth histogram
  std::vector<double> histogram_smoothed( (int)max_value + 1, 0 );
  for ( int i = 4; i < (int)max_value + 1 - 3; ++i )
  {
    std::vector<double>::iterator begin = histogram.begin() + i - 3;
    std::vector<double>::iterator end   = histogram.begin() + i + 3 + 1; // add 1 for end position
    histogram_smoothed[i] = floor( (std::accumulate( begin, end, 0.0 ) / 7) + 0.5 );
    //histogram[i] = floor( (std::accumulate( begin, end, 0.0 ) / 7) + 0.5 );
  }
  histogram.erase(histogram.begin(), histogram.end());

  // start from histogram tail, look for first bin containing more than the cutoff number of voxels,
  // call this cutoff_value_1
  int cutoff_value_1 = -1;
  for ( int i = histogram_smoothed.size() - 1; i >= 0; --i )
  {
    if ( histogram_smoothed[i] > cutoff_count )
    {
      std::cout << "Starting normalization value = " << i << " (" << (float)1/i << ")" << std::endl;
      cutoff_value_1 = i;
      break;
    }
  }

  // shift cutoff_value_1 away from the tail a small amount (25%)
  // start from the shifted cutoff_value_1, look for first bin containing less than the cutoff number of voxels,
  // call this cutoff_value_2
  int cutoff_value_2 = int( cutoff_value_1 - cutoff_value_1 * 0.25 );
  for ( int i = cutoff_value_2; i < histogram_smoothed.size(); ++i )
  {
    if ( histogram_smoothed[i] <= cutoff_count )
    {
      std::cout << "Final normalization value = " << i  << " (" << (float)1/i << ")" << std::endl;
      cutoff_value_2 = i;
      break;
    }
  }

  return cutoff_value_2;
}

//===================================================================================
// MAIN
//===================================================================================
int main( int argc, char *argv[] )
{
	//=========================
  // Verify Input
  //=========================
  if ( argc == 1 )
    Usage( argv );

  Parameters parameters;
  GetCommandLineParameters( argc, argv, parameters );
  
  if ( !parameters.IsInitialized() )
    return 1;

	//=========================
  // Read Input Images
  //=========================
	const int Dim = 3;
	typedef double PixelType;
	typedef itk::Image< PixelType, Dim > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

  ReaderType::Pointer t2 = ReaderType::New();
  t2->SetFileName( parameters.t2_filename );

	ReaderType::Pointer pd = ReaderType::New();
  pd->SetFileName( parameters.pd_filename );

  ReaderType::Pointer avg = ReaderType::New();
  avg->SetFileName( parameters.avg_filename );

  if ( !egitk::ReadImage< PixelType, Dim >( t2  ) || 
       !egitk::ReadImage< PixelType, Dim >( pd  ) ||
       !egitk::ReadImage< PixelType, Dim >( avg ) )
    exit(1);

	//=========================
  // Calculate HFB mask
  //=========================

  // Find normalization factors cx (T2) and cy (PD)
  std::cout << std::endl << "Searching for T2 normalization value..." << std::endl;
  int cx = GetCutoff< PixelType, Dim>( t2->GetOutput(), parameters.cutoff_value ); 

  std::cout << std::endl << "Searching for PD normalization value..." << std::endl;
  int cy = GetCutoff< PixelType, Dim>( pd->GetOutput(), parameters.cutoff_value );

  // Normalize images
  itk::ImageRegionIterator< ImageType > t2_it( t2->GetOutput(), t2->GetOutput()->GetLargestPossibleRegion() );
  for ( t2_it.GoToBegin(); !t2_it.IsAtEnd(); ++t2_it )
  {
    t2_it.Set( t2_it.Get() / cx );
  }

  itk::ImageRegionIterator< ImageType > pd_it( pd->GetOutput(), pd->GetOutput()->GetLargestPossibleRegion() );
  for ( pd_it.GoToBegin(); !pd_it.IsAtEnd(); ++pd_it )
  {
    pd_it.Set( pd_it.Get() / cy );
  }

  // Smooth avg mask
  itk::ImageRegionIterator< ImageType > avg_it( avg->GetOutput(), avg->GetOutput()->GetLargestPossibleRegion() );
  for ( avg_it.GoToBegin(); !avg_it.IsAtEnd(); ++avg_it )
  {
    avg_it.Set( avg_it.Get() * 1000 );
  }

  
  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > BlurType;
  BlurType::Pointer blur = BlurType::New();
  blur->SetInput( avg->GetOutput() ); 
  blur->SetNormalizeAcrossScale( false );
  blur->SetSigma( 2 );
  blur->SetDirection( 0 );
  blur->Update();
  blur->SetInput( blur->GetOutput() );
  blur->SetDirection( 1 );
  blur->Update();
  blur->SetInput( blur->GetOutput() );
  blur->SetDirection( 2 );
  blur->Update();

  itk::ImageRegionIterator< ImageType > blur_it( blur->GetOutput(), blur->GetOutput()->GetLargestPossibleRegion() );
  for ( blur_it.GoToBegin(); !blur_it.IsAtEnd(); ++blur_it )
  {
    blur_it.Set( std::floor( blur_it.Get() * 2 + 0.5 ) );
  }

  // Start creation of HfB
  ImageType::Pointer hfb = ImageType::New();
  hfb = egitk::CreateEmptyImage< PixelType, Dim >( t2->GetOutput() );
  itk::ImageRegionIterator< ImageType > hfb_it( hfb, hfb->GetLargestPossibleRegion() );

  for ( hfb_it.GoToBegin(), blur_it.GoToBegin(), t2_it.GoToBegin(), pd_it.GoToBegin();
       !hfb_it.IsAtEnd();
       ++blur_it, ++hfb_it, ++t2_it, ++pd_it )
  {
    if ( blur_it.Get() > 1995 ) 
    {
      hfb_it.Set( 1000 );
    }
    if ( t2_it.Get() > parameters.n_t2_threshold && pd_it.Get() > parameters.n_pd_threshold && blur_it.Get() > 50 )
    {
      hfb_it.Set( 1000 );
    }
  }

  // Smooth again

  blur->SetInput( hfb );
  blur->SetNormalizeAcrossScale( false );
  blur->SetSigma( 2.0 );
  blur->SetDirection( 0 );
  blur->Update();
  blur->SetInput( blur->GetOutput() );
  blur->SetDirection( 1 );
  blur->Update();
  blur->SetInput( blur->GetOutput() );
  blur->SetDirection( 2 );
  blur->Update();

  for ( blur_it.GoToBegin(); !blur_it.IsAtEnd(); ++blur_it )
  {
    //blur_it.Set( std::floor( blur_it.Get() * 2 + 0.5 ) );
 
    if ( std::floor( blur_it.Get() * 2 + 0.5 )  > 1000  )
    {
      blur_it.Set( 1 );
    }
    else 
    {
      blur_it.Set( 0 );
    }

  }

  //=========================
  // Write Output Image
  //=========================

  typedef itk::ImageFileWriter< ImageType > WriterType;

  WriterType::Pointer writer = WriterType::New();
  egitk::WriteImage< PixelType, unsigned char, 3 >( parameters.out_filename.c_str(),  blur->GetOutput() );
	return 0;

}

