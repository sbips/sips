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

// Last Modified: 10/15/2008

#include "itkImage.h"

#include "../Common/itk_imgops_util_templates.h"
#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgmorph_util_templates.h"

#include <vector>

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkNumericTraits.h"
#include "itkNeighborhoodOperator.h"

// Global functions.

int usage( char* argv[] )
{
  std::cout << std::endl;
	std::cout << "Usage: " << argv[0] << std::endl;
  std::cout << "       <flhyper.img> <wmt.img> <wmt_threshold> <outfn_prefix>" << std::endl;
  std::cout << "       " << std::endl;
	std::cout << "  <flhyper.img>: binary image of segmented FLAIR hyperintensities" << std::endl;
	std::cout << "  <wmt.img>: coregistered MNI white matter template" << std::endl;
	std::cout << "  <wmt_threshold>: threshold for white matter to apply to the template" << std::endl;
	std::cout << "  <outfn_prefix>:  output filename prefix" <<std::endl <<std::endl;
	std::cout << "Notes:" <<std::endl;
	std::cout << "  - Rejects any hyperintensity not connected to the thresholded" <<std::endl;
	std::cout << "      white matter template (using 6-face 3D connectivity)" << std::endl;
	std::cout << "  - Output image contains \"1\"s for accepted and \"2\"s for rejected lesions"  <<std::endl;
	std::cout << "  - Must supply the .img extension for image files";
  std::cout << std::endl;
  return 1;
}

// Global program.

int main( int argc, char* argv[] )
{
	if ( argc != 5 )
    return usage( argv );
  
  typedef double PixelType;
  const int dims = 3;
  typedef itk::Image< PixelType, dims > ImageType;

	std::cout<< "Reading input volumes..." << std::endl;
	
	// Load in white matter template (wmt) and lesion (flhyper) images.
  ImageType::Pointer flhyper = ImageType::New();
  ImageType::Pointer wmt = ImageType::New();
	flhyper = io_utils::ReadItkVol< PixelType, dims >( argv[1] );
  wmt = io_utils::ReadItkVol< PixelType, dims >( argv[2] );
  
	// Get threshold for white matter template.
	int WM_threshold = atoi( argv[3] );

	std::cout << "Thresholding white matter template...";
	
	// Threshold wmt image.
	wmt = itk_iops::ThresholdImageToBinary< PixelType, PixelType, dims >( wmt, 1, 0, WM_threshold, itk::NumericTraits< PixelType >::max() );

	std::cout << "Preparing to reject hyperintensities using the template..." << std::endl;

	// Create Structuring element ("1"s for voxels face-connected to the center of a 3x3x3 cube).
	int sevals[] = {0,0,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0};
	typedef itk::Neighborhood< PixelType, dims > NeighborhoodType;
	NeighborhoodType se;
	se.SetRadius( 1 );
	int i=0;
	for ( itk::Neighborhood< PixelType,dims>::Iterator it = se.Begin(); it != se.End(); ++it, ++i )
		*it = sevals[i];

	// Dilate wmt using the above structuring element.
	//   This will define the search space for 3D connected hyperintensities.
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, NeighborhoodType > DilateFilterType;
	DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
	dilateFilter->SetInput( wmt );
	dilateFilter->SetKernel( se );
	dilateFilter->SetDilateValue( 1 );
	dilateFilter->Update();
	
	ImageType::Pointer wmt_dilate = ImageType::New();
  wmt_dilate=dilateFilter->GetOutput();
	wmt_dilate->DisconnectPipeline();

	// Assign a unique label to each connected hyperintensity on the flhyper image.
	ImageType::Pointer flhyper_lab = ImageType::New();
	flhyper_lab = itk_iops::RelabelConnectedComponents< PixelType, PixelType, dims > ( flhyper, 1 );
    flhyper_lab = itk_imath::ImageFindValueAndSet< PixelType, PixelType, dims >( flhyper_lab, flhyper, 0, 0 );
		 		
		//  Create a lookup vector for each labelled hyperintensity. 
	int num_lesions = itk_imath::ImageMax<PixelType, dims>( flhyper_lab );
  std::vector< int >* lookup = new std::vector< int >( num_lesions + 1 );

	// Make a copy of the labelled fl_hyper image; Zero any voxels on the fl_hyper;
	//  image that correspond to zero voxels on the dilated thresholded wmt image.
	ImageType::Pointer keep = ImageType::New();
	keep = itk_iops::DuplicateImage< PixelType, dims >( flhyper_lab );
	keep = itk_imath::ImageFindValueAndSet< PixelType, PixelType, dims >( keep, wmt_dilate, 0, 0 );

	std::cout << "First pass over image..." << std::endl;
	
	// First pass over image to find lesion labels to accept.
	typedef itk::ImageRegionIterator< ImageType > ItType;
	ItType it2( keep, keep->GetLargestPossibleRegion() );
  for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2  )
	{
		// Store the labels of any hyperintensity not rejected after
		//  masking with the dilated thresholded wmt image in the lookup vector.
		int value = it2.Get();
		if ( value != 0 ) 
			lookup->at( value ) = value;
	}

	std::cout << "Second pass over image..." << std::endl;
	// Second pass over image to remove hyperintensities not accepted in the first pass.
	ItType it3( flhyper_lab, flhyper_lab->GetLargestPossibleRegion() );
	for ( it3.GoToBegin(); !it3.IsAtEnd(); ++it3 )
	{
		int value = it3.Get();
		if ( lookup->at( value ) == 0 && it3.Get() != 0 )
			it3.Set( 2 ); // reject
		if ( lookup->at( value ) != 0 )
			it3.Set( 1 ); // accept
	}

	std::cout << "Done..." <<std::endl;

  // Save output.
	
	std::cout << "Writing output volume..." << std::endl;
  std::string outfn;
  outfn = argv[4]; outfn += "_flwmt_lesions.img";
  io_utils::WriteItkVol< PixelType, unsigned char, dims >( outfn.c_str(), flhyper_lab );


	return 0;
}
