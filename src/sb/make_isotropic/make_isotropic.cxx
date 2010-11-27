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

#include "../Common/itk_imgio_util_templates.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"

int main( int argc, char *argv[] )
{

	// Check input arguments.
	if ( argc < 3 ) 
	{
		std::cout << std::endl;
		std::cout << "Usage: " << std::endl;
		std::cout << "  " << argv[0]   << " imgIN imgOUT" << std::endl;
		std::cout << std::endl;
		std::cout << "Notes:" << std::endl;
		std::cout << "  -file extension required (.img/.nii/.nii.gz)" << std::endl;
		std::cout << std::endl;
		return 1;
	}

	// Read input.
	const unsigned int dim = 3;
	typedef float PixelType;
	typedef itk::Image<PixelType, dim> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	ReaderType::Pointer inputReader = ReaderType::New();
	inputReader->SetFileName( argv[1] );

	if ( !egitk::ReadImage< PixelType, dim >( inputReader ) )
		exit(1);

	ImageType::Pointer inputImage = ImageType::New();
	inputImage = inputReader->GetOutput();

	// Get input voxel sizes and dimensions.
	ImageType::SpacingType inSize = inputImage->GetSpacing();
	ImageType::SizeType inSpacing = inputImage->GetLargestPossibleRegion().GetSize();

	// Get isotropic voxel size.
	float isotropic_voxel_size = inSize[0];

	for (int i=1; i < dim; ++i)
	{
		if (inSize[i] - isotropic_voxel_size < 0.001)
		{
			isotropic_voxel_size = inSize[i];
		}
	}

	// Calculate new image dimensions.
	ImageType::SizeType outSize;
	for (int i=0; i < dim; ++i )
	{
		outSize[i] = std::floor(( inSize[i] * (inSpacing[i]-1) / isotropic_voxel_size ) + 1);
	}

	// Fix to ensure lobmask command receives input with x dimension = y dimension
	if ( outSize[0] > outSize[1] ) 
	{
		outSize[1] = outSize[0];
	}
	else if ( outSize[1] > outSize[0] )
	{
		outSize[0] = outSize[1];
	}

	// Set new image isotropic spacing.
	ImageType::SpacingType outSpacing;
	outSpacing[0] = isotropic_voxel_size; 
	outSpacing[1] = isotropic_voxel_size; 
	outSpacing[2] = isotropic_voxel_size;

	// Resample input image.
	typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleType;
	
	ResampleType::Pointer resampler = ResampleType::New();

	resampler->SetInput( inputImage );
	resampler->SetOutputSpacing( outSpacing );
	resampler->SetSize( outSize );
	resampler->SetOutputDirection( inputImage->GetDirection() );
	resampler->Update();

	// Write output image.
	egitk::WriteImageSamePixelType< PixelType, dim > ( argv[2], inputReader, resampler->GetOutput() );

  return 0;

}

