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

#include "itkMatrix.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "../Common/itk_imgio_util_templates.h"


// - Inverts a 4x4 rotation matrix obtained using oblique module in Analyze for use with SABRE.

// - Input Matrix = [  R   R   R   0
//                     R   R   R   0
//                     R   R   R   0
//                    -Tx -Ty -Tz  1 ];

// - Output Matrix = [ R   R   R   Tx
//                     R   R   R   Ty
//                     R   R   R   Tz
//                     0   0   0   1 ];

int main( int argc, char *argv[] )
{
	//----------------------------------------------------
	// Check input
	//----------------------------------------------------

	if ( argc < 5 ) 
	{
		std::cout << std::endl;
		std::cout << "Usage: " << std::endl;
		std::cout << "  " << argv[0] << " isotropicT1 T1 T1toACPC.mat prefix" << std::endl;
		std::cout << std::endl;
		std::cout << "Input:" << std::endl;
		std::cout << "  isotropicT1:  output from make_isotropic or zoomer" << std::endl;
		std::cout << "  T1:           T1 as acquired" << std::endl;
		std::cout << "  T1toACPC.mat: T1 to ACPC rotation matrix obtained in Analyze" << std::endl;
		std::cout << "  prefix:       prefix to be combined with _T1_to_ACPC.mat & _ACPC_to_T1.mat" << std::endl;
		std::cout << std::endl;
		std::cout << "Notes:" << std::endl;
		std::cout << "  -file extension required (.img/.nii/.nii.gz)" << std::endl;
		std::cout << std::endl;
		return 1;
	}

	//----------------------------------------------------
	// Load input images and matrix
	//----------------------------------------------------

	const unsigned int dim = 3;
	typedef float PixelType;
	typedef itk::Image<PixelType, dim> ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	// Load input and "output" (reference) images.
	ReaderType::Pointer inputReader = ReaderType::New();
	ReaderType::Pointer outputReader = ReaderType::New();
	inputReader->SetFileName( argv[1] );
	outputReader->SetFileName( argv[2] );
	egitk::ReadImage< PixelType, dim >( inputReader );
	egitk::ReadImage< PixelType, dim >( outputReader );

	ImageType::Pointer inputImage = ImageType::New();
	ImageType::Pointer outputImage = ImageType::New();
	inputImage = inputReader->GetOutput();
	outputImage = outputReader->GetOutput();
 
	// Get input and output voxel sizes.
	ImageType::SpacingType insize = inputImage->GetSpacing();
	ImageType::SpacingType outsize = outputImage->GetSpacing();
	
	// Load analyze matrix to invert.
	itk::Matrix<float, 4, 4> e;

	std::ifstream matrix( argv[3] );

	for (int i=0; i < 4; i++) {
		for (int j=0; j < 4; j++) {
			// load transposed
			matrix >> e[j][i];
		}
	}
	// multiple translations by -1
	e[0][3] *= -1; 	e[1][3] *= -1; 	e[2][3] *= -1;

	//----------------------------------------------------
	// Write ACPC_to_T1 matrix
	//----------------------------------------------------

	std::string outputACPCtoT1_fn = argv[4];
	outputACPCtoT1_fn += "_acpc_to_T1.mat";
	std::ofstream outputACPCtoT1( outputACPCtoT1_fn.c_str() );
	outputACPCtoT1 << e;

	//----------------------------------------------------
	// Invert matrix
  //----------------------------------------------------
	e = e.GetInverse();
	
	//----------------------------------------------------
	// Set tolerance for voxel size differences
  //----------------------------------------------------

	// Voxel size differences are considered equal unless greater than this value:
	const float voxel_size_diff_tolerance = 0.01;
	
	//----------------------------------------------------
	// Adjust matrix for differences in input voxel sizes
	//----------------------------------------------------

	float pixel_size_i = insize[0];
	if (insize[1] < pixel_size_i) { pixel_size_i = insize[1]; }
	if (insize[2] < pixel_size_i) { pixel_size_i = insize[2]; }
	
	if ( fabs(insize[0] - pixel_size_i) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying x pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][0] /= (insize[0] / pixel_size_i);
		}
	}
	if ( fabs(insize[1] - pixel_size_i) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying y pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][1] /= (insize[1] / pixel_size_i);
		}
	}
	if ( fabs(insize[2] - pixel_size_i) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying z pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][2] /= (insize[2] / pixel_size_i);
		}
	}

	//----------------------------------------------------
	// Adjust matrix for differences in output voxel sizes
	//----------------------------------------------------

	float pixel_size_o = outsize[0];
	if (outsize[1] < pixel_size_o) { pixel_size_o = outsize[1]; }
	if (outsize[2] < pixel_size_o) { pixel_size_o = outsize[2]; }

	if ( fabs(outsize[0] - pixel_size_o) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying x pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[0][i] /= (outsize[0] / pixel_size_o);
		}
	}
	if ( fabs(outsize[1] - pixel_size_o) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying y pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[1][i] /= (outsize[1] / pixel_size_o);
		}
	}
		if ( fabs(outsize[2] - pixel_size_o) > voxel_size_diff_tolerance )
	{
		//std::cout << "Modifying z pixel sizes" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[2][i] /= (outsize[2] / pixel_size_o);
		}
	}

	//----------------------------------------------------
	// Write T1_to_ACPC matrix
	//----------------------------------------------------

	std::string outputT1toACPC_fn = argv[4];
	outputT1toACPC_fn += "_T1_to_acpc.mat";
	std::ofstream outputT1toACPC( outputT1toACPC_fn.c_str() );
	outputT1toACPC << e;

  return 0;

}

