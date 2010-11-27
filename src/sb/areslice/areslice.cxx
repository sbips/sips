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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#include "../Common/itk_imgio_util_templates.h"
#include "itkMatrix.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageSliceIteratorWithIndex.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Notes:
// - Reslices input image to space defined by a reference image and a rotation matrix.
// - Rotation matrix specifies a voxel-to-voxel mapping (consistent with AIR).
// - Intended to be used in SABRE to rotate images to and from ACPC/T1 space.
//      To use in more general situations, the size and spacing of the input and reference
//      images should be added to the rotation matrix (i.e. size and spacing info added
//      to the output from the atransform command).  Currently, atransform creates a rotation
//      matrix using an input and output image, and for areslice to reslice the image correctly,
//      it must receive as input the rotation matrix, as well as the input and output images passed
//      to atransform in order to generate the rotation matrix.  Alternatively, using a rotation 
//      matrix that specifies a world-to-world mapping would advoid this problem altogether.


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
void
Usage( char* argv[] )
{
	std::cout << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << "  " << argv[0] << " imgIN imgREF INtoREF.mat imgOUT [options]" << std::endl;
	std::cout << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "  -zpad <int>: number of z slices to add to the reference volume" << std::endl;
	std::cout << "  -nn:         nearest neighbor interpolation (default: trilinear)" << std::endl;
	std::cout << std::endl;
	std::cout << "Notes:" << std::endl;
	std::cout << "  -file extension required (.img/.nii/.nii.gz)" << std::endl;
	std::cout << "  -reslices an image according to a voxel-to-voxel rotation matrix" << std::endl;
	std::cout << std::endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

template< typename TPixelType, int TDims >
itk::Matrix< float, 4, 4 >
ModifyVoxelSizes( typename itk::Matrix<float,4,4> e, 
								  typename itk::Image< TPixelType, TDims>::SpacingType inSpacing,
									typename itk::Image< TPixelType, TDims>::SpacingType refSpacing)
{

	// Voxel size differences are considered equal unless greater than this value:
	const float voxel_spacing_diff_tolerance = 0.01;
	
	float pixel_spacing_o = refSpacing[0];
	if (refSpacing[1] < pixel_spacing_o) { pixel_spacing_o = refSpacing[1]; }
	if (refSpacing[2] < pixel_spacing_o) { pixel_spacing_o = refSpacing[2]; }

	if ( fabs(refSpacing[0] - pixel_spacing_o) > voxel_spacing_diff_tolerance )
	{
		//std::cout << "Modifying x pixel Spacings" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][0] *= (refSpacing[0] / pixel_spacing_o);
		}
	}
	if ( fabs(refSpacing[1] - pixel_spacing_o) > voxel_spacing_diff_tolerance )
	{
		//std::cout << "Modifying y pixel Spacings" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][1] *= (refSpacing[1] / pixel_spacing_o);
		}
	}
		if ( fabs(refSpacing[2] - pixel_spacing_o) > voxel_spacing_diff_tolerance )
	{
		//std::cout << "Modifying z pixel Spacings" << std::endl;
		for (int i = 0; i < 4; ++i)
		{
			e[i][2] *= (refSpacing[2] / pixel_spacing_o);
		}
	}

		return e;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// Nearest Neighbor Interpolation
template< typename TPixelType, int TDims >
typename itk::Image< TPixelType, TDims >::Pointer
ResampleImageNN( typename itk::Image< TPixelType, TDims >::Pointer input,
								 typename itk::Image< TPixelType, TDims >::Pointer reference,
								 typename itk::Matrix<float, 4, 4> e,
								 int zpad
)
{
	typedef itk::Image< TPixelType, TDims > ImageType;

	typename ImageType::SizeType inSize        = input->GetLargestPossibleRegion().GetSize();
	typename ImageType::SizeType refSize       = reference->GetLargestPossibleRegion().GetSize();
	typename ImageType::SpacingType refSpacing = reference->GetSpacing();

	refSize[2] += zpad;
	typename ImageType::IndexType start;
	start[0] = 0; start[1] = 0; start[2] = 0;
	typename ImageType::RegionType outputRegion;
	outputRegion.SetSize( refSize );
	outputRegion.SetIndex( start );

	typename ImageType::Pointer output = ImageType::New();
	output->SetSpacing( refSpacing );
	output->SetDirection( input->GetDirection() );
	output->SetOrigin( input->GetOrigin() );
	output->SetRegions( outputRegion );
	output->Allocate();
	output->FillBuffer( 0 );

		typename ImageType::IndexType inputPixelIndex; 
		typename ImageType::IndexType outputPixelIndex;


		float x_k = e(0,3);
		float y_k = e(1,3);
		float z_k = e(2,3);
		
		std::cout << "Resampling Image...  ";
		for (int k = 0; k < refSize[2]; ++k )
		{
			float x_j = x_k; 
			float y_j = y_k;
			float z_j = z_k;

			for (int j = 0; j < refSize[1]; ++j )
			{
			  float x_i = x_j;
        float y_i = y_j;
        float z_i = z_j;
			 
				for (int i = 0; i < refSize[0]; ++i )
				{
					  float x_p = x_i;
            float y_p = y_i;
            float z_p = z_i;

						if ( x_p >= 0 && x_p < inSize[0] - 1 )
						{
							if ( y_p >=0 && y_p < inSize[1] - 1 )
							{
								if ( z_p >= 0 && z_p < inSize[2] - 1 )
								{
									inputPixelIndex[0]  = floor( x_p + 0.5 );
									inputPixelIndex[1]  = floor( y_p + 0.5 );
									inputPixelIndex[2]  = floor( z_p + 0.5 );
									outputPixelIndex[0] = i; 
									outputPixelIndex[1] = j; 
									outputPixelIndex[2] = k;
									output->SetPixel( outputPixelIndex, input->GetPixel(inputPixelIndex) );
								}
							}
						}
					x_i += e(0,0);
					y_i += e(1,0);
					z_i += e(2,0);
				}
			  x_j += e(0,1);
				y_j += e(1,1);
				z_j += e(2,1);
			}
			x_k += e(0,2);
			y_k += e(1,2);
			z_k += e(2,2);
		}
		std::cout << "Finished." << std::endl;

	return output;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Trilinear Interpolation
template< typename TPixelType, int TDims >
typename itk::Image< TPixelType, TDims >::Pointer
ResampleImageTLIN( typename itk::Image< TPixelType, TDims >::Pointer input,
								 typename itk::Image< TPixelType, TDims >::Pointer reference,
								 typename itk::Matrix<float, 4, 4> e,
								 int zpad
)
{
	typedef itk::Image< TPixelType, TDims > ImageType;

	typename ImageType::SizeType inSize        = input->GetLargestPossibleRegion().GetSize();
	typename ImageType::SizeType refSize       = reference->GetLargestPossibleRegion().GetSize();
	typename ImageType::SpacingType refSpacing = reference->GetSpacing();

	refSize[2] += zpad;
	typename ImageType::IndexType start;
	start[0] = 0; start[1] = 0; start[2] = 0;
	typename ImageType::RegionType outputRegion;
	outputRegion.SetSize( refSize );
	outputRegion.SetIndex( start );

	typename ImageType::Pointer output = ImageType::New();
	output->SetSpacing( refSpacing );
	output->SetDirection( input->GetDirection() );
	output->SetOrigin( input->GetOrigin() );
	output->SetRegions( outputRegion );
	output->Allocate();
	output->FillBuffer( 0 );



		float x_k = e(0,3);
		float y_k = e(1,3);
		float z_k = e(2,3);
		
		std::cout << "Resampling Image...  ";
		for (int k = 0; k < refSize[2]; ++k )
		{
			float x_j = x_k; 
			float y_j = y_k;
			float z_j = z_k;

			for (int j = 0; j < refSize[1]; ++j )
			{
			  float x_i = x_j;
        float y_i = y_j;
        float z_i = z_j;
			 
				for (int i = 0; i < refSize[0]; ++i )
				{
					  float x_p = x_i;
            float y_p = y_i;
            float z_p = z_i;

						if ( x_p >= 0 && x_p < inSize[0] - 1 )
						{
							if ( y_p >=0 && y_p < inSize[1] - 1 )
							{
								if ( z_p >= 0 && z_p < inSize[2] - 1 )
								{

					 				float dx = x_p - floor(x_p);
									float dy = y_p - floor(y_p);
									float dz = z_p - floor(z_p);
								
									float x000, x001, x010, x011, x100, x101, x110, x111;
									float y000, y001, y010, y011, y100, y101, y110, y111;
									float z000, z001, z010, z011, z100, z101, z110, z111;
									float v000, v001, v010, v011, v100, v101, v110, v111;

									x000 = floor(x_p);
									x100 = x000 + 1;

									if ( floor(x_p) < 0 )
									{
										x000 = 1;
										x100 = x000;
									}
									else if ( floor(x_p) > inSize[0] - 2  )
									{
										x000 = inSize[0] - 1;
										x100 = x000;
									}

									x010 = x000;
									x001 = x000;
									x011 = x000;

									x110 = x100;
									x101 = x100;
									x111 = x100;

									y000 = floor(y_p);
									y010 = y000 + 1;

									if ( floor(y_p) < 0 )
									{
										y000 = 1;
										y100 = y000;
									}
									else if ( floor(y_p) > inSize[1] - 2 )
									{
										y000 = inSize[1] - 1;
										y010 = y000;
									}
									
									y100 = y000;
									y001 = y000;
									y101 = y000;

									y110 = y010;
									y011 = y010;
									y111 = y010;

									z000 = floor(z_p);
									z001 = z000 + 1;

									if ( floor(z_p) < 0 )
									{
										z000 = 1;
										z001 = z000;
									}
									else if ( floor(z_p) > inSize[2] - 2 )
									{
										z000 = inSize[2] - 1;
										z001 = z000;
									}

									z100 = z000;
									z010 = z000;
									z110 = z000;

									z101 = z001;
									z011 = z001;
									z111 = z001;

									x010 = x000;
									x001 = x000;
									x011 = x000;

									x110 = x100;
									x101 = x100;
									x111 = x100;
				
									typename ImageType::IndexType pixelIndex;
									
									pixelIndex[0] = x000; pixelIndex[1] = y000; pixelIndex[2] = z000;
									//std::cout << pixelIndex << std::endl;
									v000 = input->GetPixel( pixelIndex );
							
									pixelIndex[0] = x010; pixelIndex[1] = y010; pixelIndex[2] = z010;
									//std::cout << pixelIndex << std::endl;
									v010 = input->GetPixel( pixelIndex );

									pixelIndex[0] = x001; pixelIndex[1] = y001; pixelIndex[2] = z001;						
									//std::cout << pixelIndex << std::endl;
									v001 = input->GetPixel( pixelIndex );

									pixelIndex[0] = x011; pixelIndex[1] = y011; pixelIndex[2] = z011;
									//std::cout << pixelIndex << std::endl;
									v011 = input->GetPixel( pixelIndex );

									pixelIndex[0] = x100; pixelIndex[1] = y100; pixelIndex[2] = z100;
									//std::cout << pixelIndex << std::endl;
									v100 = input->GetPixel( pixelIndex );
									
									pixelIndex[0] = x110; pixelIndex[1] = y110; pixelIndex[2] = z110;
									//std::cout << pixelIndex << std::endl;
									v110 = input->GetPixel( pixelIndex );
									
									pixelIndex[0] = x101; pixelIndex[1] = y101; pixelIndex[2] = z101;
									//std::cout << pixelIndex << std::endl;
									v101 = input->GetPixel( pixelIndex );

									pixelIndex[0] = x111; pixelIndex[1] = y111; pixelIndex[2] = z111;
									//std::cout << pixelIndex << std::endl;
									v111 = input->GetPixel( pixelIndex );

									float value = 0;
									value += v000 * (1-dx) * (1-dy) * (1-dz);
									value += v010 * (1-dx) *  dy    * (1-dz);
									value += v001 * (1-dx) * (1-dy) *   dz;
									value += v011 * (1-dx) *  dy    *   dz;
									
									value += v100 *  dx    * (1-dy) * (1-dz);
									value += v110 *  dx    *  dy    * (1-dz);
									value += v101 *  dx    * (1-dy) *  dz;
									value += v111 *  dx    *  dy    *  dz;
									
									pixelIndex[0] = i; pixelIndex[1] = j; pixelIndex[2] = k;
									output->SetPixel( pixelIndex, value );	

								}
							}
						}
					x_i += e(0,0);
					y_i += e(1,0);
					z_i += e(2,0);
				}
			  x_j += e(0,1);
				y_j += e(1,1);
				z_j += e(2,1);
			}
			x_k += e(0,2);
			y_k += e(1,2);
			z_k += e(2,2);
		}
		std::cout << "Finished." << std::endl;

	return output;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
	// Check arguments
	
	if ( argc < 5 ) 
	{
		Usage( argv );
		exit( 1 );
	}

	// Check for z padding
	
	int zpad = 0;
	bool nn =  false;
	for ( int i=0; i < argc; ++i )
	{
		std::string cl(argv[i]);
		if ( cl == "-zpad" && i != argc - 1)
		{
			zpad = std::atoi(argv[i+1]);
			std::cout << "Padding z dimension by " << zpad << std::endl;
		}
		if ( cl == "-nn" )
		{
		  nn = true;
			std::cout << "Using nearest neighbor interpolation" << std::endl;
		}

	}

	// Define image types and load input and reference images

	const unsigned int dim = 3;
	typedef float PixelType;
	typedef itk::Image<PixelType, dim> ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;

	ReaderType::Pointer inputReader     =  ReaderType::New();
	ReaderType::Pointer referenceReader =  ReaderType::New();

	inputReader->SetFileName( argv[1] );
	referenceReader->SetFileName( argv[2] );

	if ( !egitk::ReadImage< PixelType, dim >( inputReader ) )     return 1;
	if ( !egitk::ReadImage< PixelType, dim >( referenceReader ) ) return 1;
	
	ImageType::Pointer input = ImageType::New();
	ImageType::Pointer reference = ImageType::New();

	input = inputReader->GetOutput();
	reference = referenceReader->GetOutput();

	ImageType::SpacingType inSpacing =  input->GetSpacing();
	ImageType::SpacingType refSpacing = reference->GetSpacing();

	// Load rotation matrix and modify voxel sizing if necessary

	itk::Matrix<float, 4, 4> e;

	std::ifstream matrix( argv[3] );

	for (int i=0; i < 4; i++) {
		for (int j=0; j < 4; j++) {
			matrix >> e[i][j];
		}
	}  

	e = ModifyVoxelSizes< PixelType, dim >(e, inSpacing, refSpacing );

	// Resample and save output

	ImageType::Pointer output = ImageType::New();

	if ( nn )
	{
		output = ResampleImageNN< PixelType, dim>( input, reference, e, zpad );
	}
	else 
	{
		output = ResampleImageTLIN< PixelType, dim>( input, reference, e, zpad );
	}

	egitk::WriteImageSamePixelType< PixelType, dim >(argv[4], inputReader, output );

	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
