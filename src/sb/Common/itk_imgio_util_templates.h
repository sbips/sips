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

// Written 2006 by Erin Gibson //
//   Last Modified: 04/05/09  //

#ifndef _itk_imgio_util_templates_h_
#define _itk_imgio_util_templates_h_

#include "itkImage.h"
#include "itkAnalyzeImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNiftiImageIO.h"

namespace egitk
{

	template< typename TInputPixelType, int TDim >
	bool
	ReadImage( typename itk::ImageFileReader< itk::Image< TInputPixelType, TDim > >::Pointer reader )
	{
    reader->GlobalWarningDisplayOff();
	  try 
  	{
 	  	reader->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error reading input!" << std::endl << err << std::endl;
			return false;
	  }
		return true;
	} // ReadImage(...)

	template< typename TInputPixelType, typename TOutputPixelType, int TDims >
	int
	WriteImage( const char* filename, typename itk::Image< TInputPixelType, 3 >::Pointer inputVol )
	{
		// Writes itk::Image< TInputPixelType, TDims > to file using datatype specified by TOutputPixelType
		// Image saved as Analyze for filenames ending in ".img", otherwise image saved as Nifti
		//!!!  Unnecessary rescaling is done when input and output datatypes are the same
		//     ( rescaling will almost always be necessary when loading data in as doubles )
			
		typedef itk::Image< TInputPixelType, TDims > InputImageType;
		typedef itk::Image< TOutputPixelType, TDims > OutputImageType;

		// Calculate min and max values in the input image. 
		
		typedef itk::MinimumMaximumImageFilter< InputImageType > MinMaxFilterType;

		typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		
		minMaxFilter->SetInput( inputVol );
		minMaxFilter->Update();

		const TInputPixelType imageMin = minMaxFilter->GetMinimum();
		const TInputPixelType imageMax = minMaxFilter->GetMaximum();
   
   
    if ( imageMax > itk::NumericTraits< TOutputPixelType >::max() || imageMin < itk::NumericTraits< TOutputPixelType >::NonpositiveMin() )
		{
      std::cerr << std::endl << "Error saving image:";
			std::cerr << "  Numerical limit for " << typeid( TOutputPixelType).name() << " exceeded."  <<std::endl;
      return false;
		}

		// Rescale datatype from double preserving the original min & max. 
		
		typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;

		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();

		rescaleFilter->SetInput( inputVol );
		rescaleFilter->SetOutputMinimum( imageMin );
		rescaleFilter->SetOutputMaximum( imageMax );
		rescaleFilter->Update();
		
		// Write out rescaled volume 
		
		typedef itk::ImageFileWriter< OutputImageType > WriterType;

		typename WriterType::Pointer writer = WriterType::New();
	

		writer->SetInput( rescaleFilter->GetOutput() );

		std::string filename_ext = filename;
		filename_ext = filename_ext.substr( filename_ext.length()-4, filename_ext.length() );
		
		if ( filename_ext == ".img" )
		{
				itk::AnalyzeImageIO::Pointer io_analyze = itk::AnalyzeImageIO::New();
    		io_analyze->SetGlobalWarningDisplay( false );
				writer->SetImageIO( io_analyze );
		}
		else
		{
			itk::NiftiImageIO::Pointer io_nifti = itk::NiftiImageIO::New();
   		io_nifti->SetGlobalWarningDisplay( false );
			writer->SetImageIO( io_nifti );
		}

		writer->SetGlobalWarningDisplay( false );
		writer->SetFileName( filename );
	  try 
  	{
 	  	writer->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error writing input!" << std::endl << err << std::endl;
			return false;
	  }

		return true;
	}; // writeImage(...)

	template< typename TPixelType, int TDims >
	bool
	WriteImageSamePixelType( const char* filename, typename itk::ImageFileReader< typename itk::Image< TPixelType, 3 > >::Pointer inputReader, typename itk::Image< TPixelType, TDims >::Pointer outputImage )
	{
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned char ) )
		WriteImage<TPixelType, unsigned char, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned int ) )
		WriteImage<TPixelType, unsigned int, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned short ) )
		WriteImage<TPixelType, unsigned short, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( unsigned long ) )
		WriteImage<TPixelType, unsigned long, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( signed char ) )
		WriteImage<TPixelType, signed char, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( signed int ) )
		WriteImage<TPixelType, signed int, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( signed short ) )
		WriteImage<TPixelType, signed short, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( signed long ) )
		WriteImage<TPixelType, signed long, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid( float ) )
		WriteImage<TPixelType, float, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid(double ) )
		WriteImage<TPixelType, double, 3>( filename, outputImage );
	if ( inputReader->GetImageIO()->GetComponentTypeInfo() == typeid(short ) )
		WriteImage<TPixelType, short, 3>( filename, outputImage );

	return true;
	}; // WriteImageSamePixelType


}

#endif
