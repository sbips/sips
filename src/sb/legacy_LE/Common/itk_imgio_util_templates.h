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
//   Last Modified: 10/12/08  //

#ifndef _itk_IO_util_templates_h_
#define _itk_IO_util_templates_h_

#include "itkImage.h"
#include "itkAnalyzeImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


namespace io_utils
{

	template< typename TPixelType, int TDim >
	typename itk::Image< TPixelType, TDim >::Pointer
	ReadItkVol( const char* filename )
	{
		// Reads image and stores it as the datatype specified by TPixelType
		//   !! IO mode determined by filename extension, no other checking done
		//   !! no checking done to determine if TPixelType is appropriate for
		//        input image ( expecting that all data will be stored as doubles
		//        for program duration and then saved using a smaller datatype )
		
		typedef itk::Image< TPixelType, TDim > ImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;

		typename ReaderType::Pointer reader = ReaderType::New();
    reader->GlobalWarningDisplayOff();
		reader->SetFileName( filename );
		reader->Update();
		
		typename ImageType::Pointer temp = ImageType::New();
		temp = reader->GetOutput();
		temp->DisconnectPipeline();

		return temp;
	} // readItkVol(...)


	template< typename TInputPixelType, typename TOutputPixelType, int TDims >
	int
	WriteItkVol( const char* filename, typename itk::Image< TInputPixelType, 3 >::Pointer inputVol )
	{
		// Writes itk::Image< TInputPixelType, 3 > to file using datatype specified by TOutputPixelType
		//!!!  No checking is done to determine whether the datatype conversion is safe
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

		if ( imageMax > itk::NumericTraits< TOutputPixelType >::max() )
		{
      std::cerr << std::endl << "Error saving image:";
			std::cerr << "  Numerical limit for " << typeid( TOutputPixelType).name() << " exceeded."  <<std::endl;
      return 1;
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
                itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
		io->SetGlobalWarningDisplay( false );

		writer->SetInput( rescaleFilter->GetOutput() );
		writer->SetImageIO( io );
		writer->SetGlobalWarningDisplay( false );

		writer->SetFileName( filename );
		writer->Update();

		return 0;
	}; // writeItkVol(...)




}

#endif
