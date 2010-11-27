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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkAnalyzeObjectMap.h"
#include "itkAnalyzeObjectLabelMapImageIO.h"
#include "itkAnalyzeObjectLabelMapImageIOFactory.h"
#include "../Common/itk_imgio_util_templates.h"

int main( int argc, char ** argv )
{

	if ( argc < 3 )
	{
		std::cout << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << "  " << argv[0] << " objIN [imgIN] imgOUT" << std::endl;
		std::cout << std::endl;
		std::cout << "Input:" << std::endl;
		std::cout << "  objIN:  Analyze object map" << std::endl;
		std::cout << "  imgOUT: Output image" << std::endl;
		std::cout << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "  imgIN:  Image on which the object map was drawn" << std::endl;
		std::cout << std::endl;
		std::cout << "Notes:" << std::endl;
		std::cout << "  -file extension required (.img/.nii/.nii.gz/.obj)" << std::endl;
		std::cout << "  -if imgIN is supplied, imgOUT will have the correct voxel dimensions" << std::endl;
		std::cout << std::endl;
		exit(1);
	}

	const int Dim = 3;
	typedef unsigned char CharPixelType;
	typedef double        DoublePixelType;
	typedef itk::Image< CharPixelType, Dim > CharImageType;
	typedef itk::Image< DoublePixelType, Dim > DoubleImageType;
	typedef itk::ImageFileReader< CharImageType > CharReaderType;
	typedef itk::ImageFileReader< DoubleImageType > DoubleReaderType;

  itk::ObjectFactoryBase::RegisterFactory( itk::AnalyzeObjectLabelMapImageIOFactory::New() );

	itk::AnalyzeObjectLabelMapImageIO::Pointer io = itk::AnalyzeObjectLabelMapImageIO::New();

	typedef itk::ImageFileReader< CharImageType > ObjectReaderType;
	ObjectReaderType::Pointer objectReader = ObjectReaderType::New();
	objectReader->SetFileName( argv[1] );
	objectReader->SetImageIO( io );
	try 
  	{
 	  	objectReader->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error reading input!" << std::endl << err << std::endl;
			exit(1);
	  }
	
	typedef itk::AnalyzeObjectMap< CharImageType > ObjectMapType;
	ObjectMapType::Pointer objmap = ObjectMapType::New();
	objmap->ImageToObjectMap( objectReader->GetOutput() );
	
	typedef itk::ImageFileWriter< CharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	
	if ( argc == 4 )
	{
	  typedef itk::ImageFileReader< DoubleImageType > ImageReaderType;
	  ImageReaderType::Pointer imageReader = ImageReaderType::New();
	  imageReader->SetFileName( argv[2] );
	try 
  	{
 	  	imageReader->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error reading input!" << std::endl << err << std::endl;
			exit(1);
	  }
	
		DoubleImageType::SpacingType imgSpacing;
		imgSpacing = imageReader->GetOutput()->GetSpacing();
		CharImageType::Pointer output = CharImageType::New();
    output = objmap;
    output->SetSpacing( imgSpacing );
	
		writer->SetFileName( argv[3] );
		writer->SetInput( objmap );
	}
	else
	{
		writer->SetFileName( argv[2] );
		writer->SetInput( objmap );
	}
	writer->SetGlobalWarningDisplay( false );
	
	try 
  	{
 	  	writer->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error reading input!" << std::endl << err << std::endl;
			exit(1);
	  }
	
}