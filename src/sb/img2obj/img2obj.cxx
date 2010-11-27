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
#include "itkCastImageFilter.h"
#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"

int main( int argc, char ** argv )
{

	if ( argc < 3 )
	{
		std::cout << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << "  " << argv[0] << " imgIN objOUT" << std::endl;
		std::cout << std::endl;
		std::cout << "Input:" << std::endl;
		std::cout << "  imgIN:  Input image" << std::endl;
		std::cout << "  objOUT: Output object map" << std::endl;
		std::cout << std::endl;
		std::cout << "Notes:" << std::endl;
		std::cout << "  -file extension required (.img/.nii/.nii.gz/.obj)" << std::endl;
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

	typedef itk::ImageFileReader< DoubleImageType > DoubleReaderType;
	DoubleReaderType::Pointer imageReader = DoubleReaderType::New();
	imageReader->SetFileName( argv[1] );
	try 
  	{
 	  	imageReader->Update();
	  }
	  catch( itk::ExceptionObject& err )
	  {
	  	std::cerr << std::endl << "Error reading input!" << std::endl << err << std::endl;
			exit(1);
	  }

	float imgmax = egitk::ImageMax< DoublePixelType, Dim >( imageReader->GetOutput() );
	if ( imgmax > 255 )
	{
		std::cout << "Error:  Object maps allow for only 255 labels"  << std::endl;
		std::cout << "Maximum value of input image = " << imgmax << std::endl;
		exit(1);
	}

	typedef itk::CastImageFilter< DoubleImageType, CharImageType > CastType;
	CastType::Pointer caster = CastType::New();
	caster->SetInput( imageReader->GetOutput() );
	caster->Update();

	itk::AnalyzeObjectLabelMapImageIO::Pointer io = itk::AnalyzeObjectLabelMapImageIO::New();
	io->SetNumberOfObjects( imgmax );
	
	typedef itk::ImageFileWriter< CharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetImageIO( io );
	writer->SetInput( caster->GetOutput() );
	writer->SetGlobalWarningDisplay( false );
	writer->SetFileName( argv[2] );
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