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

//Last Modified: 10/12/2008

#ifndef _utils_h_
#define _utils_h_


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMinimumMaximumImageCalculator.h"


namespace utils
{

	typedef itk::Image< short, 3 > ImageType3D16;
	typedef itk::Image< short, 2 > ImageType2D16;

	std::vector< double > GetVectorSliceFromITKvol ( ImageType3D16::Pointer input, int slicenum )
	{
		typedef itk::ExtractImageFilter< ImageType3D16, ImageType2D16 > ExtractFilterType;

		ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
		
		extractFilter->SetInput( input );
	
		ImageType3D16::IndexType  extractIndex = {{0, 0, slicenum}};

		ImageType3D16::SizeType   extractSize;
		extractSize = input->GetLargestPossibleRegion().GetSize();
		extractSize[2] = 0;

		ImageType3D16::RegionType extractRegion;
		extractRegion.SetSize( extractSize );
		extractRegion.SetIndex( extractIndex );
		
		extractFilter->SetExtractionRegion( extractRegion );
		extractFilter->UpdateLargestPossibleRegion();

		ImageType2D16::Pointer inputSlice = ImageType2D16::New();

		inputSlice = extractFilter->GetOutput();
		//inputSlice->DisconnectPipeline();
		//extractFilter->ResetPipeline();

		const int xdim = input->GetLargestPossibleRegion().GetSize().operator [](0);
		const int ydim = input->GetLargestPossibleRegion().GetSize().operator [](1);
		
		std::vector< double > currentSliceVector( xdim * ydim );

			
		typedef itk::ImageRegionIterator< ImageType2D16 > it2D16Type;

		it2D16Type it( inputSlice, inputSlice->GetLargestPossibleRegion() );
		
		for (int i = 0 ; !it.IsAtEnd(); ++it, ++i )
			currentSliceVector[i] = it.Get();	

		return currentSliceVector;
	}

	void SaveVectorSliceToITKvol ( std::vector< double >* U, ImageType3D16::Pointer output, int slicenum )
	{

		
		ImageType3D16::IndexType  extractIndex = {{0, 0, slicenum}};

		ImageType3D16::SizeType   extractSize;
		extractSize = output->GetLargestPossibleRegion().GetSize();
		extractSize[2] = 1;

		ImageType3D16::RegionType copyRegion;
		copyRegion.SetSize( extractSize );
		copyRegion.SetIndex( extractIndex );

		typedef itk::ImageRegionIterator< ImageType3D16 > it3D16Type;

		it3D16Type it( output, copyRegion );
		
		for (int i = 0 ; !it.IsAtEnd(); ++it, ++i )
			it.Set( U->operator [](i) * 1000 );	

		return;

	}

	bool IsEmptySliceVector( const std::vector< double >& slice )
	{
    int voxels = 0;
		for ( int i = 0 ; i < slice.size() ; ++ i ) 
			if ( slice[i] != 0 )
      {
        ++voxels;
        if ( voxels > 25 )
				  return false;
      }
		return true;
	}



	ImageType3D16::Pointer ItkReadInput( std::string infn  )
	{
		typedef itk::ImageFileReader< ImageType3D16 > ReaderType3D16;

		ReaderType3D16::Pointer reader3D16 = ReaderType3D16::New();
		reader3D16->SetFileName( infn.c_str() );
		reader3D16->Update();

		ImageType3D16::Pointer input = ImageType3D16::New();
		input = reader3D16->GetOutput();
		
		return input;
	}

	void ItkSaveOutput( const std::string& infn, ImageType3D16::Pointer output )
	{
		typedef itk::ImageFileWriter< ImageType3D16 > WriterType3D16;

		WriterType3D16::Pointer writer3D16 = WriterType3D16::New();

		writer3D16->SetInput( output );
		writer3D16->SetFileName( infn.c_str() );
		writer3D16->Update();


	}

}

#endif
