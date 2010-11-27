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

#ifndef _itk_imgmath_util_templates_h_
#define _itk_imgmath_util_templates_h_


#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageIOBase.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkAnalyzeImageIO.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"


namespace egitk
{
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  CreateEmptyImage( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
    // Create an empty image with the same dimensions as the input image.
  		
		typename itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer empty = ImageType::New();
   
    empty->SetLargestPossibleRegion( input->GetLargestPossibleRegion() );
    empty->SetBufferedRegion( input->GetLargestPossibleRegion() );
    empty->SetRequestedRegion( input->GetLargestPossibleRegion() );
    empty->SetSpacing( input->GetSpacing() );
    
		/*
		itk::Matrix<double,3,3> cosines;
    cosines[0][0]= 1;
    cosines[0][1]= 0;
    cosines[0][2]= 0;
    cosines[1][0]= 0;
    cosines[1][1]=-1;
    cosines[1][2]= 0;
    cosines[2][0]= 0;
    cosines[2][1]= 0;
    cosines[2][2]= 1;
	  */

  empty->SetDirection( input->GetDirection() );
  empty->Allocate();
  typedef itk::ImageRegionIterator< ImageType > ItType;

  ItType it( empty, empty->GetLargestPossibleRegion() ); 

  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	  it.Set( 0 );

    return empty;
  } // CreateEmptyImage(...)


  template < typename TPixelType, int TDims >
  TPixelType
  ImageMax( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
		// Find minimum value in image.
  
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxFilterType;

    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput( input );
    minMaxFilter->Update();

    return minMaxFilter->GetMaximum();
  } // ImageMax(...)


  template < typename TPixelType, int TDims >
  TPixelType
  ImageMin( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
    // Find maximum value in image.
  
		typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxFilterType;

    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput( input );
    minMaxFilter->Update();

    return minMaxFilter->GetMinimum();
  } // ImageMin(...)


  template < typename TPixelType, int TDims >
  int
  ImageCount( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
	 // Count non-zero pixels in an image.
 
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    int count = 0;

    for ( input_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it )
        if ( input_it.Get() != 0  ) 
          ++count;
    
    return count;

  }// ImageCount(...)

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                    //// Operations with Constants /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageUThr( typename itk::Image< TPixelType, TDims >::Pointer input,
             float Uthr)
  {
    // Zero voxels below Uthr.

    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::ImageRegionIterator< ImageType > ItType;

		typename ImageType::Pointer result = ImageType::New();

		result = CreateEmptyImage< TPixelType, TDims>( input );

    ItType input_it(  input, input->GetLargestPossibleRegion() );
    ItType result_it( result, input->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it, ++result_it )
        if ( input_it.Get() > Uthr  ) 
          result_it.Set( 0 );
				else
					result_it.Set( input_it.Get() );
    
    return result;

  }// ImageUThr(...)


	template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageThr( typename itk::Image< TPixelType, TDims >::Pointer input,
             float thr)
  {
		// Zero voxels above thr.
  
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::ImageRegionIterator< ImageType > ItType;

		typename ImageType::Pointer result = ImageType::New();

		result = CreateEmptyImage< TPixelType, TDims>( input );

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, input->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it, ++result_it )
        if ( input_it.Get() < thr  ) 
          result_it.Set( 0 );
				else
					result_it.Set( input_it.Get() );
    
    return result;

  }// ImageThr(...)

  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer 
  ImageAddConstant( typename itk::Image< TPixelType, TDims>::Pointer input,
                    TPixelType constant )
  {
		// Add a constant value to an image.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin();
         !input_it.IsAtEnd(); ++input_it, ++result_it )
      result_it.Set( input_it.Get() + constant );

    return result;
  } // ImageAddConstant( ... )


  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer 
  ImageSubConstant( typename itk::Image< TPixelType, TDims>::Pointer input,
                    TPixelType constant )
  {
    // Subtract a constant from an image.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin(); 
         !input_it.IsAtEnd(); ++input_it, ++result_it )
      result_it.Set( input_it.Get() - constant );

    return result;
  } // ImageSubConstant( ... )
  
  
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMulConstant( typename itk::Image< TPixelType, TDims >::Pointer input,
                    TPixelType constant )
  {
		// Multiply and image by a constant.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin(); 
         !input_it.IsAtEnd(); ++input_it, ++result_it )
      result_it.Set( input_it.Get() * constant );

    return result;
  }  // ImageMulConstant(...)


  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImagePower( typename itk::Image< TPixelType, TDims >::Pointer input,
                    TPixelType constant )
  {
		// Raise image to a constant power.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input_it.GoToBegin(), result_it.GoToBegin(); 
          !input_it.IsAtEnd(); ++input_it, ++result_it )
      result_it.Set( std::pow( input_it.Get(), constant ) );

    return result;
  }  // ImagePower(...)

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                 //// Operations with One Image /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                 //// Operations with Two Images /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
  template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMinImage( typename itk::Image< TPixelType, TDims >::Pointer input1,
	               typename itk::Image< TPixelType, TDims >::Pointer input2)
	{	
		typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input1 );
		
    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input1_it( input1, input1->GetLargestPossibleRegion() );
		ItType input2_it( input2, input2->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input1_it.GoToBegin(), input2_it.GoToBegin(), result_it.GoToBegin(); 
          !input1_it.IsAtEnd(), !input2_it.IsAtEnd();
					++input1_it, ++input2_it, ++result_it )
					result_it.Set( std::min( input1_it.Get(), input2_it.Get() ) );

    return result;
	} // ImageMinImage

  template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMaxImage( typename itk::Image< TPixelType, TDims >::Pointer input1,
	               typename itk::Image< TPixelType, TDims >::Pointer input2)
	{	
		typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims>( input1 );
		
    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input1_it( input1, input1->GetLargestPossibleRegion() );
		ItType input2_it( input2, input2->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );

    for ( input1_it.GoToBegin(), input2_it.GoToBegin(), result_it.GoToBegin(); 
          !input1_it.IsAtEnd(), !input2_it.IsAtEnd();
					++input1_it, ++input2_it, ++result_it )
					result_it.Set( std::max( input1_it.Get(), input2_it.Get() ) );

    return result;
	} // ImageMaxImage

  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageAddImage( typename itk::Image< TPixelType, TDims >::Pointer input1, 
                 typename itk::Image< TPixelType, TDims >::Pointer input2 )
  {
		// Add two images.
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typedef itk::AddImageFilter< ImageType, ImageType, ImageType > AddFilterType;

    typename AddFilterType::Pointer addFilter = AddFilterType::New();

    addFilter->SetInput1( input1 );
    addFilter->SetInput2( input2 );
    addFilter->Update();

    return addFilter->GetOutput();
  } // ImageAddImage(...)

  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageSubImage( typename itk::Image< TPixelType, TDims >::Pointer input1, 
                 typename itk::Image< TPixelType, TDims >::Pointer input2 )
  {
		// Subtract one image from another.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractFilterType;

    typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();

    subFilter->SetInput1( input1 );
    subFilter->SetInput2( input2 );
    subFilter->Update();

    return subFilter->GetOutput();
  } // ImageSubImage(...)


  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMulImage( typename itk::Image< TPixelType, TDims >::Pointer input1,
                 typename itk::Image< TPixelType, TDims >::Pointer input2 )
  {
		//  Multiply one image by another.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typedef itk::MultiplyImageFilter< ImageType, ImageType, ImageType > MultiplyFilterType;

    typename MultiplyFilterType::Pointer mulFilter = MultiplyFilterType::New();

    mulFilter->SetInput1( input1 );
    mulFilter->SetInput2( input2 );
    mulFilter->Update();

    return mulFilter->GetOutput();
  } // ImageMulImage(...)


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                   //// Operations with a Mask /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageFindAndSet( typename itk::Image< TPixelType, TDims >::Pointer input,
                   typename itk::Image< TPixelType, TDims >::Pointer mask,
                   TPixelType replaceValue )
  {
		// Replace any voxel in the input image with a new value,
    // iff the voxel is non-zero in the mask image.
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims >( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );
    ItType mask_it( mask, mask->GetLargestPossibleRegion() );

    for ( mask_it.GoToBegin(), input_it.GoToBegin(), result_it.GoToBegin();
         !input_it.IsAtEnd(); ++input_it, ++result_it, ++mask_it )
      if ( mask_it.Get() != 0 )
        result_it.Set( replaceValue );
      else
        result_it.Set( input_it.Get() );

    return result;
  } // ImageFindAndSet(...)


  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMask( typename itk::Image< TPixelType, TDims >::Pointer input,
             typename itk::Image< TPixelType, TDims >::Pointer mask )
  {
		// Replace any voxel in the input image with a new value,
    // iff the voxel is non-zero in the mask image.

    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims >( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );
    ItType mask_it( mask, mask->GetLargestPossibleRegion() );

    for ( mask_it.GoToBegin(), input_it.GoToBegin(), result_it.GoToBegin();
         !input_it.IsAtEnd(); ++input_it, ++result_it, ++mask_it )
      if ( mask_it.Get() == 0 )
        result_it.Set( 0 );
      else
        result_it.Set( input_it.Get() );

    return result;
  } // ImageMask(...)

  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageFindValueAndSet( typename itk::Image< TPixelType, TDims >::Pointer input,
                        typename itk::Image< TPixelType, TDims >::Pointer mask,
                        TPixelType findValue,
                        TPixelType replaceValue )
  {
		// Replace any voxel in the input image with the replaceValue,
    // iff the voxel has the findValue in the mask image.
  
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TPixelType, TDims >( input );

    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    ItType result_it( result, result->GetLargestPossibleRegion() );
    ItType mask_it( mask, mask->GetLargestPossibleRegion() );

    for ( mask_it.GoToBegin(), input_it.GoToBegin(), result_it.GoToBegin();
         !input_it.IsAtEnd(); ++input_it, ++result_it, ++mask_it )
      if ( mask_it.Get() == findValue )
        result_it.Set( replaceValue );
      else
        result_it.Set( input_it.Get() );

    return result;
  } // ImageFindValueAndSet(...)

template< typename TPixelType, int TDims >
typename itk::Image< TPixelType, TDims >::Pointer
KeepOnAConnectedToB( typename itk::Image< TPixelType, TDims >::Pointer inputA,
                     typename itk::Image< TPixelType, TDims >::Pointer inputB,
                     double accept_value, double reject_value,
                     int min_object_size, int max_object_size )
  {
    typedef  itk::Image< TPixelType, TDims > ImageType;
    typedef  itk::Image< double, TDims > OutputImageType;
	  	std::cout << "Preparing to find connected objects..." << std::endl;
  
	// Create structuring element ("1"s for voxels face-connected to the center of a 3x3x3 cube).
	int sevals[] = {0,0,0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,0,1,0,0,0,0};
	typedef typename itk::Neighborhood< TPixelType, TDims > NeighborhoodType;
	NeighborhoodType se;
	se.SetRadius( 1 );
	int i=0;
	for ( typename itk::Neighborhood< TPixelType,TDims>::Iterator it = se.Begin(); it != se.End(); ++it, ++i )
		*it = sevals[i];

	// Dilate inputB (search space for connected objects)
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, NeighborhoodType > DilateFilterType;
  typename DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
	dilateFilter->SetInput( inputB );
	dilateFilter->SetKernel( se );
	dilateFilter->SetDilateValue( 1 );
	dilateFilter->Update();
	typename ImageType::Pointer inputB_dilate = ImageType::New();
	inputB_dilate = dilateFilter->GetOutput();
	inputB_dilate->DisconnectPipeline();

  // Assign a unique label to each connected object
	typename ImageType::Pointer inputA_lab = ImageType::New();
  
  typedef itk::ScalarConnectedComponentImageFilter< ImageType, OutputImageType > ConnectedFilterType;
  typename ConnectedFilterType::Pointer connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput( inputA );
  connectedFilter->SetBackgroundValue( 0 );
  connectedFilter->FullyConnectedOn();
  
  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelFilterType;
  typename RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  relabelFilter->SetMinimumObjectSize( min_object_size );
  relabelFilter->SetInput( connectedFilter->GetOutput() );
  relabelFilter->Update();

  // Remove background object (typically, this will be the largest object and so could be removed later
  //                           with the removal of large objects, but to be safe, it is removed explicitly here)
  inputA_lab = egitk::ImageFindValueAndSet< TPixelType, TDims >( relabelFilter->GetOutput(), inputA, 0, 0 );
  
	//  Create a lookup vector for each labeled object
	int num_objects = egitk::ImageMax< TPixelType, TDims >( inputA_lab );
  std::vector< int > lookup( num_objects + 1 );
 
	// Make a copy of the labeled image
  typedef itk::ImageRegionIterator< ImageType > ItType;
  typename ImageType::Pointer keep = egitk::CreateEmptyImage< TPixelType, TDims >( inputA_lab );
  ItType it0(  inputA_lab, inputA_lab->GetLargestPossibleRegion() );
  ItType it1( keep, keep->GetLargestPossibleRegion() );
  for ( it0.GoToBegin(), it1.GoToBegin(); !it0.IsAtEnd(); ++it0, ++it1 )
  {
    it1.Set( it0.Get() );
  }

  // Zero any voxels on the labelled image that correspond to zero voxels on the dilated image.
  keep = egitk::ImageFindValueAndSet< TPixelType, TDims >( keep, inputB_dilate, 0, 0 );
  
	std::cout << "First pass over image..." << std::endl;
	
	// First pass over image to find object labels to accept.
  ItType it2( keep, keep->GetLargestPossibleRegion() );
  for ( it2.GoToBegin(); !it2.IsAtEnd(); ++it2  )
	{
		// Store the labels of any hyperintensity not rejected after
		//  masking with the dilated thresholded wmt image in the lookup vector.
		int value = it2.Get();
	  if ( value != 0 ) 
      lookup.at( value ) = value; 
	}

  for ( int i=0; i < num_objects; ++i )
  {
    if ( relabelFilter->GetSizeOfObjectsInPixels()[i] > max_object_size )
    {
      std::cout << "Removing object #" << i + 1 << " with " << relabelFilter->GetSizeOfObjectsInPixels()[i] << " voxels" << std::endl;
      lookup.at( i + 1 ) = 0;
    }
  }

	std::cout << "Second pass over image..." << std::endl;
	// Second pass over image to remove objects not accepted in the first pass.
	ItType it3( inputA_lab, inputA_lab->GetLargestPossibleRegion() );
	for ( it3.GoToBegin(); !it3.IsAtEnd(); ++it3 )
	{
		int value = it3.Get();
    if ( value != 0 )
    {
      if ( lookup[ value ] == 0 )
			  it3.Set( reject_value ); // reject
      if ( lookup[ value ] != 0 )
        it3.Set( accept_value );  // accept
    }

	}
  return inputA_lab;
	  
  } //KeepOnAConnectedToB
} // namespace egitk


#endif
