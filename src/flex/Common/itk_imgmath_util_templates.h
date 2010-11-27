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


namespace itk_imath
{
  // Create an empty image with the same dimensions as the input image.
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  CreateEmptyImage( typename itk::Image< TPixelType, TDims >::Pointer input)
  {    typename itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typename ImageType::Pointer empty = ImageType::New();
   
    empty->SetLargestPossibleRegion( input->GetLargestPossibleRegion() );
    empty->SetBufferedRegion( input->GetLargestPossibleRegion() );
    empty->SetRequestedRegion( input->GetLargestPossibleRegion() );
    empty->SetSpacing( input->GetSpacing() );
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
	empty->SetDirection( cosines );
	empty->Allocate();

    return empty;
  } // CreateEmptyImage(...)


  // Find minimum value in image.
  template < typename TPixelType, int TDims >
  TPixelType
  ImageMax( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxFilterType;

    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput( input );
    minMaxFilter->Update();

    return minMaxFilter->GetMaximum();
  } // ImageMax(...)


  // Find maximum value in image.
  template < typename TPixelType, int TDims >
  TPixelType
  ImageMin( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::MinimumMaximumImageFilter< ImageType > MinMaxFilterType;

    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

    minMaxFilter->SetInput( input );
    minMaxFilter->Update();

    return minMaxFilter->GetMinimum();
  } // ImageMin(...)


  // Count non-zero pixels in an image.
  template < typename TPixelType, int TDims >
  int
  ImageCount( typename itk::Image< TPixelType, TDims >::Pointer input)
  {
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

  // Zero voxels below Uthr.
  template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageUThr( typename itk::Image< TPixelType, TDims >::Pointer input,
             float Uthr)
  {
    typedef itk::Image< TPixelType, TDims > ImageType;
    typedef itk::ImageRegionIterator< ImageType > ItType;

    ItType input_it( input, input->GetLargestPossibleRegion() );
    int count = 0;

    for ( input_it.GoToBegin(); !input_it.IsAtEnd(); ++input_it )
        if ( input_it.Get() > Uthr  ) 
          input_it.Set( 0 );
    
    return input;

  }// ImageUthr(...)


  // Add a constant value to an image.
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer 
  ImageAddConstant( typename itk::Image< TPixelType, TDims>::Pointer input,
                    TPixelType constant )
  {
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


  // Subtract a constant from an image.
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer 
  ImageSubConstant( typename itk::Image< TPixelType, TDims>::Pointer input,
                    TPixelType constant )
  {
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
  

  // Multiply and image by a constant.
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImageMulConstant( typename itk::Image< TPixelType, TDims >::Pointer input,
                    TPixelType constant )
  {
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


  // Raise image to a constant power.
  template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ImagePowConstant( typename itk::Image< TPixelType, TDims >::Pointer input,
                    TPixelType constant )
  {
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
  }  // ImagePowConstant(...)

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                 //// Operations with Two Images /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

  // Add two images.
  template< typename TInputPixelType1, typename TInputPixelType2, typename TOutputPixelType, int TDims >
  typename itk::Image< TOutputPixelType, TDims >::Pointer
  ImageAddImage( typename itk::Image< TInputPixelType1, TDims >::Pointer input1, 
                 typename itk::Image< TInputPixelType2, TDims >::Pointer input2 )
  {
    typedef typename itk::Image< TInputPixelType1, TDims > InputImageType1;
    typedef typename itk::Image< TInputPixelType2, TDims > InputImageType2;
    typedef typename itk::Image< TOutputPixelType, TDims > OutputImageType;
    typedef itk::AddImageFilter< InputImageType1, InputImageType2, OutputImageType > AddFilterType;

    typename AddFilterType::Pointer addFilter = AddFilterType::New();

    addFilter->SetInput1( input1 );
    addFilter->SetInput2( input2 );
    addFilter->Update();

    return addFilter->GetOutput();
  } // ImageAddImage(...)


  // Subtract one image from another.
  template< typename TInputPixelType1, typename TInputPixelType2, typename TOutputPixelType, int TDims >
  typename itk::Image< TOutputPixelType, TDims >::Pointer
  ImageSubImage( typename itk::Image< TInputPixelType1, TDims >::Pointer input1, 
                 typename itk::Image< TInputPixelType2, TDims >::Pointer input2 )
  {
    typedef typename itk::Image< TInputPixelType1, TDims > InputImageType1;
    typedef typename itk::Image< TInputPixelType2, TDims > InputImageType2;
    typedef typename itk::Image< TOutputPixelType, TDims > OutputImageType;

    typedef itk::SubtractImageFilter< InputImageType1, InputImageType2, OutputImageType > SubtractFilterType;

    typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();

    subFilter->SetInput1( input1 );
    subFilter->SetInput2( input2 );
    subFilter->Update();

    return subFilter->GetOutput();
  } // ImageSubImage(...)


  //  Multiply one image by another.
  template< typename TInputPixelType1, typename TInputPixelType2, 
            typename TOutputPixelType, int TDims >
  typename itk::Image< TOutputPixelType, TDims >::Pointer
  ImageMulImage( typename itk::Image< TInputPixelType1, TDims >::Pointer input1,
                 typename itk::Image< TInputPixelType2, TDims >::Pointer input2 )
  {
    typedef typename itk::Image< TInputPixelType1, TDims > InputImageType1;
    typedef typename itk::Image< TInputPixelType2, TDims > InputImageType2;
    typedef typename itk::Image< TOutputPixelType, TDims > OutputImageType;
    typedef itk::MultiplyImageFilter< InputImageType1, InputImageType2, OutputImageType >
                                                                      MultiplyFilterType;

    typename MultiplyFilterType::Pointer mulFilter = MultiplyFilterType::New();

    mulFilter->SetInput1( input1 );
    mulFilter->SetInput2( input2 );
    mulFilter->Update();

    return mulFilter->GetOutput();
  } // ImageMulImage(...)


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
                   //// Operations with a Mask /////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//


  // Replace any voxel in the input image with a new value,
  // iff the voxel is non-zero in the mask image.
  template< typename TInputPixelType, typename TMaskPixelType, int TDims >
  typename itk::Image< TInputPixelType, TDims >::Pointer
  ImageFindAndSet( typename itk::Image< TInputPixelType, TDims >::Pointer input,
                   typename itk::Image< TMaskPixelType, TDims >::Pointer mask,
                   TInputPixelType replaceValue )
  {
    typedef typename itk::Image< TInputPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TInputPixelType, TDims >( input );

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


    // Replace any voxel in the input image with a new value,
  // iff the voxel is non-zero in the mask image.
  template< typename TInputPixelType, typename TMaskPixelType, int TDims >
  typename itk::Image< TInputPixelType, TDims >::Pointer
  ImageMask( typename itk::Image< TInputPixelType, TDims >::Pointer input,
                   typename itk::Image< TMaskPixelType, TDims >::Pointer mask )
  {
    typedef typename itk::Image< TInputPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TInputPixelType, TDims >( input );

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

  // Replace any voxel in the input image with the replaceValue,
  // iff the voxel has the findValue in the mask image.
  template< typename TInputPixelType, typename TMaskPixelType, int TDims >
  typename itk::Image< TInputPixelType, TDims >::Pointer
  ImageFindValueAndSet( typename itk::Image< TInputPixelType, TDims >::Pointer input,
                        typename itk::Image< TMaskPixelType, TDims >::Pointer mask,
                        TMaskPixelType findValue,
                        TInputPixelType replaceValue )
  {
    typedef typename itk::Image< TInputPixelType, TDims > ImageType;

    typename ImageType::Pointer result = ImageType::New();

    result = CreateEmptyImage< TInputPixelType, TDims >( input );

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


} // namespace itk_imath


#endif
