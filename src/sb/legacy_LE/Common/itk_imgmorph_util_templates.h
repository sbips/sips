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

#ifndef _itk_imgmorph_templates_h_
#define _itk_imgmorph_templates_h_


#include "itkImage.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"


#include "itk_imgops_util_templates.h" // itk_iops namespace

namespace itk_imorph
{
   template< typename TPixelType, int TDims >
  typename itk::BinaryBallStructuringElement< TPixelType, TDims>
  CreateBallSE( TPixelType radius )
  {
      typedef typename itk::BinaryBallStructuringElement< TPixelType, TDims > BinaryBallStructSEType;
      BinaryBallStructSEType binaryBallSE;
      binaryBallSE.SetRadius( radius );
      binaryBallSE.CreateStructuringElement();
      return binaryBallSE;
  } // CreateBallSE(...)
  
   template< typename TPixelType, int TDims >
  typename itk::BinaryCrossStructuringElement< TPixelType, TDims>
  CreateCrossSE( TPixelType radius )
  {
      typedef typename itk::BinaryCrossStructuringElement< TPixelType, TDims > BinaryCrossStructSEType;
      BinaryCrossStructSEType binaryCrossSE;
      binaryCrossSE.SetRadius( radius );
      binaryCrossSE.CreateStructuringElement();
      return binaryCrossSE;
  } // CreateCrossSE(...)

  // 3D Dilation.

    template< typename TInputPixelType, typename TOutputPixelType, 
              int TDims, typename TStructElementType >
  typename itk::Image< TOutputPixelType, TDims >::Pointer
  DilateBinaryImage( typename itk::Image< TInputPixelType, TDims >::Pointer input, TStructElementType structElement )
  {
    typedef itk::Image< TInputPixelType, TDims > InputImageType;
    typedef itk::Image< TOutputPixelType, TDims > OutputImageType;
    typedef itk::BinaryDilateImageFilter< InputImageType, OutputImageType, TStructElementType > BinaryDilateFilterType;

    typename BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();

    dilateFilter->SetInput( input );
    dilateFilter->SetKernel( structElement );
    dilateFilter->SetDilateValue( 1 );
    dilateFilter->Update();

    std::cout << std::endl << dilateFilter->GetNameOfClass() << ": Processing..." << std::endl;

    return dilateFilter->GetOutput();
  } // DilateBinaryImage(...)



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 // 2D processing of 3D volume
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   template< typename TInputImageType, typename TOutputImageType, typename TStructElementType2D >
  typename TOutputImageType::Pointer
  DilateBinaryVol2D( typename TInputImageType::Pointer input, TStructElementType2D structElement )
  {

    typedef typename TInputImageType::PixelType InputPixelType;
    typedef typename TOutputImageType::PixelType OutputPixelType;

    typedef typename itk::Image< InputPixelType, 2 >InputImageType2D;
    typedef typename itk::Image< OutputPixelType, 2 >OutputImageType2D;

    typedef typename itk::BinaryDilateImageFilter< InputImageType2D, OutputImageType2D, 
                                               TStructElementType2D > BinaryDilateFilterType2D;

    typename BinaryDilateFilterType2D::Pointer dilateFilter = BinaryDilateFilterType2D::New();
    dilateFilter->SetKernel( structElement );
    dilateFilter->SetForegroundValue( 1 );
    dilateFilter->SetBackgroundValue( 0 );
    dilateFilter->SetDilateValue( 1 );

    return itk_iops::FilterVolAsSlice< InputPixelType, OutputPixelType >( input, dilateFilter );
     
    } // DilateBinaryVol2D(...)


   template< typename TInputImageType, typename TOutputImageType, typename TStructElementType2D >
  typename TOutputImageType::Pointer
  ErodeBinaryVol2D( typename TInputImageType::Pointer input, TStructElementType2D structElement )
  {
    typedef typename TInputImageType::PixelType InputPixelType;
    typedef typename TOutputImageType::PixelType OutputPixelType;

    typedef typename itk::Image< InputPixelType, 2 >InputImageType2D;
    typedef typename itk::Image< OutputPixelType, 2 >OutputImageType2D;

    typedef itk::BinaryErodeImageFilter< InputImageType2D, OutputImageType2D, TStructElementType2D > BinaryErodeFilterType2D;

    typename BinaryErodeFilterType2D::Pointer erodeFilter = BinaryErodeFilterType2D::New();
    erodeFilter->SetKernel( structElement );
    erodeFilter->SetForegroundValue( 1 );
    erodeFilter->SetBackgroundValue( 0 );
    erodeFilter->SetErodeValue( 1 );

    return itk_iops::FilterVolAsSlice< InputPixelType, OutputPixelType >( input, erodeFilter );
     
  } // ErodeBinaryVol2D(...)


   template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  GetObjectPerimeter( typename itk::Image< TPixelType, TDims>* input )
  {
   
    typedef typename itk::Image< TPixelType, TDims > ImageType;
    typedef typename itk::Image< TPixelType, 2 > ImageType2D;
   
    typename ImageType::Pointer brain = ImageType::New();
    
    brain = itk_iops::ThresholdImageToBinary< TPixelType, TPixelType, TDims>( input, 1, 0, 1, 7 );
    brain = itk_iops::FillHolesInBinaryVol2D< TPixelType, TPixelType >( brain );

    typename ImageType::Pointer temp = ImageType::New();
    temp = itk_iops::DuplicateImage< TPixelType, TDims>( &*brain );
  
    temp = itk_imorph::ErodeBinaryVol2D< ImageType, ImageType >( temp, itk_imorph::CreateBallSE< double, 2 >( 1 ) );
    temp = itk_imorph::ErodeBinaryVol2D< ImageType, ImageType >( temp, itk_imorph::CreateBallSE< double, 2 >( 1 ) );

    typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType > SubtractFilterType;

    typename SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();

    subtractFilter->SetInput1( brain );
    subtractFilter->SetInput2( temp );
    subtractFilter->Update();

    return subtractFilter->GetOutput();
  } // GetObjectPerimeter(...)


} // namespace itk_imorph
#endif
