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

#ifndef _itk_imgop_util_templates_h_
#define _itk_imgop_util_templates_h_

#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryMedianImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkThresholdImageFilter.h"

#include <iostream>


namespace itk_iops // itk image operations = itk_iops
{

    template< typename TInputPixelType, typename TOutputPixelType, typename TFilter >
  typename itk::Image< TOutputPixelType, 3 >::Pointer
  FilterVolAsSlice( typename itk::Image< TInputPixelType, 3>::Pointer inputVol, TFilter inputFilter )
	{
	
    typedef itk::Image< TInputPixelType, 3 >InputImageType;
    typedef itk::Image< TOutputPixelType, 3 >OutputImageType;
		
    typedef itk::Image< TInputPixelType, 2>InputImageType2D;
    typedef itk::Image< TOutputPixelType, 2 >OutputImageType2D;
		
    // Create output image with the same dimensions and spacing as input volume
    typename OutputImageType::Pointer outputVol = OutputImageType::New();
    outputVol->SetLargestPossibleRegion( inputVol->GetLargestPossibleRegion() );
    outputVol->SetBufferedRegion( inputVol->GetLargestPossibleRegion() );
    outputVol->SetRequestedRegion( inputVol->GetLargestPossibleRegion() );
    outputVol->SetSpacing( inputVol->GetSpacing() );
    outputVol->Allocate();

		// Count slices in volume
		int numslices = inputVol->GetLargestPossibleRegion().GetSize().operator []( 2 );
		
		// Filter all slices 
		for ( int slicenum = 0; slicenum < numslices; ++slicenum )  // filter all slices: loop start
		{
			/////////////////////////////////////////////////////////////
			///// Extract single slice from input volume and filter /////
			/////////////////////////////////////////////////////////////

			typedef itk::ExtractImageFilter< InputImageType, InputImageType2D > ExtractFilterType;

			typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

			typename InputImageType::IndexType extractStart = {{0, 0, slicenum}};
	
			typename InputImageType::SizeType extractSize = inputVol->GetLargestPossibleRegion().GetSize();
			extractSize[2] = 0; //(size = 0: signals a conversion from 3D to 2D for extactFilter)
			
			typename InputImageType::RegionType extractRegion;
			extractRegion.SetSize( extractSize );
			extractRegion.SetIndex( extractStart );
			
			extractFilter->SetInput( inputVol );
			extractFilter->SetExtractionRegion( extractRegion );
			extractFilter->Update();

			inputFilter->SetInput( extractFilter->GetOutput() );
			inputFilter->Update();

			typename InputImageType2D::Pointer tempslice = InputImageType2D::New();
			tempslice = inputFilter->GetOutput();

			// Output progress message
      std::cout << inputFilter->GetNameOfClass() <<": Processing slice " << slicenum + 1 << std::endl;
		

			/////////////////////////////
			///// Skip blank slices /////
			/////////////////////////////
	/*		
			typedef itk::MinimumMaximumImageFilter< InputImageType2D > MinMaxFilterType;
	
			MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();

			minMaxFilter->SetInput( tempslice );
			minMaxFilter->Update();

			if ( minMaxFilter->GetMaximum() < 1 )  continue;
	*/					

			//////////////////////////////////////////////////////
			///// Copy filtered single slice to input volume /////
			//////////////////////////////////////////////////////

      // Define insert region on output vol where the filtered 2D should be copied to

      typename OutputImageType::IndexType insertStart = {{0, 0, slicenum}};
	
			typename OutputImageType::SizeType insertSize = inputVol->GetLargestPossibleRegion().GetSize();
			insertSize[2] = 1;
			
			typename OutputImageType::RegionType insertRegion;
			insertRegion.SetSize( insertSize );
			insertRegion.SetIndex( insertStart );
			

      // Copy data
			typedef itk::ImageRegionIterator< OutputImageType > It3dType;
			typedef itk::ImageRegionConstIterator< InputImageType2D > It2dType;

			It3dType it3d( outputVol, insertRegion );
			It2dType it2d( tempslice, tempslice->GetLargestPossibleRegion() );

      for ( it2d.GoToBegin(), it3d.GoToBegin(); !it2d.IsAtEnd(); ++it3d, ++it2d )
        it3d.Set( it2d.Get() );  
        
		}  // filter all slices: loop end

		return outputVol;
	}; // FilterVolAsSlice(...)


    template < typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  ThresholdImage( typename itk::Image< TPixelType, TDims >::Pointer input, 
                      TPixelType lowerThresh,
                      TPixelType upperThresh )
  {
    typedef itk::Image< TPixelType, TDims > ImageType;
    
    typedef itk::ThresholdImageFilter< ImageType > ThresholdFilterType;

    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();

    thresholdFilter->SetInput( input );
    thresholdFilter->SetLower( lowerThresh );
    thresholdFilter->SetUpper( upperThresh );
    thresholdFilter->Update();

    return thresholdFilter->GetOutput();
  }
  


	template< typename TInputPixelType, typename TOutputPixelType, int TDims >
  typename itk::Image< TOutputPixelType, TDims>::Pointer
  ThresholdImageToBinary( typename itk::Image< TInputPixelType, TDims>::Pointer input,
                        int insideValue,  int outsideValue,
                        TInputPixelType lowerThresh, 
                        TInputPixelType upperThresh )
	{
    typedef itk::Image< TInputPixelType, TDims > InputImageType;
    typedef itk::Image< TOutputPixelType, TDims > OutputImageType;
    
    typedef itk::BinaryThresholdImageFilter< InputImageType, OutputImageType > BinaryFilterType;

		typename BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();

		binaryFilter->SetInput ( input );
		binaryFilter->SetInsideValue( insideValue );
		binaryFilter->SetOutsideValue( outsideValue );
		binaryFilter->SetLowerThreshold( lowerThresh );
		binaryFilter->SetUpperThreshold( upperThresh );
		binaryFilter->Update();

    std::cout << std::endl << binaryFilter->GetNameOfClass() << ": Processing..." << std::endl;

		return binaryFilter->GetOutput();
	}  // ThresholdImageToBinary(...)


  // 2D slices of 3D volume processed seperately
  
  template< typename InputPixelType, typename OutputPixelType >
  typename itk::Image< InputPixelType, 3>::Pointer
  FillHolesInBinaryVol2D( typename itk::Image< InputPixelType, 3>::Pointer input )
  {
    //////////////////////////////////////////////////////////////////////////
    //  Fills holes in a binary image, slice by slice.
    //    Holes are defined as any area in the image that can not be reached
    //    by flood-filling from the far edge of the image
    //////////////////////////////////////////////////////////////////////////

    typedef itk::Image< InputPixelType, 3 > InputImageType;
    typedef itk::Image< OutputPixelType, 3 > OutputImageType;

    typedef itk::Image< InputPixelType, 2 > InputImageType2D;
    typedef itk::Image< OutputPixelType, 2 > OutputImageType2D;

    typedef itk::ConnectedThresholdImageFilter< InputImageType2D, OutputImageType2D > ConnectedFilterType;
  
    itk::Index< 2 > index = {{0, 0}};  // set start seed as the top corner of every slice

    typename ConnectedFilterType::Pointer connectedFilter = ConnectedFilterType::New();
    connectedFilter->SetSeed( index );
    connectedFilter->SetLower( 0 );
    connectedFilter->SetUpper( 0 );
    connectedFilter->SetReplaceValue( 1 );

    // Flood fill each slice from top corner:
    //   Resulting image will have have 1s where flood operation took place,
    //   and 0s for both the object and any holes in the object

//!!!!!!!!!!! create new image here !!!!!!!!!!!!!!!!!!!!!!!!
    input = FilterVolAsSlice<InputPixelType, InputPixelType>( input, connectedFilter ); 

    // Return binary image from the 0s in the flood filled image
    return ThresholdImageToBinary< InputPixelType, OutputPixelType, 3 >( input, 1, 0, 0, 0 );

  } // FillHolesInBinaryVol2D(...)




    template< typename TPixelType, int TDims >
  typename itk::Image< TPixelType, TDims >::Pointer
  DuplicateImage( typename itk::Image< TPixelType, TDims>*  input )
  {
    typedef typename itk::Image< TPixelType, TDims > ImageType;

    typedef itk::ImageDuplicator< ImageType > DuplicatorType;

    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();

    duplicator->SetInputImage( input );
    duplicator->Update();

    return duplicator->GetOutput();
  } // DupliateImage(...)



    template< typename TInputPixelType, typename TOutputPixelType, int TDims >
  typename itk::Image< TOutputPixelType, TDims >::Pointer
  RelabelConnectedComponents( typename itk::Image< TInputPixelType, TDims>::Pointer input, int minObjectSize )
  {
    typedef itk::Image< TInputPixelType, TDims > InputImageType;
    typedef itk::Image< double, TDims > OutputImageType;


  
  typedef itk::ScalarConnectedComponentImageFilter< InputImageType, OutputImageType > ConnectedFilterType;
  typename ConnectedFilterType::Pointer connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput( input );
  connectedFilter->SetBackgroundValue( 0 );
  connectedFilter->SetFullyConnected(0);
  
  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelFilterType;
  typename RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  relabelFilter->SetMinimumObjectSize( minObjectSize );
  relabelFilter->SetInput( connectedFilter->GetOutput() );
  relabelFilter->Update();

    return relabelFilter->GetOutput(); 

  } // RelabelConnectedComponents(...)


  template< typename TInputPixelType, typename TOutputPixelType >
  typename itk::Image< TOutputPixelType, 3 >::Pointer
  MedianFilterVol2D( typename itk::Image< TInputPixelType, 3 >::Pointer inputVol, int radius )
  {
    typedef itk::Image< TInputPixelType, 3 > InputImageType;
    typedef itk::Image< TOutputPixelType, 3> OutputImageType;
    
    typedef itk::Image< TInputPixelType, 2 > InputImageType2D;
    typedef itk::Image< TOutputPixelType, 2> OutputImageType2D;
    
    typedef itk::BinaryMedianImageFilter< InputImageType2D, OutputImageType2D > MedianFilterType;

    typename MedianFilterType::Pointer medianFilter = MedianFilterType::New();

    typename InputImageType2D::SizeType neighRadius = {{radius,radius}};

    medianFilter->SetRadius( neighRadius );
    medianFilter->SetForegroundValue( 1 );
    medianFilter->SetBackgroundValue( 0 );

    
    return FilterVolAsSlice< TInputPixelType, TOutputPixelType >( inputVol, medianFilter );

  } //MedianFilterVol2D(...)




} // namespace itk_iops



#endif
