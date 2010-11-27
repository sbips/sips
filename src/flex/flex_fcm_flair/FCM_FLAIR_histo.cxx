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

#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"

#include "FCM_FLAIR_histo.h"

itk::Image< double, 2>::Pointer
fl_hist::ExtractSliceFromVol( itk::Image< double, 3 >::Pointer input, int slicenum )
{
  typedef itk::Image< double, 3 > ImageType;
  typedef itk::Image< double, 2 > ImageType2D;

  typedef itk::ExtractImageFilter< ImageType, ImageType2D > ExtractFilterType;

  ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

  extractFilter->SetInput( input );

  ImageType::IndexType extractStart = {{0, 0, slicenum}};
	
	ImageType::SizeType extractSize = input->GetLargestPossibleRegion().GetSize();
  extractSize[2] = 0; //(size = 0: signals a conversion from 3D to 2D for extactFilter)
			
	ImageType::RegionType extractRegion;
	extractRegion.SetSize( extractSize );
	extractRegion.SetIndex( extractStart );
			
	extractFilter->SetInput( input );
	extractFilter->SetExtractionRegion( extractRegion );
	extractFilter->Update();

  return extractFilter->GetOutput();
}


std::vector<int>
fl_hist::histo( itk::Image< double, 3>::Pointer flInput, itk::Image<double, 3>::Pointer flfcInput )
{
  int numslices = flInput->GetLargestPossibleRegion().GetSize( 2 );
  
  std::vector< int > result;//( numslices, 0 );
  
  typedef itk::Image< double, 2 > ImageType2D;
  
  typedef itk::Statistics::ScalarImageToHistogramGenerator<
           itk::Image< double, 2> > HistogramGeneratorType;

  HistogramGeneratorType::Pointer flHistoGen =
                       HistogramGeneratorType::New();
  HistogramGeneratorType::Pointer flfcHistoGen =
                       HistogramGeneratorType::New();

  ImageType2D::Pointer fl = ImageType2D::New();
  ImageType2D::Pointer flfc = ImageType2D::New();

for ( int slicenum = 0; slicenum < numslices; ++slicenum )
  {

    fl = ExtractSliceFromVol( flInput, slicenum );
    flfc = ExtractSliceFromVol( flfcInput, slicenum );

    if ( itk_imath::ImageMax< double, 2 >( fl ) == 0 ) 
    {
      result.push_back( 0 );
      continue;
    }

    double flMax = itk_imath::ImageMax< double, 2 >( fl );
    

    flHistoGen->SetInput( fl );
    flHistoGen->SetNumberOfBins( flMax + 1 );
    flHistoGen->SetMarginalScale( 1 );
    
    flHistoGen->Compute();


    flfcHistoGen->SetInput( flfc );
    flfcHistoGen->SetNumberOfBins( flMax + 1 );
    flfcHistoGen->SetMarginalScale( 1 );
    
    flfcHistoGen->Compute();
    

    typedef HistogramGeneratorType::HistogramType HistogramType;
    
    const HistogramType* flHisto = flHistoGen->GetOutput();
    const HistogramType* flfcHisto = flfcHistoGen->GetOutput();

    HistogramType::ConstIterator flIt = flHisto->Begin();
    HistogramType::ConstIterator flItEnd = flHisto->End();
    HistogramType::ConstIterator flfcIt = flfcHisto->Begin();
    HistogramType::ConstIterator flfcItEnd = flfcHisto->End();

    std::vector< double > flHistoVector;
    std::vector< double > flfcHistoVector;

   
    while ( flIt != flItEnd )
    {
     // std::cout << flIt.GetFrequency() << " " <<flfcIt.GetFrequency() << std::endl;
      flHistoVector.push_back( flIt.GetFrequency() );
      flfcHistoVector.push_back( flfcIt.GetFrequency() );
      ++flIt;  ++ flfcIt;
    }
   

    unsigned int bin = flMax;
    for ( int i = flMax; i >= 0; --i )
    {
      if ( int(flfcHistoVector[i]) != int(flHistoVector[i])  ) break;
      --bin;
      
    }
    
    result.push_back( bin );
  } // for 


  return result;

} //histo


itk::Image< double, 3 >::Pointer
fl_hist::GetFlairLesions( itk::Image< double, 3>::Pointer fl, std::vector< int >& flThresh )
{
  
  typedef itk::Image< double, 3 > ImageType;

    typedef itk::ImageRegionIterator< ImageType > ItType;

    std::vector<int>::iterator it;
    int slicenum = 0;
    for ( it = flThresh.begin(); it != flThresh.end(); ++it )
    {

      ImageType::IndexType start = {{0, 0, slicenum}};

      ImageType::SizeType size = fl->GetLargestPossibleRegion().GetSize();
	    size[2] = 1; 

		  ImageType::RegionType region;
		  region.SetSize( size );
		  region.SetIndex( start );
  
      fl->SetLargestPossibleRegion( region );

      ItType flIt( fl, fl->GetLargestPossibleRegion() );
      

      int thresh = *it;
      std ::cout << "Lesion threshold for slice " << slicenum+1 << ": " << thresh << std::endl;
      
      for ( flIt.GoToBegin(); !flIt.IsAtEnd(); ++flIt )
        if ( flIt.Get() < thresh )
          flIt.Set( 0 );

      ++slicenum;

    }

    return fl;
}
