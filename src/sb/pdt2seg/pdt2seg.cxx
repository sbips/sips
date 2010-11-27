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

#include <itkImage.h>
#include <itkImageFileReader.h>

#include "Common/itk_imgio_util_templates.h"
#include "Common/itk_imgmath_util_templates.h"
#include "egRunningStats.h"

#include <numeric>
#include <string>
#include <vector>
#include <iostream>

/*************************************************************************/
/** GLOBAL FUNCTIONS. */
/*************************************************************************/

int Usage( char* argv[] )
{
  std::cout << std::endl << "PD/T2 Hyperintensity Segmentation Tool" << std::endl << std::endl;
  std::cout << "Usage" << std::endl;
  std::cout << argv[0] << std::endl << "   <PD> <T2> <PDt> <T2t> <win> <OUT>" << std::endl << std::endl;;
  std::cout << " <PD>:  (PD, masked, intensity inhomogeneity corrected, extension req'd)" << std::endl;
  std::cout << " <T2>:  (T2, masked, intensity inhomogeneity corrected, extension req'd)" << std::endl;
  std::cout << " <PDt>: (PD threshold, floating point value)" << std::endl;
  std::cout << " <T2t>: (T2 threshold, floating point value)" << std::endl;
  std::cout << " <win>: (Size of window used to smooth the histogram, integer value)" << std::endl;
  std::cout << " <OUT>: (Output image filename, suffix will be added, but extension req'd)" << std::endl;
  std::cout << std::endl;
  std::cout << "Output" << std::endl;
  std::cout << " <OUT_pdseg>:   (PD segmentation)" << std::endl; 
  std::cout << " <OUT_pdt2seg>: (Combined PD/T2 segmentation)" << std::endl; 
  std::cout << " <OUT_t2seg>:   (T2 segmentation)" << std::endl;
  std::cout << std::endl;
  std::cout << "Algortihm" << std::endl;
  std::cout << " -PD lesion threshold = PD mean NWM + (PD adj. max - PD adj. min)*PDt" << std::endl;
  std::cout << " -T2 lesion threshold = T2 histo max index + (T2 FWHM - T2 histo max index)*T2t" << std::endl;
  std::cout << std::endl;
  return 1;
}

template < typename TPixelType, int TDim >
std::vector< float > GetSegParameters( typename itk::Image< TPixelType, TDim >::Pointer pd, 
                                                             typename itk::Image< TPixelType, TDim >::Pointer t2,
                                                             float histothresh, int histowindow )
{
  typedef typename itk::Image< TPixelType, TDim > ImageType;

  /** Create PD and T2 histograms. */
  float t2max = egitk::ImageMax< TPixelType, TDim >( t2 );
  float pdmax = egitk::ImageMax< TPixelType, TDim >( pd );
  float t2count = 0;
  float pdcount = 0;

  std::vector< float > t2histo( t2max + 1, 0 );
  std::vector< float > pdhisto( pdmax + 1, 0 );

  typedef itk::ImageRegionIterator< ImageType > ItType; 
  
  ItType t2It( t2, t2->GetLargestPossibleRegion() );
  ItType pdIt( pd, pd->GetLargestPossibleRegion() );

  for ( t2It.GoToBegin(), pdIt.GoToBegin(); !t2It.IsAtEnd(); ++t2It, ++pdIt )
  {
    if ( t2It.Get() != 0 )
    {
      t2histo[ t2It.Get() ]++;
      ++t2count;
    }
    if ( pdIt.Get() != 0 )
    {
      pdhisto[ pdIt.Get() ]++;
      ++pdcount;
    }
  }

  
  /** Trim PD and T2 histograms.
      Remove large and small outlier values from the histogram.
      (by looking for bins that contain less than some small proportion of 
       total brain voxels.) */
  std::vector< unsigned int> t2trim;
  for ( unsigned int i = 0; i < t2histo.size(); ++i )
  {
    if ( t2histo[i] > t2count*histothresh )
    {
      t2trim.push_back( i );
    }
  }
 
  std::vector< unsigned int> pdtrim;
  for ( unsigned int i = 0; i < pdhisto.size(); ++i )
  {
    if ( pdhisto[i] > pdcount*histothresh )
    {
      pdtrim.push_back( i );
    }
  } 

  /** Smooth T2 histogram. */
  std::vector< float > t2histoSmoothed( t2max + 1, 0 );

  for ( unsigned int i = 1 + histowindow;  i < t2histo.size()-1 - histowindow; ++i )
  {
    std::vector< float >::iterator start = t2histo.begin() + i - histowindow;
    std::vector< float >::iterator end =   t2histo.begin() + i + histowindow + 1; // 1 past the desired end
    t2histoSmoothed[i] = std::accumulate( start, end, 0.0 );
    t2histoSmoothed[i] = t2histoSmoothed[i] / ( histowindow*2 + 1 );
  }

  /** Find index and value of bin containing the largest number of voxels. */
  unsigned histoMaxIndex = 0;
  float    histoMaxValue = 0;
  for ( unsigned int i = 1; i < t2histoSmoothed.size()-1; ++i )
  {
    if ( t2histoSmoothed[i] > histoMaxValue )
    {
      histoMaxIndex = i;
      histoMaxValue = t2histoSmoothed[i];
    }
  }

  std::vector<float> segParameters;
  segParameters.push_back( histoMaxIndex );
  segParameters.push_back( histoMaxValue );

  /** Find the right FWHM. */
  unsigned int t2cut = 0;
  for ( unsigned int i = histoMaxIndex; i < t2histoSmoothed.size()-1; ++i )
  {
    if ( t2histoSmoothed[i] < histoMaxValue * 0.5 )
    {
      t2cut = i;
      break;
    }
  }

  /** Find small intensity window around the peak of the smoothed histogram.
      This window will be used to estimate the mean of normal appearing WM
      on the PD image. */
  unsigned int t2cut4pd = 0;
  for ( unsigned int i = histoMaxIndex; i < t2histoSmoothed.size()-1; ++i )
  {
    if ( t2histoSmoothed[i] < histoMaxValue * 0.995 )
    {
      t2cut4pd = i;
      break;
    }
  }

  t2cut4pd = t2cut4pd - histoMaxIndex;


  /** Find mean of normal appearing WM on the PD. */
  eg::RunningStats runningStats;

  unsigned int pdminTrimmed = pdtrim[0];
  unsigned int pdmaxTrimmed =  pdtrim[ pdtrim.size()-1 ];
  for ( t2It.GoToBegin(), pdIt.GoToBegin(); !t2It.IsAtEnd(); ++t2It, ++pdIt )
  {
    unsigned int lcut = histoMaxIndex - t2cut4pd;
    unsigned int ucut = histoMaxIndex + t2cut4pd;

    if ( t2It.Get() > lcut && t2It.Get() < ucut )
    {
      if ( pdIt.Get() > pdmaxTrimmed || pdIt.Get() < pdminTrimmed )
        continue;
      runningStats.Push( pdIt.Get() );
    } 
  }

  float pdmean = runningStats.Mean();

  segParameters.push_back( t2cut );
  segParameters.push_back( t2cut4pd );
  segParameters.push_back( pdminTrimmed );
  segParameters.push_back( pdmaxTrimmed );
  segParameters.push_back( pdmean );

  /** Return histogram parameters required for segmentation. */
  return segParameters;
}

/*************************************************************************/
/** MAIN FUNCTION. */
/*************************************************************************/

int main( int argc, char* argv[] )
{
  /*************************************************************************/
  /** Standard typedefs. */
  /*************************************************************************/
  const int dim = 3;
	typedef double PixelType;

	typedef itk::Image< PixelType, dim > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

  /*************************************************************************/
  /** Organize command line input. */
  /*************************************************************************/

  if ( argc != 7 )
  {
   return Usage( argv );
  }

  std::string pdfn  = argv[1];
  std::string t2fn  = argv[2];
  std::string outfn = argv[6];

  float pdthresh    = atof(argv[3]);
  float t2thresh    = atof(argv[4]);
  int   histowindow = atoi(argv[5]);

 float histothresh = 0.00001;

  std::cout << std::endl;
  std::cout << "PD image:             " << pdfn << std::endl;
  std::cout << "T2 image:             " << t2fn << std::endl;
  std::cout << "PD threshold:         " << pdthresh << std::endl;
  std::cout << "T2 threshold:         " << t2thresh << std::endl;
  std::cout << "Histogram smoothing:  " << histowindow << std::endl;
  std::cout << "Output image:         " << outfn << std::endl;
  std::cout << std::endl;

  /** Check for filename extension. */
  int pos = outfn.find('.');

  if ( pos == std::string::npos ) 
  {
    std::cout<< "Error, must supply file extension for output." << std::endl; 
    return 1;
  }

  /** Create output filename prefix and image type suffix. */
  std::string outpfx( outfn, 0, pos );
  std::string outtype( outfn, pos, outfn.size() );

  /*************************************************************************/
  /** Read Input. */
  /*************************************************************************/
  
  ReaderType::Pointer pdReader = ReaderType::New();
  pdReader->SetFileName( pdfn.c_str() );
  if ( !egitk::ReadImage<PixelType, dim>( pdReader ) ) { exit(1);	}

	ReaderType::Pointer t2Reader = ReaderType::New();
  t2Reader->SetFileName( t2fn.c_str() );
  if ( !egitk::ReadImage<PixelType, dim>( t2Reader ) ) { exit(1);	}


  /*************************************************************************/
  /** Calculate Segmenation Parameters. */
  /*************************************************************************/

  std::vector<float> segParameters = GetSegParameters<PixelType, dim>( pdReader->GetOutput(), 
                                                                     t2Reader->GetOutput(), 
                                                                     histothresh,
                                                                     histowindow           );

  float histoMaxIndex  = segParameters[0];
  float histoMaxValue  = segParameters[1];
  float histoFWHM      = segParameters[2];
  float nwmWindowForPD = segParameters[3];
  float pdMin          = segParameters[4];
  float pdMax          = segParameters[5];
  float pdMean         = segParameters[6];

  std::cout << "Smoothed T2 histogram max index: " << histoMaxIndex<< std::endl;
  std::cout << "Smoothed T2 histogram max value: " << histoMaxValue << std::endl;
  std::cout << "Smoothed T2 FWHM:                " << histoFWHM << std::endl;
  std::cout << "Distance from peak for PD NAWM:  +/- " << nwmWindowForPD << std::endl;
  std::cout << "PD adjusted minimum:             " << pdMin << std::endl;
  std::cout << "PD adjusted maximum:             " << pdMax << std::endl;
  std::cout << "PD mean NWM:                     " << pdMean << std::endl;

  float t2LesionThreshold = histoMaxIndex + ( histoFWHM - histoMaxIndex ) * t2thresh;
  float pdLesionThreshold = pdMean + ( pdMax - pdMin ) * pdthresh;

  std::cout << "PD lesion threshold:             " << pdLesionThreshold << std::endl;
  std::cout << "T2 lesion threshold:             " << t2LesionThreshold << std::endl;
  
  /*************************************************************************/
  /** Segment And Save. */
  /*************************************************************************/

 typedef itk::ImageRegionIterator< ImageType > ItType; 
 
  /** Segment T2. */
 ItType t2It( t2Reader->GetOutput(), t2Reader->GetOutput()->GetLargestPossibleRegion() );
  
  for ( t2It.GoToBegin(); !t2It.IsAtEnd(); ++t2It )
  {
    if ( t2It.Get() <  t2LesionThreshold )
    {
      t2It.Set( 0.0 );
    }
    else
    {
      t2It.Set( 1.0 );
    }
  }

  /** Save T2 segmentation. */
  std::string out;
  out = outpfx + "_t2seg" + outtype;
  std::cout << std::endl << "Saving T2 segmentation...";
  egitk::WriteImage< PixelType, char, dim >( out.c_str(), t2Reader->GetOutput() );
  std::cout << "  Done."  << std::endl;
  
  /** Segment PD. */
  ItType pdIt( pdReader->GetOutput(), pdReader->GetOutput()->GetLargestPossibleRegion() );
    
  for ( pdIt.GoToBegin(); !pdIt.IsAtEnd(); ++pdIt )
  {
    if ( pdIt.Get() < pdLesionThreshold )
    {
      pdIt.Set( 0.0 );
    }
    else
    {
      pdIt.Set( 1.0 );
    }
  }

  /** Save PD segmentation. */
  out = outpfx + "_pdseg" + outtype;
  std::cout << std::endl << "Saving PD segmentation...";
  egitk::WriteImage< PixelType, char, dim >( out.c_str(), pdReader->GetOutput() );
  std::cout << "  Done."  << std::endl;

  /** Combine PD and T2 segmentations. */

  for ( t2It.GoToBegin(), pdIt.GoToBegin(); !t2It.IsAtEnd(); ++t2It, ++pdIt )
  {
    if ( t2It.Get() > 0 && pdIt.Get() > 0  )
    {
      pdIt.Set( 1.0 );
    } 
    else 
    {
      pdIt.Set( 0.0 );
    }
  }

  /** Save combined PD and T2 segmentation. */
  out = outpfx + "_pdt2seg" + outtype;
  std::cout << std::endl << "Saving combined PD/T2 segmentation...";
  egitk::WriteImage< PixelType, char, dim >( out.c_str(), pdReader->GetOutput() );
  std::cout << "  Done."  << std::endl;

  return 0;

}