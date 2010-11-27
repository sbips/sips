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

#include "itkImage.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"

#include <vector>



namespace fl_hist
{

  itk::Image< double, 2>::Pointer
  ExtractSliceFromVol( itk::Image< double, 3 >::Pointer input, int slicenum );


  std::vector<int>
  histo( itk::Image< double, 3>::Pointer flInput, itk::Image<double, 3>::Pointer flfcInput );


  itk::Image< double, 3 >::Pointer
  GetFlairLesions( itk::Image< double, 3>::Pointer fl, std::vector< int >& flThresh );

}
