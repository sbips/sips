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
#include "itkOrientImageFilter.h"
#include "itkAnalyzeImageIO.h"
#include "itkSpatialOrientation.h"
#include "itkIOCommon.h"
#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkNormalizeImageFilter.h"
#include <typeinfo>
#include "../Common/itk_imgmath_util_templates.h"
#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgops_util_templates.h"



int main( int argc, char* argv[] )
{

  if ( argc != 3 ) return 1;

  itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();

  typedef itk::Image< double, 3 > ImageType;  
  typedef itk::ImageFileReader< ImageType > ReaderType;
    
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetImageIO( io );
  reader->SetFileName( argv[1] );
  reader->Update();
  ImageType::Pointer input = reader->GetOutput();
 
  typedef itk::NormalizeImageFilter< ImageType, ImageType > NormalizeFilterType;
  NormalizeFilterType::Pointer normalize = NormalizeFilterType::New();

  normalize->SetInput( input );
  normalize->Update();




  input = normalize->GetOutput();

  io_utils::WriteItkVol< double, float, 3 > ( argv[2], input );
 
 /*

  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned char ) )
  {
    io_utils::WriteItkVol< double, unsigned char, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( char ) )
  {   
    io_utils::WriteItkVol< double, char, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned short ) )
  {   
    io_utils::WriteItkVol< double, unsigned short, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid (  short ) )
  {   
    io_utils::WriteItkVol< double, short, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned long ) )
  {
    io_utils::WriteItkVol< double, unsigned long, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( long ) )
  {   
    io_utils::WriteItkVol< double, long, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( float ) )
  {   
    io_utils::WriteItkVol< double, float, 3 > ( argv[2], input );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( double ) )
  {   
   io_utils::WriteItkVol< double, double, 3 > ( argv[2], input );
  }
  else 
  
*/
  return 0;
}
