#include "itkImage.h"
#include "itkOrientImageFilter.h"
#include "itkAnalyzeImageIO.h"
#include "itkSpatialOrientation.h"
#include "itkIOCommon.h"
#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
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
 
  itk::SpatialOrientation::ValidCoordinateOrientationFlags fileOrientation;
  itk::ExposeMetaData< itk::SpatialOrientation::ValidCoordinateOrientationFlags >
    ( input->GetMetaDataDictionary(), itk::ITK_CoordinateOrientation, fileOrientation );

  typedef itk::OrientImageFilter< ImageType, ImageType > OrientFilterType;

  OrientFilterType::Pointer orienter = OrientFilterType::New();

  orienter->SetGivenCoordinateOrientation( fileOrientation );
  orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI);
  orienter->SetInput( input );
  orienter->Update();

  input = orienter->GetOutput();
 
 

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
  

  return 0;
}