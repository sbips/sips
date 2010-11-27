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
//   Last Modified: 10/12/2008  //

//**********************************************************************************************************

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImageIOBase.h"
#include "itkAnalyzeImageIO.h"

#include "../Common/itk_imgmath_util_templates.h"



//**********************************************************************************************************
// Templated Diffusion Filter Pipeline.
template < typename InputImageType, typename OutputImageType >
int FilterImage( typename itk::Image< double, 3 >::Pointer input, char* argv[] )
{
  // Filter input image, convert to output image type, and write output //

  // Standard typedefs.
  typedef typename itk::GradientAnisotropicDiffusionImageFilter
                   < InputImageType, InputImageType >  FilterType;
  typedef typename itk::RescaleIntensityImageFilter
                   < InputImageType, OutputImageType > RescaleFilterType;
  typedef typename itk::ImageFileWriter
                   < OutputImageType >                 WriterType;

  // Get filter parameters from command line arguments.
  const unsigned int numberOfIterations = atof( argv[3] ); // default 5
  const double       timeStep = atof ( argv[4] ); // default 0.0625;
  const double       conductance = atof ( argv[5] ); //1.95;

  // Setup diffusion filter.
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( input );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  filter->Update();



  // Setup rescale filter.
  // Ensure that the output image has the same max and min as the input image
  //    or else the entire allowable range of the output datatype will be used.
  typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum( itk_imath::ImageMin< double, 3 >( input ) );
  rescaler->SetOutputMaximum( itk_imath::ImageMax< double, 3 >( input ) );
  rescaler->SetInput( filter->GetOutput() );
  rescaler->Update();


  itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();

  // Setup image file writer.
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetImageIO( io );
  writer->SetInput( rescaler->GetOutput() );

  // Write image output.
  try 
  { 
	  writer->Update(); 
  } 
  catch( itk::ExceptionObject& err ) 
  { 
	  std::cerr << std::endl << "Exception caught!!" << std::endl; 
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }	
  
  return EXIT_SUCCESS;
}
//**********************************************************************************************************
// Display command line usage information.
int Usage( char* argv[] )
{
	  std::cerr << std::endl;
    std::cerr << "Usage: " << std::endl;
	  std::cerr << argv[0] << std::endl;
    std::cerr << "  <in.img> <out.img> <num_iterations> <time_step> <conductance>" << std::endl;
	  std::cerr << std::endl;
    std::cerr << "  Applies 3D anisotropic diffusion filter to an image." << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Must supply the .img extension for both input and output filenames." << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Typical values: num_iterations = 5" << std::endl;
    std::cerr << "                  time_step = 0.0625" << std::endl;
    std::cerr << "                  conductance= 1.95" << std::endl;
  	
	  return EXIT_FAILURE;

}
//**********************************************************************************************************
// Run program.
int main( int argc, char* argv[] )
{
  // Output usage message if fewer than 5 arguments are provided.
  if( argc < 6 ) 
  {
    return Usage( argv );
  }
  
  // Standard typedefs.
  typedef itk::Image< double, 3 >                InputImageType;
  typedef itk::Image< float, 3 >                 ProcImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;

  // Setup image file reader.
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] ); 
  
  try 
  { 
	  reader->Update(); 
  } 
  catch( itk::ExceptionObject& err ) 
  { 
	  std::cerr << std::endl << "Exception caught!!" << std::endl; 
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }	
  
  // Apply diffusion filtering on input image.
  // Output image will have same type as input image.
  //   ( Query reader's ImageIOBase object about the actual pixel type contained in the image
  //     and then call the appropriate templated version of FilterImage function ) 
  

  // Save output.

  if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned char ) )
  {
    return FilterImage< InputImageType, itk::Image< unsigned char, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( char ) )
  {   
    return FilterImage< InputImageType, itk::Image<  char, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned short ) )
  {   
    return FilterImage< InputImageType, itk::Image< unsigned short, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid (  short ) )
  {   
    return FilterImage< InputImageType, itk::Image<  short, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( unsigned long ) )
  {
    return FilterImage< InputImageType, itk::Image< unsigned long, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( long ) )
  {   
    return FilterImage< InputImageType, itk::Image< long, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( float ) )
  {   
    return FilterImage< InputImageType, itk::Image< float, 3> >( reader->GetOutput(), argv );
  }
  else if ( reader->GetImageIO()->GetComponentTypeInfo() == typeid ( double ) )
  {   
    return FilterImage< InputImageType, itk::Image< double, 3> >( reader->GetOutput(), argv );
  }
  else 
  {
    std::cerr << "Error: Input image has an unrecognized datatype" << std::endl;
    return EXIT_FAILURE;
  }

}
  
