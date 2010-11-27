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
#include <string>
#include <iostream>

#include "../Common/itk_imgio_util_templates.h"
#include "../Common/itk_imgmath_util_templates.h"

#include "itkImage.h"
#include "itkImageFileReader.h"

void Usage( char* argv[] )
{
	std::cout << "Usage:" << std::endl;
	std::cout << "  " << argv[0] << " imgIN1 [operations] imgOUT [options]" << std::endl;
	std::cout << std::endl;
	std::cout << "Operations:" << std::endl;
	std::cout << "  -add  <imgIN2>:  add two images" << std::endl;
	std::cout << "  -sub  <imgIN2>:  subtract two images"  << std::endl;
	std::cout << "  -mul  <imgIN2>:  multiply two images"  << std::endl;
	std::cout << "  -div  <imgIN2>:  divide two images" << std::endl;
	std::cout << "  -mas  <imgIN2>:  set zero valued voxels in imgIN2 to zero on imgIN1" << std::endl;
	std::cout << "  -min  <imgIN2>:  find minimum value at each voxel" << std::endl;
	std::cout << "  -max  <imgIN2>:  find maximum value at each voxel" << std::endl;
	std::cout << "  -thr  <float>:   zero voxels above floating point threshold" << std::endl;
	std::cout << "  -uthr <float>:   zero voxels below floating point threshold" << std::endl;
	std::cout << "  -bin:            set non-zero voxels equal to 1" << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "  -float:          write output as floats" << std::endl;
  std::cout << "  -ushort:         write output as unsigned short integers" << std::endl;
  std::cout << "  -short:          write output as signed short integers" << std::endl;
  std::cout << "  -uchar:           write output as unsigned char" << std::endl;
  std::cout << "  -char:           write output as signed char" << std::endl;	
	std::cout << std::endl;
	std::cout << "Notes:" << std::endl;
	std::cout << "  -file extension required (.img/.nii/.nii.gz)" << std::endl;
	std::cout << std::endl;
	exit(1);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] )
{
  int add     = 0;
	int sub     = 0;
	int mul     = 0;
	int div     = 0;
	int mas     = 0;
	int thr     = 0;
	int uthr    = 0;
	int flt     = 0;
  int uchr     = 0;
  int chr     = 0;	
  int ushrt   = 0;
  int shrt    = 0;
	int bin     = 0;
	int max     = 0;
	int min     = 0;
	float val   = 0;
	int opIndex = 0;
	int opType  = 0;
	int opTotal = 0;
	std::string img1;
	std::string img2;
	std::string out;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for ( int i=1; i<argc; ++i )
	{
		std::string cl( argv[i] );

		if ( cl == "-add" )
		  { ++add;   opType = 2;  opIndex = i;	++opTotal;  }
		if ( cl == "-sub" )
		  { ++sub;   opType = 2;  opIndex = i;  ++opTotal;  }
		if ( cl == "-mul" )
		  { ++mul;   opType = 2;  opIndex = i;  ++opTotal;  }
		if ( cl == "-div" )
		  {	++div;   opType = 2;  opIndex = i;  ++opTotal;  }
		if ( cl == "-min" )
		  { ++ min;  opType = 2;  opIndex = i;  ++opTotal; }
		if ( cl == "-max" )
		  { ++max; opType = 2;    opIndex = i;  ++opTotal; }
		if ( cl == "-mas" )
		  {	++mas;   opType = 2;  opIndex = i;  ++opTotal;  }
		if ( cl == "-thr" )
		  {	++thr;   opType = 3;  opIndex = i;  ++opTotal;  }
		if ( cl == "-uthr" )
		  { ++uthr;  opType = 3;  opIndex = i;  ++opTotal;  }
		if ( cl == "-bin" )
		  { ++bin;   opType = 1;  opIndex = i;  ++opTotal;  }

		if ( opType == 1 )
		{
			if ( opIndex + 1 > argc - 1 || opIndex== 0 )
			{
				std::cout << "input error" << std::endl;
				exit(1);
			}
			img1 = argv[opIndex-1];
		  out  = argv[opIndex+1];	
		}

		if ( opType == 2 )
		{
			if ( opIndex + 2 > argc - 1 || opIndex == 0 )
			{
				std::cout << "input error" << std::endl;
				exit(1);
			}
		  img1 = argv[opIndex-1];
		  img2 = argv[opIndex+1];
		  out  = argv[opIndex+2];
		}

		if ( opType == 3 )
		{
			if ( opIndex + 2 > argc - 1 || opIndex == 0 )
			{
				std::cout << "input error" << std::endl;
				exit(1);
			}
		  img1 = argv[opIndex-1];
		  val  = atof( argv[opIndex+1] );
		  out  = argv[opIndex+2];
		}

		if ( cl == "-float" )
			flt = 1;
		if ( cl == "-char" )
			chr = 1;
		if ( cl == "-uchar" )
			uchr = 1;
	       if ( cl == "-short" )
			shrt = 1;
		if ( cl == "-ushort" )
			ushrt = 1;
		}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( opTotal != 1 )
	{
		Usage( argv );
	}

	if ( opType == 1 )
	{
		if ( argc < 3 ) 
		  Usage( argv );
	}

	else if ( opType == 2 || opType == 3 )
	{
    if ( argc < 5 )
			Usage( argv );
	}

	const int dim = 3;
	typedef double PixelType;

	typedef itk::Image< PixelType, dim > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	ReaderType::Pointer img1Reader = ReaderType::New();
	ReaderType::Pointer img2Reader = ReaderType::New();

	img1Reader->SetFileName( img1.c_str() );
	if ( !egitk::ReadImage<PixelType, dim>( img1Reader ) )
	{
		exit(1);
	}

	if ( opType == 2 )
	{
	  img2Reader->SetFileName( img2.c_str() );
		if ( !egitk::ReadImage<PixelType, dim>( img2Reader ) )
		{
			exit(1);
		}
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ImageType::Pointer output = ImageType::New();

	if ( add )
	  output = egitk::ImageAddImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( sub )
		output = egitk::ImageSubImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( mul )
		output = egitk::ImageMulImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( div )
		output = egitk::ImageSubImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( min )
		output = egitk::ImageMinImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( max )
		output = egitk::ImageMaxImage< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( mas )
		output = egitk::ImageMask< PixelType, dim >( img1Reader->GetOutput(), img2Reader->GetOutput() );
	else if ( bin )
		output = egitk::ImageFindAndSet< PixelType, dim >(img1Reader->GetOutput(), img1Reader->GetOutput(), 1 );
	else if ( uthr )
		output = egitk::ImageUThr< PixelType, dim >( img1Reader->GetOutput(), val);
	else if ( thr )
		output = egitk::ImageThr< PixelType, dim >( img1Reader->GetOutput(), val);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( flt )
  {
    std::cout << "Writin output as float" << std::endl;
		egitk::WriteImage< PixelType, float, dim >( out.c_str(), output );
  }
  else if ( shrt )
  {
    std::cout << "Writing output as short"  << std::endl;
     egitk::WriteImage< PixelType, short, dim >( out.c_str(), output );
  }
  else if ( ushrt )
  {
    std::cout << "Writing output as unsigned short" << std::endl;
     egitk::WriteImage< PixelType, unsigned short, dim >( out.c_str(), output );
  }
  else if ( chr )
  {
    std::cout << "Writing output as signed char"  << std::endl;
     egitk::WriteImage< PixelType, char, dim >( out.c_str(), output );
  }
  else if ( uchr )
  {
    std::cout << "Writing output as unsigned char"  << std::endl;
     egitk::WriteImage< PixelType, unsigned char, dim >( out.c_str(), output );
  }  
	else
  {
		egitk::WriteImageSamePixelType< PixelType, dim >( out.c_str(), img1Reader, output );
  }

  return 0;

}

