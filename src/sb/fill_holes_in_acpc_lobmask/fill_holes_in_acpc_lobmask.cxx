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

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "../Common/itk_imgio_util_templates.h"
#include <vector>

void Usage( char* argv[] )
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << argv[0] << " lobmaskACPC lobmaskACPCfilled" << std::endl;
  std::cout << std::endl;
  std::cout << "Notes:" << std::endl;
  std::cout << "-file extension required (.img/.nii/.nii.gz)" << std::endl;
  std::cout << std::endl;
  exit(1);
}

int main( int argc, char* argv[] )
{

  if ( argc != 3 ) Usage( argv );
  
  typedef double PixelType;
  const int Dim = 3;
  typedef itk::Image< PixelType, Dim > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  // Read input.
  if ( !egitk::ReadImage< PixelType, Dim >( reader ) ) { std::cout << "Error Reading File."  << std::endl; return 1; };

  ImageType::Pointer input = ImageType::New();
  input = reader->GetOutput();

  ImageType::SizeType inSize = input->GetLargestPossibleRegion().GetSize();

  // Set repeat flag 
  // (handles the possibility that there may be two separate holes on one side of the image)
  int repeat = 1;
 
  while ( repeat > 0 )
  {
  	  
    for ( int z=0; z <inSize[2]; ++ z )
    {

      //=========================================================================
      // RIGHT SIDE
      //==========================================================================
      ImageType::IndexType inIndex; 
      inIndex[0] = 0;
      inIndex[1] = 0;
      inIndex[2] = z;
      std::vector< int > startY_indices;
      
      // 0 0 0 0 X X X X ....
      // ~
      // (find the (0,Y) locations that are zero)
      for ( int i = inSize[2] - 1; i >= 0; --i )
      {
        inIndex[1] = i;
        if ( input->GetPixel( inIndex ) == 0 )
        // hole found
        {
          // check that hole is continuous
          if ( startY_indices.size() != 0 )
	        {
	          if ( abs( i - startY_indices.back() ) > 1 )
	          {
              // then we have two holes on one side of the image and the process should be repeated
	            ++repeat; break;
	           }
	        }
         // add hole Y location to vector
	       startY_indices.push_back(i);
        }
      }

      // if no holes found, skip slice
      if ( startY_indices.size() == 0 ) continue;
    
      // 0 0 0 0 X X X X ....
      //       ~
      // find the last zero X location for each of the zero (0,Y) starting locations
      std::vector< int > end2X_indices;
      if ( startY_indices.size() != 0 )
      {
          for ( int i = 0; i < startY_indices.size(); ++i )
          {
            inIndex[1] = startY_indices[i];
            inIndex[0] = 1;
            while( input->GetPixel( inIndex ) == 0 )
            {
              inIndex[0] = inIndex[0] + 1;
            }
            end2X_indices.push_back(inIndex[0] - 1 );       
          }
      }

      // find the ROI values in the (x,y) locations just above the first (0,Y) line of zeros
      // (this line will contain either 1 or 2 ROI values)
      std::vector< int > fill_colour;
      int region_change_marker = -1;
      inIndex[0] = 0;
      inIndex[1] = startY_indices[0]+1;
      fill_colour.push_back( input->GetPixel( inIndex ) );
      for ( int i = 0; i < end2X_indices[0]; ++i )
      {
        if ( input->GetPixel( inIndex ) != fill_colour[0] )
        {
          fill_colour.push_back( input->GetPixel( inIndex ) );
          region_change_marker = inIndex[0] - 1;
          break;
        }
        ++inIndex[0];
      }

      if ( fill_colour.size() > 1 && region_change_marker >= end2X_indices.back() ) 
      { 
        std::cout << std::endl << "Warning: unusual hole found on right side of slice " << z+1 << ", trying to fill, but check the result. " << std::endl;
        // fill hole with the first fill_colour, ignore any other fill_colours found
	fill_colour.resize(1);
      }
      
      
      // if just one ROI value found, fill all zeros voxels with that same value
      if ( fill_colour.size() == 1 ) // fill with one ROI value
      {
        // fill
        for (int i=0; i< startY_indices.size(); ++i)
        {
          inIndex[1] = startY_indices[i];
          for (int x=0; x <= end2X_indices[i]; ++x )
          {
            inIndex[0] = x;
            input->SetPixel( inIndex, fill_colour[0] );
          }
        }
      }
      else // fill with two ROI values, along the appropriate diagonal
      {
        //if ( region_change_marker >= end2X_indices.back() ) { std::cout << "Warning: unknown hole type found and not filled. " << std::endl; continue; }
        
        // calculate the new X position of the 2nd ROI value for each (0,Y) zero position
        int gradient = std::abs( region_change_marker - end2X_indices.back() );
        
        std::vector< int > start2_indices;
        for (int i = 0; i < startY_indices.size(); ++i )
        {
          start2_indices.push_back( region_change_marker + gradient/startY_indices.size()*(i+1) );
        }

        // fill
        for ( int i = 0; i < startY_indices.size(); ++i )
        {
          inIndex[1] = startY_indices[i];
          for ( int x = 0; x <= end2X_indices[i]; ++x )
          {
            inIndex[0] = x;
            if ( x < start2_indices[i] )
              input->SetPixel( inIndex, fill_colour[0] );
            else
              input->SetPixel( inIndex, fill_colour[1] );
           }
        }
      }

      
    
    } // end slice loop for right side

    
   for ( int z=0; z <inSize[2]; ++ z )
   { 
    //=========================================================================
    // LEFT SIDE
    //==========================================================================

      ImageType::IndexType inIndex; 
      inIndex[0] = inSize[0] - 1;
      inIndex[1] = 0;
      inIndex[2] = z;
      std::vector< int > startY_indices;
   
      // ... X X X X 0 0 0 0 
      //                   ~
      // (find the (inSize[0],Y) locations that are zero)
      for ( int i = inSize[2]-1; i >= 0; --i )
      {
        inIndex[1] = i;
        if ( input->GetPixel( inIndex ) == 0 )
        {
	        if ( startY_indices.size() != 0 )
	        {
	          if ( abs( i - startY_indices.back() ) > 1)
	          {
	            ++repeat; break;
	          }
	        }
	       startY_indices.push_back(i);
        }
      }

      // if no holes found, skip slice
      if ( startY_indices.size() == 0 ) continue;
    
      // 0 0 0 0 X X X X ....
      //       ~
      // find the last zero X location for each of the zero (inSize[1]-1,Y) starting locations
      std::vector< int > end2X_indices;
      for ( int i = 0; i < startY_indices.size(); ++i )
      {
        inIndex[0] = inSize[0] - 1;
        inIndex[1] = startY_indices[i];
        while( input->GetPixel( inIndex ) == 0 )
        {
          inIndex[0] = inIndex[0] - 1;
        }
        end2X_indices.push_back(inIndex[0] + 1 );       
      }

      // find the ROI values in the (x,y) locations just above the first (0,Y) line of zeros
      // (this line will contain either 1 or 2 ROI values)
      std::vector< int > fill_colour;
      int region_change_marker = -1;
      inIndex[0] = inSize[0] - 1;
      inIndex[1] = startY_indices[0] + 1;
      fill_colour.push_back( input->GetPixel( inIndex ) );
      for ( int i = inSize[0]-1; i > end2X_indices[0]; --i )
      {
        if ( input->GetPixel( inIndex ) != fill_colour[0] )
        {
          fill_colour.push_back( input->GetPixel( inIndex ) );
          region_change_marker = inIndex[0] + 1;
          break;
        }
        --inIndex[0];
      }
      
      if ( fill_colour.size() > 1 && region_change_marker <= end2X_indices.back() ) 
      { 
        std::cout << std::endl << "Warning: unusual hole found on left side of slice " << z+1 <<", trying to fill, but check the result. " << std::endl;
	fill_colour.resize(1);
      }
      
      // if just one ROI value found, fill all zeros voxels with that same value
      if ( fill_colour.size() == 1 ) // fill with one ROI value
      {
        // fill
        for ( int i = 0; i< startY_indices.size(); ++i)
        {
          inIndex[1] = startY_indices[i];
          for ( int x = inSize[0] - 1; x >= end2X_indices[i]; --x )
          {
            inIndex[0] = x;
            input->SetPixel( inIndex, fill_colour[0] );
          }
        }
      } 
      else // fill with two ROI values, along the appropriate diagonal
      {
        //if ( region_change_marker <= end2X_indices.back() ) { std::cout << "Warning: unknown hole type found and not filled. " << std::endl; continue; }
        
        // calculate the new X position of the 2nd ROI value for each (inSize[0]-1,Y) zero position
        int gradient = std::abs( region_change_marker - end2X_indices.back() );
        
        std::vector< int > start2_indices;
        for (int i = 0; i < startY_indices.size(); ++i )
        {
          start2_indices.push_back( region_change_marker - gradient/startY_indices.size()*(i+1) );
        }

        // fill
        for ( int i = 0; i < startY_indices.size(); ++i )
        {
          inIndex[1] = startY_indices[i];
          for ( int x = inSize[0]-1; x >= end2X_indices[i]; --x )
          {
            inIndex[0] = x;
            if ( x < start2_indices[i] )
              input->SetPixel( inIndex, fill_colour[1] );
            else
              input->SetPixel( inIndex, fill_colour[0] );
           }
        }
      }
      
   
   } // end z slice loop left side
   --repeat;

 } // end while loop
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  egitk::WriteImage< PixelType, unsigned char, 3>( argv[2], input );

  
  return 0;
}