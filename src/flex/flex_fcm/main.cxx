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

//   Written 2006 by Erin Gibson     //
// Last Modified: 10/12/2008

// Applies Fuzzy-C Means (FCM) algorithm to 2D slices of a 3D image volume.

#include "FLEX_fcm.h"
#include "utils.h"

#include <vector>
#include <string>
#include <sstream>
#include <time.h>
#include <iostream>


int UsageMessage(char* argv[] )
{
	std::cerr << std::endl;
  std::cerr << "Usage:" << argv[0] << std::endl;
  std::cerr << "        <in.img>" << std::endl;
  std::cerr << "        <cluster_n>" << std::endl;;
  std::cerr << "        <exponent>" << std::endl;;
  std::cerr << "        <max_iter>" << std::endl;;
  std::cerr << "        <min_improv>" << std::endl;;
  std::cerr << "Typical values: exponent = 2" << std::endl;
  std::cerr << "                max_iter = 100" << std::endl;
  std::cerr << "                min_improv = 0.00001" << std::endl;
  std::cerr << std::endl;
	return 1;
}



int main( int argc, char* argv[] )
{
 
	if ( argc != 6 ) 
    return UsageMessage( argv );

  // Get input parameters.
	const std::string infn  = argv[1];
	const int    cluster_n  = atoi( argv[2] ); 
	const double expo       = atof( argv[3] );  
	const int    max_iter   = atoi( argv[4] ); 
	const double min_improv = atof( argv[5] ); 
	
	// Read image file.
	utils::ImageType3D16::Pointer input = utils::ItkReadInput( infn );
	
	// Create vector of pointers to images with a length equal to the cluster_n.
  //   This vector will be used to hold the cluster results.
  //   Each pointer is intialized by re-reading the input image data.
	std::vector< utils::ImageType3D16::Pointer > allOutput( cluster_n );
	for ( int i = 0; i < cluster_n; ++i )
		allOutput[i] = utils::ItkReadInput( infn );
	
  // Determine xyz dimensions of input image.
	const int xdim = input->GetLargestPossibleRegion().GetSize(0);
	const int ydim = input->GetLargestPossibleRegion().GetSize(1);
	const int zdim = input->GetLargestPossibleRegion().GetSize(2);

	//  Create a vector which will hold a single slice of data.
	std::vector< double > currentSliceVector(xdim * ydim);

  // Start FCM.
	for (int i = 0; i < zdim; ++i )
	{
		// Get image data for the current slice.
    currentSliceVector = utils::GetVectorSliceFromITKvol( input, i );

    // Output slice progress message.
		std::cout << "Fuzzy clustering slice:  " << i + 1 << std::endl;
		
    // Skip empty slices (or slices that have only a few non-zero voxels)
		if ( ! utils::IsEmptySliceVector( currentSliceVector) )
		{	
      // Run FCM on current slice.
			FCM fcm( currentSliceVector, cluster_n, expo, max_iter, min_improv );
			fcm.RunFcm();

      // Copy FCM results to allOutput vector.
			for ( int j = 0; j < cluster_n; ++j )
				utils::SaveVectorSliceToITKvol( fcm.GetU( j + 1 ), allOutput[j], i );
		}
	}



	// Save output with "_fcX_mcY_exX.img" appended to input prefix
  //   where X = number of clusters, Y = membership class center ranking
  //   (i.e. m1 == bg ) and Z = exponent used
  std::string outfn(infn);
  outfn.erase( outfn.length() - 4 );
  for ( int i = 0; i < cluster_n; ++i )
	{  
    std::ostringstream os;
    os << outfn << "_fc" << cluster_n << "_mc" << ( i + 1 ) << "_ex" << expo << ".img";    

		utils::ItkSaveOutput( os.str(), allOutput[ i ] );
	}

  return 0;

}



