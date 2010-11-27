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

// To get Thalamus.
// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o getGBT getBGT.cpp
//By: Dr. Azhar Quddus, July 14, 2003.

//#include <assert.h>
//#include <iomanip.h>
//#include <vector>
//#include <algorithm>
//#include "globals.h"
//#include "grid.h"
//#include "make_grid.h"
//#include "bresenhem.h"
#include "../../../include_from_cortex/analyze_utils.h"

using namespace nk;

typedef unsigned char uchar;

uchar  BACKGROUND = 0;


//////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc != 3 )
      {
	cout <<"Usage: getBGT <LobMaskFile> <OutputFile>"<<endl;
	cout<<"To extract only Basal Ganglia and Thalamus (left-12: right-25)."<<endl;
	cout<<"By: Dr. Azhar Quddus"<<endl;
	exit(1);
      }

    info in_info,out_info;
    string fname(argv[1]);

    uchar *ip;
    ip = read_analyze<uchar>(fname, in_info);

    if(in_info.datatype != 2)
      {
	cout<<"DataType not 8-bit/pp!.. exiting..."<<endl;
	return -1;
      }
   long vol = in_info.sz[0] * in_info.sz[1] * in_info.sz[2];

   uchar *mask = new uchar[vol];

   if(!mask) {
       cout << "Memory allocation failed. Exiting..." << endl;
       exit(1);
   }
   fill(mask, mask + vol, BACKGROUND);

   for(int i=0; i<vol;i++)
     {
       int pix=(int) ip[i];

       if(pix == 12 || pix == 13)
	 mask[i]=(uchar)12;

       if(pix == 25 || pix == 26)
	 mask[i]=(uchar)25;
     }

    // output
   out_info=in_info;
   out_info.max = 25;
   out_info.min = 0;
   out_info.datatype = 2;
   cout << "Writing: "<< argv[2]<<"...";
   write_analyze(string(argv[2]), (char *)mask,  out_info);
   cout <<"Done!"<<endl;

    delete [] mask;
    delete [] ip;

    return 0;
}  // end of main


