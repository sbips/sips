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

/////////////////////////////////////////////////////////////
//
// written by Natasa Kovacevic, 2000
//
////////////////////////////////////////////////////////////
//#include "/net/axon/axon-local/natasa/analyze_utils/analyze_utils.h"
#include "../../../include_from_cortex/analyze_utils.h"
#include <iomanip.h>

// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o count16 count16.cpp

using namespace nk;

int main(int argc, char **argv)
{
    if( argc < 2 )
        {
            cout <<"Usage: count16 <fileroot> - no .img or .hdr extensions.\n\n"
                 <<"Example: count16 garv_T1seg \n";
            return(-1);
        }
     // read input image
    info hdr;
    unsigned short *Img;
    Img = read_analyze<unsigned short>(argv[1], hdr);
    
    int histo[65536];
    
    for(int i=0; i < 65536; i++) histo[i]=0;

    int volume = hdr.sz[0]*hdr.sz[1]*hdr.sz[2];
    for(int i=0; i < volume; i++) ++histo[ Img[i] ];

    delete [] Img;
    
    // output voxel size
    float voxel_vol = hdr.voxsz[0]*hdr.voxsz[1]*hdr.voxsz[2];
    cout <<"---------------------------------------------------------------\n";
    cout << "        count            volume\n";
    cout <<"---------------------------------------------------------------\n";
    int brain_count=0;  // to hold total number of non-zero voxels
    cout.setf(ios::fixed, ios::floatfield);
    for(int i=0; i<65536; i++) {
        if(i) brain_count += histo[i];
        if( histo[i] != 0 ){
            cout << setw(5) << setiosflags(ios::right) << i
             << setw(10) << setiosflags(ios::right)
             << histo[i] << "\t" << setw(15) << setiosflags(ios::right)
             << voxel_vol * histo[i] << endl;
        }
    }
    cout <<"---------------------------------------------------------------\n";
    cout << "total count of non-zero voxels: "
         << setw(10) << setiosflags(ios::right)
         << brain_count << "\t" << setw(15) << setiosflags(ios::right)
         << voxel_vol * brain_count << endl;
    cout <<"---------------------------------------------------------------\n";
    return 0;
}
