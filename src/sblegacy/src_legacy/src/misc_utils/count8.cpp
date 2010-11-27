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

//#include "/net/axon/axon-local/natasa/analyze_utils/analyze_utils.h"
#include "../../../include_from_cortex/analyze_utils.h"
#include <iomanip.h>

// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o count8 count8.cpp 

using namespace nk;

int main(int argc, char **argv)
{
    if( argc < 2 )
        {
            cout <<"
Usage: count8 <fileroot> - no .img or .hdr extensions.
             [-e] - if just one long line of volumes is required
    \n\n"
                 <<"Example: count garv_T1_seg \n";
            exit(1);
        }
           

    int sz[3];         
    float voxsize[3];
    
     // read input image
    unsigned char *Img;
    string fname = argv[1];
    info hdr_info;
    Img = read_analyze<unsigned char>(fname, hdr_info);
    
    int volume = hdr_info.sz[0]*hdr_info.sz[1]*hdr_info.sz[2];
    float voxel_vol = hdr_info.voxsz[0]*hdr_info.voxsz[1]*hdr_info.voxsz[2];

    int histo[256];
    for(int i=0; i<256; i++) histo[i]=0;
    for(int i=0; i<volume; i++) ++histo[ Img[i] ];
    
    CommandLine cl(argc, argv);
    
    int nozero_count=0;  // to hold total number of non-zero voxels
    int no2_count=0;     // to hold total number of non-2 voxels
    
    if(cl.find("-e")) { // if excel format for output wanted
        cout << fname << "\t";
        cout.setf(ios::fixed, ios::floatfield);
        for(int i=0; i<256; i++) {
            if(i) nozero_count += histo[i];
            if( histo[i] != 0 ){
                cout << setw(15) << setiosflags(ios::right)
                     << voxel_vol * histo[i] << ",";
            }
        }
        cout << setw(15) << setiosflags(ios::right)
             << voxel_vol * nozero_count << endl;
    }

    else {
        cout <<"---------------------------------------------------------------\n";
        cout << "        count            volume\n";
        cout <<"---------------------------------------------------------------\n";
        cout.setf(ios::fixed, ios::floatfield);
        for(int i=0; i<256; i++) {
            if(i) nozero_count += histo[i];
            if(i != 2) no2_count += histo[i];
            if( histo[i] != 0 ){
                cout << setw(3) << setiosflags(ios::right) << i
                     << setw(10) << setiosflags(ios::right)
                     << histo[i] << "\t" << setw(15) << setiosflags(ios::right)
                     << voxel_vol * histo[i] << endl;
            }
        }
        cout <<"---------------------------------------------------------------\n";
        cout << "total count of non-zero voxels: "
             << setw(10) << setiosflags(ios::right)
             << nozero_count << "\t" << setw(15) << setiosflags(ios::right)
             << voxel_vol * nozero_count << endl;
        cout <<"---------------------------------------------------------------\n";
        
        cout << "total count of non-2 voxels:    "
             << setw(10) << setiosflags(ios::right)
             << no2_count << "\t" << setw(15) << setiosflags(ios::right)
             << voxel_vol * no2_count << endl;
        cout <<"---------------------------------------------------------------\n";
    }
    cout << flush;
    delete [] Img;
    return 0;
}
