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

#include "../../../include_from_cortex/analyze_utils.h"
#include <iomanip.h>
//#include "/home/sela/codeLib/c++/gs-2.0/include/gs/command_line.h"

// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o stats2 stats2.cpp

using namespace nk;

typedef unsigned char uchar;

uchar  BACKGROUND = 0;
#include "../../../include_from_cortex/lobe_struct2.h"

int main(int argc, char **argv)
{
    if( argc != 3 )
        {
            cout <<"Usage: stats2 <lobmask file> <seg file>"  << endl;
            exit(1);
        }


    int sz[3];
    float voxsize[3];

    string lobmask_fname=argv[1];
    string seg_fname=argv[2];
    info lobmask_info, seg_info;
    unsigned char *lobmask=0, *seg=0;

    // read all relevant images

    lobmask = read_analyze<unsigned char>(lobmask_fname, lobmask_info);
    seg= read_analyze<unsigned char>(seg_fname, seg_info);

    // check for size
    if( !equal(lobmask_info.sz, lobmask_info.sz+3, seg_info.sz))
    {
         cout << "stats2: Error, size mismatch!  Exiting ... " << endl;
         exit(1);
    }
   
    long volume = lobmask_info.sz[0] * lobmask_info.sz[1] * lobmask_info.sz[2];
    double vox_vol = lobmask_info.voxsz[0] * lobmask_info.voxsz[1] * lobmask_info.voxsz[2];
   
    //Forming Lobes
    lobe_struct lobe[256];
    for(int i = 0; i < 256; ++i) lobe[i].set(i);
  
    uchar *lob_ptr = lobmask, *seg_ptr = seg;
    int MAXLOBE=0, MINLOBE=255;
    for(int t=0; t < volume; ++t) {
        lobe[*lob_ptr].update_counts(seg_ptr);
        ++seg_ptr;
        if(*lob_ptr > MAXLOBE) MAXLOBE = *lob_ptr;
        if(*lob_ptr < MINLOBE) MINLOBE = *lob_ptr;
        ++lob_ptr;
    }

    // output statistics
    cout <<"      tot lob vol        Vcsf           Scsf         gray          white" << endl;
    //    cout << seg_fname << endl;
    for(int i = MINLOBE; i<= MAXLOBE; ++i) 
            lobe[i].print_excel_stats(vox_vol);
   
    //Delete Pointers
    delete [] lobmask;
    delete [] seg;
    return 0;
}
