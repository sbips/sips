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

// compile with:
// g++ -O2 -o lobmask lobmask.cpp 


#include <assert.h>
#include <iomanip.h>
#include <vector>
#include <algorithm>
#include "globals.h"
#include "grid.h"
#include "make_grid.h"
#include "bresenhem.h"
#include "tracings.cpp"
#include "floodfill.cpp"
#include "basal_ganglia.cpp"
#include "orbital_frontal.cpp"


void copy_trace_mask(uchar *trace_mask, uchar *mask, int xstride,
                     int beg, int med, int end, uchar lobe_offset);
void fill_top_and_bottom_slices(uchar *mask);
void fill_front_and_back_slices(uchar *mask);
void fill_right_and_left_slices(uchar *mask);
//////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc < 6 ) print_usage_msg();
    
    info hdr_info;
    string fname(argv[1]);
    read_analyze_header(fname + ".hdr", hdr_info);

    if(hdr_info.datatype == 2) make_grid<uchar>(argc, argv);
    else {
        if(hdr_info.datatype == 4) make_grid<short>(argc, argv);
        else {
            cout << "Input file has wrong datatype:"
                 <<" only uchar and short are allowed. Exiting... " << endl;
            exit(1);
        }
    }
  
   uchar *mask = new uchar[volume];
   if(!mask) {
       cout << "Memory allocation failed. Exiting..." << endl;
       exit(1);
   }
   fill(mask, mask + volume, (uchar)BACKGROUND);
  // make right and left mask edges filled with gridline
   for(int z=0; z < sz[2]; ++z) 
       for(int y=0; y < sz[1]; ++y)
           mask[z*area + y*sz[0]] = mask[z*area + y*sz[0] + sz[0]-1] =GRIDLINE;
       
   uchar *left_trace_mask, *right_trace_mask ;
   load_tracings(argv[2], left_trace_mask, right_trace_mask);
  
   
   connect_trancings(left_trace_mask, 0, lpron);
   fill_trace_mask(left_trace_mask, 0, 0);
   recolor_lines(left_trace_mask, 0, 0);
   
   connect_trancings(right_trace_mask, LINE_OFFSET, rpron);
   fill_trace_mask(right_trace_mask, LINE_OFFSET, LOBE_OFFSET);
   recolor_lines(right_trace_mask, LINE_OFFSET, LOBE_OFFSET);
   
   CommandLine cl(argc, argv);
   if(cl.find("-debug")) {
       info info_temp;
       info_temp.max = 255;
       info_temp.datatype = 2;
       info_temp.min = 0;
       info_temp.sz[0]=sz[1]; info_temp.sz[1]=sz[2]; info_temp.sz[2]=1;
       cout << "Writing left trace mask" << endl;
       write_analyze("left", left_trace_mask,  info_temp);
       cout << "Writing right trace mask" << endl;
       write_analyze("right", right_trace_mask,  info_temp);
   }
 
    // now copy trace_masks all the way in
    copy_trace_mask(left_trace_mask, mask, -1,
                     sz[0]-2, grid.sagital[5], grid.sagital[4], 0);
    copy_trace_mask(right_trace_mask, mask, 1,
                     1, grid.sagital[3], grid.sagital[4]+1, LOBE_OFFSET);

    // fill OBF lobes
    orbital_frontal(mask);
   
    // fill out basal ganglia
    basal_ganglia(mask);

    // reassign C lines
     uchar *ptr = mask;
     int z;
     for(z=0; z < volume; ++z) {
         if( *ptr == LC_LINE || *ptr == RC_LINE ) *ptr = *(ptr-sz[0]);
         ++ptr;
     }

    // now that all floodfills are done all edge slices are grey
    fill_top_and_bottom_slices(mask);
    fill_front_and_back_slices(mask);
    fill_right_and_left_slices(mask);

    // if there are any weird (not any lobe type) values
    // issue a warning
    ptr = mask;
    uchar min = 255, max = 0;
    for(z=0; z < volume; ++z) {
         if( *ptr > max ) max = *ptr;
         if( *ptr < min ) min = *ptr;
         ++ptr;
     }
    if( min != LSUPF || max != RMIF )
        cout << "Warning: lobar mask has holes." << endl;
        
    // output
    hdr_info.max = 255;
    hdr_info.min = 0;
    hdr_info.datatype = 2;
    cout << "Writing " << argv[3] << endl;
    write_analyze(string(argv[3]), mask,  hdr_info);
    
    
    delete [] mask;
    delete [] left_trace_mask;
    delete [] right_trace_mask;
    return 0;
}  // end of main

//////////////////////////////////////////////////////////////////////
void copy_trace_mask(uchar *trace_mask, uchar *mask, int xstride,
                     int beg, int med, int end, uchar lobe_offset) {
    int x, y, z;
    uchar val;
    for(z=0; z < sz[2]; ++z)
        for(y=0; y < sz[1]; ++y) {
            
            for( x=beg; x != med; x += xstride)
                mask[z*area + y*sz[0] + x] =
                    trace_mask[z*sz[1] + y];

            for(x=med; x != end; x += xstride){
                val =  trace_mask[z*sz[1] + y];
                switch( val - lobe_offset ) {
                case LSUPF:
                    mask[z*area + y*sz[0] + x] = LMSF + lobe_offset;
                    break;
                case LIF:
                    mask[z*area + y*sz[0] + x] = LMIF + lobe_offset;
                    break;
                default:
                    mask[z*area + y*sz[0] + x] = val;
                    break;
                }
            }
        }
}
////////////////////////////////////////////////////////////////////
void fill_top_and_bottom_slices(uchar *mask) {
    uchar *bot_ptr = mask, *top_ptr = mask + area * (sz[2]-1);
    for(int t=0; t < area; ++t) {
        *bot_ptr = *(bot_ptr + area);
        *top_ptr = *(top_ptr - area);
        ++bot_ptr;
        ++top_ptr;
    }
}
//////////////////////////////////////////////////////////////////////
void fill_front_and_back_slices(uchar *mask) {
    int y2 = sz[1]-1;
    int x, z, t;
    for(z=0; z < sz[2]; ++z)
        for(x=0; x < sz[0]; ++x) {
            t = z * area + x;
            mask[t] = mask[t + sz[0]];
            mask[t + y2*sz[0]] = mask[t + (y2-1)*sz[0]];
        }
}
//////////////////////////////////////////////////////////////////////
void fill_right_and_left_slices(uchar *mask) {
    int x2 = sz[0]-1;
    int x, y, z, t;
    for(z=0; z < sz[2]; ++z)
        for(y=0; y < sz[1]; ++y) {
            t = z * area + y * sz[0];
            mask[t] = mask[t + 1];
            mask[t + x2] = mask[t + x2-1];
        }
}
