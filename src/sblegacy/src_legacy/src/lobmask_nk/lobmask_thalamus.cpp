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
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o lobmask_thalamus lobmask_thalamus.cpp
//By: Dr. Azhar Quddus, June 4, 2003.

#include <assert.h>
#include <iomanip.h>
#include <vector>
#include <algorithm>
#include "make_grid.h"
#include "../../../include_from_cortex/grid2.h"
#include "bresenhem.h"
#include "globals.h"
#include "floodfill.cpp"


void copy_trace_mask(uchar *trace_mask, uchar *mask, int xstride,
                     int beg, int med, int end, uchar lobe_offset);
void copy_trace2(uchar *trace_mask, uchar *mask, int beg, int end);
void fill_top_and_bottom_slices(uchar *mask);
void fill_front_and_back_slices(uchar *mask);
void fill_right_and_left_slices(uchar *mask);
//////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc != 5 )
      {
	cout <<"Usage: lobmask_thalamus ErodeFile(for size info) OutputFile -g GridFile"<<endl;
	cout<<"To extract Thalamus. By: Dr. Azhar Quddus"<<endl;
	exit(1);
      }

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

   //Make Blank Axial slice
   int axial_area=sz[0] * sz[1];
   uchar *axial_plane=new uchar[axial_area];
   if(!axial_plane)
     {
       cout << "Failed to allocate memory for trace masks. Exiting..." << endl;
       exit(1);
     }
   fill(axial_plane, axial_plane + axial_area, BACKGROUND);

   //Now define region points [0]->x & [1]->y
   int A[2],B[2],C[2],D[2],E[2],F[2],G[2],H[2],I[2],J[2];

   A[0]=grid.m;
   A[1]=grid.pc;

   B[0]=(int)((float)grid.m-(2.0/3.0)*((float)grid.m-(float)grid.sagital[3]));
   B[1]=grid.coronal[7];

   C[0]=midpt(grid.sagital[2],grid.sagital[3]);
   C[1]=B[1];

   D[0]=C[0];
   D[1]=A[1];

   E[0]=B[0];
   E[1]=(int)((float)grid.ac-(1.0/4.0)*((float)grid.ac-(float)grid.pc));

   F[0]=A[0];
   F[1]=E[1];

   G[0]=(int)((float)grid.m+(2.0/3.0)*((float)grid.sagital[5]-(float)grid.m));
   G[1]=F[1];

   H[0]=midpt(grid.sagital[5],grid.sagital[6]);
   H[1]=A[1];

   I[0]=H[0];
   I[1]=B[1];

   J[0]=G[0];
   J[1]=I[1];

   //Draw lines Right side & fill
   int dim[2];
   dim[0]=sz[0];
   dim[1]=sz[1];

   uchar val=(uchar)1;

   draw_line(A[0],A[1],B[0],B[1],axial_plane,dim,val+1);
   draw_line(B[0],B[1],C[0],C[1],axial_plane,dim,val+1);
   draw_line(C[0],C[1],D[0],D[1],axial_plane,dim,val+1);
   draw_line(D[0],D[1],E[0],E[1],axial_plane,dim,val+1);
   draw_line(E[0],E[1],F[0],F[1],axial_plane,dim,val+1);
   draw_line(F[0],F[1],A[0],A[1],axial_plane,dim,val+1);

   //Flood Fill
   int center[2];
   center[0]=(A[0]+B[0]+C[0]+D[0]+E[0]+F[0]+A[0])/7;
   center[1]=(A[1]+B[1]+C[1]+D[1]+E[1]+F[1]+A[1])/7;
   //   cout <<"center: "<<center[0]<<"-"<<center[1]<<endl;
   {
   int dirs[] = {1, sz[1], -1, -sz[1]};

   uchar * mseed=axial_plane+(center[1])*sz[0]+center[0];

   Floodfill f_l(mseed,dirs,4,acceptable,val+1);
   }
   //Draw lines Left side & Fill
   dim[0]=sz[0];
   dim[1]=sz[1];

   draw_line(A[0],A[1],F[0],F[1],axial_plane,dim,val+3);
   draw_line(F[0],F[1],G[0],G[1],axial_plane,dim,val+3);
   draw_line(G[0]+1,G[1],H[0],H[1],axial_plane,dim,val+3);
   draw_line(H[0],H[1],I[0],I[1],axial_plane,dim,val+3);
   draw_line(I[0],I[1],J[0],J[1],axial_plane,dim,val+3);
   draw_line(J[0]+1,J[1],A[0],A[1],axial_plane,dim,val+3);
   //Flood Fill
   //int center[2];
   center[0]=(A[0]+F[0]+G[0]+H[0]+I[0]+J[0]+A[0])/7;
   center[1]=(A[1]+F[1]+G[1]+H[1]+I[1]+J[1]+A[1])/7;
   //   cout <<"center: "<<center[0]<<"-"<<center[1]<<endl;
   {
   int dirs[] = {1, sz[1], -1, -sz[1]};

   uchar * mseed=axial_plane+(center[1])*sz[0]+center[0];

   Floodfill f_l(mseed,dirs,4,acceptable,val+3);
   }

   //Cut Thalamus interior and posterior
   int TH=grid.pc+(int)(((double)grid.ac-(double)grid.pc)/3.0);

   for(int r=0; r<sz[1];r++)
     for(int c=0;c<sz[0];c++)
       {
	 if(axial_plane[r*sz[0]+c] && r>TH)
	   axial_plane[r*sz[0]+c] = axial_plane[r*sz[0]+c] + (uchar)5;
       }

    // now copy trace_masks all the way in
   int lastSlice=(int)((float)grid.axial[6]-(1.0/4.0)*((float)grid.axial[6]-(float)grid.axial[5]));
    copy_trace2(axial_plane, mask, grid.ap, lastSlice);

    // output
    hdr_info.max = 9;
    hdr_info.min = 0;
    hdr_info.datatype = 2;
    cout << "Writing: "<< argv[2]<<"...";
    write_analyze(string(argv[2]), mask,  hdr_info);
    cout <<"Done!"<<endl;
    
    delete [] mask;
    delete [] axial_plane;
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

//////////////////////////////////////////////////////////////////////
//Axial slice copy
void copy_trace2(uchar *trace_mask, uchar *mask, int beg, int end) 
{
    int x, y, z;
    uchar *buf=new uchar[sz[0]*sz[1]];
    memcpy(buf,trace_mask,sz[0]*sz[1]);

    //copy
    for(z=beg; z <= end; ++z){
      memcpy((mask+z*area),buf,sz[0]*sz[1]);
    }
    
    delete [] buf;
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
