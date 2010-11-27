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
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o lobmask_cingulate lobmask_cingulate.cpp
//By: Dr. Azhar Quddus, Feb 28, 2003.

#include <assert.h>
#include <iomanip.h>
#include <vector>
#include <algorithm>
#include "../../../include_from_cortex/make_grid.h"
#include "../../../include_from_cortex/grid2.h"
#include "../../../include_from_cortex/bresenhem.h"
#include "../../../include_from_cortex/globals.h"
#include "../../../include_from_cortex/floodfill.cpp"
#define FALSE 0
#define TRUE 1

void copy_trace_mask(uchar *trace_mask, uchar *mask, int xstride,
                     int beg, int med, int end, uchar lobe_offset);
void copy_trace2(uchar *trace_mask, uchar *mask, int beg, int end, uchar color_offset);
void fill_top_and_bottom_slices(uchar *mask);
void fill_front_and_back_slices(uchar *mask);
void fill_right_and_left_slices(uchar *mask);

//////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc != 5 )
      {
	cout <<"Usage: lobmask_cingulate Object_File(after obj2img) OutputFile -g GridFile"<<endl<<"To make cigulate mask By: Dr. Azhar Quddus"<<endl;
	exit(1);
      }

    info hdr_info;
    string fname(argv[1]);
    uchar *img=read_analyze<uchar>(fname, hdr_info);

    if(hdr_info.datatype == 2) 
      make_grid<uchar>(argc, argv);
    else
      {
            cout << "Input file has wrong datatype: only 8 bpp is allowed. Exiting... " << endl;
            exit(1);
      }
   
   uchar *mask = new uchar[volume];
   if(!mask) {
       cout << "Memory allocation failed. Exiting..." << endl;
       exit(1);
   }
   fill(mask, mask + volume, (uchar)BACKGROUND);
    
   //Get First Sagital slice
   int sagital_area=sz[1] * sz[2];
   uchar *sagital_plane=new uchar[sagital_area];
   if(!sagital_plane) 
     {
       cout << "Failed to allocate memory for trace masks. Exiting..." << endl;
       exit(1);
     }
   fill(sagital_plane, sagital_plane + sagital_area, BACKGROUND);

   int Slice1=0;
   int Slice2=0;


   int L_curv_row=0;
   int L_curv_col=0;
   int R_curv_row=0;
   int R_curv_col=0;

   int *curv_tbl=new int[sz[2]];

   //For Side ONE
   //============
   fill(curv_tbl,curv_tbl+sz[2],0);

   //Find the slice1
   for(int x=0; x < sz[0]; ++x)
     {
       for(int y=0; y < sz[1]; ++y)
	 for(int z=0; z < sz[2]; ++z)
	   sagital_plane[z*sz[1] + y]=img[z*area+y*sz[0]+x];

       bool flag=0;

       for(int p=0;p< sagital_area;p++)
	 if(sagital_plane[p]==(uchar)1)
	   {
	       flag=1;
	       Slice1=x;
	       break;
	   }
      
       if(flag)
	 break;

     }
       
   cout <<"Left Mask at: "<<Slice1+1<<endl;

   //Scan for curv
   for(int z=0; z < sz[2]; z++)
     {
       for(int y=0; y < sz[1]-1; y++)
	 { 
	   int loc=z*sz[1]+y;
	   if(sagital_plane[loc]==0 && sagital_plane[loc+1]>0)
	     curv_tbl[z]=y+1;
	 }
     }

   //Find curv_col
   for(int z=0; z < sz[2]; z++)
     if(curv_tbl[z]>L_curv_col)
       L_curv_col=curv_tbl[z];

   //Find curv_row
   for(int z=0; z < sz[2]; z++)
     if(curv_tbl[z]==L_curv_col)
       L_curv_row=z;

   cout<<"Left Mask: curv_row="<<L_curv_row<<":: curv_col="<<L_curv_col<<endl;

   for(int y=0; y < sz[1]; ++y)
     for(int z=0; z < sz[2]; ++z)
       {
	 int pix=(int)sagital_plane[z*sz[1]+y];

	 if(pix > 0 && y < grid.pc)
	   sagital_plane[z*sz[1]+y]=(uchar)1;

	 if(pix > 0 && y >= grid.pc && y < grid.ac)
	   sagital_plane[z*sz[1]+y]=(uchar)2;

	 if(pix > 0 && y >= grid.ac && y < L_curv_col && z > L_curv_row)
	   sagital_plane[z*sz[1]+y]=(uchar)3;

	 if(pix > 0 && y >= L_curv_col)
	   sagital_plane[z*sz[1]+y]=(uchar)4;

	 if(pix > 0 && y >= grid.ac && y < L_curv_col && z <= L_curv_row)
	   sagital_plane[z*sz[1]+y]=(uchar)5;
       }

    // now copy trace_masks all the way in
   int startS=grid.sagital[3]+(int)((1.0/5.0)*(double)(grid.m-grid.sagital[3]));
   copy_trace2(sagital_plane, mask, startS, grid.m,0);

   //For Side TWO
   //============
   fill(curv_tbl,curv_tbl+sz[2],0);
   fill(sagital_plane, sagital_plane + sagital_area, BACKGROUND);

   //Find the slice1
   for(int x=0; x < sz[0]; ++x)
     {
       for(int y=0; y < sz[1]; ++y)
	 for(int z=0; z < sz[2]; ++z)
	   sagital_plane[z*sz[1] + y]=img[z*area+y*sz[0]+x];

       bool flag=0;

       for(int p=0;p< sagital_area;p++)
	 if(sagital_plane[p]==(uchar)2)
	   {
	       flag=1;
	       Slice2=x;
	       break;
	   }
      
       if(flag)
	 break;

     }
       
   cout <<"Right Mask at: "<<Slice2+1<<endl;

 //Scan for curv
   for(int z=0; z < sz[2]; z++)
     {
       for(int y=0; y < sz[1]-1; y++)
	 { 
	   int loc=z*sz[1]+y;
	   if(sagital_plane[loc]==0 && sagital_plane[loc+1]>0)
	     curv_tbl[z]=y+1;
	 }
     }

   //Find curv_col
   for(int z=0; z < sz[2]; z++)
     if(curv_tbl[z]>R_curv_col)
       R_curv_col=curv_tbl[z];

   //Find curv_row
   for(int z=0; z < sz[2]; z++)
     if(curv_tbl[z]==R_curv_col)
       R_curv_row=z;

   cout<<"Right Mask: curv_row="<<R_curv_row<<":: curv_col="<<R_curv_col<<endl;

   for(int y=0; y < sz[1]; ++y)
     for(int z=0; z < sz[2]; ++z)
       {
	 int pix=(int)sagital_plane[z*sz[1]+y];

	 if(pix > 0 && y < grid.pc)
	   sagital_plane[z*sz[1]+y]=(uchar)1;

	 if(pix > 0 && y >= grid.pc && y < grid.ac)
	   sagital_plane[z*sz[1]+y]=(uchar)2;

	 if(pix > 0 && y >= grid.ac && y < R_curv_col && z > R_curv_row)
	   sagital_plane[z*sz[1]+y]=(uchar)3;

	 if(pix > 0 && y >= R_curv_col)
	   sagital_plane[z*sz[1]+y]=(uchar)4;

	 if(pix > 0 && y >= grid.ac && y < R_curv_col && z <= R_curv_row)
	   sagital_plane[z*sz[1]+y]=(uchar)5;
       }

    // now copy trace_masks all the way in
   int endS=grid.sagital[5]-(int)((1.0/5.0)*(double)(grid.sagital[5]-grid.m));
  
   copy_trace2(sagital_plane, mask, grid.m, endS,1);
    
    // Write output
    hdr_info.max = 10;
    hdr_info.min = 0;
    hdr_info.datatype = 2;
    cout << "Writing: "<< argv[2]<<"...";
    write_analyze(string(argv[2]), mask,  hdr_info);
    cout <<"Done!"<<endl;
    
    delete [] mask;
    delete [] sagital_plane;
    delete [] curv_tbl;

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
void copy_trace2(uchar *trace_mask, uchar *mask, int beg, int end,uchar color_offset) 
{
    int x, y, z;
    uchar *buf=new uchar[sz[1]*sz[2]];
    memcpy(buf,trace_mask,sz[1]*sz[2]);

    //Put Color
    for(z=0; z < sz[2]; ++z)
        for(y=0; y < sz[1]; ++y)
	  {
	    if(buf[z*sz[1] + y])
	      buf[z*sz[1] + y]=buf[z*sz[1] + y]+color_offset;
	  }

    //copy
    for(z=0; z < sz[2]; ++z)
        for(y=0; y < sz[1]; ++y) 
	  {
            
            for( x=beg; x != end; x ++)
                mask[z*area + y*sz[0] + x] = buf[z*sz[1] + y];
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
