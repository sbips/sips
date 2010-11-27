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
// /srutility/gnu/gcc-2.95.2/bin/g++ -O2 -o lobmaskCF2 lobmaskCF2.cpp
//By: Dr. Azhar Quddus, August 27, 2003.

#include <assert.h>
#include <iomanip.h>
#include <vector>
#include <algorithm>
#include "globals.h"
#include "grid.h"
#include "make_grid.h"
#include "bresenhem.h"
//#include "tracings.cpp"
#include "floodfill.cpp"
//#include "basal_ganglia.cpp"
//#include "orbital_frontal.cpp"


void copy_trace_mask(uchar *trace_mask, uchar *mask, int xstride,
                     int beg, int med, int end, uchar lobe_offset);
void copy_trace2(uchar *trace_mask, uchar *mask, int beg, int end);
void fill_top_and_bottom_slices(uchar *mask);
void fill_front_and_back_slices(uchar *mask);
void fill_right_and_left_slices(uchar *mask);
void get_Axial_Slice(uchar *slice, uchar *vol,int N);
void set_Axial_Slice(uchar *slice, uchar *vol,int N);
bool Blank_Slice(uchar *slice, int iarea);
int MaxWidth(uchar *slice,int r, int c);
bool FindLine(uchar *slice,int R, int C, int r);
//////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc != 5 )
{      
	cout <<"Usage: lobmaskCF2 <BasalGangliaThalamus_File(BGT)> <OutputFile> -g GridFile"<<endl;
	cout<<"To generate mask for Collenergic Fibres. By: Dr. Azhar Quddus"<<endl;
	exit(1);
      }
    
    info hdr_info;
    string fname(argv[1]);
//    string cfname(argv[2]);

    //Read Header & Make Grid info
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
   
    int xgrid=grid.sagital[2]-grid.sagital[1];

    //Read Image files
    uchar *img=read_analyze<uchar>(fname, hdr_info);
 //   uchar *cing=read_analyze<uchar>(cfname, hdr_info);

    //Make Blank Volume
   uchar *mask = new uchar[volume];
   if(!mask) {
       cout << "Memory allocation failed. Exiting..." << endl;
       exit(1);
   }
   fill(mask, mask + volume, (uchar)BACKGROUND);

   //Copy Singulate into output mask
//   memcpy(mask,cing,volume);

    //Make Axial slice
   int axial_area=sz[0] * sz[1];
   uchar *axial_in=new uchar[axial_area];
   uchar *axial_out=new uchar[axial_area];

   if(!axial_in || !axial_out)
     {
       cout << "Failed to allocate memory for trace masks. Exiting..." << endl;
       exit(1);
     }

   fill(axial_in, axial_in + axial_area, BACKGROUND);

   //Find Strip end (from top)
   int oldWidth=0;
   int Width=0;
   int slice_break=0;
   int breakSlicesfromTop=0;

   for(int i=sz[2]-1;i>=0;i--)
     {
       get_Axial_Slice(axial_in,img,i);
       if(!Blank_Slice(axial_in,axial_area))
	 {
	   Width=MaxWidth(axial_in,sz[1], sz[0]);
	   breakSlicesfromTop++;
	 }

       if(oldWidth==0)
	 oldWidth=Width;

       if(oldWidth!=0 && Width!=0 && oldWidth!=Width)
	 {
	   slice_break=i+1;
	   break;
	 }
     }
   cout<<"Slice Break 1:"<<slice_break<<endl;

   int stSlices=3 * breakSlicesfromTop;

   //Find full curve end (from top)
   int slice_break2=0;
   for(int i=slice_break;i>=0;i--)
     {
       get_Axial_Slice(axial_in,img,i);
       bool line=true;
       int r=grid.coronal[9];

       if(!FindLine(axial_in,sz[1],sz[0],r))
	 {
	   slice_break2=i+1;
	   break;
	 }
     }
   cout<<"Slice Break 2:"<<slice_break2<<endl;


   cout<<"Valid Slices:";

//Processing Half-Curved regions...
   for(int i=0;i<slice_break2;i++)
     {
       int LastRow=0;

       get_Axial_Slice(axial_in,img,i);
       get_Axial_Slice(axial_out,mask,i);

       if(!Blank_Slice(axial_in,axial_area))
	 {cout<<i<<" " ;
	   for(int r=0;r<sz[1];r++)
	     for(int c=0;c<sz[0];c++)
	       {
		 int pos=r*sz[0]+c;
		 if(axial_in[pos]==0 && axial_in[pos+1]>0)
		   {
		     int C=c+(int)((double)xgrid/4.0);
		     for(int lc=C;lc > (C-xgrid/2);lc--)
		       axial_out[r*sz[0]+lc]=12;
		   }

		 if(axial_in[pos]==0 && axial_in[pos-1]>0)
		   {
		     int C=c-(int)((double)xgrid/4.0);
		     for(int lc=C;lc < (C+xgrid/2);lc++)
		       axial_out[r*sz[0]+lc]=14;
		   }
	       }
	 }

       set_Axial_Slice(axial_out, mask,i);
     }

   //Processing Full-Curved regions...
   for(int i=slice_break2;i<slice_break;i++)
     {
       get_Axial_Slice(axial_in,img,i);
       get_Axial_Slice(axial_out,mask,i);

       if(!Blank_Slice(axial_in,axial_area))
	 {cout<<i<<" " ;
	   for(int r=0;r<sz[1];r++)
	     for(int c=0;c<sz[0];c++)
	       {
		 int pos=r*sz[0]+c;
		 if(axial_in[pos]==0 && axial_in[pos+1]>0)
		   {
int C=c+3;
		     for(int lc=C;lc > (C-xgrid/2);lc--)
		       axial_out[r*sz[0]+lc]=12;
		   }

		 if(axial_in[pos]==0 && axial_in[pos-1]>0)
		   {
int C=c-3;
		     for(int lc=C;lc < (C+xgrid/2);lc++)
		       axial_out[r*sz[0]+lc]=14;
		   }
	       }
	 }

       set_Axial_Slice(axial_out, mask,i);
     }
   
   //Generating straight regions...
/*
   int LTR=grid.coronal[6];
   int LBR=grid.coronal[13];
   int LRC=(int)((double)grid.sagital[3]-(double)xgrid/2.0);
   int LLC=LRC-(int)((double)xgrid/2.0);
LRC=LRC+3;
LLC=LLC+3;

   int RTR=LTR;
   int RBR=LBR;
   int RLC=(int)((double)grid.sagital[5]+(double)xgrid/2.0);
   int RRC=RLC+(int)((double)xgrid/2.0);
RLC=RLC-3;
RRC=RRC-3;

   for(int i=slice_break;i < (slice_break + stSlices); i++)
     {
       get_Axial_Slice(axial_in,img,i);

       cout<<i<<" ";

       get_Axial_Slice(axial_out,mask,i);

       for(int r=LTR;r<=LBR;r++)
	 for(int c=LLC;c<=LRC;c++)
	   axial_out[r*sz[1]+c]=12;

       for(int r=RTR;r<=RBR;r++)
	 for(int c=RLC;c<=RRC;c++)
	   axial_out[r*sz[1]+c]=14;

       set_Axial_Slice(axial_out, mask,i);
     }
*/     
    // output
    hdr_info.max = 15;
    hdr_info.min = 0;
    hdr_info.datatype = 2;
    cout<<endl;
    cout << "Writing: "<< argv[2]<<"...";
    write_analyze(string(argv[2]), mask,  hdr_info);
    cout <<"Done!"<<endl;
    
    delete [] mask;
    delete [] img;
//    delete [] cing;
    delete [] axial_in;
    delete [] axial_out;
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

///////////////////////////////////////////////////////////////////
//Get Axial Slice
 void get_Axial_Slice(uchar *slice, uchar *vol,int N)
   {
     int iarea=sz[0]*sz[1];
     memcpy(slice,(vol+N*area),iarea);
   }

///////////////////////////////////////////////////////////////////
//Set Axial Slice
 void set_Axial_Slice(uchar *slice, uchar *vol,int N)
   {
     int iarea=sz[0]*sz[1];
     memcpy((vol+N*area),slice,iarea);
   }

///////////////////////////////////////////////////////////////////
//Blank Axial Slice
 bool Blank_Slice(uchar *slice, int iarea)
   {
     bool flag=true;
     for(int i=0;i<iarea;i++)
       {
	 if(slice[i]!=0)
	     flag=false;
       }

     return flag;
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

//////////////////////////////////////////////////////////////////////
int MaxWidth(uchar *slice,int r, int c)
{

  int MxWidth=0;
  for(int i=0;i<r;i++)
    {
      int width=0;
      for(int j=0;j<c;j++)
	{
	  if(slice[i*c+j]>0)
	    width++;
	}

      if(width>MxWidth)
	MxWidth=width;
    }

  return MxWidth;
}

///////////////////////////////////////////////////////////////////
bool FindLine(uchar *slice,int R, int C, int r)
{

  bool flag=false;

  for(int j=0;j<C;j++)
    if(slice[r*C+j]>0)
      {
	flag=true;
	break;
      }

  return flag;
}
