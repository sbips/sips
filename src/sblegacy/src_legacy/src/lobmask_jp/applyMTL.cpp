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
#include <string>
#include <iostream>
#include <iomanip.h>
//#include "/home/sela/codeLib/c++/gs-2.0/include/gs/command_line.h"

//By: Dr. Azhar Quddus
// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o applyMTL applyMTL.cpp

using namespace nk;

typedef unsigned char uchar;

uchar  BACKGROUND = 0;

#define W 3
#define G 4
#define SC 5
#define VC 7

int main(int argc, char *argv[])
{
  double RMTL=0;
  double LMTL=0;
  int segFlag=0;

    if( argc < 4 || argc > 5)
        {
            cout <<"Usage: applyMTL <lobmask file> <MTL mask file> <Segmentation File(Optional)> <OutPut File>"  << endl;
            exit(1);
        }

    string lobmask_fname=argv[1];
    string mtlmask_fname=argv[2];
    string segmask_fname,output_fname;

    if(argc==5)
      {
          segmask_fname=argv[3];
          output_fname=argv[4];
	  segFlag=1;
      }
    else
      {
	output_fname=argv[3];
      }

    info lobmask_info, mtlmask_info, segmask_info, output_info;

    uchar *lobmask=0, *mtlmask=0, *segmask=0, *output=0;

    // read all relevant images

    lobmask = read_analyze<uchar>(lobmask_fname, lobmask_info);
    mtlmask = read_analyze<uchar>(mtlmask_fname, mtlmask_info);
    if(segFlag)
	segmask= read_analyze<uchar>(segmask_fname, segmask_info);

    // check for size
    if( !equal(lobmask_info.sz, lobmask_info.sz+3,mtlmask_info.sz))
    {
         cout << "MTLMask: Size Mismatch!  Exiting ... " << endl;
        exit(1);
    }

   if(segFlag && !equal(lobmask_info.sz, lobmask_info.sz+3,segmask_info.sz))
    {
         cout << "SegMask: Size Mismatch!  Exiting ... " << endl;
        exit(1);
    }

    long vol = lobmask_info.sz[0] * lobmask_info.sz[1] * lobmask_info.sz[2];
    double vox_vol = lobmask_info.voxsz[0] * lobmask_info.voxsz[1] * lobmask_info.voxsz[2];

    //copy
    output=new uchar[vol];
    memcpy(output,lobmask,vol);

    //Do the job
    if(!segFlag)
      {
        for(int i = 0; i < vol; i++)
        {
	  int val=(int)*(lobmask+i);

	  if (*(mtlmask+i)>0)
	    {
	      *(output+i)=0;
	      if(val>=3 && val<=15)
	        LMTL++;
	      if(val>=16 && val<=28)
	        RMTL++;
	    }
        }
         //Write Output
         cout <<"Left MTL Vol :"<<(double)(LMTL*vox_vol)<< endl;
         cout <<"Right MTL Vol:"<<(double)(RMTL*vox_vol)<< endl;
      }
    else
      {
	double lw=0,lg=0,lsc=0,lvc=0,rw=0,rg=0,rsc=0,rvc=0;
        for(int i = 0; i < vol; i++)
        {
	  int val=(int)*(lobmask+i);

	  if (*(mtlmask+i)>0)
	    {
	      *(output+i)=0;
	      if(val>=3 && val<=15)
		{
		  LMTL++;
		  if(*(segmask+i)==W) lw++;
		  if(*(segmask+i)==G) lg++;
		  if(*(segmask+i)==SC) lsc++;
		  if(*(segmask+i)==VC) lvc++;
		}
	      if(val>=16 && val<=28)
		{
		  RMTL++;
		  if(*(segmask+i)==W) rw++;
		  if(*(segmask+i)==G) rg++;
		  if(*(segmask+i)==SC) rsc++;
		  if(*(segmask+i)==VC) rvc++;
		}
	    }
        }
         //Write Output
	//	cout <<endl<<"Left MTL Vol :"<<(double)(LMTL*vox_vol)<< endl;
	//cout <<"Right MTL Vol:"<<(double)(RMTL*vox_vol)<< endl<<endl;
	//cout <<"Left MTL:"<<endl<<"========="<<endl;

	 cout <<endl<<"Left/Right"<<setw(20);
	 cout <<"Total_LobMaskVol."<<setw(15);
	 cout <<"White"<<setw(15);
	 cout <<"Gray"<<setw(15);
	 cout <<"SCSF"<<setw(15);
	 cout <<"VCSF"<<setw(20);
	 cout <<"Total_SegmVol."<< endl;

	 cout <<"Left-MTL:"<< setw(15);
	 cout <<(double)(LMTL*vox_vol)<<setw(20);
	 cout <<(double)(lw*vox_vol)<< setw(15);
	 cout <<(double)(lg*vox_vol)<< setw(15);
	 cout <<(double)(lsc*vox_vol)<< setw(15);
	 cout <<(double)(lvc*vox_vol)<< setw(20);
	 cout <<(double)(lw*vox_vol+lg*vox_vol+lsc*vox_vol+lvc*vox_vol)<<endl;

	 cout <<"Right-MTL:"<< setw(15);
	 cout <<(double)(RMTL*vox_vol)<<setw(20);
	 cout <<(double)(rw*vox_vol)<< setw(15);
	 cout <<(double)(rg*vox_vol)<< setw(15);
	 cout <<(double)(rsc*vox_vol)<< setw(15);
	 cout <<(double)(rvc*vox_vol)<< setw(20);
	 cout <<(double)(rw*vox_vol+rg*vox_vol+rsc*vox_vol+rvc*vox_vol)<<endl<<endl;
      }

    //Write Output File
    output_info=lobmask_info;
    write_analyze(output_fname, (char*)output, output_info);

    delete [] lobmask;
    delete [] mtlmask;
    if(segFlag) delete [] segmask;
    delete [] output;

    return 0;
}
