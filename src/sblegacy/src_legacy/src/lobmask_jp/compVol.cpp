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

#include "/net/cortex/csf/pub/include/analyze_utils.h"
#include <string>
#include <iostream>
#include <iomanip.h>

//#include "/home/sela/codeLib/c++/gs-2.0/include/gs/command_line.h"

// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o compVol compVol.cpp

using namespace nk;

//#include "lobe_struct.h"
typedef unsigned char uchar;
#define LOBS 29
 
int main(int argc, char **argv)
{
  double LobC[LOBS];
  double TotVol=0;

  for (int i=0;i<LOBS;i++)
    LobC[i]=0;

  string Lab[LOBS]={
/*0*/ "BACKGROUND",
/*1*/ " GRIDLINE",
/*2*/ "HIGHLIGHT",
/*3*/ "LSUPF",
/*4*/ "LIF",
/*5*/ "LOBF",
/*6*/ "LMOBF",
/*7*/ "LSP",
/*8*/ "LIP",
/*9*/ "LO",
/*10*/  "LAT",
/*11*/  "LPT",
/*12*/  "LABGT",
/*13*/  "LPBGT",
/*14*/  "LMSF",
/*15*/  "LMIF",
/*16*/  "RSUPF",
/*17*/  "RIF",
/*18*/  "ROBF",
/*19*/  "RMOBF",
/*20*/  "RSP",
/*21*/  "RIP",
/*22*/  "RO",
/*23*/  "RAT",
/*24*/  "RPT",
/*25*/  "RABGT",
/*26*/  "RPBGT",
/*27*/  "RMSF",
/*28*/  "RMIF"
};

    if( argc != 2 )
        {
            cout <<"
Usage: compVol <lobmask file>"  << endl;
            exit(1);
        }
           
    int sz[3];         
    float voxsize[3];

    string lobmask_fname=argv[1];
    info lobmask_info;
    uchar *lobmask=0;
    
    // read all relevant image
    
    lobmask = read_analyze<uchar>(lobmask_fname, lobmask_info);
    
    long volume = lobmask_info.sz[0] * lobmask_info.sz[1] * lobmask_info.sz[2];
    double vox_vol = lobmask_info.voxsz[0] * lobmask_info.voxsz[1] * lobmask_info.voxsz[2];
   
    for(int i=0;i<volume;i++)
      {
	int val=*(lobmask+i);
      LobC[val]++;
      }

    //ToTal Vol
    int j=0;

    for(j=3;j<LOBS;j++)
      TotVol= TotVol+LobC[j];

    TotVol= TotVol*vox_vol;

    cout <<"Report:"<< endl<<"======="<<endl;
    cout <<"Total Volume"<<"= "<< TotVol<< endl;

    for(j=0;j<LOBS;j++)
      cout <<Lab[j]<< "    "<<(double)(LobC[j]*vox_vol)<< endl;
   
    delete [] lobmask;
   
    return 0;
}
