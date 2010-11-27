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

//By: Dr. Azhar Quddus
// To compile:
// g++ -O2 -o cleanEdgeSeg cleanEdgeSeg.cpp /net/cortex/csf/pub/include/MyVolume.cpp \
//	/net/cortex/csf/pub/include/MyImage.cpp /net/cortex/csf/pub/include/connexe.c
// use at least gcc-3.2


#include<iostream>
#include <string>
//#include "MyImage.h"
#include "../../../include_from_cortex/MyVolume.h"

using namespace std;

/////////////Local Min/Max Functions (templated) FOR DEBUGGING /////////////
//****Max******
template <class T>
T bufMax(T *mybuf,int len)
{
  T mmax=(T)0;
  for(int i=0;i<len;i++)
    {
      if(mybuf[i]> mmax)
	mmax=mybuf[i];
    }
  return mmax;
}

//****Min*******
template <class T>
T bufMin(T *mybuf,int len)
{
  T mmin=(T)255;
  for(int i=0;i<len;i++)
    {
      if(mybuf[i] < mmin)
	mmin=mybuf[i];
    }
  return mmin;
}

///////////////////////////////////////

int main(int argc, char *argv[])
{
  if(argc!=5)
    {
      cout <<"Usage:   cleanEdgeSeg <Segm. binary Volume> <T1Segm> <Output> <VCSFvalue>"<< endl;
      cout<<"By: Dr. A. Quddus"<< endl;
      return -1;
    }

  string inFname=argv[1];
  string T1SegFname=argv[2];
  string outFname=argv[3];
  int vcsf_value = atoi( argv[4] );

  CLVolume inV,T1V,outV;

  if(!T1V.read_analyze(T1SegFname) || !inV.read_analyze(inFname))
    {
      cout <<"Couldnt Read input Files! Exiting..."<<endl;
      return -1;
    }

  if(!T1V.CheckVolSize(inV))
    {
      cout <<"Size Mismatch... exiting!"<<endl;
      return -1;
    }
  
  //for Debugging only
  //inV.DispSlice(80);

  ushort N=inV.NumberOfSlices();
  cout<<"N="<<N<<endl;
  
  outV=T1V;

  //Prepare output Header
  hdr myhdr=outV.GetHeader();
  myhdr.Imax=1;
  myhdr.Imin=0;

  //do the stuff
  CLImage imgIn(256,256,BPP8);
  CLImage tmp1(256,256,BPP8);
  CLImage tmp2(256,256,BPP8);
  CLImage imgOut(256,256,BPP8);

  for(int i=1;i<=N;i++)
    {
      cout <<"Processing Slice# "<<i<<"...";

      if(!inV.GetImagePlane(imgIn,i))
	{
	  cout <<"Couldnt get input slice#:"<<i<<endl;
	  break;
	}
  
      if(imgIn.IsBlank())
	{
	  cout<<"Blank Slice ??"<<endl;
	  imgOut=imgIn;
	  outV.PutImagePlane(imgOut,i);
	  continue;
	}

      if(!T1V.GetImagePlane(tmp1,i))
	{
	  cout <<"Couldnt get T1 slice#:"<<i<<endl;
	  return -1;
	}

      //binarize tmp1
      int xl=0;
      int yl=0;

      tmp1.GetImgSize(&xl,&yl);

      uchar *binT1=new uchar[xl*yl];
      if(!tmp1.ImageBuff(binT1))
	{
	  cout <<"T1 Buff Read failed!"<<endl;
	  return 1;
	}
    
      for(int p=0;p<(xl*yl);p++)
	{
	  if(binT1[p])
	    binT1[p]=(uchar)1;
	}

      //Binary erosion
      tmp2=tmp1;
      tmp2.SetImage(binT1,xl,yl,BPP8);

      if(!tmp2.BinErode(tmp1,vcsf_value))
	{
	  cout <<"BinErode Failed !..."<<endl;
	  continue;
	}

      if(!tmp1.BinErode(tmp2,vcsf_value))
	{
	  cout <<"BinErode Failed !..."<<endl;
	  continue;
	}

      if(!tmp2.BinErode(tmp1,vcsf_value))
	{
	  cout <<"BinErode Failed !..."<<endl;
	  continue;
	}
      
      //Image substract
      uchar *imgSub=new uchar [xl*yl];
      memset(imgSub,0,(xl*yl));

      uchar *imgErode=new uchar[xl*yl];

      if(!tmp1.ImageBuff(imgErode))
	{
	  cout<<"erode buff read failed!"<<endl;
	  return 1;
	}

      for(int p=0;p<(xl*yl);p++)
      	imgSub[p] =(uchar)((ushort)binT1[p] - (ushort)imgErode[p]);

      //Labelling input Image (binary)
      CLImage labImg(xl,yl,BPP16);     //make 16-bit Labeled Image (all zeros)

      int M=imgIn.Label_Connected_Regions(labImg)-1;
    
      cout << "No of Labels= " << M << "\n";

      //printf("No of Labels= %d\n",M);

      ushort *labp=new ushort[xl*yl];
      if(!labImg.ImageBuff16(labp))
	{
	  cout<<"cant get labeled buff!"<<endl;
	  return -1;
	}

      //Generate index table for the objects at the edge
      bool *labTbl=new bool[M];
      for(int idx=0;idx < M;idx++)   //clear table
	labTbl[idx]=(bool)0;

      for(int idx=0; idx < (xl*yl);idx++)   //fill table
	{
	  if(labp[idx]>(ushort)0 && imgSub[idx]>(uchar)0)
	    labTbl[labp[idx]-1]=(bool)1;
	}
      //cout<<"Max= "<<(int)bufMax(labTbl,M)<<"Min= "<<(int)bufMin(labTbl,M)<<endl;
      //cleaning
      uchar *pp=new uchar[xl*yl];

      if(!imgIn.ImageBuff(pp))
	{
	  cout<<"couldnt get buffer!"<<endl;
	  return -1;
	}
      //cout<<"Max= "<<(int)bufMax(imgSub,xl*yl)<<"Min= "<<(int)bufMin(imgSub,xl*yl)<<endl;
      for(int tidx=0;tidx < M; tidx++)
	{
	  if(labTbl[tidx])
	    {
	      for(int idx=0; idx < (xl*yl);idx++)
		{
		  if((int)labp[idx]==(tidx+1))
		    {//cout<<"Modifying here..."<<endl;
		    pp[idx]=(uchar)0;
		    }
		}
	    }
	}

      if(!imgOut.SetImage(pp,xl,yl,BPP8))
	{
	  cout <<"SetImage for output buf failed!"<<endl;
	  return -1;
	}

      if(!outV.PutImagePlane(imgOut,i))
      	{
      	  cout <<"Couldnt Prepare output Volume- Slice#"<<i<<endl;
      	  break;
      	}

      delete [] labTbl;
      delete [] labp;
      delete [] pp;
      delete [] binT1;
      delete [] imgErode;
      delete [] imgSub;
    }

  if(!outV.SetHeader(myhdr))
    {
      cout<<"Output header writing failed!"<<endl;
      return -1;
    }

  if(!outV.write_analyze(outFname))
    {
      cout<<"Output file writing failed!"<<endl;
      return -1;
    }

  cout <<"Done."<<endl;
  exit( 1 );
  return 0;
}
