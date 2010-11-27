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

// This Class was written by : Dr. Azhar Quddus
// Nov 21, 2002.
// All rights reserved. Dept. of Cognitive Neurology, Sunnybrook & Womens'
// Toronto, Canada.

//Compiled with: g++ -c -I. MyVolume.cpp -o MyVolume.o
//Making static library with : ar -rc libmyVolume.a MyVolume.o
//Ranlib: ranlib libmyVolume.a

#include <iostream>    //For memcpy()
#include "dbh.h"       // Analyze Read/Write structures
#include "MyVolume.h"
#include "connexe.h"   //3-D Connected-Component-Labeling
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdio.h>

using namespace std;

CLVolume::CLVolume()
{
  uchar * mVBuff=0;
  double mVxSize[3]={0.0, 0.0, 0.0};
  ushort mVolSize[3]={0,0,0};
  int mImax=0;
  int mImin=0;
  short mDataType=0;
}

CLVolume::CLVolume(const CLVolume &rhs)
{
  long sz=rhs.mVolSize[0]*rhs.mVolSize[1]*rhs.mVolSize[2]*rhs.mDataType/2;
  delete [] mVBuff;    //delete old buff
  mVBuff=new uchar[sz];
  memcpy(mVBuff,rhs.mVBuff,sz);

  memcpy(mVxSize,rhs.mVxSize,3*sizeof(double));
  memcpy(mVolSize,rhs.mVolSize,3*sizeof(ushort));
  mImax=rhs.mImax;
  mImin=rhs.mImin;
  mDataType=rhs.mDataType;
}

CLVolume::~CLVolume()
{
  if(mVBuff)
    {
      delete [] mVBuff;
      mVBuff=0;
    }
}

bool CLVolume::read_analyze_header(const string &hdr_fname, CLVolume &vol)
{
  bool flag=FALSE;
    ifstream in(hdr_fname.c_str());
    if(!in)
      {
	flag=FALSE; 
	return flag;
      }

    in.seekg(42);
    in.read(reinterpret_cast<char *>(vol.mVolSize), 3*sizeof(short));
    in.seekg(80);
    in.read(reinterpret_cast<char *>(vol.mVxSize), 3*sizeof(float));
    in.seekg(70);
    in.read(reinterpret_cast<char *>(&(vol.mDataType)), sizeof(short));
    in.seekg(140);
    in.read(reinterpret_cast<char *>(&(vol.mImax)), sizeof(int));
    in.read(reinterpret_cast<char *>(&(vol.mImin)), sizeof(int));

    flag=TRUE;
    return flag;
}

bool CLVolume::read_analyze(const string &froot)
{
    bool flag=FALSE;
    //    CLVolume Tinfo;

    flag=read_analyze_header(froot + ".hdr", *this);
    if(flag==FALSE)
      return flag; 

    //*this = Tinfo;
     
    // make sure that image is either 8 or 16 bit,
    if(!(this->mDataType == 2 || this->mDataType == 4)) 
      {
	flag=FALSE;
	return flag;
      }
       
    long volume = this->mVolSize[0] * this->mVolSize[1] * this->mVolSize[2];
    
    this->mVBuff  = new uchar[volume*this->mDataType/2];

    if( !this->mVBuff ) 
      {
	flag=FALSE;
	return flag;
      }
    
    ifstream in((froot + ".img").c_str());
    if( !in ) 
      {
	flag=FALSE;
	return flag;
      }
    
    // Read properly
    in.read((char *)this->mVBuff, volume * this->mDataType/2);

    flag=TRUE;
    return flag;
  
}

bool CLVolume:: make_hdr(const string &hdr_fname, const CLVolume &vol)
{
    int i;
    bool flag=FALSE;

    struct dsr hdr;
    FILE *fp;
    
    memset(&hdr,0, sizeof(struct dsr));
    for(i=0;i<8;i++) hdr.dime.pixdim[i]=0.0;
   
    hdr.dime.vox_offset = 0.0;
    hdr.dime.roi_scale   = 1.0;
    hdr.dime.funused1    = 0.0;
    hdr.dime.funused2    = 0.0;
    hdr.dime.cal_max     = 0.0;
    hdr.dime.cal_min     = 0.0;
  
   
    hdr.dime.datatype = vol.mDataType;
    hdr.dime.bitpix = vol.mDataType * 4;
   
 
    if((fp=fopen(hdr_fname.c_str(),"w"))==0)
      {
	flag=FALSE;
	return flag;
      }

    hdr.dime.dim[0] = 4;  /* all Analyze images are taken as 4 dimensional */
    hdr.hk.regular = 'r';
    hdr.hk.sizeof_hdr = sizeof(struct dsr);

    hdr.dime.dim[1] = vol.mVolSize[0];  /* slice width  in pixels */
    hdr.dime.dim[2] = vol.mVolSize[1];  /* slice height in pixels */
    hdr.dime.dim[3] = vol.mVolSize[2];  /* volume depth in slices */
    hdr.dime.dim[4] = 1;  /* number of volumes per file */

    hdr.dime.glmax  = vol.mImax;  /* maximum voxel value  */
    hdr.dime.glmin  = vol.mImin;  /* minimum voxel value */

    hdr.dime.pixdim[1] = vol.mVxSize[0]; /* voxel width in mm per pixel */ 
    hdr.dime.pixdim[2] = vol.mVxSize[1]; /* voxel height in mm per pixel */
    hdr.dime.pixdim[3] = vol.mVxSize[2]; /* voxel depth in mm per pixel*/

/*   Assume zero offset in .img file, byte at which pixel
       data starts in the image file */

    hdr.dime.vox_offset = 0.0; 
    
/*   Planar Orientation;    */
/*   Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
     Movie flag ON:  3 = transverse, 4 = coronal, 5 = sagittal  */  

    hdr.hist.orient     = 0;  
    
/*   up to 3 characters for the voxels units label; i.e. 
        mm., um., cm.               */

    strcpy(hdr.dime.vox_units," ");
   
/*   up to 7 characters for the calibration units label; i.e. HU */

    strcpy(hdr.dime.cal_units," ");  
    
/*     Calibration maximum and minimum values;  
       values of 0.0 for both fields imply that no 
    calibration max and min values are used    */

    hdr.dime.cal_max = 0.0; 
    hdr.dime.cal_min = 0.0;

    fwrite(&hdr,sizeof(struct dsr),1,fp);
    fclose(fp);

    flag=TRUE;
    return flag;
}

bool CLVolume::write_analyze(const string &froot)
{
    bool flag=FALSE;
    make_hdr(froot + ".hdr", *this);
    
    ofstream out((froot + ".img").c_str());
    if( !out ) 
      {
	flag=FALSE;
	return flag;
      }

    int Tvol=this->mVolSize[0] * this->mVolSize[1] * this->mVolSize[2];

    out.write((char *)this->mVBuff, Tvol* this->mDataType/2); 

    flag=TRUE;
    return flag;    
}

bool CLVolume:: GetImagePlane(CLImage &img, int Slice)
{
  bool flag=FALSE;

  if(!mVBuff)
    return flag;

  int r=mVolSize[0];
  int c=mVolSize[1];
  int N=mVolSize[2];

  uchar *lbuff=new uchar[r*c*mDataType/2];

  if(!lbuff)
    return flag;

  long idx=r*c*(Slice-1)*mDataType/2;
  memcpy(lbuff,(mVBuff+idx),r*c*mDataType/2);

  if(!img.SetImage(lbuff,r,c,mDataType))
    return flag;

  delete [] lbuff;
 
  flag=TRUE;

  return flag;
}

bool CLVolume:: PutImagePlane(CLImage &img, int Slice)
{
  bool flag=FALSE;

  int r=mVolSize[0];
  int c=mVolSize[1];
  int N=mVolSize[2];

  if(!mVBuff)
    return flag;

  long idx=r*c*(Slice-1)*mDataType/2;

  uchar *buf=new uchar[r*c*mDataType/2];

  if(!img.ImageBuff(buf))
    return flag;

  memcpy((mVBuff+idx),buf,r*c*mDataType/2);

  delete [] buf;

  flag=TRUE;
  
  return flag;
}

bool CLVolume::DoMask(const CLVolume LabV,CLVolume *out,const int val)
{
  bool flag=FALSE;
 
  int lab=0;
  
  if(!CheckVolSize(LabV))
    return flag;

  long vol=this->mVolSize[0]*this->mVolSize[1]*this->mVolSize[2];
 
  *out=*this;

  int lf=LabV.mDataType/2;
  int tf=this->mDataType/2;
  int of=out->mDataType/2;

  uchar *op=out->mVBuff;
  uchar *lp=LabV.mVBuff;
 
  for(int i=0;i<vol;i++)
    {
      if(lf>1)
	memcpy(&lab,(lp+i*lf),lf);
      else
	lab=(int)*(lp+i);
      
      if(lab==val)
	memset((op+i*of),0,of);

      lab=0;     //clear lab... important!
    }

  flag=TRUE;
  return flag;
}

bool CLVolume::CheckVolSize(const CLVolume &vol)
{
  bool flag=FALSE;
  if(this->mVolSize[0]!=vol.mVolSize[0] || this->mVolSize[1]!=vol.mVolSize[1] || this->mVolSize[2]!=vol.mVolSize[2])
    return flag;

flag=TRUE;
return flag;

}

ushort CLVolume::NumberOfSlices()
{
  if(mVBuff)
    return mVolSize[2];
  else
    return 0;
}

hdr CLVolume::GetHeader()
{
  hdr myhdr;    //make new header

  myhdr.VxSize[0]=mVxSize[0];
  myhdr.VxSize[1]=mVxSize[1];
  myhdr.VxSize[2]=mVxSize[2];

  myhdr.VolSize[0]=mVolSize[0];
  myhdr.VolSize[1]=mVolSize[1];
  myhdr.VolSize[2]=mVolSize[2];

  myhdr.Imax=mImax;
  myhdr.Imin=mImin;

  myhdr.DataType=mDataType;

  return myhdr;
}

bool CLVolume::SetHeader(const hdr &myhdr)
{
  bool flag=TRUE;

  this->mVxSize[0]=myhdr.VxSize[0];
  this->mVxSize[1]=myhdr.VxSize[1];
  this->mVxSize[2]=myhdr.VxSize[2];

  this->mVolSize[0]=myhdr.VolSize[0];
  this->mVolSize[1]=myhdr.VolSize[1];
  this->mVolSize[2]=myhdr.VolSize[2];

  this->mImax=myhdr.Imax;
  this->mImin=myhdr.Imin;

  this->mDataType=myhdr.DataType;

  return flag;
}

bool CLVolume::GetVolBuff(uchar *buff)
{
  bool flag=FALSE;
  long vol=mVolSize[0]*mVolSize[1]*mVolSize[2];

  if(!buff || !mVBuff)
    return flag;

  memcpy(buff,mVBuff,vol*mDataType/2);

  flag=TRUE;
  return flag;
}

bool CLVolume::SetVolBuff(const uchar *buff,long sz)
{
  bool flag=FALSE;

  if(!buff)
    return flag;

  if(sz!=(mVolSize[0]*mVolSize[1]*mVolSize[2]*mDataType/2))
    return flag;

  /*
  if(!mVBuff || mVBuff==0)
    {
      cout<<"mVBuff not present."<<endl;
      mVBuff=new uchar[sz];
    }
  */

  mVBuff=new uchar[sz];
  memcpy(mVBuff,buff,sz);

  flag=TRUE;
  return flag;
}

int CLVolume::Label3D(CLVolume &out)
{
  if(mDataType!=2)
    return -1;

  int dimV[3];
  dimV[0]=this->mVolSize[0];
  dimV[1]=this->mVolSize[1];
  dimV[2]=this->mVolSize[2];

  //Prepare output Volume
  out.mVxSize[0]=this->mVxSize[0];
  out.mVxSize[1]=this->mVxSize[1];
  out.mVxSize[2]=this->mVxSize[2];

  out.mVolSize[0]=this->mVolSize[0];
  out.mVolSize[1]=this->mVolSize[1];
  out.mVolSize[2]=this->mVolSize[2];

  out.mDataType=(ushort)4;

  long vol=dimV[0]*dimV[1]*dimV[2];
  ushort *op=new ushort [vol];

  Connexe_SetConnectivity((int)6);
  Connexe_SetMinimumSizeOfComponents((int)3);   //discard less than 3 voxels
  Connexe_SetMaximumNumberOfComponents((int)65536); //keep 65536 components
  int nLab=CountConnectedComponents(this->mVBuff,UCHAR,op,USHORT,dimV);

  out.mImin=0;
  out.mImax=nLab;

  out.mVBuff=new uchar[out.mDataType/2*vol];

  memcpy(out.mVBuff,op,sizeof(ushort)*vol);

  delete [] op;
  return nLab;
}

bool CLVolume::DispSlice(int n)
{
  bool flag=FALSE;
  int xl=(int)this->mVolSize[0];
  int yl=(int)this->mVolSize[1];
  short typ=(short)this->mDataType;

  CLImage img(xl,yl,typ);

  if(!this->GetImagePlane(img,n))
    {
      cout <<"Couldnt read Slice!"<<endl;
      return flag;
    }

  if(!img.DispImg())
    {
      cout <<"Display Failed!"<<endl;
      return flag;
    }

  flag=TRUE;
  return flag;
}
