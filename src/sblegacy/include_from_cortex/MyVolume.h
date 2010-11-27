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

#ifndef _VOLUME_H_
#define _VOLUME_H_

#include <string>
#include "MyImage.h"

using namespace std;

#define TRUE 1
#define FALSE 0

typedef unsigned char uchar;
typedef unsigned short ushort;

struct hdr
{ float VxSize[3];     // voxel size in mm
  ushort VolSize[3];    // image dimensions
  int Imax;             // absolute max and min of image data [min max]
  int Imin;     
  short DataType;
};

class CLVolume
{
 private:
  uchar *mVBuff;
  float mVxSize[3];     // voxel size in mm
  ushort mVolSize[3];    // image dimensions
  int mImax;             // absolute max and min of image data [min max]
  int mImin;     
  short mDataType;

  bool read_analyze_header(const string &hdr_fname, CLVolume &vol);
  bool make_hdr(const string &hdr_fname, const CLVolume &vol);

 public:
  CLVolume();        //Beware, will not allocate any memory!
  CLVolume(const CLVolume &);   //Copy constructor.. will allocate memory.
  ~ CLVolume();

  // bool VolumeBuff(uchar &buf);  //Will allocate memory& give you buff
  bool read_analyze(const string &froot);
  bool write_analyze(const string &froot);
 
  bool CheckVolSize(const CLVolume &vol);
  ushort NumberOfSlices();
 
  //Volume Processing ...
  bool DoMask(const CLVolume LabV, CLVolume *out,const int val); //Mask val values from this volume (does memory allocation
  bool GetImagePlane(CLImage &img,int Slice);
  bool PutImagePlane(CLImage &img, int Slice); 
  hdr GetHeader();
  bool SetHeader(const hdr &myhdr);
  bool GetVolBuff(uchar *buff);  //will NOT alloc memory & copy private buff.
  bool SetVolBuff(const uchar *buff,long sz); //will initialize private buff if needed.
  int Label3D(CLVolume &out);     //3-D Labelling of binary volume
  bool DispSlice(int n);       //for debugging only
};

#endif
