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
// Nov 20, 2002.
// All rights reserved. Dept. of Cognitive Neurology, Sunnybrook & Womens'
// Toronto, Canada.

//Compiled with: g++ -c -I. MyImage.cpp -o MyImage.o
//Making static library with : ar -rc libmyImage.a MyImage.o
//Ranlib: ranlib libmyImage.a

#include <iostream>    //For memcpy()
#include <fstream>     //For file I/O
#include "MyList.h"    //Used in Labelling Algorithm
#include "MyImage.h"

using namespace std;

#define XLEN 256
#define YLEN 256

uchar  dialate_table[256] = {
	1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1,
	2, 0, 2, 2, 2, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2,
	2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 2, 2,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 3, 3,
	2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 2, 2,
	3, 0, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
	3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 3, 3,
	4, 0, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	3, 0, 3, 3, 3, 0, 3, 3, 0, 0, 3, 3, 0, 0, 3, 3,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 4, 4,
	3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 3, 3,
	4, 0, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 4, 4,
	5, 0, 5, 5, 5, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
};

CLImage::CLImage()
{
  int mXlen=0;
  int mYlen=0;

  uchar *mBuff=0;

  short mType=0;
}

CLImage::CLImage(int Xl, int Yl,short typ)
{
  long len=Xl*Yl;
  
  //delete [] mBuff;
  //mBuff=0;

  mBuff=new uchar[len*typ/2];
 
  memset(mBuff,0,len*typ/2);

  mXlen=Xl;
  mYlen=Yl;
  mType=typ;
}

CLImage::~CLImage()
{
  if(mBuff)
    {
      delete [] mBuff;
      mBuff=0;
    }
}

CLImage::CLImage(const CLImage &rhs)
{
  mXlen=rhs.mXlen;
  mYlen=rhs.mYlen;

  long sz=mXlen*mYlen*mType/2;
  delete [] mBuff;     //delete old buf if any
  mBuff=new uchar[sz];
  memcpy(mBuff,rhs.mBuff,sz);

  mType=rhs.mType;
}

bool CLImage::SetImage(uchar *buf,int Xl, int Yl,short typ)
{
  bool flag=FALSE;

  long len=Xl*Yl;
  
  if(!mBuff)
    return flag;

  if(!buf)
   return flag;
 
  memcpy(this->mBuff,buf,len*typ/2);

  mXlen=Xl;
  mYlen=Yl;
  mType=typ;

  flag=TRUE;
  return flag;
}

bool CLImage::IsImage()
{
  bool flag;
  if(!this->mBuff)
    flag=FALSE;
  else
    flag=TRUE;

  return flag;
}

bool CLImage::IsBlank()
{
  bool flag=TRUE;

  // if(!IsImage())
  //  return flag;

  int len=this->mXlen*this->mYlen;

  ushort *cp=new ushort;

  for(int i=0;i<len;i++)
    {
      if(mType>2)
	memcpy(cp,(mBuff+i*mType/2),mType/2);
      else
	*cp=(ushort)*(mBuff+i);

      if(*cp)
	{
	  flag=FALSE;
	  break;
	}
    }
 
  delete cp;

  return flag;
}

bool CLImage::GetImgSize(int *xlen,int *ylen)
{
  bool flag=FALSE;
  if(!this->IsImage())
    return flag;

  *xlen=mXlen;
  *ylen=mYlen;

  flag=TRUE;
  return flag;
}

void CLImage::FreeImage()
{
  if (mBuff)
    {
      delete [] mBuff;
      mBuff=0;
    }
}

bool CLImage::ImageBuff16(ushort *buf)
{
  bool flag=FALSE;

  if(!mBuff || !buf)
    return flag;

  long len=mXlen*mYlen;

  short typ=mType;

  if(typ==BPP16)
    memcpy(buf,mBuff,len*2);    //Reading 16-bit data into 16bit buff
  else
    {
      for(int i=0;i<len;i++)
	buf[i]=(ushort)mBuff[i];  // Reading 8-bit data into 16bit buff
    }

  flag=TRUE;
  return flag;
}

bool CLImage::ImageBuff(uchar *buf)
{
  bool flag=FALSE;

  if(!mBuff || !buf)
    return flag;

  long len=mXlen*mYlen;
  memcpy(buf,mBuff,len*mType/2);

  flag=TRUE;
  return flag;
}

//Image Processing Functions (Public)
bool CLImage::GetBox(ushort *bbuf,int Xbl,int Ybl,int Cr,int Cc)
{
  bool flag=FALSE;

  if(Cc-Xbl<0 || Cr-Ybl<0)
    return flag;

  ushort *buf=new ushort[this->mXlen*this->mYlen];

  if(!this->ImageBuff16(buf))
    return flag;

  int count=0;
  int jpx=(Xbl-1)/2;
  int jpy=(Ybl-1)/2;
  for(int i=(Cr-jpy);i<(Cr-jpy+Ybl);i++)
    {
      for(int j=(Cc-jpx);j<(Cc-jpx+Xbl);j++)
      { 
	*(bbuf+count)=*(buf+i*this->mXlen+j);
	count++;
      }
    }

  delete [] buf;

  flag=TRUE;
  return flag;
}

bool CLImage::Local_Threshold(CLImage &out,int Xbl,int Ybl,double percent)
{
  bool flag=FALSE;
  if(!this->IsImage())
    return flag;

  int len=Xbl*Ybl;
  ushort *bbuf=new ushort[len];

  int xlen=this->mXlen;
  int ylen=this->mYlen;

  if(!(Xbl%2) || !(Ybl%2))
    {
      cout <<"Block Size must be an odd number"<<endl;
      return flag;
    }

  int jpx=(Xbl-1)/2;
  int jpy=(Ybl-1)/2;
  for(int i=jpx;i<(xlen-jpx);i++)
    {
      for(int j=jpy;j<(ylen-jpy);j++)
	{
	  if(GetBox(bbuf,Xbl,Ybl,i,j))
	    {
	      int sum=0;
	      int mx=0;
	      double avg;
	      int Vcount=0;
	      for(int k=0;k<len;k++)
		{
		  if(bbuf[k]>0)
		    {
		      Vcount++;
		      sum=sum+(int)bbuf[k];
		      if(bbuf[k]>mx)
			mx=bbuf[k];
		    }
		}
	      avg=(double)sum/Vcount;

	      if(Vcount<(int)(len/2))
		continue;
	      else
		{
		  int val=*(this->mBuff+i*xlen+j);

		  if(val>(int)((mx-avg)*percent))
		    *(out.mBuff+i*xlen+j)=1;
		}
	    }
	}
    }

  flag=TRUE;
  return flag;
}

bool CLImage::BinErode(CLImage &out,int order)
{

        bool flag=FALSE;

        if(this->GetType()!=BPP8 || out.GetType()!=BPP8)
           return flag;

        uchar *ap, *rp;
	int xl, yl, i, j, ylm, xlm, xlp;

	this->GetImgSize(&xl,&yl);

	ap = new uchar [xl*yl];

 
	if(!this->ImageBuff(ap))
	  return flag;

	rp = new uchar [xl *yl];

	if(!out.ImageBuff(rp))
	  return flag;

	uchar *rpp=rp;

	ylm = yl - 1;
	xlm = xl - 1;
	xlp = xl + 1;

	for (i = 0; i < xl; i++)
		rp[i] = 0;

	for (j = 1, ap += xl, rp += xl; j < ylm; j++, rp++, ap++)
	{
		rp[0] = 0;
		rp[xlm] = 0;

		for(i = 1, rp++, ap++; i<xlm; i++, rp++, ap++)
		{
			*rp = ((*ap!=0) && ((ap[-1] + ap[1] + ap[xlm]
				+ ap[xlp] + ap[xl] + ap[-xlm] + ap[-xlp]
				+ ap[-xl]) > order)) ? 1 : 0;
		}
	}

	for (i = 0; i < xl; i++)
		rp[i]=0;

	if(!out.SetImage(rpp,xl,yl,2))
	  cout<<"BInErode():SetImage failed!"<<endl;
        //delete [] ap;
	//delete [] rp;
        flag=TRUE;
	return flag;
}

bool CLImage::BinDilate(CLImage &out,int order)
{
        bool flag=FALSE;

        if(this->GetType()!=BPP8 || out.GetType()!=BPP8)
           return flag;

	uchar *ap, *rp;
	int xl, yl, i, j, ylm, xlm, xlp;

	this->GetImgSize(&xl,&yl);

	ap = new uchar [xl*yl];
	if(!this->ImageBuff(ap))
	  return flag;

	rp = new uchar [xl*yl];
	if(!out.ImageBuff(rp))
	  return flag;

	uchar *rpp=rp;

	ylm = yl - 1;
	xlm = xl - 1;
	xlp = xl + 1;

	for (i = 0; i < xl; i++)
		rp[i] = ap[i];

	for (j = 1, ap += xl, rp += xl; j < ylm; j++, rp++, ap++)
	{
		rp[0] = ap[0];
		rp[xlm] = ap[xlm];

		for (i = 1, rp++, ap++; i < xlm; i++, rp++, ap++)
		{
			*rp = ((*ap != 0) || ((ap[-1] + ap[1] + ap[xlm]
				+ ap[xlp] + ap[xl] + ap[-xlm] + ap[-xlp]
				+ ap[-xl]) >= order)) ? 1 : 0;
		}
	}

	for (i = 0; i < xl; i++)
		rp[i] = ap[i];

	if(!out.SetImage(rpp,xl,yl,BPP8))
	  cout <<"BinDilate():SetImage Failed!"<<endl;

	delete [] ap;
	delete [] rp;

	flag=TRUE;
	return flag;
}

bool CLImage::BinDilateNJ(CLImage &out,int order)
{
        bool flag=FALSE;

        if(this->GetType()!=BPP8 || out.GetType()!=BPP8)
           return flag;

	uchar *ap, *rp, tot;
	int xl, yl, i, j, ylm, xlm, xlp, bin_adr;
	this->GetImgSize(&xl,&yl);

	ap = new uchar[xl*yl];
	if(!this->ImageBuff(ap))
	  return flag;

	rp = new uchar [xl*yl];
	if(!out.ImageBuff(rp))
	  return flag;

	uchar *rpp=rp;

	ylm = yl - 1;
	xlm = xl - 1;
	xlp = xl + 1;

        for (i = 0; i < xl; i++)
		rp[i] = ap[i];

	for (j = 1, ap += xl, rp += xl; j < ylm; j++, rp++, ap++)
	{
		rp[0] = ap[1];
		rp[xlm] = ap[xlm-1];

		for (i = 1, rp++, ap++; i < xlm; i++, rp++, ap++)
		{
			if (!*ap)
			{
				bin_adr = (((((((((((((ap[xlp] << 1) | ap[xl]) << 1)
					| ap[xlm]) << 1) | ap[1]) << 1) | rp[-1]) << 1)
					| rp[-xlm]) << 1) | rp[-xl]) << 1) | rp[-xlp];

				*rp = ((tot = dialate_table[bin_adr]) != 0) ?
					(((tot + ap[-1] + ap[-xlm] + ap[-xl] + ap[-xlp])
					> order) ? 1 : 0) : 0;
			}
			else 
				*rp = 1;
		}
	}

	for (i = 0; i < xl; i++)
		rp[i] = ap[i];

	if(!out.SetImage(rpp,xl,yl,BPP8))
	  cout<<"BinDilateNJ():SetImage failed!"<<endl;

	delete [] ap;
	delete [] rp;

	flag=TRUE;
	return flag;
}


/********************************************************************
 *								    *
 *  Binary images labeling algorithm:                               *
 *	    puts a unigue label on each object within binary ARRAY a *
 *
 *  - assign numbers in an increasing order to adjacent pixels      *
 *  - choose the smallest number among connected pixels             *
 *  - Link list of all labels connected labels will be in seperate  *
 *      lists  (see routine link for list format)		    *
 *  - count the number of connected regions (label_index)           *
 *  - assign a distinct label to each connected region (1,2,...)    *
 *								    *
 *			Array a is the binary selection array	    *
 *			Array r is the resulting label array  (ushort) *
 *								    *
 *			Returns (1 + maximum label #) 		    *
 *					  (-1) if error	            *
 *								    *
 *By: Dr. Azhar Quddus Nov 20,2002. SWCHSC, ON, Canada.             *
 ********************************************************************/
int CLImage::Label_Connected_Regions(CLImage &out)
{
  if(this->GetType()!=BPP8 || out.GetType()!=BPP16)
    return -1;

        int x,y,xl,yl,xlp,xlm,xlength,ylength;
	this->GetImgSize(&xl,&yl);

	uchar *ap=new uchar[xl*yl];
	if(!this->ImageBuff(ap))
	  return -1;

	ushort *rp=new ushort[xl*yl];
	if(!out.ImageBuff16(rp))
	  return -1;

	struct list *labels;
        unsigned short *over,*overlay, current_label,label_index, next, n, f;
	
	xlength = xl+2;
	ylength = yl+1;
	xlm = xl+1;
	xlp = xl+3;

	ushort *rpp=rp;

/*
	Allocate a array of zeros for temp. label storage.
	Note the array is extended	by 2 columns and 1 row
		so that special edge conditions can be ignored
*/
   if((over=(unsigned short *)calloc(xlength*ylength,sizeof(unsigned short)))==NULL){
			return -1;
			}
  
/*
   Assign initial labels
	The label assigned is the the lowest non zero label of the neigbourhood.
	If there are no labels then assign the next available label.
	Note that we only have to check 4 out of 8 nearest neigbours as the rest
		have not yet been assigned a label.
*/ 
   for (y=0, next = 1, overlay=&over[xlm]; y<yl; y++)
    for (x=0,overlay+=2; x<xl; x++, overlay++, ap++)
		if(ap[0]){				  // if source image is none zero then label
			*overlay = next;
			if((overlay[-xlength])&&(*overlay > overlay[-xlength]))
				*overlay = overlay[-xlength];
			if((overlay[-xlp])&&(*overlay  > overlay[-xlp]))
				*overlay = overlay[-xlp];
			if((overlay[-xlm])&&(*overlay > overlay[-xlm]))
				*overlay = overlay[-xlm];
			if((overlay[-1])&&(*overlay > overlay[-1]))
				*overlay = overlay[-1];
			if(next == *overlay)++next;
			}

//		Allocate zeroed storage for linked lists
   if((labels=(struct list *)calloc(next,sizeof(struct list)))==NULL){
			return -1;
			}
//	 For all element assign root and tail equal to element address
	for (x=1; x<next; x++)labels[x].tail=labels[x].root=x;
/*
	Do one pass to create link lists	by looking at nearest neigbours
	for lists with a different root.

	Note that we only have to check 4 out of 8 nearest neigbours as the rest
		will be checked later in the same pass.

	This step could be speeded up by pre allocateing link list array and
		combining this pass with the first pass.
*/
   for (y=0, overlay = &over[xlm]; y<yl; y++)
    for (x=0, overlay+=2; x<xl; x++, overlay++)
		if(f=overlay[0]){
			if((n=overlay[-xlength])&&(f != n)
				&&	((f=labels[f].root)!=(n=labels[n].root)))  link(labels,f,n);
			if((n=overlay[-xlm])&&(f != n)
				&& ((f=labels[f].root)!=(n=labels[n].root))) link(labels,f,n);
			if((n=overlay[-xlp])&&(f != n)
				&& ((f=labels[f].root)!=(n=labels[n].root))) link(labels,f,n);
			if((n=overlay[-1])&&(f != n)
				&& ((f=labels[f].root)!=(n=labels[n].root))) link(labels,f,n);
			}


//	Relabel the root to remove gaps 	place the new label in .next
	for (x=1,label_index=1, current_label=0; x<next; x++)
		if(labels[x].root>current_label){
			current_label = labels[x].root;
			labels[x].next = label_index++;
			}

//	Relabel the root of all labels to the new labels
//   labels[].root becomes a direct addressable translation table
	for (x=1; x<next; x++)labels[x].root = labels[labels[x].root].next;

//	Final pass to relabel and place the new labels in the result array
   for (y=0, overlay=&over[xlm]; y<yl; y++)
    for (x=0,overlay+=2; x<xl; x++, overlay++, rp++)
		*rp = labels[*overlay].root;
   
   out.SetImage((uchar *)rpp,xl,yl,BPP16);

   free(over);
   free(labels);

  // delete [] ap;
   //delete [] rp;

   return label_index;
}

bool CLImage::DispImg()
{
  bool flag=FALSE;

  ofstream ofile ("_TMP_.tmp",ofstream::out | ofstream::binary);

  if(!ofile.is_open())
    {
      cout<<"Cant open tmp file!"<<endl;
      return flag;
    }

  if(!this->IsImage())
    {
      cout<<"No buffer!"<<endl;
      return flag;
    }

  //Save Data
  //ofile << (ushort)this->mXlen <<(ushort)this->mYlen <<(ushort)this->mType;
  //long pos=ofile.tellp();
  //ofile.seekp(pos);
  ofile.write((char *)this->mBuff,(this->mXlen * this->mYlen * this->mType/2));

  ofile.close();
  system("debugImg");

  flag=TRUE;
  return flag;
}
