// This Class was written by : Dr. Azhar Quddus
// Nov 20, 2002.
// All rights reserved. Dept. of Cognitive Neurology, Sunnybrook & Womens'
// Toronto, Canada.

#ifndef _IMAGE_H_
#define _IMAGE_H_

#define BPP8 2
#define BPP16 4

#define TRUE 1
#define FALSE 0

typedef unsigned char uchar;
typedef unsigned short ushort;

class CLImage
{
 private:
  uchar *mBuff;
  int mXlen;
  int mYlen;
  short mType;   //2 or 4 for 8-bit/16-bit images, respectively

 public:
  //General Functions
  CLImage();        //This will not Allocate any memory... just sets zeros
  CLImage(int xl, int yl, short typ);   //Will allocate memory
  CLImage(const CLImage &);     //Copy constructor... will allocate memory.
 
  bool SetImage(uchar *buf,int xl, int yl,short typ);
  ~CLImage();

  bool IsImage();   //This should be called after constructors
  bool IsBlank();   //To check if it is blank (all zeros)
  void FreeImage();
  short GetType() {return mType;}
  bool GetImgSize(int *xlen, int *ylen);
  //  bool IsBinary(CLImage img);   //time comsuming... will do later
  bool ImageBuff16(ushort *buf);   //will not allocate memory
  bool ImageBuff(uchar *buf);      //will not allocate memory

  //Image processing functions
  bool GetBox(ushort *bbuf,int Xbl,int Ybl,int Cx,int Cy); //will give Xbl*Ybl(odd) ushorts in Bbuff, Cx,Cy is the center of the box
  bool BinErode(CLImage &out,int order);
  bool BinDilate(CLImage &out,int order);
  bool BinDilateNJ(CLImage &out, int order);
  int Label_Connected_Regions(CLImage &out);
  bool Local_Threshold(CLImage &out,int Xbl,int Ybl,double percent); //Thrsholding is done on the percentage of (max-mean in the box)

  //For debugging only
  bool DispImg();
};


#endif

