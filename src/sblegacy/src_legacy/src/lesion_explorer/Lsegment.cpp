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

////////////////////////////////////////////////////////
//
// Copyright by Azhar Quddus, December 2002
//
////////////////////////////////////////////////////////

// To compile:
// g++ -O2 -o Lsegment Lsegment.cpp /net/cortex/csf/pub/include/MyVolume.cpp \
//	/net/cortex/csf/pub/include/MyImage.cpp /net/cortex/csf/pub/include/connexe.c
// use at least gcc-3.2

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <math.h>
#include "../../../include_from_cortex/MyVolume.h"

#define maxval(x, y) ((x)>(y) ? (x) : (y))
#define minval(x, y) ((x)<(y) ? (x) : (y))
#define round(x)  ((int)(floor((x) + 0.5)))

unsigned char BGD=0;

using namespace std;

/////////////////// forward declarations ///////////////////////////////////
template <class T>
void find_borders(const T* Img, const int *sz, int *B)
{
    int area = sz[0] * sz[1], volume = area * sz[2];
    int i, j, k;
   
    // find lowest non-zero z-level
    i=0;
   
    while( Img[i]==(T)0 && i<volume) i++;
    B[4] = i / area; 
    
    // find highest non-zero z-level
    i=volume-1;
    while( Img[i] == (T)0 && i>=0) i--;
    B[5] = i / area; 
    
    // find lowest non-zero x-level
    bool found=false;
    for(i=0; i < sz[0] && !found; i++)
        for(j=0; j < sz[1] && !found; j++)
            for(k=0; k < sz[2] && !found; k++) 
                if( Img[k*area + j*sz[0] + i] ) found=true;
    B[0] = maxval(i-1, 0);

   // find highest non-zero x-level
    found=false;
    for(i=sz[0]-1; i>=0 && !found; i--)
        for(j=0; j < sz[1] && !found; j++)
            for(k=0; k < sz[2] && !found; k++) 
                if( Img[k*area + j*sz[0] + i] ) found=true;
    B[1] = minval(i+1, sz[0]-1);

    // find lowest non-zero y-level
    found=false;
    for(j=0; j < sz[1] && !found; j++)
        for(i=0; i < sz[0] && !found; i++)
                  for(k=0; k < sz[2] && !found; k++)
                       if( Img[k*area + j*sz[0] + i] ) found=true;
    B[2] = maxval(j-1, 0);

    // find highest non-zero y-level
    found=false;
    for(j=sz[1]-1; j >=0 && !found; j--)
        for(i=0; i < sz[0] && !found; i++)
                  for(k=0; k < sz[2] && !found; k++)
                       if( Img[k*area + j*sz[0] + i] ) found=true;
    B[3] = minval(j+1, sz[1]-1);

    cout << "Borders :"
         << "  x_range=[" << B[0] << ":" << B[1] << "]"
         << "  y_range=[" << B[2] << ":" << B[3] << "]"
         << "  z_range=[" << B[4] << ":" << B[5] << "]"
         << endl;
    
    // check that borders in x, y direction are not same as image edges
    if(B[0]==0 || B[1]==sz[0]-1 ||
       B[2]==0 || B[3]==sz[1]-1) {
        cout << "Image borders suspicious. Check your mask." << endl;
        cout << "Algorithm will proceed, results unpredictable. " << endl;
    }
}

/////////////////////////////////////////////////////////////////////////////
template <class T>
int num_of_nonzero_voxels(const T* Img, const int *B, const int *sz)
{
    int count=0;
    int x, y, z, area = sz[0]*sz[1];
    for(z=B[4]; z <= B[5]; ++z)
            for(y=B[2]; y <= B[3]; ++y)
                for(x=B[0]; x <= B[1]; ++x) {
                    if( Img[ z * area + y * sz[0] + x] ) ++count;
                }
    return count;
}

///////////////////////////////////////////////////////////////////////////
template <class T>
inline void segment(const T* Img, unsigned char *Seg,const int *B, const int *sz, double cutoff )
    // Perform segmentation within borders defined by borders B
{  
        int x, y, z, t, area = sz[0]*sz[1];

	//Get Mean & Max value
	double pmean=0.0;
	int pmax=0;
	int psum=0;
	int pcount=0;

	for(z=B[4]; z <= B[5]; ++z)
            for(y=B[2]; y <= B[3]; ++y)
                for(x=B[0]; x <= B[1]; ++x) 
		  {
                    t = z * area + y * sz[0] + x;
		    //cout <<"t= "<<t<<endl;
		    int pix=(int)Img[t];

		    //cout<<x<<"-"<<y<<"-"<<z<<"-- pix= "<< pix <<endl;

		    if(pix)
		      {
			++pcount;
			psum=psum+pix;

			if(pix > pmax) pmax=pix;
		      }
                       
		  }

	if(pcount<1)
	  pmean=0;
	else
	  pmean=psum/pcount;

	int th=round(pmean+(pmax-pmean)*cutoff);

	//Do the Thresholding
	if(th>0)
	  {
	    for(z=B[4]; z <= B[5]; ++z)
            for(y=B[2]; y <= B[3]; ++y)
                for(x=B[0]; x <= B[1]; ++x) 
		  {
                    t = z * area + y * sz[0] + x;
	     
		    if( (int)Img[t] > th )  Seg[t] =(uchar)1;
                       
		  }
	  }
}


//////////////////////PROTOTYPES***/////////////////////////////////
class CommandLine;
void print_usage_message(const char *command);

void parse_command_line(int argc, char **argv, string &in_fname, string &out_fname, double &cutoff, int *step);

////////////////////******MAIN-FUNCTION****/////////////////////
int main(int argc, char **argv)
{
    if(argc < 3 ){
        print_usage_message(argv[0]);
        exit(1);
    }
    
    string in_fname, out_fname;
    double cutoff;

    int step[3], increment[3];

    // increments
    increment[0] = 2; increment[1] = 2; increment[2] = 1;
   
    CLVolume inV, outV;

    //parse command-line
    parse_command_line(argc, argv, in_fname, out_fname, cutoff,step);

    //Read input file
    if(!inV.read_analyze(in_fname))
    {
      cout <<"Couldnt Read Input File! Exiting..."<<endl;
      return -1;
    }
   
    //Get Header info
    hdr myhdr=inV.GetHeader();

    int sz[3];
    sz[0]=myhdr.VolSize[0];
    sz[1]=myhdr.VolSize[1];
    sz[2]=myhdr.VolSize[2];

    //get Volume buffer
    long volume=myhdr.VolSize[0]*myhdr.VolSize[1]*myhdr.VolSize[2];

    uchar *Img=new uchar[volume*myhdr.DataType/2];

    if(!inV.GetVolBuff(Img))
      {
	cout << "Read Buffer Failed ! exiting..." <<endl;
	return -1;
      }
   
    // Creat 8-bit seg_image of the same size as in_img and fill it with BGD
    uchar *Seg = new uchar [volume];
    if(!Seg){
        cout << "Memory allocation error. Exiting ... " << endl;
        exit(1);
    }
    fill(Seg, Seg + volume, BGD);

    //prepare/set header for seg_image
    hdr ohdr=myhdr;
    ohdr.Imax=1;
    ohdr.Imin=0;
    ohdr.DataType=2;
    outV.SetHeader(ohdr);
 
    // find borders of brain
    int GB[6];

    ushort *shortImg;
    if(myhdr.DataType==4)
      {
	shortImg=(ushort *)Img;
	find_borders(shortImg, sz, GB);
      }
    else
      find_borders(Img, sz, GB);

    //Here is the meat...
    cout << "Performing localized segmentations. Local box size ["
              << step[0] << ", " << step[1] << ", " << step[2] << "]"<< endl;
   
    int minVB=step[0]*step[1]*step[2]/7;

    int shift[3];
    for(int t=0; t < 3; t++)  shift[t] = step[t]/3;

    // LB=borders of local box, CLB=borders of the core of local box
    int LB[6], CLB[6];
        
    int t, index[3];
        
    int boxnumber=0;
    for(index[2]=GB[4]; index[2] < GB[5] - shift[2] + 1; index[2] += shift[2])
        for(index[1]=GB[2]; index[1] < GB[3] - shift[1] + 1; index[1] += shift[1])
            for(index[0]=GB[0]; index[0] < GB[1] - shift[0] + 1; index[0] += shift[0]) 
	      {
		boxnumber++;
                    for(t=0; t < 3; ++t) 
		      {
                        // lower local box border
                        CLB[2*t] = LB[2*t] = index[t];
                        // upper local box border
                        CLB[2*t+1] = LB[2*t+1] = minval(LB[2*t] + step[t] - 1, GB[2*t+1]);
		      }
		   
		    //Number of Non-Zero Voxels
		    if(myhdr.DataType==4)
		      {
			while( num_of_nonzero_voxels(shortImg, LB, sz) < minVB)
			  {
			    //extend local box size by increment
			    for(t=0; t < 3; ++t)
			      {
				LB[2*t] = maxval(LB[2*t] - increment[t], GB[2*t]);
				LB[2*t+1] = minval(LB[2*t+1] + increment[t], GB[2*t+1]);
			      }
			  }
		      }
		    else
		      {
			while( num_of_nonzero_voxels(Img, LB, sz) < minVB )
			  {
			    //extend local box size by increment
			    for(t=0; t < 3; ++t)
			      {
				LB[2*t] = maxval(LB[2*t] - increment[t], GB[2*t]);
				LB[2*t+1] = minval(LB[2*t+1] + increment[t], GB[2*t+1]);
			      }
			  }
		      }

                    // Reduce box to its middle core for all but edge boxes
                    for(t=0; t < 3; ++t)
		      {
                        if(CLB[2*t]==GB[2*t]) // inital box was on the lower edge
                            CLB[2*t+1] = GB[2*t] + 2*shift[t];
                        else if(CLB[2*t+1]==GB[2*t+1])// inital box was on the upper edge
                            CLB[2*t] += shift[t];
                        else      // inital box was away from the edges
			  {
                            CLB[2*t] += shift[t]; 
                            CLB[2*t+1] = CLB[2*t] + shift[t];
			  }
                        
		      }

                    // Now do segmentation on the core box
		    if (myhdr.DataType==4)
		      segment(shortImg, Seg, CLB, sz, cutoff);
		    else
		      segment(Img, Seg, CLB, sz, cutoff);

		    // cout<<"reached here"<<endl;
	      }

    //Save output File
    if(!outV.SetHeader(ohdr))
      {
	cout <<"Header substitution failed!... exiting.."<<endl;
	return -1;
      }

    if(!outV.SetVolBuff(Seg,volume))
      {
	cout <<"Buffer subtitution failed! ... exiting.."<<endl;
	return -1;
      }

 if(!outV.write_analyze(out_fname))
   {
     cout <<"Error writing output file!... exiting.."<<endl;
     return -1;
   }

    delete [] Img;
    delete [] Seg;

    cout <<"Done."<<endl;

    return 0;
}  // end of main

///////////////////////////////////////////////////////////////////////////
class CommandLine // CommandLine class written by Gal Sela
{
  public:
    CommandLine(int argc, char **argv)
      : _argc(argc), _argv(argv), _current(0) { }

    CommandLine &find(const char *tag, int offset = 0) {
      _current = _argc;
      for(int j=0; j<_argc; j++){
        if(!strncmp(_argv[j], tag, strlen(tag))){
          if(j+offset < _argc) _current = j+offset;
          break;
        }
      }
      return(*this);
    }
    
    operator int() { return(_current < _argc); }
    operator const char *() { return(_argv[_current]); }

    CommandLine &operator++() {
      ++_current;
      if(_current > _argc) _current = _argc;
      return(*this);
    }
  
    CommandLine operator++(int) {
      CommandLine tmp = *this;
      ++*this;
      return(tmp);
    }
 
    char * operator*() {
      return(_argv[_current]);
    }

    char * operator[](int offset) {
      return(_argv[_current+offset]);
    }    

    int current() const { return(_current); }

  private:
    int _argc;
    char **_argv;
    int _current;
};
////////////////////////////////////////////////////////////////////////////
void print_usage_message(const char *command) {

    cout << "Usage: " << command << "\n" << endl;
    cout << "<PD file>  PD image(all non-brain matter voxels set to 0, 16bit)\n" ;
    cout << "<T2 File>  T2 image(all non-brain matter voxels set to 0, 16bit)\n" ;
    cout << "\n<output file>  segmented Output (8bit, [0 1])\n" ;
    cout << "\n[-x <int>] x-size of local box (default 96)\n" ;
    cout << "[-y <int>] y-size of local box (default 96)\n" ;
    cout << "[-z <int>] z-size of local box (default 72)\n" ;
    cout << "\n[-c <double>] Cut-Off for hyperintensity (default 0.25)\n" ;
    cout << "\nExample1: segment inFile outFile 256 256 200 -x 42 -y 42 -z 42\n" ;
    cout << "          Perform segmentation in local boxes of size 42x42x42.\n" ;
    cout << "          when local box has less than 10000 brain voxels, its size is\n" <<endl;

}
///////////////////////////////////////////////////////////////////////
void parse_command_line(int argc, char **argv, string &in_fname, string &out_fname,double &cutoff,int *step)
{
    CommandLine cl(argc, argv);

    in_fname = string(argv[1]); out_fname = string(argv[2]);
    cout << "Input image: " << in_fname << endl;
    cout << "Segmented image: " << out_fname << endl;
    
    cutoff = 0.25;
    if(cl.find("-c")) cutoff = atof(cl[1]);

    step[0] = 96; step[1] = 96; step[2] = 72;
    if(cl.find("-x")) step[0] = atoi(cl[1]); step[0] = (step[0]/3)*3;
    if(cl.find("-y")) step[1] = atoi(cl[1]); step[1] = (step[1]/3)*3;
    if(cl.find("-z")) step[2] = atoi(cl[1]); step[2] = (step[2]/3)*3;  
}

