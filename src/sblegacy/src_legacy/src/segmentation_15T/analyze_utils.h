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

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include <stdio.h>
#include <math.h>
#include "../../../include_from_cortex/dbh.h"
#include "../../../include_from_cortex/conv.h"


#ifndef ANALYZE_UTILITES_H
#define ANALYZE_UTILITES_H

namespace nk{
#define NAMELENGTH 32
#define HEADERLENGTH 120

#define absval(x) ((x) > 0 ? (x) : (-(x)) )
#define maxval(x, y) ( (x) > (y) ? (x) : (y) )
#define minval(x, y) ( (x) < (y) ? (x) : (y) )
#define round(x)  ((int)(floor((x) + 0.5)))

using namespace std;

typedef vector<string> ROINameVector;

template <class T>
T *load_analyze_object_map(char *fname, ROINameVector &roi_names, int *sz,
                           short &numobj)
{
    int dim[5];//dimension parameters: version #,width,height,depth,
                          //# of objects
    int i;
    char *name;//name of objects
    char *header;//other header info for each object

    ifstream obj_file(fname);
    if(!obj_file){
        cout << "File " << fname << " not found. Exiting ..." << endl;
        exit (1);
    }
    obj_file.read(reinterpret_cast<char *>(dim), 5*sizeof(int));

    // check if byte-swap needed
    if (dim[4] > 10000)  // if # of objects is absurdly large than swap
    	for (int temp = 0; temp < 5; temp++)
    		dim[temp] = iconv(dim[temp]);

    numobj = dim[4];

    //read in object names and headers
    name = new char[NAMELENGTH];
    header = new char [HEADERLENGTH];

    for (i=0;i<dim[4];i++) {
        obj_file.read(name, NAMELENGTH*sizeof(char));
        string temp_name(name);
        roi_names.push_back(temp_name);
        obj_file.read(header, HEADERLENGTH*sizeof(char));
    }
    delete [] name;
    delete [] header;

    sz[0] = dim[1];
    sz[1] = dim[2];
    sz[2] = dim[3];
    int volume = sz[0] * sz[1] * sz[2];

    //get object map
    T *roi_mask = new T [volume];
    if (!roi_mask) {
        cout<<"Failed to allocate enough memory!\n";
        exit(0);
    }
    fill(roi_mask, roi_mask + volume, (unsigned char)0);

    unsigned char num;
    unsigned char val;
    int prev,total = 0;

    T *ptr = roi_mask;
    while (total < volume) {
        prev = total;
        obj_file.read(reinterpret_cast<char *>(&num),sizeof(char));
        obj_file.read(reinterpret_cast<char *>(&val), sizeof(char));
        while (total < prev + (int)num) {
            *ptr++= (T)val;
            total++;
        }
    }
    //return object map
    return(roi_mask);
}
/*
//////////////////////////////////////////////////////////////////////////////
//loads an analyze object map
template <class T>
T *load_analyze_object_map(char *fname, ROINameVector &roi_names, int *sz,
                           short &numobj)
{
    int dim[5];//dimension parameters: version #,width,height,depth,
                          //# of objects
    int i;
    char *name;//name of objects
    char *header;//other header info for each object
    
                                                                
    ifstream obj_file(fname);
    if(!obj_file){
        cout << "File " << fname << " not found. Exiting ..." << endl;
        exit (1);
    }
    obj_file.read(dim,5*sizeof(int));
    numobj = dim[4];

    //read in object names and headers
          name = new char[NAMELENGTH];    
          header = new char [HEADERLENGTH];
    for (i=0;i<dim[4];i++) {
        obj_file.read(name,NAMELENGTH*sizeof(char));
        roi_names.push_back(name);
        obj_file.read(header,HEADERLENGTH*sizeof(char));
    }
    delete [] name;
    delete [] header;
    
    sz[0] = dim[1];
    sz[1] = dim[2];
    sz[2] = dim[3];
    int volume = sz[0] * sz[1] * sz[2];
    
    //get object map
    T *roi_mask = new T [volume];
    if (!roi_mask) {
        cout<<"Failed to allocate enough memory!\n";
        exit(0);
    }
    fill(roi_mask, roi_mask + volume, (unsigned char)0);
    
    unsigned char num;
    unsigned char val;
    int prev,total = 0;

    T *ptr = roi_mask;
    while (total < volume) {
        prev = total;
        obj_file.read(&num,sizeof(char));
        obj_file.read(&val,sizeof(char));
        while (total < prev + (int)num) {
            *ptr++= (T)val;
            total++;
        }
    }
    //return object map
    return(roi_mask);
}
*/


////////////////////////////////////////////////////////////////////////

struct info{
    short sz[3];       // image dimensions
    float voxsz[3];    // voxel size in mm
    short datatype;    // alowable values: 0, 1, 2, 4, 16, 32, 64, 128, 255
                       // although here we consider only 2 and 4, corresponding
                       // to unsigned char and short
    int max, min;      // absolute max and min of image data
};
///////////////////////////////////////////////////////////////////////////
bool read_analyze_header(const string &hdr_fname, info &hdr_info)
{
    ifstream in(hdr_fname.c_str());
    if(!in) {
        cout << "Couldn't open " << hdr_fname << " for input."
             << " Exiting..." << endl;
        exit (1);
    }
    bool swap = false;
    short datatype;
    in.seekg(70);
    in.read(reinterpret_cast<char *>(&datatype), sizeof(short));
    if(datatype != 2 && datatype != 4 && datatype !=16) {
      // assuming that 2,4, and 16 are the only choices
      // must swap bytes
      swap = true;
    }

    in.seekg(42);
    in.read(reinterpret_cast<char *>(hdr_info.sz), 3*sizeof(short));
    in.seekg(80);
    in.read(reinterpret_cast<char *>(hdr_info.voxsz), 3*sizeof(float));
    in.seekg(70);
    in.read(reinterpret_cast<char *>(&(hdr_info.datatype)), sizeof(short));
    in.seekg(140);
    in.read(reinterpret_cast<char *>(&(hdr_info.max)), sizeof(int));
    in.read(reinterpret_cast<char *>(&(hdr_info.min)), sizeof(int));
    
    if(swap) {
      for(int i=0; i < 3; ++i) {
	hdr_info.sz[i] = siconv(hdr_info.sz[i]);
	hdr_info.voxsz[i] = fconv(hdr_info.voxsz[i]);
      }
      hdr_info.datatype = siconv(hdr_info.datatype);
      hdr_info.max = iconv(hdr_info.max);
      hdr_info.min = iconv(hdr_info.min);
    }
    return swap;
}
///////////////////////////////////////////////////////////////////////////
template <class T>
T *read_analyze(const string &froot, info &hdr_info)
{
    info Tinfo;
    bool swap = read_analyze_header(froot + ".hdr", Tinfo);
    if( &hdr_info ) hdr_info = Tinfo;
    
    // make sure that image is either 8 or 16 bit,
    if(!(hdr_info.datatype == 2 || hdr_info.datatype == 4 || hdr_info.datatype == 16)) {
      cout << "Image datatype in " << froot << " is neither 8 nor 16 bits. For now we do not support anything else. So, exiting..." << endl;
      exit (1);
    }

    long volume = Tinfo.sz[0] * Tinfo.sz[1] * Tinfo.sz[2];

    T *img = new T[volume];
    if( !img ) {
        cout << "Memory allocation error in read_analyze("
             << froot << "). Exiting...." << endl;
        exit (1);
    }

    ifstream in((froot + ".img").c_str());
    if( !in ) {
        cout << "Can't open " << froot + ".img" << " for input. Exiting..."
            << endl;
        exit (1);
    }

    // if template type agrees with datatype in the image

    if((sizeof(T) == Tinfo.datatype/2) || (sizeof(T) == Tinfo.datatype/4)) {
        in.read(reinterpret_cast<char *>(img), volume * sizeof(T));
	if(Tinfo.datatype == 4 && swap) {
	  T *ptr = img;
	  for (int i=0; i < volume; i++, ptr++) *ptr = siconv(*ptr);
	}
	else if(Tinfo.datatype == 16 && swap) {
	  T *ptr = img;
	  for (int i=0; i < volume; i++, ptr++) *ptr = (T)fconv(*ptr);
	}
        return(img);
    }

    else {
      cout << "read_analyze error: wrong template. Exiting... "<< endl;
      return 0;
    }
}
//////////////////////////////////////////////////////////////////////////////
void make_hdr(const string &hdr_fname, const info &hdr_info)
    {
    int i;
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


    hdr.dime.datatype = hdr_info.datatype;
    hdr.dime.bitpix = hdr_info.datatype * 4;


    if((fp=fopen(hdr_fname.c_str(),"w"))==0)
        {
        printf("unable to create: %s\n", hdr_fname.c_str());
        exit(0);
    }

    hdr.dime.dim[0] = 4;  /* all Analyze images are taken as 4 dimensional */
    hdr.hk.regular = 'r';
    hdr.hk.sizeof_hdr = sizeof(struct dsr);

    hdr.dime.dim[1] = hdr_info.sz[0];  /* slice width  in pixels */
    hdr.dime.dim[2] = hdr_info.sz[1];  /* slice height in pixels */
    hdr.dime.dim[3] = hdr_info.sz[2];  /* volume depth in slices */
    hdr.dime.dim[4] = 1;  /* number of volumes per file */

    hdr.dime.glmax  = hdr_info.max;  /* maximum voxel value  */
    hdr.dime.glmin  = hdr_info.min;  /* minimum voxel value */
{
    hdr.dime.pixdim[1] = hdr_info.voxsz[0]; /* voxel width in mm per pixel */

        hdr.dime.pixdim[2] = hdr_info.voxsz[1]; /* voxel height in mm per pixel */
        hdr.dime.pixdim[3] = hdr_info.voxsz[2]; /* voxel depth in mm per pixel
        */
}
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
}
////////////////////////////////////////////////////////////////////////////
void write_analyze(const string &froot, void *img, info &hdr_info)
{
    make_hdr(froot + ".hdr", hdr_info);
    
    ofstream out((froot + ".img").c_str());
    if( !out ) {
        cout << "Can't open " << froot + ".img" << " for output. Exiting..."
             << endl;
        exit (1);
    }
 
    out.write(reinterpret_cast<const char *>(img), hdr_info.sz[0] * hdr_info.sz[1] * hdr_info.sz[2] * hdr_info.datatype/2);
}

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//       Strips off the substring before or after the first occurrence of  //
//       a character and returns the new string.                           //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
char *stripstr(char *roiname,char letter,bool before)
{
    char *outname,*name,*tmpname;

    tmpname = roiname;
    outname = new char[strlen(roiname)];
    name =outname;
    if (before) { 
        while ((*roiname) && (*roiname != letter))
            roiname++;
        if (*roiname++)
            while (*outname++=*roiname++)
                ;
        else
            while (*outname++=*tmpname++)
                ;
    }
    else { //strips off the characters after and including letter
        while ((*roiname != letter) && (*roiname)) 
            *outname++ = *roiname++;
        *outname = '\0';
    }
    return(name);
}

/////////////////////////////////////////////////////////////////
//  CommandLine class (written by Gal Sela)
/////////////////////////////////////////////////////////////////
class CommandLine
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
    
    //operator int() { return(_current < _argc); }
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

///////////////////////////////////////////////////////////////////
// added for viz purposes
template <class T>
unsigned char *read_analyze_slice(const std::string &froot, int ortodir, int slice, int *dim, int scalemax=0)
{
  unsigned char *display_img_slice=NULL;
  info Tinfo;
  bool swap = read_analyze_header(froot + ".hdr", Tinfo);
  // make sure that image is either 8 or 16 bit,
  if(!(Tinfo.datatype == 2 || Tinfo.datatype == 4)) {
    cout << "Image datatype in " << froot << " is neither 8 nor 16 bits. For now we do not support anything else. So, exiting..." << endl;
    exit (1);
  }

  if(Tinfo.sz[2]==1) {
    cout << "read_anlyze_slice: this is 2D data, ignoring ortodir." << endl;
  }
  else {
    if(ortodir < 0 || ortodir > 2) {
      cout<< "read_anlyze_slice: ortodir has impossible value " << ortodir << ". defaulting to ortodir=2 (z-axis)." << endl;
      ortodir=2;
    }
  }
  int axis[3], sz[3];
  for(int i=0; i < 3; i++) {
    sz[i]=Tinfo.sz[i];
    axis[i]=(ortodir+i+1)%3;
    dim[i]=Tinfo.sz[axis[i]];
  }
  if(slice<0 || slice>=sz[axis[2]]) {
     cout<< "read_anlyze_slice: impossible slice number" << slice << ". defaulting to slice 0." << endl;
     slice=0;
  }

  long area=dim[0]*dim[1], orig_area = sz[0]*sz[1];
  T *img_slice = new T [area];
  if( !img_slice ) {
    cout << "Memory allocation error in read_analyze_slice("
         << froot << "). Exiting...." << endl;
    exit (1);
  }
  ifstream in((froot + ".img").c_str());
  if( !in ) {
    cout << "Can't open " << froot + ".img" << " for input. Exiting..."
         << endl;
    exit (1);
  }

  // if template type agrees with datatype in the image
  if(sizeof(T) == Tinfo.datatype/2) {
    int i, j, *a, *b, *c;
    switch(ortodir) {
    case 0: a=&slice; b=&i; c=&j; break;
    case 1: a=&j; b=&slice; c=&i; break;
    case 2: a=&i; b=&j; c=&slice; break;
    }
    T *ptr = img_slice;
    for(j=0; j < dim[1]; j++)
      for(i=0; i < dim[0]; i++) {
        in.seekg(sizeof(T)*(*c*orig_area + *b*sz[0] + *a), ios::beg);
        in.read(reinterpret_cast<char *>(ptr), sizeof(T));
        ptr++;
      }
    if(Tinfo.datatype == 4 && swap) {
      ptr = img_slice;
      for (int i=0; i < dim[0]*dim[1]; i++, ptr++) *ptr = siconv(*ptr);
    }
    
    int min, max;
    if(scalemax==0) {
      // if scalemax not supplied find real min, max and rescale [min,max]->[0,255]
      min=65535; max=0;
      ptr=img_slice;
      for(int t=0; t < area; t++, ptr++){
	if(*ptr < min) min = *ptr;
	if(*ptr > max) max = *ptr;
      }
      cout << "true min max " << min << " " << max << endl;
    }
    else {
      // rescale [0,scalemax]->[0,255]
      min=0; max = scalemax;
    }
    display_img_slice = new unsigned char [area];
    unsigned char *displ_ptr=display_img_slice;
    ptr = img_slice;
    for(int t=0; t < area; t++, ptr++, displ_ptr++){
      *displ_ptr = static_cast< unsigned char>(round(255.0*(*ptr - min)/(max-min)));
    }
    delete [] img_slice;
  }
  else {
    cout << "read_analyze_slice error: wrong template. Exiting... "<< endl;
  }
  return display_img_slice;
}
///////////////////////////////////////////////////////////////////

} // end name space nk
#endif //ANALYZE_UTILITES_H

    

