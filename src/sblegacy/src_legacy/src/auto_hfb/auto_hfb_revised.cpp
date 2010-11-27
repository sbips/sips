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

// compile with:
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o auto_hfb auto_hfb.cpp 

#include "../../include/analyze_utils.h"
#include <stack>
#include <numeric>

using namespace nk;

//////////////////////// GLOBAL VARIABLES //////////////////////////////////
unsigned short *T2_img, *P_img;
unsigned char *seg_img;
int sz[3], area, volume;
double MAGIC_NUMBER = 0.4;  
int cx, cy;                           // histogram dimensions (scan cutoffs)
vector <int> excludes;

////////////////////////////// FUNCTIONS ///////////////////////////////////
void print_usage_msg();
int get_cutoff(const unsigned short *Img, int cutoff_count, int start);
                // return intensity range of Img
void segment(int start_index, const int *dim, double magic_number);
                // function which does initial segmentation into g, w, csf
void Floodfill_2D();
                // turn brain voxels to pink - all that are
                // 2Dconnected to centers 
void Find_unclassified_voxels(); // 2D fuv
////////////////////////////////  MAIN //////////////////////////////////////

int main(int argc, char ** argv)
{
    if(argc < 2 ){
        print_usage_msg();
        exit(1);
    }
    
    info hdr_info;
    string froot(argv[1]);
    T2_img = read_analyze<unsigned short>(froot + "_T2",
                                                          hdr_info);
    P_img = read_analyze<unsigned short>(froot + "_P",
                                                          hdr_info);
    cout << "Input images: " << froot + "_T2  " << froot + "_P" << endl;
    
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
    cout << "size " << sz[0] << ", " << sz[1] << "," << sz[2] << endl;
    area = sz[0]*sz[1];
    volume = area * sz[2];
    
    seg_img = new unsigned char [volume];
    assert(seg_img);
    fill(seg_img, seg_img + volume, 0);

    CommandLine cl(argc, argv);

    int cutoff_count = 50;
    int start = 500;
    if(cl.find("-start")) start = atoi(cl[1]);
    
    cx = get_cutoff(T2_img, cutoff_count, start);
    cy = get_cutoff(P_img, cutoff_count, start);
    cout << "cx, cy ; " << cx << ", " << cy << endl;
    
    if(cl.find("-c")) MAGIC_NUMBER = atof(cl[1]);
   
    // do segmentation slice by slice
    int dim[3];
    dim[0] = sz[0];    dim[1] = sz[1];
    dim[2] = 1;
    double delta = 0.2/sz[2];
    for( int k = 0; k < sz[2]; k++) {
        segment(area * k, dim, (0.8 + k*delta)*MAGIC_NUMBER) ;
    } 
    
    double fix = 1.0;
    if(cl.find("-r")) {    // for required slices change magic number and
                           // redo segmentation
         
        fix = atof(cl[1]);
        cout << "Redoing segmentation with fix = " << fix
             << " for slices " << flush;

        int k;
        int i = 2;
        dim[2] = 1;
        while( cl[i] != NULL ) {
            k = atoi(cl[i]);
            cout << k << " " << flush;
            segment(area * k, dim, (0.8 + k*delta) * fix * MAGIC_NUMBER);
            ++i;
        }
    }

    Floodfill_2D();
    
    cout << "\nFinding 2D unclassified voxels";
    if(cl.find("-x")) {
        cout << ", excluding slices ";
        int i = 1;
        while( cl[i] != NULL ) {
            excludes.push_back(atoi(cl[i]));
            cout << atoi(cl[i]) << " ";
            ++i;
        }
    }
    cout << endl;
    
    Find_unclassified_voxels();

    hdr_info.min = 0;
    hdr_info.max = 255;
    hdr_info.datatype = 2;
    if(cl.find("-suf")) froot = froot + cl[1];
    else froot = froot + "_HfB";
    cout << "Writing " << froot  << endl;
    write_analyze(froot , seg_img,  hdr_info);
    
    delete [] T2_img;
    delete [] P_img;
    delete [] seg_img;
    return 0;
}  // end of main

////////////////////////////////////////////////////////////////////////////

int get_cutoff(const unsigned short *Img, int cutoff_count, int start)
{
    vector<int> Full_histo(65536, 0); // initalize histogram with 0's
    
     // calculate full histogram
    for(int t=0; t < volume; t++) {
        if(Img[t] >= 65536) ++(Full_histo[65535]);
        else ++Full_histo[Img[t]];
    }

    // smooth it out
    for(int i = start; i < 2400; i++) {
        Full_histo[i] = round(accumulate(Full_histo.begin()+i-3,
                                         Full_histo.begin()+ i+ 4, 0.0) / 7.0) ;
    }
  
    // find right cutoff point
    int i = start;
    while( Full_histo[i] > cutoff_count ) ++i;
    
    return i;
} // end of get_cutoff

////////////////////////////////////////////////////////////////////////////

void segment(int start_index, const int *dim, double magic_number) {
    
    int t, i, j, k;
    double x, y, z;
    for(k=0; k < dim[2]; k++) {
        for(j=0; j < dim[1]; j++) {
            for(i=0; i < dim[0]; i++) {
                t = start_index + k*area + j*sz[0] + i;
                x =( (double)T2_img[t] )/ (double)cx;
                y = ( (double)P_img[t] )/ (double)cy;
                z = x*x + y*y;
                if( z > magic_number ) {
                    if ( z > 2.7 * magic_number ) seg_img[t] = (unsigned char)5;
                    else seg_img[t] = (unsigned char)3;  
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////

bool acceptable_fuv_2D(const unsigned char *pos){
    if (*pos == 0 ) return true;
    else return false;
}

////////////////////////////////////////////////////////////////////////////

inline int forward_dir_2D_fuv(int index) { // it is assumed that pos is acceptable
    
    int x, y, z;
    z = index / area;
    y = (index - z * area) / sz[0];
    x = index - z * area - y*sz[0];
    
    unsigned char *pos = seg_img + index;
    
    if(x < sz[0]-1) if(acceptable_fuv_2D(pos+1)) return 1;
    if(x > 0) if(acceptable_fuv_2D(pos-1)) return (-1);
    if(y < sz[1]-1) if(acceptable_fuv_2D(pos+sz[0])) return sz[0];
    if(y > 0) if(acceptable_fuv_2D(pos-sz[0])) return (-sz[0]);

    return 0;
}
////////////////////////////////////////////////////////////////////////////
void floodfill_2D_fuv(const int &seed) {
    // do floodfilling in each slice: given slice find all 0 voxels
    // that are 2D connected to given seed in a  slice and turn them
    // to background

    int t, dir;
   
    stack<int> stack_2D;
    stack_2D.push(seed); // set a seed
            
    while( !stack_2D.empty() ) {
        t = stack_2D.top();
        seg_img[t] = 2;
                
        dir = forward_dir_2D_fuv(t);
        if( dir != 0 ) stack_2D.push(t + dir);
        else stack_2D.pop();   // poping from stack means that dead end has been reached
    }
}
   
////////////////////////////////////////////////////////////////////////////

inline bool acceptable_fuv_2D_01(const unsigned char *pos){
    if (*pos == 0 || *pos == 1) return true;
    else return false;
}

////////////////////////////////////////////////////////////////////////////

inline int forward_dir_2D_fuv_01(int index) {
// it is assumed that voxel at position index is acceptable
    
    int x, y, z;
    z = index / area;
    y = (index - z * area) / sz[0];
    x = index - z * area - y*sz[0];
    
    unsigned char *pos = seg_img + index;
    
    if(x < sz[0]-1) if(acceptable_fuv_2D_01(pos+1)) return 1;
    if(x > 0) if(acceptable_fuv_2D_01(pos-1)) return (-1);
    if(y < sz[1]-1) if(acceptable_fuv_2D_01(pos+sz[0])) return sz[0];
    if(y > 0) if(acceptable_fuv_2D_01(pos-sz[0])) return (-sz[0]);

    return 0;
}
////////////////////////////////////////////////////////////////////////////

void floodfill_2D_fuv_01(const int &seed) {
    // do floodfilling in each slice: given slice find all 0 voxels
    // that are 2D connected to given seed in a  slice and turn them
    // to background (for turning large areas of yellow voxels to black)

    int t, dir;
   
    stack<int> stack_2D;
    stack_2D.push(seed); // set a seed
            
    while( !stack_2D.empty() ) {
        t = stack_2D.top();
        seg_img[t] = 2;
                
        dir = forward_dir_2D_fuv_01(t);
        if( dir != 0 ) stack_2D.push(t + dir);
        else stack_2D.pop();   // poping from stack means that dead end has been reached
    }
}
   
////////////////////////////////////////////////////////////////////////////

void Find_unclassified_voxels(){

    int i, t, count, min_count = 600, dir;
    unsigned char *ptr;
    for(int k=0; k < sz[2]; k++) {
        
        // if k not excluded by the user
        if( find(excludes.begin(), excludes.end(), k) == excludes.end() ) {
            
            // go through all 4 corners of each slice
            floodfill_2D_fuv(k*area);
            floodfill_2D_fuv(k*area + sz[0] -1);
            floodfill_2D_fuv(k*area + sz[0] * (sz[1]-1));
            floodfill_2D_fuv(k*area + sz[0]*sz[1] -1);
            
            // now check if there is some large area of unclassifieds left, if yes, plant a seed there
            // and turn it to background as well
            for(i=0; i < area; i++) {
                if( seg_img[k*area + i] == (unsigned char)0 ) {  // look for unclassified

                    // start counting from this spot to see how much unclassified area can we get
                    count = 0;
                    stack<int> stack_2D;
                    stack_2D.push(k*area + i); // set a seed
            
                    while( (!stack_2D.empty()) && count < min_count ) {
                        t = stack_2D.top();
                        seg_img[t] = 1;
                        
                        dir = forward_dir_2D_fuv(t);
                        if( dir )  stack_2D.push(t + dir);
                        else { stack_2D.pop();   count++;}
                    }

                    // if large area encountered, plant a seed and turn to background
                    if( count >= min_count - 1 ) floodfill_2D_fuv_01(k*area + i);
                }
            }

            // turn all remaining unclassified (value 0 ) to value 6
            ptr = seg_img + k * area;
            for(i=0; i < area; i++) {
                if( *ptr < 2 ) *ptr = 7;
                ++ptr;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////

inline bool acceptable_2D(const unsigned char *pos){
    if ( *pos == 5 || *pos == 3 ) return true;
    else return false;
}
///////////////////////////////////////////////////////////////////////////

inline int forward_dir_2D(int index) { // it is assumed that pos is acceptable
    
    int x, y, z;
    z = index / area;
    y = (index - z * area) / sz[0];
    x = index - z * area - y*sz[0];
    
    unsigned char *pos = seg_img + index;
    
    if(x < sz[0]-1) if(acceptable_2D(pos+1)) return 1;
    if(x > 0) if(acceptable_2D(pos-1)) return (-1);
    if(y < sz[1]-1) if(acceptable_2D(pos+sz[0])) return sz[0];
    if(y > 0) if(acceptable_2D(pos-sz[0])) return (-sz[0]);

    return 0;
}
////////////////////////////////////////////////////////////////////////////

void Floodfill_2D() {  // creating pink voxels
    
    int t, dir;
    int middle = sz[0] * round(sz[1]/2.3)  + sz[0]/2;  // lower middle of slice
    for(int k=0; k < sz[2]; k++){
            stack<int> stack_2D;
            stack_2D.push( k * area + middle ); // seed in the middle of k-th slice
            
            while( !stack_2D.empty() ) {
                t = stack_2D.top();
               seg_img[t] = 6;
                
                dir = forward_dir_2D(t);
                if( dir != 0 ) stack_2D.push(t + dir);
                else stack_2D.pop();   // poping from stack means that dead end has been reached
            }
    }
    
}
////////////////////////////////////////////////////////////////////////////
    
void print_usage_msg() {
    cout << "
Usage: auto_hfb

       <T2/PD file root> - it is assumed that they both have the same root

       [-suf] <suffix> to add suffix to the input file, by default  _HfB

       [-fuv3D]  use this switch if you want 3D find-unclassified-voxels procedure
      
       [-x] <s1 s2 s3 ...> exclude slices s1, s2, s3, .. from fuv2D procedure

       [-c] <cutoff> use this to change general cutoff, by default it is set to 0.4 (must have
             a very good reason to change this)

       [-r <percent value> <s1 s2 ...>] use this switch to
            forse redoing segmentation on slices s1, s2, ... by using specified percent
            value of the standard background-brain cutoff

Example 1: auto_hfb bunm92

will read bunm92_P and bunm92_T2 file and output bunm92_HfB file.Output file should be further edited using Pandora so that all that should be considered brain has the highest value within the image.

Example 2: auto_hfb bunm92  -r 0.9 0 1 2 3

will redo segmentation on slices 0, 1, 2, 3 by reducing the standard cutoff to its 90% value
(consequently, there will be more brain and less background on these slices).
This option is to be used typically for bottom of cerebellum if it comes out sprinkled with
background. For consistency purposes, do not use this option on more than 5 slices.

Example2: auto_hfb bunm92 -suf _hfb  -x 0 1 2 3 4

will output bunm_hfb file. Slices 0, 1, 2, 3 and 4 will be excluded in the
2D floodfill operation.
" ;
 
}
