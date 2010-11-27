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
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o erode erode.cpp

// Peel one layer of voxels from the segmented image, and then remove CSF
// i.e. all voxels with value = 5. Use this to mask T1 image.

#include "../../../include_from_cortex/analyze_utils.h"
#include <stack>

//////////////////////// GLOBAL VARIABLES //////////////////////////////////
unsigned char *seg_img;
int sz[3], area, volume;
unsigned char FINAL_VALUE = 253;
////////////////////////////// FUNCTIONS ///////////////////////////////////
void find_borders(const unsigned char *img, int *B);
void peel_one_layer_and_remove_csf(unsigned char *s, int *B);
int find_seed(int *B);
void Floodfill(int seed);
////////////////////////////////  MAIN //////////////////////////////////////

using namespace nk;

int main(int argc, char ** argv)
{
    if(argc < 4 ){
        cout <<"
Usage: erode <T1fname>, name of masked T1 image
             <segfname>, name of the T1 segmented image in the same space
             <outfname>,  name for the output eroded image

Example: erode xxx_sinc_acpc_T1_masked xxx_sinc_acpc_T1_fseg xxx_acpc_erode

 "<< endl;
        
        exit(1);
    }
   
    info hdr_info;
    seg_img = read_analyze<unsigned char>(argv[2], hdr_info);
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
   
    
    area = sz[0]*sz[1];
    volume = area * sz[2];
    info hdr_info1;
    unsigned char *T1_img = read_analyze<unsigned char>(argv[1], hdr_info1);
    
    if( !equal(sz, sz + 3, hdr_info1.sz) ) {
        cout << argv[1]<< " and " << argv[2]
             << " do not have the same size. Exiting ..." << endl;
        exit(1);
    }
    if( !equal(hdr_info.voxsz, hdr_info.voxsz + 3, hdr_info1.voxsz) ) {
        cout << argv[1] << " and " << argv[2]
             << " do not have the same voxel size. Exiting ..." << endl;
        exit(1);
    }
    
    int B[6];  // borders of seg image
    find_borders(seg_img, B);

    // peel one layer of voxels around segmented image and remove all csf
    peel_one_layer_and_remove_csf(seg_img, B);
    
    int seed = find_seed(B);
    cout << "seed " << seed << endl;
    
    // do 3D floodfill starting from seed and fill all encountered voxels
    // fith FINAL_VALUE -> all those voxels are later kept in T1 image,
    // and all the rest are masked out
    Floodfill(seed);

    // mask out any voxels in T1 image that are not equal to FINLA_VALUE
    // in seg_img
    unsigned char *ptr1 = T1_img, *ptr2 = seg_img;
    for(int t=0; t < volume; ++t) {
        if(*ptr2 != FINAL_VALUE) *ptr1 = 0;
        ++ptr1; ++ptr2;
    }
    
    cout << "Writing " << argv[3]  << endl;
    write_analyze(argv[3] , T1_img, hdr_info);
    
    delete [] seg_img;
    delete [] T1_img;
    return 0;
}  // end of main
///////////////////////////////////////////////////////////////////////////
void find_borders(const unsigned char *Img, int *B)
{
    int area = sz[0] * sz[1], volume = area * sz[2];
    int i, j, k;

    // find lowest non-zero z-level
    i=0;
    while( Img[i]==0 && i<volume) i++;
    B[4] = i / area; 
    
    // find highest non-zero z-level
    i=volume-1;
    while( Img[i] == (unsigned char)0 && i>=0) i--;
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
}
////////////////////////////////////////////////////////////////////////////
void peel_one_layer_and_remove_csf(unsigned char *s, int *B) {
    unsigned char *sp;
    int x, y, z;
    for(z=B[4]; z <= B[5]; ++z)
        for(y=B[2]; y <= B[3]; ++y)
            for(x=B[0]; x <= B[1]; ++x) {
                sp = s + z*area + y*sz[0] + x ;
                // if not on the edge or if it is CSF (sCSF) turn to 1
                if(*sp > 1 &&
                   !(*(sp-1)!=0 && *(sp+1)!=0 && *(sp-sz[0])!=0 &&
                     *(sp+sz[0])!=0 && *(sp+area)!=0 && *(sp-area)!=0))
                    *sp = 1;
                if(*sp == 5) *sp = 1;
            }
}

///////////////////////////////////////////////////////////////////////////
int find_seed(int *B) { // fand any non-0, non-1, valued voxel near the middle of
    // the image
    int t = ((B[5]-B[4])/2)*area + ((B[3]-B[2])/2)*sz[0], seed=0;
    int x = (B[1]-B[0])/2;
    bool found = false;
    while( !found && x <= B[1] ) {
        if(seg_img[t + x] > 1) {
            found = true;
            seed = t + x;
            break;
        }
        else ++x;
    }
    if(seed == 0) {
        cout << "Couldn't find seed, exiting " << endl;
        exit(1);
    }
    else return seed;
}
////////////////////////////////////////////////////////////////////////////
bool acceptable(const unsigned char *pos){
    if ( *pos > 1 && *pos != FINAL_VALUE) return true;
    else return false;
}

////////////////////////////////////////////////////////////////////////////
int forward_dir(int index) { // it is assumed that pos is acceptable
    
    int x, y, z;
    z = index / area;
    y = (index - z * area) / sz[0];
    x = index - z * area - y*sz[0];
    
    unsigned char *pos = seg_img + index;
    
    if(x < sz[0]-1) if(acceptable(pos+1)) return 1;
    if(x > 0) if(acceptable(pos-1)) return (-1);
    if(y < sz[1]-1) if(acceptable(pos+sz[0])) return sz[0];
    if(y > 0) if(acceptable(pos-sz[0])) return (-sz[0]);
    if(z < sz[2]-1) if(acceptable(pos+area)) return area;
    if(z > 0) if(acceptable(pos-area)) return (-area);

    return 0;
}
////////////////////////////////////////////////////////////////////////////
void Floodfill(int seed) {  // creating pink voxels
    
    int t, dir;
    stack<int> stk;
    stk.push( seed ); 
    
    while( !stk.empty() ) {
        t = stk.top();
        seg_img[t] = FINAL_VALUE;
        
        dir = forward_dir(t);
        if( dir != 0 ) stk.push(t + dir);
        else stk.pop();   
    }
}
////////////////////////////////////////////////////////////////////////////
