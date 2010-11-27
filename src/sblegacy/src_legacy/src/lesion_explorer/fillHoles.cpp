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
// /home/brainlab/gcc-2.95.3/bin/g++ -O2 -o fillHoles fillHoles.cpp

#include "../../../include_from_cortex/analyze_utils.h"
#include <stack>

//////////////////////// GLOBAL VARIABLES //////////////////////////////////
unsigned char *seg_img;
int sz[3], area, volume;
unsigned char FINAL_VALUE = 253;
////////////////////////////// FUNCTIONS ///////////////////////////////////
void Floodfill(int seed);
////////////////////////////////  MAIN //////////////////////////////////////

using namespace nk;

int main(int argc, char ** argv)
{
    if(argc != 3 ){
        cout <<"Usage: fillHoles <inFname (Binary volume)> <outFname>"<< endl;
        
        exit(1);
    }
   
    info hdr_info;
    seg_img = read_analyze<unsigned char>(argv[1], hdr_info);
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
   
    area = sz[0]*sz[1];
    volume = area * sz[2];
    
    //Define seed at the start of second slice (Just an educated guess!)
    int seed=256*256+1;
    
    // do 3D floodfill starting from seed and fill all encountered voxels
    // fith FINAL_VALUE -> all those voxels are later kept in T1 image,
    // and all the rest are masked out
    Floodfill(seed);

    //Prepare output
    unsigned char *ptr1 = new unsigned char[volume];
    memset(ptr1,0,volume);
    unsigned char *out_img=ptr1;

    unsigned char *ptr2 = seg_img;
    for(int t=0; t < volume; ++t) {
        if(*ptr2 != FINAL_VALUE) *ptr1 = 1;
        ++ptr1; ++ptr2;
    }

    write_analyze(argv[2] , out_img, hdr_info);
    
    delete [] seg_img;
    delete [] ptr1;

    return 0;
}  // end of main

////////////////////////////////////////////////////////////////////////////
bool acceptable(const unsigned char *pos){
    if ( *pos == (unsigned char)0  && *pos != FINAL_VALUE) return true;
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
