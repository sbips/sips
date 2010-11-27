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

#ifndef EDGE_H
#define EDGE_H

#include "globals.h"

/* find first  brain voxel going along line starting at (x,y,z) in
   given direction -> return edge found this way (slice #) */
template <class T>
int find_edge(T *img, int x, int y, int z, dir D) {
    // first determine starting point for edge, and his increase and stride
    // then determine the other two strides: stride1, stride2 
    int edge, stride, increase, stride1, stride2, s1, s2;
    int i, j;
    switch(D.d) {
    case ( XX): edge = x; increase = 1; stride = 1;
               i = y; j= z; stride1 = sz[0]; stride2 = area;
               break;
    case (_XX): edge = x; increase = -1;  stride = 1;
               i = y; j= z; stride1 = sz[0]; stride2 = area;
               break; 
    case ( YY):edge = y; increase = 1; stride = sz[0];
               i = x; j= z; stride1 = 1; stride2 = area;
               break;
    case (_YY): edge = y; increase = -1; stride = sz[0];
               i = x; j= z; stride1 = 1; stride2 = area;
               break;
    case ( ZZ): edge = z;increase = 1; stride = area;
               i = x; j= y; stride1 = 1; stride2 = sz[0];
               break;
    case (_ZZ): edge = z; increase = -1; stride = area;
               i = x; j= y; stride1 = 1; stride2 = sz[0];
               break;
    default: break;
    }
   
    bool found = false;
    while( D.min<= edge  && edge <= D.max ) {
         if( img[edge*stride + i*stride1 + j*stride2] ) {
             found = true;
             return edge;
         }
         if(found) break;
         else edge += increase;
    }
    return edge;
}
/* /////////////////////////////////////////////////////////////////////// */

template <class T>
int find_first_slice(T *img, int e, dir D) {
    /* starting at edge move in and find first slice which is not all 0's */
    int stride, increase, stride1, stride2, s1, s2;
    int edge = e, beg1, end1, beg2, end2, i, j;
    switch(D.d) {
    case ( XX): increase = 1; stride = 1;
        beg1 = 0; end1 = sz[1]; stride1 = sz[0];
        beg2 = 0; end2 = sz[2]; stride2 = area;
    	cout << "beg1: " << beg1 << " " << "beg2: " << beg2 << endl;
        break;
    case (_XX): increase = -1;  stride = 1;
        beg1 = 0; end1 = sz[1]; stride1 = sz[0];
        beg2 = 0; end2 = sz[2]; stride2 = area;
        break; 
    case ( YY): increase = 1; stride = sz[0];
        beg1 = 0; end1 = sz[0]; stride1 = 1;
        beg2 = 0; end2 = sz[2]; stride2 = area;
        break;
    case (_YY): increase = -1; stride = sz[0];
        beg1 = 0; end1 = sz[0]; stride1 = 1;
        beg2 = 0; end2 = sz[2]; stride2 = area;
        break;
    case ( ZZ): increase = 1; stride = area;
        beg1 = 0; end1 = sz[0]; stride1 = 1;
        beg2 = 0; end2 = sz[1]; stride2 = sz[0];
        break;
    case (_ZZ): increase = -1; stride = area;
        beg1 = 0; end1 = sz[0]; stride1 = 1;
        beg2 = 0; end2 = sz[1]; stride2 = sz[0];
        break;
    default: break;
    }
    bool found = false;
    int ret = D.min;
    while( D.min<= edge  && edge <= D.max && !found) {
        for(i=beg1; i < end1 && !found; ++i)
            for(j=beg2; j < end2 && !found; ++j){
                if( img[edge*stride + i*stride1 + j*stride2] ) {
                    found = true;
                    ret = edge;
                }
            }
        edge += increase;
    }
    return ret;
}
//////////////////////////////////////////////////////////////////////////
/* start going along line starting at (x,y,z) in given direction (from D.min to D.max). Return the first slice where square centered around line contains non-zero voxels ->edge */
template <class T>
int find_first_square(T *img, int x, int y, int z, dir D) {
    // first determine starting point for edge, and his increase and stride
    // then determine the other two strides: stride1, stride2 
    int edge, stride, increase, stride1, stride2, s1, s2;
    int i, j;
    switch(D.d) {
    case ( XX): edge = x; increase = 1; stride = 1;
               i = y; j= z; stride1 = sz[0]; stride2 = area;
               break;
    case (_XX): edge = x; increase = -1;  stride = 1;
               i = y; j= z; stride1 = sz[0]; stride2 = area;
               break; 
    case ( YY):edge = y; increase = 1; stride = sz[0];
               i = x; j= z; stride1 = 1; stride2 = area;
               break;
    case (_YY): edge = y; increase = -1; stride = sz[0];
               i = x; j= z; stride1 = 1; stride2 = area;
               break;
    case ( ZZ): edge = z;increase = 1; stride = area;
               i = x; j= y; stride1 = 1; stride2 = sz[0];
               break;
    case (_ZZ): edge = z; increase = -1; stride = area;
               i = x; j= y; stride1 = 1; stride2 = sz[0];
               break;
    default: break;
    }
   
    bool found = false;
    int ii, jj;
    int ret = D.min;
    while( D.min<= edge  && edge <= D.max && !found ) {
        for(ii=i-10; ii <= i+10 && !found; ++ii)
            for(jj=j-10; jj <= j+10 && !found; ++jj) {
                if( img[edge*stride + ii*stride1 + jj*stride2] ) {
                    found = true;
                    ret = edge;
                }
            }
       edge += increase;
    }
    return ret;
}
#endif
