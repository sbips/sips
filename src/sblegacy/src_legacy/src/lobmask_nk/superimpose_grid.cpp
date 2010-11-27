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
// g++  -o superimpose_grid superimpose_grid.cpp /home/nkovacev/analyze_utils_new/analyze_utils.cpp

#include "analyze_utils.h"
#include "command_line.h"
#include <assert.h>
#include <iomanip.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream.h>
#include "globals.h"
#include "grid.h"
#include "edge.h"

using namespace nk;

void print_usage_msg() {
      cout << "Usage: superimpose_grid
       <image fname>, 8 or 16 bit image (eroded or masked)
       <output fname>, image with gridlines = 0
       -ac <int>, ac coronal slice #
       -pc <int>, pc coronal slice #
       -ap <int>, ap axial slice #
       -m <int>, midsagital slice #
       [-right <int>], sagital slice # of right edge
       [-left <int>], sagital slice # of left edge
       [-at <int>], coronal slice # of the anterior edge
       [-pt <int>], coronal slice # of the posterior edge
       [-b <int>], axial slice # of bottom
       [-t <int>], axial slice # of top" << endl ;
      exit(1);
}

template <class T>
void superimpose_grid(int argc, char **argv);

int main(int argc, char ** argv)
{
    if(argc < 7) print_usage_msg();
 
    info hdr_info;
    string fname(argv[1]);
    read_analyze_header(fname + ".hdr", hdr_info);
    
    if(hdr_info.datatype == 2) superimpose_grid<uchar>(argc, argv);
    else  superimpose_grid<short>(argc, argv);
    
    return 0;
}  // end of main
////////////////////////////////////////////////////////////////////////////
template <class T>
void superimpose_grid(int argc, char **argv) {
    int p, pc, ac, a; /* special coronal planes */
    int b, ap, t;     /* special axial planes */
    int l, r, m;      /* special sagital planes: left , right, medial */
    
    gs::CommandLine cl(argc, argv);
    
    if(cl.find("-ac")) ac = atoi(cl[1]);
    else print_usage_msg(); 
    if(cl.find("-pc")) pc = atoi(cl[1]);
    else print_usage_msg(); 
    if(cl.find("-ap")) ap = atoi(cl[1]);
    else print_usage_msg(); 
    if(cl.find("-m")) m = atoi(cl[1]);
    else print_usage_msg(); 
 
    info hdr_info;
    T *img = read_analyze<T>(string(argv[1]), hdr_info);
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
    // cout << "size " << sz[0] << ", " << sz[1] << "," << sz[2] << endl;
    area = sz[0] * sz[1];
    volume = area * sz[2];
    
   
    /* (m,pc,ap) = 3D-coord of pc point */
    if(cl.find("-right")) r = atoi(cl[1]);
    else r = find_edge<T>(img, 0, pc, ap, dir(XX, 0, m));
    if(cl.find("-left")) l = atoi(cl[1]);
    else l = find_edge<T>(img, sz[0]-1, pc, ap, dir(_XX, m, sz[0]-1)); 
    cout << "r, l " << r << ", " << l<< endl;
    
    if(cl.find("-pt")) p = atoi(cl[1]);
    else p = find_edge<T>(img, m+10, 0, ap, dir(YY, 0, pc));
    if(cl.find("-at")) a = atoi(cl[1]);
    else a = find_edge<T>(img, m+10, sz[1]-1, ap, dir(_YY, ac,sz[1]-1));
    cout << "p, a " << p << ", " << a<< endl;

    if(cl.find("-b")) b = atoi(cl[1]);
    else b = find_edge<T>(img, midpt(r, m), ac, 0, dir(ZZ, 0, ap));
    if(cl.find("-t")) t = atoi(cl[1]);
    else t = find_edge<T>(img, m, pc, sz[2]-1, dir(_ZZ, ap, sz[2]-1));
    cout << "b, t " << b << ", " << t<< endl;
    
    grid.set(p, pc, ac, a, b, ap, t, r, m, l);
    grid.print_points();
    uchar *mask = new uchar [volume];
    assert(mask);
    
    uchar *mask_ptr;
    T *img_ptr;
    int x, y, z, x1, y1, z1, x2, y2, z2, index;
    
    // highlight pc - 4, pc line
    y1 = grid.coronal[4]; y2 = grid.coronal[8];
    for(z=0; z < sz[2]; ++z)
        for(x=0; x < sz[0]; ++x) {
            index = z*area + x;
            img[index + y1*sz[0]] = img[index + y2*sz[0]] = 0;
        }

    // highlight ap line
    z1 = grid.axial[4];
    for(y=0; y < sz[1]; ++y)
        for(x=0; x < sz[0]; ++x){
            img[z1*area + y*sz[0] + x] = 0;
        }
    
    cout << "Writing " << argv[2]<< endl;
    write_analyze(string(argv[2]), img,  hdr_info);
    delete [] img;
    delete [] mask;
} 
