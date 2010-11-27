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

#ifndef MAKE_GRID_H
#define MAKE_GRID_H

#include "globals.h"
#include "edge.h"
#include "../../../include_from_cortex/analyze_utils.h"
#include <string>
#include <iostream.h>
#include <fstream.h>

void print_usage_msg() {
      cout << "Usage: lobmask
       <eroded image fname>
       <object filename>
       <maskfilename>
       [-g <gridfname>], name for the grid text file, by default grid.txt
       -ac <int>, ac coronal slice #
       -pc <int>, pc coronal slice #
       -ap <int>, ap axial slice #
       -m <int>, midsagital slice #
       -rpron <int>, right pron point
       -lpron <int>, left pron point
       [-right <int>], sagital slice # of right edge
       [-left <int>], sagital slice # of left edge
       [-at <int>], coronal slice # of the anterior edge
       [-pt <int>], coronal slice # of the posterior edge
       [-b <int>], axial slice # of bottom
       [-t <int>], axial slice # of top " << endl ;
      exit(1);
}

using namespace nk;

template <class T>
void make_grid(int argc, char **argv) {
    
    int p, pc, ac, a; /* special coronal planes */
    int b, ap, t;     /* special axial planes */
    int l, r, m;      /* special sagital planes: left , right, medial */
    p=pc=ac=a=b=ap=t=l=r=m=rpron=lpron=-1;
    
    CommandLine cl(argc, argv);
    
    string gridfname("grid.txt");
    if(cl.find("-g"))   gridfname = cl[1];
    
    ifstream readgrid(gridfname.c_str());
    if(readgrid) {
        string landmark;
        while(!readgrid.eof()) {
            readgrid >> landmark;
            if(landmark == "-ac") { readgrid >> ac; continue;}
            if(landmark == "-pc") { readgrid >> pc; continue;}
            if(landmark == "-ap") { readgrid >> ap; continue;}
            if(landmark == "-m" ) { readgrid >> m; continue;}
            if(landmark == "-rpron") { readgrid >> rpron; continue;}
            if(landmark == "-lpron") { readgrid >> lpron; continue;}
            
            if(landmark == "-pt") { readgrid >> p; continue;}
            if(landmark == "-at") { readgrid >> a; continue;}
            if(landmark == "-b") { readgrid >> b; continue;}
            if(landmark == "-t") { readgrid >> t; continue;}
            if(landmark == "-right") { readgrid >> r; continue;}
            if(landmark == "-left") { readgrid >> l; continue;}
        }
        readgrid.close();
    }
    /* if grid file was not found than look for command line inputs */
    else {
        if(cl.find("-ac")) ac = atoi(cl[1]);
        if(cl.find("-pc")) pc = atoi(cl[1]);
        if(cl.find("-ap")) ap = atoi(cl[1]);
        if(cl.find("-m")) m = atoi(cl[1]);
        if(cl.find("-lpron")) lpron = atoi(cl[1]);
        if(cl.find("-rpron")) rpron = atoi(cl[1]);
        
        if(cl.find("-right")) r = atoi(cl[1]);
        if(cl.find("-left")) l = atoi(cl[1])-1;
        if(cl.find("-pt")) p = atoi(cl[1])-1;
        if(cl.find("-at")) a = atoi(cl[1])-1;
        if(cl.find("-b")) b = atoi(cl[1])-1;
        if(cl.find("-t")) t = atoi(cl[1])-1;
        
    }
    
    /* check that mandatory landmarks have been read */
    if( ac<0 || pc<0 || ap<0 || m<0 || rpron<0 || lpron<0 ) {
        cout << "Error in reading mandatory landmarks. Exiting .." << endl;
        exit(1);
    }

    /* sustract 1 from all landmarks, since they were acquired in Analyze */
    --ac; --pc; --ap; --m; --lpron; --rpron; --r; --l; --p; --a; --b; --t;
    
    info hdr_info;
    T *img = read_analyze<T>(string(argv[1]), hdr_info);
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
    area = sz[0] * sz[1];
    volume = area * sz[2];
    
    /* calculate optional landmarks if they were not provided */
    /* (m,pc,ap) = 3D-coord of pc point */
    if(r < 0) r = find_first_slice<T>(img, 0, dir(XX, 0, m));
    if(l < 0) l = find_first_slice<T>(img, sz[0]-1, dir(_XX, m, sz[0]-1));
    if(p < 0) p = find_first_slice<T>(img, 0, dir(YY, 0, pc));
    if(a < 0) a = find_first_slice<T>(img, sz[1]-1, dir(_YY, ac, sz[1]-1));
    if(t < 0) t = find_first_slice<T>(img, sz[2]-1, dir(_ZZ, ap, sz[2]-1));
    if(b < 0) b = find_first_square<T>(img, midpt(r, m), ac, 0, dir(ZZ, 0, ap));

    /* output full grid file */
    ofstream writegrid(gridfname.c_str());
    if(!writegrid) {
        cout << "Couldn't open " << gridfname << " for output. Permissions?" << endl;
    }
    writegrid << "-pt " << p+1 << " -pc " << pc+1 << " -ac " << ac+1
              << " -at " << a+1
        << " -b " << b+1 << " -ap " << ap+1 << " -t " << t+1
        << " -right " << r+1 << " -m " << m+1 << " -left " << l+1
        << " -rpron " << rpron+1 << " -lpron " << lpron+1 << endl;
    
    grid.set(p, pc, ac, a, b, ap, t, r, m, l);
    grid.print_points();

    delete [] img;

}


#endif
