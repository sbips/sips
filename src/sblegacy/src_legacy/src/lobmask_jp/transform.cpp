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
// /srutility/gnu/gcc-2.95.2/bin/c++ -O2 -o transform transform.cpp

// Note: mat saved in Analyze actually represents acpc_to_T1.air
// so that air file is created first and then invert_air is applied

extern "C" {
#include "/net/cortex/csf/pub/include/write_air16.c"
}
#include "/net/cortex/csf/pub/include/analyze_utils.h"
#include <stdio.h>
#include <stdlib.h>

using namespace nk;

const char *AIR_PATH = "/srutility/AIR5.2.5";

void read_matrix(double e[4][4], const char *matfname) {
   ifstream mat(matfname);
   if(!mat) { cout << "Matrix file not found. Exit.. " << endl; exit(1);}

   int i, j;
   for(j=0; j < 4; j++){
        for(i=0; i < 4; i++) {
            mat >> e[j][i];
        }
   }
   // multiply shifts by -1, don't know why
   e[3][0] *= -1.0; e[3][1] *= -1.0; e[3][2] *= -1.0;
}

int main(int argc, char **argv) {
    if(argc < 5) {
        cout<< "
Usage : transform <matrix_fname>
                  <air_fname>
                  <T1_fname>
                  <T1_fc_fname>

Example: transform xxx.mat xxx_T1_to_acpc.air xxx_T1 xxx_T1_fc "
            << endl;
        exit(1);
    }
    string mat_fname = string(argv[1]);
    string air_fname = string(argv[2]);
    string standard_fname = string(argv[3]); // standard for acpc_to_T1
    string reslice_fname = string(argv[4]);  // reslice for acpc_to_T1
    string invertair_fname = "acpc_to_T1.air";
    
    info standard_info;
    read_analyze_header(standard_fname + ".hdr", standard_info);
    info reslice_info;
    read_analyze_header(reslice_fname +".hdr", reslice_info);

    // make inverted air file first
    air16 Air;
    air16 *air = &Air;
    strncpy(air->s_file, standard_fname.c_str() , 127);
    strncpy(air->r_file, reslice_fname.c_str() , 127);
    air->s.bits=8;
    air->s.x_dim=standard_info.sz[0];
    air->s.y_dim=standard_info.sz[1];
    air->s.z_dim=standard_info.sz[2];
    air->s.x_size=standard_info.voxsz[0];
    air->s.y_size=standard_info.voxsz[1];
    air->s.z_size=standard_info.voxsz[2];
    air->r.bits=8;
    air->r.x_dim=reslice_info.sz[0];
    air->r.y_dim=reslice_info.sz[1];
    air->r.z_dim=reslice_info.sz[2] + 30;
    air->r.x_size=reslice_info.voxsz[0];
    air->r.y_size=reslice_info.voxsz[1];
    air->r.z_size=reslice_info.voxsz[2];
    air->s_hash=1;
    air->r_hash=2;
    air->s_volume=0;
    air->r_volume=0;

    // read .mat file into air's matrix (interpret it as acpc_to_T1)
    double e[4][4];
    read_matrix(e, mat_fname.c_str());
    
    // write air's matrix as air file 
    write_air16(invertair_fname.c_str(), 1, e, 1, air);

    char cmd[1000];
    // get T1_to_acpc air file 
    sprintf(cmd, "invert_air %s %s y ", invertair_fname.c_str(),
            air_fname.c_str() );
    system(cmd);

    // remove inverted air file
    sprintf(cmd, "rm -f %s ", invertair_fname.c_str());
    system(cmd);
    return 0;
}

   
