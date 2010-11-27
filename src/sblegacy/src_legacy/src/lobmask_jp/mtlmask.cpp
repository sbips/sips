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
// /home/brainlab/gcc-2.95.3/bin/g++ -o mtlmask mtlmask.cpp

#include "../../../include_from_cortex/grid3.h"

using namespace nk;

void print_grid_error_message() {
    cout << "Misssing grid landmarks... Exiting..." << endl;
    exit(1);
}

Grid grid;
int sz[3], area, volume;
enum Value {BACKGROUND=0, LMTL=40, RMTL=41};

int main(int argc, char ** argv) {

    if(argc < 4) {
        cout << "
Usage: mtlmask <grid fname> <acpc fname> <outmask fname>

Example: mtlmask xxx_grid.txt xxx_acpc xxx_mtlmask
"
             << endl;
        exit(1);
    }

    
    ifstream ingrid(argv[1]);
    vector<string> item;
    string temp;

    while(!ingrid.eof()) {
        ingrid >> temp;
        item.push_back(temp);
        cout << item.back() << "\t" ;
    }

    int p, pc, ac, a; /* special coronal planes */
    int b, ap, t;     /* special axial planes */
    int l, r, m;      /* special sagital planes: left , right, medial */
    
    vector<string>::iterator it;

    if((it = find(item.begin(), item.end(), "-ac")) != item.end()) {
        ++it;
        ac = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
    if((it = find(item.begin(), item.end(), "-pc")) != item.end()) {
        ++it;
        pc = atoi((*it).c_str())-1;
    }
    else print_grid_error_message(); 
    if((it = find(item.begin(), item.end(), "-ap")) != item.end()) {
        ++it;
        ap = atoi((*it).c_str())-1;
    }
    else print_grid_error_message(); 
    if((it = find(item.begin(), item.end(), "-m")) != item.end()) {
        ++it;
        m = atoi((*it).c_str())-1;
    }
    else print_grid_error_message(); 
    if((it = find(item.begin(), item.end(), "-right")) != item.end()) {
        ++it;
        r = atoi((*it).c_str())-1;
    }
    if((it = find(item.begin(), item.end(), "-left")) != item.end()) {
        ++it;
        l = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
    if((it = find(item.begin(), item.end(), "-pt")) != item.end()) {
        ++it;
        p = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
    if((it = find(item.begin(), item.end(), "-at")) != item.end()) {
        ++it;
        a = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
    if((it = find(item.begin(), item.end(), "-b")) != item.end()) {
        ++it;
        b = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
    if((it = find(item.begin(), item.end(), "-t")) != item.end()) {
        ++it;
        t = atoi((*it).c_str())-1;
    }
    else print_grid_error_message();
  
    grid.set(p, pc, ac, a, b, ap, t, r, m, l);
    // grid.print_points();

    info hdr_info;
    read_analyze_header(string(argv[2]) + ".hdr", hdr_info);
    copy(hdr_info.sz, hdr_info.sz + 3, sz);
    area = sz[0] * sz[1];
    volume = area * sz[2];

    unsigned char *mask = new unsigned char [volume];
    if(!mask) {
        cout << "memory couldn't be allocated. Exit " << endl;
        exit(1);
    }
    fill(mask, mask+volume, BACKGROUND);
    // mtlboxes;
    int ytop, ybot, ztop, zbot, xr, xl, x, y, z;

    xr = grid.sagital[2];
    xl = grid.sagital[6];
    zbot = 0;
    
    ytop = grid.inverse_talairach_coord1(7);
    ybot = grid.inverse_talairach_coord1(0);
    ztop = grid.inverse_talairach_coord2(-15);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }

    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-7);
    ztop = grid.inverse_talairach_coord2(-10);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }
    
    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-16);
    ztop = grid.inverse_talairach_coord2(-10);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }
    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-24);
    ztop = grid.inverse_talairach_coord2(-7);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }
    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-28);
    ztop = grid.inverse_talairach_coord2(-7);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }
    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-30);
    ztop = grid.inverse_talairach_coord2(-4);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y > ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }
    ytop = ybot;
    ybot = grid.inverse_talairach_coord1(-42);
    ztop = grid.inverse_talairach_coord2(0);
    for(z=ztop; z >= zbot; --z) {
        for(y=ytop; y >= ybot; --y){
            for(x=xr; x <=m; ++x) mask[z*area + y*sz[0] +x] = RMTL;
            for(x=m+1; x <= xl; ++x) mask[z*area + y*sz[0] +x] = LMTL;
        }
    }

    cout << "\nWriting " << argv[3] << endl;
    write_analyze(argv[3], mask,  hdr_info);

    delete [] mask;
    return 0;
}
