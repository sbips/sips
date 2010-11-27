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

#include "globals.h"


int find_SF(uchar *plane, int xlev, uchar sf_val) {
    int ret = -1;
    bool found = false;
    for(int y=grid.a; y >= grid.pc && !found; --y){
        if(plane[y*sz[0] + xlev] == sf_val){
            found = true;
            ret = y;
        }
    }
    return ret;
}

void draw_first_SF(uchar *plane, int &lyprev, int &ryprev) {
    int x, y, ly, ry, t;
    ly = find_SF(plane, sz[0]-2, LSF_LINE);
    ry = find_SF(plane, 1, RSF_LINE);
    if(ly == -1) ly = lyprev;
    if(ry == -1) ry = ryprev;
    for(x=sz[0]-2; x > grid.m; --x){
        t = ly*sz[0] + x;
        if(plane[t] == ORB_LINE) {
            break;
        }
        else  plane[t] = ORB_LINE;
    }
    for(x=1; x <= grid.m; ++x){
        t = ry*sz[0] + x;
        if(plane[t] == ORB_LINE) {
            break;
        }
        else plane[t] = ORB_LINE;
    }
    lyprev = ly; ryprev = ry;
}

void outline(uchar *mask, int topz, int botz, int *xx, int *yy) {
    // draw outline of orbital frontal lobes in the same way for
    // all axial slices (zbot, ztop]
    uchar *plane;
    int x, y, z, ly, ry, t;
    for(z = topz; z > botz; --z){
        plane = mask + z*area;
        draw_line(xx[0], yy[0], xx[1], yy[1], plane, sz, ORB_LINE);
        draw_line(xx[1], yy[1], xx[2], yy[2], plane, sz, ORB_LINE);
        draw_line(xx[2], yy[2], xx[3], yy[3], plane, sz, ORB_LINE);
        draw_line(xx[3], yy[3], xx[4], yy[4], plane, sz, ORB_LINE);
        draw_line(xx[4], yy[4], xx[5], yy[5], plane, sz, ORB_LINE);
    }
}
///////////////////////////////////////////////////////////////////////

void orbital_frontal(uchar *mask) {
    // temporarily build midline sagital wall from ap-2.5 to ap level
    int x = grid.m, y, z, t; 
    for(z=grid.ap ; z > midpt(grid.axial[2], grid.axial[1]); --z)
        for(y=grid.pc; y < sz[1]-1; ++y)
            mask[z*area + y*sz[0] + x] = M_LINE;
    // this midsagital wall is needed for floodfilling left and rigt
    // OBF lobes

   
    unsigned char *plane;
    int xx[6], yy[6], ztop, zbot;
    int lyprev, ryprev;
    
    // make sure we have a good start for lyprev/ryprev on ap level
    // this may happen if C and SF cross exactely at ap level so
    // that there is no sf on that level
    int lyprevtemp, ryprevtemp;
    z = grid.ap;
    ryprevtemp = find_SF(mask + z*area, 1, RSF_LINE);
    lyprevtemp = find_SF(mask + z*area, sz[0]-2, LSF_LINE);
    
    ryprev = maxval(ryprevtemp, grid.ac);
    lyprev = maxval(lyprevtemp, grid.ac);
    
    // top level: [ap, ap - 0.5)
    ztop = grid.axial[4];  zbot = midpt(grid.axial[3], grid.axial[4]);
    xx[0] = grid.sagital[1];  yy[0] = grid.coronal[15];   
    xx[1] = grid.sagital[1];  yy[1] =grid.coronal[13];
    xx[2] = midpt(grid.sagital[3],grid.m);  yy[2] = grid.ac;
    xx[3] = midpt(grid.m, grid.sagital[5]);  yy[3] = yy[2];
    xx[4] = grid.sagital[7];  yy[4] = yy[1];
    xx[5] = grid.sagital[7];  yy[5] = yy[0];  
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }
    
    // next level: [ap - 0.5, ap - 1)
    ztop =  zbot;   zbot = grid.axial[3];
    xx[0] = grid.sagital[1];  yy[0] = midpt(grid.coronal[14],grid.coronal[15]);
    xx[1] = grid.sagital[2];  yy[1] = grid.coronal[12];
    xx[2] = grid.sagital[3];  yy[2] = grid.ac;
    xx[3] = grid.sagital[5];  yy[3] = yy[2];
    xx[4] = grid.sagital[6];  yy[4] = yy[1];
    xx[5] = grid.sagital[7];  yy[5] = yy[0];      
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }
    
    // next level: [ap - 1, ap - 1.25) OKOKOKOKOKOKOKOK
    ztop =  zbot;   zbot = midpt(ztop, midpt(grid.axial[2], grid.axial[3])); 
    xx[0] = grid.sagital[1];  yy[0] = grid.coronal[16];        
    xx[1] = grid.sagital[2];  yy[1] = grid.coronal[12];
    xx[2] = grid.sagital[3];  yy[2] = grid.ac;
    xx[3] = grid.sagital[5];  yy[3] = yy[2];
    xx[4] = grid.sagital[6];  yy[4] = yy[1];
    xx[5] = grid.sagital[7];  yy[5] = yy[0];
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }

    // next level: [ap - 1.25, ap - 1.5) OKOKOKOKOKOKOKOK
    ztop =  zbot;   zbot = midpt(grid.axial[3], grid.axial[2]); 
    xx[0] = grid.sagital[1] + round((grid.sagital[2]-grid.sagital[1])/3);
    yy[0] = grid.coronal[15];        
    xx[1] = grid.sagital[2];  yy[1] = grid.coronal[12];
    xx[2] = grid.sagital[3];  yy[2] = grid.ac;
    xx[3] = grid.sagital[5];  yy[3] = yy[2];
    xx[4] = grid.sagital[6];  yy[4] = yy[1];
    xx[5] = grid.sagital[7] - round((grid.sagital[7]-grid.sagital[6])/3);
    yy[5] = yy[0];
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }
    
    // next level: [ap - 1.5, ap - 1.75) OKOKOKOKOKOKOKOK
    ztop =  zbot;   zbot = midpt(ztop, grid.axial[2]);
    xx[0] = grid.sagital[2];  yy[0] = grid.coronal[17];        
    xx[1] = grid.sagital[2];  yy[1] = midpt(grid.coronal[12], grid.coronal[13]);
    xx[2] = midpt(grid.sagital[3], grid.m);  yy[2] = grid.ac;
    xx[3] = midpt(grid.m, grid.sagital[5]);  yy[3] = yy[2];
    xx[4] = grid.sagital[6];  yy[4] = yy[1];
    xx[5] = grid.sagital[6];  yy[5] = yy[0]; 
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }
    
    // next level: [ap - 1.75, ap - 2) OKOKOKOKOKOKOKOK
    ztop =  zbot;   zbot = grid.axial[2]; 
    xx[0] = grid.sagital[1];  yy[0] = grid.coronal[18];        
    xx[1] = midpt(grid.sagital[2], grid.sagital[3]);
    yy[1] = midpt(grid.coronal[12], grid.coronal[13]);
    xx[2] = midpt(grid.sagital[3], grid.m);  yy[2] = grid.ac;
    xx[3] = midpt(grid.m, grid.sagital[5]);  yy[3] = yy[2];
    xx[4] = midpt(grid.sagital[5], grid.sagital[6]);
    yy[4] = yy[1];
    xx[5] = grid.sagital[7];  yy[5] = yy[0]; 
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_first_SF(plane, lyprev, ryprev);
    }
    
    // last level: [ap - 2, ap - 2.5) OKOKOKOKOKOKOKOKOKOKOK
    ztop = zbot;  zbot = midpt(grid.axial[2], grid.axial[1]);
    xx[0] = grid.sagital[2]; yy[0] = grid.a;
    xx[1] = grid.sagital[3]; yy[1] = grid.coronal[13];
    xx[2] = midpt(grid.sagital[3], grid.m);  yy[2] = grid.ac;
    xx[3] = midpt(grid.m, grid.sagital[5]);  yy[3] = yy[2];
    xx[4] = grid.sagital[5];  yy[4] = yy[1];
    xx[5] = grid.sagital[6];  yy[5] = yy[0];        
    outline(mask, ztop, zbot, xx, yy);
    for(z=ztop; z > zbot; --z){
        plane = mask + z*area;
        draw_line(1, yy[0], xx[0], yy[0], plane, sz, ORB_LINE);
        draw_line(xx[5], yy[5], sz[0]-2, yy[5], plane, sz,ORB_LINE);
    }
    
    int dirs[] = {1, sz[1], -1, -sz[1]};
    bool found;
    for(z=grid.axial[4]; z > zbot; --z) {
        plane = mask + z*area;
        // fill with ROBF 
        Floodfill robf_f(plane + (sz[1]-2)*sz[0] + grid.r,
                          dirs, 4, acceptable_OBF, ROBF);
        // fill with LOBF 
        Floodfill lobf_f(plane + (sz[1]-2)*sz[0] + grid.l,
                         dirs, 4, acceptable_OBF, LOBF);
        // reassign RMOBF lobe
        for(x=grid.sagital[3]; x <= grid.m; ++x) {
            for(y=sz[1]-2; y > grid.pc; --y){
                t = y*sz[0] + x;
                if(plane[t] == ROBF) plane[t] = RMOBF;
            }
        }
        // reassign LMOBF lobe
        for(x=grid.m + 1; x <= grid.sagital[5]; ++x) {
            for(y=sz[1]-2; y > grid.pc; --y){
                t = y*sz[0] + x;
                if(plane[t] == LOBF ) plane[t] = LMOBF;
            }
        }
        // reassign ORB  line
        for(x=1; x < sz[0]-1; ++x) {
            for(y=grid.pc; y <= grid.a; ++y){
                t = y*sz[0] + x;
                if(plane[t] == ORB_LINE)
                    plane[t] = plane[(sz[1]-2)*sz[0] + x] ;
            }
        }
        // reassign midsagital wall
        for(y=sz[1]-2; y >= grid.pc; --y){
            plane[y*sz[0] + grid.m]=plane[y*sz[0]+grid.m-1];
        }
    }
    // reassign SF lines
    for(z=1; z < grid.t; ++z) {
        for(x=1; x < sz[0]-1; ++x) {
            for(y=grid.p; y <= grid.a; ++y){
                t = z*area +  y*sz[0] + x;
                if(mask[t] == LSF_LINE || mask[t] == RSF_LINE)
                    mask[t] = mask[(z-1)*area + y*sz[0] + x] ;
            }
        }
    }
    
}
