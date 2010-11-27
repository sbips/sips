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

void basal_ganglia(uchar *mask) {
    int x, y, z, xpt[9], ypt[9], Zbot, Ztop, zbot, ztop;
    // temporarily build midline sagital wall for all levels of BGT
    x = grid.m;
    Zbot = midpt(grid.axial[3], grid.axial[4]);
    Ztop = grid.axial[6]; 
    for(z=Zbot; z <= Ztop; ++z)
        for(y= midpt(grid.coronal[7], grid.coronal[8]); y <= grid.coronal[13] ; ++y)
            mask[z*area + y*sz[0] + x] = M_LINE;
    // this midsagital wall is needed for floodfilling left and rigt
    // BGT lobes
    
    // basal ganglia is divided into 3 z levels
    int i;
    uchar *plane;
    int dirs[] = {1, sz[1], -1, -sz[1]};
   
     // level [ap - 0.5, ap - 0.25]
    zbot = Zbot;  ztop = midpt(zbot, grid.ap);
    xpt[0] = grid.sagital[3] - ((grid.sagital[3]- grid.sagital[2])/3);
    ypt[0] = grid.coronal[13];
    xpt[1] = grid.sagital[5] + ((grid.sagital[6]- grid.sagital[5])/3);
    ypt[1] = ypt[0];
    xpt[2] =  grid.sagital[6] - ((grid.sagital[6]- grid.sagital[5])/3);
    ypt[2] = grid.ac;
    xpt[3] = xpt[2];  ypt[3] = grid.coronal[10];
    xpt[4] = grid.sagital[2] +  ((grid.sagital[3]- grid.sagital[2])/3);
    ypt[4] = ypt[3];
    xpt[5] = xpt[4]; ypt[5] = ypt[2];
    xpt[6]= xpt[0];  ypt[6] = ypt[0];   // close the curve
    for(z=zbot; z <= ztop; ++z) {
        plane = mask + z*area;
        for(i=0; i < 6; ++i)
            draw_line(xpt[i], ypt[i], xpt[i+1], ypt[i+1], plane, sz, BASG_LINE);
        Floodfill LABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LABGT);
        Floodfill RABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RABGT);

        Floodfill LPBGT_f(plane + (ypt[3] + 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LPBGT);
        Floodfill RPBGT_f(plane + (ypt[3] + 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RPBGT);
    }
    // level (ap-0.25, ap + 1.5]
    zbot = ztop + 1;  ztop = midpt(grid.axial[5], grid.axial[6]);
    xpt[0] = grid.sagital[3] - ((grid.sagital[3]- grid.sagital[2])/3);
    ypt[0] = grid.coronal[13];
    xpt[1] = grid.sagital[5] + ((grid.sagital[6]- grid.sagital[5])/3);
    ypt[1] = ypt[0];
    xpt[2] =  grid.sagital[6] - ((grid.sagital[6]- grid.sagital[5])/3);
    ypt[2] = grid.ac;
    xpt[3] = xpt[2];  ypt[3] = grid.coronal[9];
    xpt[4] = xpt[1];  ypt[4] = midpt(grid.coronal[7], grid.coronal[8]);
    xpt[5] = xpt[0];  ypt[5] = ypt[4];
    xpt[6] = grid.sagital[2] +  ((grid.sagital[3]- grid.sagital[2])/3);
    ypt[6] = ypt[3];
    xpt[7] = xpt[6];  ypt[7] = ypt[2];
    xpt[8] = xpt[0];  ypt[8] = ypt[0]; // close the curve
    for(z=zbot; z <= ztop; ++z) {
        plane = mask + z*area;
        for(i=0; i < 8; ++i)
            draw_line(xpt[i], ypt[i], xpt[i+1], ypt[i+1], plane, sz, BASG_LINE);
        Floodfill LABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LABGT);
        Floodfill RABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RABGT);
        Floodfill LPBGT_f(plane + (ypt[4] + 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LPBGT);
        Floodfill RPBGT_f(plane + (ypt[4] + 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RPBGT);
    }                   
    // level (ap + 1.5, ap + 2]
    zbot = ztop + 1; ztop = grid.axial[6];
    xpt[0] = grid.sagital[3];  ypt[0] = grid.coronal[13];
    xpt[1] = grid.sagital[5];  ypt[1] = ypt[0];
    xpt[2] = xpt[1];  ypt[2] = midpt(grid.coronal[7], grid.coronal[8]);
    xpt[3] = xpt[0];  ypt[3] = ypt[2];
    xpt[4] = xpt[0];  ypt[4] = ypt[0]; // close the curve
    for(z=zbot; z <= ztop; ++z) {
        plane = mask + z*area;
        for(i=0; i < 4; ++i)
            draw_line(xpt[i], ypt[i], xpt[i+1], ypt[i+1], plane, sz, BASG_LINE);
        Floodfill LABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LABGT);
        Floodfill RABGT_f(plane + (ypt[0] - 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RABGT);
        Floodfill LPBGT_f(plane + (ypt[2] + 1)*sz[0] + (grid.m + 1),
                         dirs, 4, acceptable_BGT, LPBGT);
        Floodfill RPBGT_f(plane + (ypt[2] + 1)*sz[0] + (grid.m - 1),
                         dirs, 4, acceptable_BGT, RPBGT);
    }

    // reassign BASG_LINE and mid sagital wall
    bool found_C;
    int t;
    for(z=Zbot-1; z <= Ztop+1; ++z){
        found_C = false;
        for(y=grid.coronal[7]; y <= grid.coronal[14] ; ++y) {
            for(x=grid.sagital[1]; x <= grid.m; ++x){
                t = z*area + y*sz[0] + x;
                if(mask[t + sz[0]] == RC_LINE) found_C = true;
                if( mask[t] == BASG_LINE){
                    if(!found_C)  mask[t] = RPBGT;
                    else mask[t] = RABGT;
                }
            }
        }
        found_C = false;
        for(y=grid.coronal[7]; y <= grid.coronal[14] ; ++y) {
            for(x=grid.sagital[7]; x > grid.m; --x){
                t = z*area + y*sz[0] + x;
                if(mask[t + sz[0]] == LC_LINE ) found_C = true;
                if( mask[t] == BASG_LINE){
                    if(!found_C)  mask[t] = LPBGT;
                    else mask[t] = LABGT;
                }
            }
        }
    }
    x = grid.m;
    for(z=Zbot; z <= Ztop; ++z)
        for(y=grid.coronal[7]; y <= grid.coronal[14] ; ++y) {
            mask[z*area + y*sz[0] + x] = mask[z*area + y*sz[0] + x - 1] ;
        }
    
  
}
