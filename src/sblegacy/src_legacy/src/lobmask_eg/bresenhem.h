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

#ifndef BRESENHEM_H
#define BRESENHEM_H

void draw_line(const int &x1, const int &y1,
               const int &x2, const int &y2,
               unsigned char *plane, int *dim,
               const unsigned char &val) {
    /* plane points to the 0,0 position of the rectangle
       dim[0] = x dimension of plane, dim[1] = y dimension of plane
       val is value to put on the line pixels
       */
    
    plane[y1*dim[0] + x1] = plane[y2*dim[0] + x2] = val; /* set endpoints */

    int x = x1, y = y1;
    int dx = x2 - x1, dy = y2 - y1;
    int xchange, ychange;
    
    if( dx < 0 ) {
        xchange = -1;
        dx = - dx;
    }
    else xchange = 1;
    
    if( dy < 0 ) {
        ychange = -1;
        dy = -dy;
    }
    else ychange = 1;

    int error = 0, i = 0, length;

    if( dx < dy ) {
        length = dy;
        while( i < length) {
            y += ychange;
            error += dx;
            if( error > dy ) {
                x += xchange;
                error -= dy;
            }
            plane[y*dim[0] + x] = val;
            ++i;
        }
    }
    else {
        length = dx;
        while( i < length) {
            x += xchange;
            error += dy;
            if(error > dx ) {
                y += ychange;
                error -= dx;
            }
            plane[y*dim[0] + x] = val;
            ++i;
        }
    }
}

#endif
