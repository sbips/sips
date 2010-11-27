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

#ifndef FLOODFILL_H
#define FLOODFILL_H

#include "globals.h"

class Floodfill {
  public:
    Floodfill(uchar *seed, int *dirs, int numdirs,
              bool (*acceptable)(const uchar *), uchar val);
    void perform();
    int forward(uchar *pt);
    
  private:
    uchar *_seed;
    int *_dirs;
    int _numdirs;
    bool (*_acceptable)(const uchar *pt);
    uchar _val;
};

/////////////////////////////////////////////////////////////////

inline bool acceptable(const uchar *pt){
    if(*pt == BACKGROUND) return true;
    else return false;
}

inline bool acceptable_OBF(const uchar *pt){
    if(*pt != ORB_LINE && *pt != GRIDLINE && *pt != ROBF && *pt !=M_LINE && *pt != LOBF) return true;
    else return false;
}

inline bool acceptable_BGT(const uchar *pt){
    if(*pt != BASG_LINE && *pt != M_LINE && *pt != LC_LINE && *pt != RC_LINE &&
       *pt != LABGT && *pt != RABGT && *pt != LPBGT && *pt != RPBGT && *pt != GRIDLINE) return true;
    else return false;
}


    
#endif
