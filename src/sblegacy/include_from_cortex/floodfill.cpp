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

#include "floodfill.h"
#include <stack>
#include <algorithm>

Floodfill::Floodfill(uchar *seed, int *dirs, int numdirs,
                     bool(*acceptable)(const uchar *), uchar val) :
    _seed(seed), _dirs(dirs), _numdirs(numdirs),
    _acceptable(acceptable), _val(val) {
    perform ();
}
/////////////////////////////////////////////////////////////////////////

void Floodfill::perform() {
  stack<uchar *> Stack;
  Stack.push(_seed);
  uchar *pt;
  int t;
  while( !Stack.empty() ) {
    pt = Stack.top();
    *pt = _val;
    t = forward(pt);
    if( t ) Stack.push(pt + t);
    else Stack.pop();
  }
}
/////////////////////////////////////////////////////////////////////////////

int Floodfill::forward(uchar *pt) {
  for(int i=0; i < _numdirs; ++i) {
      if( _acceptable(pt + _dirs[i]) ) return (_dirs[i]);
  }
  return 0;
}
