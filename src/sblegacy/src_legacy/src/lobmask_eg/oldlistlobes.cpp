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

#include <iostream.h>

int main() {
    cout  << "
0       BACKGROUND
1       GRIDLINE
2       HIGHLIGHT

3       LSUPF
4       LIF
5       LOBF
6       LMOBF
7       LSP
8       LIP
9       LO
10      LT
11      LABGT
12      LPBGT
13      LMSF
14      LMIF

15      RSUPF
16      RIF
17      ROBF
18      RMOBF
19      RSP
20      RIP
21      RO
22      RT
23      RABGT
24      RPBGT
25      RMSF
26      RMIF

"
          << endl;
    return 0;
}
