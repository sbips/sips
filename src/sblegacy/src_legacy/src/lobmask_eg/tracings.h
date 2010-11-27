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

#ifndef TRACINGS_H
#define TRACINGS_H

#include "globals.h"

using namespace nk;

typedef vector<uchar> ROIValues;
typedef vector<string> ROINames;

void load_tracings(char * objfname, 
		   uchar *&left_trace_mask, 
		   uchar *&right_trace_mask);
// load object map and reassign tracings new values

void connect_trancings(uchar *trace_mask, uchar offset, int pron);
// connect tarcings and draw horizonatal lines

void fill_trace_mask(uchar *trace_mask, uchar line_offset, uchar lobe_offset);
// floodfill all areas

void recolor_lines(uchar *trace_mask, int line_offset, int lobe_offset);

void connect_C(uchar trace_val, uchar mark_val, uchar *plane);

void connect_OP(uchar trace_val, uchar *plane,
		const int &y1, int &z1);
// on entrance y1 is pc-4
// on return (y1,z1) is the lower endpt of OP_LINE

void connect_SF(uchar trace_val, uchar *plane,
		int &y2, int &z2);
// on return (y2,z2) is the upper endpt of SF_LINE

void draw_horizontal_line(uchar *plane, int z, vector <uchar> stop_values,
			   uchar val);
// draw horizonatal line at level z starting from right side until
// one of the stop values is encountered
// give value val to the horizontal line

void draw_ATPT_line(uchar *plane, uchar stop_val, uchar val);
// draw line that divides anterior and posterior temporal lobes
#endif
