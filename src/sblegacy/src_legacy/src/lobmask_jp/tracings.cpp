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

#include "tracings.h"
#include "grid.h"
#include "bresenhem.h"
#include "floodfill.h"
#include "/net/cortex/csf/pub/include/conv.h"
#include <iostream.h>
#include <string.h>
#include <vector.h>
#include <algorithm>

void load_tracings(char *objfname,
		   uchar *&left_trace_mask,
		   uchar *&right_trace_mask) {
  int num_of_tracings = 8;
  // set up expected trace names
  char *ntemp[] = { "lc", "lsc", "lsf", "lop",
		   "rc", "rsc", "rsf", "rop" };
  ROINames trace_names(ntemp, ntemp + num_of_tracings);
  // set up new trace values
  Value vtemp[] = { LC_LINE, LSC_MARK, LSF_LINE, LOP_LINE,
		    RC_LINE, RSC_MARK, RSF_LINE, ROP_LINE };
  ROIValues trace_values(vtemp, vtemp + num_of_tracings);

  ifstream obj_file(objfname);
  if(!obj_file) {
    cout << "Can't find object file " << objfname << ". Exiting..." << endl;
    exit(1);
  }
  int numobj;
  int dim[3];
  bool swap = false;
  obj_file.seekg(sizeof(int));
  obj_file.read(&dim, 3*sizeof(int));


  if(!(sz[0] == dim[0] && sz[1] == dim[1] && sz[2] == dim[2])){
      swap = true;	// try byte-swapping
      for (int temp = 0; temp < 3; temp++)
            dim[temp] = iconv(dim[temp]);
	  
      if(!(sz[0] == dim[0] && sz[1] == dim[1] && sz[2] == dim[2])){
            swap = false;
            cout << "Object map created on a file of different dimensions."
                 << " Exiting... " << endl;
            exit(1);
      }
  }
  obj_file.read(&numobj, sizeof(int));
  if (swap)
  	numobj = iconv(numobj);

  // read in actual object names and values from given object file
  ROINames roi_names;
  ROIValues roi_values;
  char name[NAMELENGTH];
  string objname;
  int loc;
  for(int val=0; val < numobj; ++val) {
    obj_file.read(name, NAMELENGTH*sizeof(char));
    objname = string(name);
    // doing this ensures that empty characters on the right are not included
    loc = objname.find(".");
    objname.assign(objname, loc + 1, objname.length() - loc - 1);
    // get rid of all stuff preceeding .
    // e.g. "10. l c" becomes "l c"
    roi_names.push_back(objname);
    roi_values.push_back(val);
    //cout << objname << "\t\t" << val << endl;
    obj_file.seekg(HEADERLENGTH, ios::cur);  // skip header
  }
  //cout << "roi_values.size() " << roi_values.size() << endl;
  // make sure that all expected trace names are present
  ROINames::iterator it = trace_names.begin();
  while(it != trace_names.end()) {
    if( find(roi_names.begin(), roi_names.end(), *it) == roi_names.end() ){
      cout << "Missing tracing : " << *it << ". Exiting..." << endl;
      exit(1);
    }
    ++it;
  }

  // load tracings as flat sagital masks
  int sagital_area = sz[1] * sz[2];
  left_trace_mask = new uchar [sagital_area];
  right_trace_mask = new uchar [sagital_area];
  if(!left_trace_mask || !right_trace_mask) {
    cout << "Failed to allocate memory for trace masks. Exiting..." << endl;
    exit(1);
  }
  fill(left_trace_mask, left_trace_mask + sagital_area, BACKGROUND);
  fill(right_trace_mask, right_trace_mask + sagital_area, BACKGROUND);

  uchar num, val, new_val, *trace_mask;
  int prev, total = 0, y, z, i, j;
  while(total < volume) {
    prev = total;
    obj_file.read(&num, sizeof(char));
    obj_file.read(&val, sizeof(char));

    if( val > 0 && val < roi_names.size() ){
        j = find(trace_names.begin(), trace_names.end(), roi_names[val]) -
         trace_names.begin();

        if( j < num_of_tracings ) {
            // object is one of the expected ones
            new_val = trace_values[j];
            if( j < num_of_tracings/2 )  trace_mask = left_trace_mask;
            // if val corresp to left -> load into left mask and give it new value
            else trace_mask = right_trace_mask;
            while( total < prev + (int)num ) {
                trace_mask[ total / sz[0] ] = new_val;
                ++total;
            }
        }
        else total += (int)num;
    }
    else total += (int)num; // skip irrelevant objects
  }
  // now all tracings have been loaded and reassigned new values
}
//////////////////////////////////////////////////////////////////////

void connect_trancings(uchar *trace_mask, uchar line_offset, int pron){
    int y1 = grid.coronal[4];
    int op_lowend_y = y1, op_lowend_z;
    connect_OP(LOP_LINE + line_offset, trace_mask, op_lowend_y, op_lowend_z);
    connect_C(LC_LINE + line_offset, LSC_MARK + line_offset, trace_mask);
    int sf_upend_y, sf_upend_z;
    connect_SF(LSF_LINE + line_offset, trace_mask, sf_upend_y, sf_upend_z);

    // connect end of SF to end of OP -> PT -line
    draw_line( op_lowend_y, op_lowend_z, sf_upend_y, sf_upend_z,
               trace_mask, sz+1, LPT_LINE + line_offset);
    // connect end of OP to pron -> OT line
    draw_line( op_lowend_y, op_lowend_z, y1, pron,
               trace_mask, sz+1, LOT_LINE + line_offset);
    // connect pron to (p, ap-3) and extend to nearest corner -> OC line
    draw_line(y1, pron, grid.coronal[0], grid.axial[1], trace_mask, sz+1,
              LOC_LINE + line_offset);
    draw_line(grid.coronal[0], grid.axial[1], 1, 1, trace_mask, sz+1,
              LOC_LINE + line_offset);

    // now draw horizontal lines at levels ap and ap+4 in both masks
    // starting from sz[1]-1
    vector<uchar> stop_values;
    stop_values.push_back(LOP_LINE + line_offset);
    draw_horizontal_line(trace_mask, grid.axial[8], stop_values, AP4_LINE);
    stop_values.erase(stop_values.begin());
    stop_values.push_back(LSF_LINE + line_offset);
    draw_horizontal_line(trace_mask, grid.axial[4], stop_values, AP_LINE);

    // divide temporal into AT and PT
    draw_ATPT_line(trace_mask, LSF_LINE + line_offset, LATPT_LINE + line_offset);
}
/////////////////////////////////////////////////////////////////////////////

void fill_trace_mask(uchar *trace_mask, uchar line_offset, uchar lobe_offset) {
  // prepare for floodfill by setting all borders to GRIDLINE
  int z, z1 = 0, z2 = sz[2]-1, y, y1 = 0, y2 = sz[1]-1;
  for(y=0; y < sz[1]; ++y)
    trace_mask[z1*sz[1] + y] = trace_mask[z2*sz[1] + y] = GRIDLINE;
  for(z=0; z < sz[2]; ++z)
    trace_mask[z*sz[1] + y1] = trace_mask[z*sz[1] + y2] = GRIDLINE;

  int dirs[] = {1, sz[1], -1, -sz[1]};
  // fill O, SP, IP,SUPF, IF, T
  Floodfill o_f(trace_mask + grid.ap*sz[1] + 1,
                dirs, 4, acceptable, LO + lobe_offset);
  Floodfill sp_f(trace_mask + (sz[2]-2)*sz[1] + grid.p,
                 dirs, 4, acceptable, LSP + lobe_offset);
  Floodfill ip_f(trace_mask + (grid.axial[8]-1)*sz[1] + grid.coronal[4],
                 dirs, 4,acceptable, LIP + lobe_offset);
  Floodfill if_f(trace_mask + grid.axial[6]*sz[1] + sz[1]-2,
                 dirs, 4, acceptable, LIF + lobe_offset);
  Floodfill supf_f(trace_mask + (sz[2]-2)*sz[1] + sz[1]-2,
		   dirs, 4, acceptable, LSUPF + lobe_offset);
  Floodfill pt_f(trace_mask + sz[1] + grid.p,
                  dirs, 4, acceptable, LPT + lobe_offset);
  // because of possibilty that C_LINE might  prevent AT filling
  // in all the way (crossing between C and SF happens below grid.ap)
  // fill AT lobe without floodfill
  int t;
  for(z=1; z <= grid.ap; ++z)
      for(y=grid.coronal[6]; y < sz[1]-1; ++y) {
          t = z*sz[1] + y;
          if(trace_mask[t] == BACKGROUND) trace_mask[t] = LAT + lobe_offset;
      }

}
////////////////////////////////////////////////////////////////////////

void recolor_lines(uchar *trace_mask, int line_offset, int lobe_offset) {
    // recolor all lines except LC and RC
  uchar *pt = trace_mask;
  int y, z;
  for(z=0; z < sz[2]; ++z)
      for(y=0; y < sz[1]; ++y){
          switch( *pt - line_offset ) {
          case(LOP_LINE):
              if( z > grid.axial[8] ) *pt = LSP + lobe_offset;
              else *pt = *(pt - 3*sz[0]);
              break;
          case(LOC_LINE):
          case(LOT_LINE): *pt = LO + lobe_offset; break;
          case(LPT_LINE): *pt = LIP + lobe_offset; break;
          case(LATPT_LINE): *pt = LPT + lobe_offset;
          default : break;
          }
          switch(*pt) {
          case GRIDLINE:
          case BACKGROUND:        break;
          case AP_LINE:
          case AP4_LINE: *pt = *(pt - sz[1]); break;
          default: break;
          }
          ++pt;
      }
}
/////////////////////////////////////////////////////////////////////////////

void connect_C(uchar trace_val, uchar mark_val, uchar *plane) {
  // find upper endpoint ->(y1,z1)
  int y1, z1;
  bool found = false;
  for(z1=sz[2]-1; z1 >=0 && !found; --z1)
      for(y1=0; y1 < sz[1] && !found; ++y1)
          if(plane[z1*sz[1] + y1] == trace_val) found = true;

  // now find center of SC mark -> (y2,z2)
  // as you do it, erase points belonging to mark, because we only need
  // its center
  int y, z, y2, z2;
  int ymin = sz[1]-1, ymax = 0, zmin = sz[2]-1, zmax = 0;
  for(z=sz[2]-1; z >= 0; --z)
      for(y=0; y < sz[1]; y++) {
          if(plane[z*sz[1] + y] == mark_val) {
              if(y < ymin) ymin = y;
              if(y > ymax) ymax =y;
              if(z < zmin) zmin =z;
              if(z > zmax) zmax = z;
              plane[z*sz[1] + y] = BACKGROUND;
          }
      }

  if( ymax !=  0 ) { 
      y2 = midpt(ymin, ymax);
      z2 = midpt(zmin, zmax);
  }
  else{
      cout << "Couldn't find mark value " << (int)mark_val
           << " . Extending endpoint of central sulcus vertically." << endl;
      y2 = y1; z2 = z1;
  }
  // connect (y1, z1) to (y2, z2) and then extend vertically
  draw_line(y1, z1, y2, z2, plane, sz + 1, trace_val);
  draw_line(y2, z2, y2, sz[2]-2, plane, sz + 1, trace_val);

  // reconnect C tracing so that if SF intersects C
  // and actual point that belongs to both is denoted as SF then
  // we reconnect two disjoint ends of C

// first find highest C point, going from bottom ->(y11, z11)
  int y11, z11;
  for(z=0; z < sz[2]-1; ++z) {
      for(y=0; y < sz[1]; ++y) {
          if(plane[z*sz[1] + y] == trace_val) {
              y11 = y;  z11 = z;
              break;
          }
      }
  }
  
  if(z11 <  z1) { // then there is an interuption in C tracing at level z11+1
      // find next lowest point ->(y22, z22) and reconnect
      int y22, z22;
      found = false;
      for(z22=z11+1; z22 < z1 && !found; ++z22)
          for(y22=0; y22 < sz[1] && !found; ++y22) 
              if(plane[z22*sz[1] + y22] == trace_val) found = true;
      
      draw_line(y11, z11, y22, z22, plane, sz+1, trace_val);
  }

  // Now find the lowest point of C and if neccessary extend verically
  // to the bottom level of basal_ganglia
  int y3, z3;
  found = false;
  for(z3=0; z3 < z1 && !found; ++z3) 
      for(y3=0; y3 < sz[1] && !found; ++y3) 
          if(plane[z3*sz[1] + y3] == trace_val) found = true;
             
  
  int bgt_bot = midpt(grid.axial[3], grid.axial[4]);
  if(z3 > bgt_bot)
     draw_line(y3, z3, y3, bgt_bot, plane, sz+1, trace_val);
}
///////////////////////////////////////////////////////////////////////////
 
void connect_OP(uchar trace_val, uchar *plane,
		const int &y1, int &z1) {
  // on entrance y1 should have value pc - 4
  // find first endpt of OP tracing that lies on pc - 4 = y1 -> (y1, z1)
  for(z1=0; z1 < sz[2] ; ++z1)
    if( plane[z1*sz[1] + y1] == trace_val) break;
  if( z1 > sz[2]-3 ) {
    cout << "OP tracing of value " << trace_val << " not crossing pc - 4 line."
	 << " Exiting..." << endl;
    exit(1);
  }
  // find upper endpoint ->(y2,z2) and connect to the nearest corner
  int y2,z2;
  bool found = false;
  for(y2=0; y2 < sz[1] && !found; ++y2) {
    for(z2=sz[2]-1; z2 >= 0 && !found; --z2) 
      if( plane[z2*sz[1] + y2] == trace_val )found = true;
  }
  if(!found) {
      cout << "Couldn't find OP tracing valued at " << trace_val
	     << ". Exiting..." << endl;
      exit(1);
  }
  draw_line(y2, z2, 1, sz[2]-2, plane, sz + 1, trace_val);
}
////////////////////////////////////////////////////////////////////////

void connect_SF(uchar trace_val, uchar *plane,
		int &y2, int &z2) {
  // find lower endpt (y1,z1) and make sure that its z coord is lower than
  // the ap level (this is neccerssary for T lobe)
  int y1, z1;
  bool found = false;
  for(y1=sz[1]-1; y1 >= 0 && !found; --y1) {
    for(z1=0; z1 < sz[2] && !found; ++z1)
      if( plane[z1*sz[1] + y1] == trace_val )found = true;
  }
  if(!found) {
    cout << "Couldn't find SF tracing valued at " << trace_val
	 << ". Exiting ..." << endl;
    exit(1);
  }
  // now find the other endpt of SF -> (y2,z2)
  found = false;
  for(y2=0; y2 < sz[1] && !found; ++y2) {
  for(z2=sz[2]-1; z2 >= 0 && !found; --z2) 
      if( plane[z2*sz[1] + y2] == trace_val ) found = true;
  }
}
//////////////////////////////////////////////////////////////////////////
  
void draw_horizontal_line(uchar *plane, int z, vector<uchar> stop_val, uchar val) {
  uchar *ptr = plane + z*sz[1] + sz[1] - 1;
  while( find(stop_val.begin(), stop_val.end(), *ptr) ==
	 stop_val.end() && ptr > plane)  *ptr--= val;
}

 //////////////////////////////////////////////////////////////////////////
void draw_ATPT_line(uchar *plane, uchar stop_val, uchar val) {
  int y1,y2,z1,z2, y3, z3;
  y1 = midpt(grid.coronal[9], grid.coronal[10]);
  z1 = grid.axial[1];
  // y2 = grid.ac;
  // z2 = grid.axial[3];
  y2 = grid.coronal[12];
  z2 = grid.ap;
  float slope = (y2 - y1)/((float)(z2 - z1));

  // find start point on bottom level
  int ybeg, zbeg;
  zbeg = grid.b;
  ybeg = round(y1 + slope * (zbeg - z1));

  // staring from (ybeg, zbeg) go up until SF is hit
  // -> get coordinates of this point ->(y3,z3)
  int y, z, t;
  bool found = false;
  for(z=grid.b; z < grid.axial[6] && !found; ++z) {
      y = round(y1 + slope * (z - z1));
      t = z * sz[1] + y;
      if(plane[z*sz[1] + y] != stop_val) { y3 = y;  z3 = z;}
      else found = true;
      if(plane[(z+1)*sz[1] + y] == stop_val ||
         plane[z*sz[1] + y+1] == stop_val) {
          found = true;
      }
  }
  if (!found) {
      cout << "Warning: SF tracing has holes?" << endl;
  }
  draw_line(ybeg, zbeg, y3, z3, plane, sz + 1, val);

  // draw vertical line from (ybeg,zbef) down
  draw_line(ybeg, zbeg, ybeg, 1, plane, sz + 1, val);
}
   
