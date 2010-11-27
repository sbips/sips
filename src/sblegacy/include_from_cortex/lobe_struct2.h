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

#ifndef LOBE_STRUCT_H
#define LOBE_STRUCT_H

#include <iostream.h>
#include <iomanip.h>
#include <algorithm>
#include <numeric>

#define MAX_LESION_VALUE 3

class lobe_struct {
  public:

    lobe_struct(): val(0), total(0) {
            fill(seg_count, seg_count + 256, 0);
	    //            fill(les_count, les_count + MAX_LESION_VALUE, 0);
    }

    void set(int v) {  val = v; }
    
    void print_stats(double vox_vol) const {
        cout << setw(3) << val << endl;
        cout.setf(ios::fixed, ios::floatfield);
        if(total) cout << "Total lobe volume : " << total * vox_vol << endl;
        cout << "Vcsf : " << seg_count[7] * vox_vol << endl;
        cout << "Scsf : " << seg_count[5] * vox_vol << endl;
        cout << "gray : "<< seg_count[4] * vox_vol << endl;
        /* cout << "white : " << (seg_count[0] + seg_count[3]) * vox_vol
             << endl;*/
        cout << "white : " << (seg_count[3] + seg_count[2]) * vox_vol << endl;
	/*   
        for(int i=1; i < MAX_LESION_VALUE; ++i) {
            cout << "lesion " << i << " : "
                 << les_count[i] * vox_vol << endl;
        }
	*/
    }
    
    void print_excel_stats(double vox_vol) const {
        cout << setw(3) << val;
        cout.setf(ios::fixed, ios::floatfield);
        if(total) cout << setw(15) << total * vox_vol;
        else cout << "               ";
        cout << setw(15) << seg_count[7] * vox_vol ;
        cout << setw(15) << seg_count[5] * vox_vol ;
        cout << setw(15) << seg_count[4] * vox_vol ;
        cout << setw(15) << (seg_count[3] + seg_count[2]) * vox_vol;
	/*
        for(int i=1; i < MAX_LESION_VALUE; ++i) {
            cout << setiosflags(ios::right) << setw(15) <<  les_count[i] * vox_vol;
        }
	*/
        cout << endl;
    }
   
    void update_counts(const uchar *seg) {
        if(seg) {
            if(*seg != BACKGROUND) {
                ++total;
                seg_count[*seg] += 1;
            }
        }
	/*
        if( lesion ) {
            if(*lesion) {
                if(*lesion < MAX_LESION_VALUE) {
                    les_count[*lesion] += 1;
                    return;
                }
                else {
                    cout << "found lesion at value >= "
                         << MAX_LESION_VALUE << " Exiting ... " << endl;
                    exit(1);
                }
            }
        }
	*/
    }
    
    int val;
    int total;
    int seg_count[256];
    //    int les_count[MAX_LESION_VALUE];
  

};
    

#endif
    
