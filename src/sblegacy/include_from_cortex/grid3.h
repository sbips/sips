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

#include "analyze_utils.h"

inline int midpt(int x, int y) { return (x + y) / 2; } /*same as floor*/

class Grid {
    public :
        Grid(): p(0), pc(0), ac(0), a(0), b(0), ap(0), t(0),
        r(0), m(0), l(0) { }
    
    void print_points() const;
                                                  
    int set(int P, int PC, int AC, int A, int B, int AP, int T,
             int R, int M, int L) {
        p=P; pc=PC; ac=AC; a=A; b=B; ap=AP; t=T; r=R; m=M; l=L; 

        if(!(p < pc && pc < ac && ac < a && b < ap && ap < t && r < m &&
             m < l) ) {
            return 0;
        }
        
        coronal[0] = p;    coronal[4] = pc;
        coronal[7] = ac;  coronal[11] = a;
        axial[0] = b;      axial[4] = ap;     axial[12] = t;
        sagital[0] = r;    sagital[4] = m;    sagital[8] = l;
        
        /* divide area between p and pc into 4 equals */
        coronal[2] = midpt(coronal[0], coronal[4]);
        coronal[1] = midpt(coronal[0], coronal[2]);
        coronal[3] = midpt(coronal[2], coronal[4]);
        
        /* divide area between pc and ac into 3 equal parts */
        int pc_ac_width = (ac - pc + 1)/3;
        coronal[5] = pc +  pc_ac_width;
        coronal[6] = ac - pc_ac_width;
        
        /* divide area between ac and a lines into 4 equal parts */
        coronal[9] = midpt(coronal[7], coronal[11]);
        coronal[8] = midpt(coronal[7], coronal[9]);
        coronal[10] = midpt(coronal[9], coronal[11]);
        
        /* divide area between b and ap lines into 4 equal parts */
        axial[2] = midpt(axial[0], axial[4]);
        axial[1] = midpt(axial[0], axial[2]);
        axial[3] = midpt(axial[2], axial[4]);
        
        /* divide area between ap and t lines into 8 equal parts */
        axial[8] = midpt(axial[4], axial[12]);
        axial[6] = midpt(axial[4], axial[8]);
        axial[5] = midpt(axial[4], axial[6]);
        axial[7] = midpt(axial[6], axial[8]);
        axial[10] = midpt(axial[8], axial[12]);
        axial[9] = midpt(axial[8], axial[10]);
        axial[11] = midpt(axial[10], axial[12]);
        
        /* divide area between r and m into 4 equal parts */
        sagital[2] = midpt(sagital[0], sagital[4]);
        sagital[1] = midpt(sagital[0], sagital[2]);
        sagital[3] = midpt(sagital[2], sagital[4]);
        
        /* divide area between m and l into 4 equal parts */
        sagital[6] = midpt(sagital[4], sagital[8]);
        sagital[5] = midpt(sagital[4], sagital[6]);
        sagital[7] = midpt(sagital[6], sagital[8]);
        
        custom_rm = 68.0/(m - r); custom_lm = 68.0/(l - m);
        custom_aca = 69.5/(a - ac);  custom_acp = 103.0/(ac - p);
        custom_apt = 74.2/(t - ap);  custom_apb = 42.5/(ap - b);

        return 1;
    }
    
    
    inline bool talairach_coord(int *c) {
        if(c[0] >= sagital[0] && c[0] <= sagital[8] &&
           c[1] >= coronal[0] && c[1] <= coronal[11] &&
           c[2] >= axial[0] && c[2] <= axial[11] ) {
            if(c[0] <= m) c[0] = round((m - c[0]) * custom_rm);
            else c[0] = round((m - c[0]) * custom_lm);
            
            if(c[1] >= ac) c[1] = round((c[1] - ac) * custom_aca);
            else c[1] = round((c[1] - ac) * custom_acp);
                
            if(c[2] >= ap) c[2] = round((c[2] - ap) * custom_apt);
            else c[2] = round((c[2] - ap) * custom_apb);
            
            return true;
        }
        else return false;
    }

    inline bool inverse_talairach_coord(int *c) {
            if(c[0] >= 0) c[0] = round(m - c[0]/custom_rm);
            else c[0] = round(m - c[0]/custom_lm);
            
            if(c[1] >= 0) c[1] = round(ac + c[1]/custom_aca);
            else c[1] = round(ac + c[1]/custom_acp);
                
            if(c[2] >= 0) c[2] = round(ap + c[2]/custom_apt);
            else c[2] = round(ap + c[2]/custom_apb);
            
            return true;
    }
    inline int inverse_talairach_coord0(int x) {
            if(x >= 0) return  round(m - x/custom_rm);
            else return round(m - x/custom_lm);
    }
    inline int inverse_talairach_coord1(int y) {
            if(y >= 0) return round(ac + y/custom_aca);
            else return round(ac + y/custom_acp);
    }
    inline int inverse_talairach_coord2(int z) {  
            if(z >= 0) return round(ap + z/custom_apt);
            else return round(ap + z/custom_apb);
    }
    
    int coronal[12], axial[13], sagital[9];
    int p, pc, ac, a; /* special coronal planes */
    int b, ap, t;     /* special axial planes */
    int l, r, m;      /* special sagital planes: left , right, medial */
    
    double custom_rm, custom_lm, custom_aca, custom_acp,
        custom_apt, custom_apb;  
};


/* ////////////////////////////////////////////////////////////////// */

void Grid::print_points() const {
        ostream_iterator<double> out(cout, ", ");
        cout << "coronal grid points: \n";
        copy(coronal, coronal + 20, out);
        cout << endl;
        cout << "axial grid points: \n";
        copy(axial, axial + 13, out);
        cout << endl;
        cout << "sagital grid points: \n";
        copy(sagital, sagital + 9, out);
        cout << endl;
}
