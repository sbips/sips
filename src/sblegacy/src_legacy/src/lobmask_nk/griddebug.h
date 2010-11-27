#ifndef GRID_H
#define GRID_H

#include <iostream.h>
#include <algorithm>

inline int midpt(int x, int y) { return (x + y) / 2; } /*same as floor*/

class Grid {
    public :
        Grid(): p(0), pc(0), ac(0), a(0), b(0), ap(0), t(0),
        r(0), m(0), l(0) { }
                                                         
        void set(int P, int PC, int AC, int A, int B, int AP, int T,
                 int R, int M, int L);
        
        void print_points() const;
    
        int coronal[20], axial[13], sagital[9];
        int p, pc, ac, a; /* special coronal planes */
        int b, ap, t;     /* special axial planes */
        int l, r, m;      /* special sagital planes: left , right, medial */
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
/* //////////////////////////////////////////////////////////////// */

void Grid::set(int P, int PC, int AC, int A, int B, int AP, int T,
               int R, int M, int L) {
    p=P; pc=PC; ac=AC; a=A; b=B; ap=AP; t=T; r=R; m=M; l=L; 
    //r=20;
    //l=160;
   //a=195;
   //p=16;
   //t=173;
    if(!(p < pc && pc < ac && ac < a && b < ap && ap < t && r < m &&
         m < l) ) {
        cout << "Arguments passed to Grid set function"
             << " do not make sense. Exiting..." << endl;
        exit(1);
    }
    std::cout << "new grid" << std::endl;
    coronal[0] = p;    coronal[8] = pc;
    coronal[11] = ac;  coronal[19] = a;
    axial[0] = b;      axial[4] = ap;     axial[12] = t;
    sagital[0] = r;    sagital[4] = m;    sagital[8] = l;
    
    /* divide area between p and pc into 8 equals */
    coronal[4] = midpt(coronal[0], coronal[8]);
    coronal[2] = midpt(coronal[0], coronal[4]);
    coronal[1] = midpt(coronal[0], coronal[2]);
    coronal[3] = midpt(coronal[2], coronal[4]);
    coronal[6] = midpt(coronal[4], coronal[8]);
    coronal[5] = midpt(coronal[4], coronal[6]);
    coronal[7] = midpt(coronal[6], coronal[8]);
    
    /* divide area between pc and ac into 3 equal parts */
    int pc_ac_width = (ac - pc + 1)/3;
    coronal[9] = pc +  pc_ac_width;
    coronal[10] = ac - pc_ac_width;
    
    /* divide area between ac and a lines into 8 equal parts */
    coronal[15] = midpt(coronal[11], coronal[19]);
    coronal[13] = midpt(coronal[11], coronal[15]);
    coronal[12] = midpt(coronal[11], coronal[13]);
    coronal[14] = midpt(coronal[13], coronal[15]);
    coronal[17] = midpt(coronal[15], coronal[19]);
    coronal[16] = midpt(coronal[15], coronal[17]);
    coronal[18] = midpt(coronal[17], coronal[19]);
    
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
}
  
    
#endif
