#ifndef GLOBALS_H
#define GLOBALS_H

#include "griddebug.h"

typedef unsigned char uchar;

enum Value {BACKGROUND, GRIDLINE, HIGHLIGHT,
	    LSUPF=3, LIF, LOBF, LMOBF, LSP, LIP, LO, LAT, LPT,
            LABGT, LPBGT, LMSF, LMIF,
	    RSUPF=16, RIF, ROBF, RMOBF, RSP, RIP, RO, RAT, RPT,
            RABGT, RPBGT, RMSF, RMIF,
	    CEREBELLUM=29,
	    LC_LINE = 100, LSC_MARK, LSF_LINE, LOP_LINE,
	    LPT_LINE, LOT_LINE, LOC_LINE, LATPT_LINE,
	    RC_LINE = 120, RSC_MARK, RSF_LINE, ROP_LINE,
	    RPT_LINE, ROT_LINE, ROC_LINE, RATPT_LINE,
	    M_LINE = 140, LM_LINE, RM_LINE,
	    ORB_LINE=150, BASG_LINE, AP_LINE, AP4_LINE};

int LOBE_OFFSET = 13, LINE_OFFSET = 20;

Grid grid;
int sz[3], area, volume;
int lpron, rpron;



enum Direction { XX, _XX, YY, _YY, ZZ, _ZZ};

class dir {
public:
    dir(unsigned char D, int m, int M): d(D), min(m), max(M) {}
    unsigned char d;
    int min;
    int max;
};

#endif
