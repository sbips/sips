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
