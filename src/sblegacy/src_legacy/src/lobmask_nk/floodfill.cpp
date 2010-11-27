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
