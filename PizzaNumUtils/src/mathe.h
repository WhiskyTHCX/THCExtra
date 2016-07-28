#ifndef MATHE_H
#define MATHE_H mathe_h

#include "cctk_Types.h"
//#include <math.h>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <algorithm>

namespace Pizza {
typedef CCTK_REAL pz_real;
typedef std::complex<pz_real> pz_cmplx;

template <class T> inline T sqr(T t) {return t*t;}
}
#endif
