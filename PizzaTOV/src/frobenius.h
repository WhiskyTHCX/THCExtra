#ifndef FROBENIUS_H
#define FROBENIUS_H

#include "polynomial.h"


namespace Pizza {
namespace TOV {

using NumUtils::polynomial_r;

///Class for solving singular ODE using Frobenius' method
class frobenius_method
{
  polynomial_r a,b;
  public:
  frobenius_method(const polynomial_r& a_, const polynomial_r& b_) : a(a_), b(b_) {}
  double root0() const;
  double root1() const;
  polynomial_r first_solution(const double root) const;
};

}
}
#endif

