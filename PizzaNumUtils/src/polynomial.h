#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include "functors.h"
#include <vector>
#include <cstring>

namespace Pizza {
namespace NumUtils {

///Class representing a polynomial, can be used also as function
class polynomial_r : public functor_r2r
{
  std::vector<double> c;
  public:
  polynomial_r() : c(1,0) {}
  polynomial_r(const std::vector<double>& c_);
  polynomial_r(const polynomial_r& p2) : c(p2.c) {}

  virtual double operator()(const double x) const;
  size_t size() const {return c.size();}
  size_t degree() const {return size() - 1;}

  const double& operator[](size_t n) const {return c.at(n);}
  double& operator[](size_t n) {return c.at(n);}
  polynomial_r deriv() const;
  void operator*=(const double a);
};

polynomial_r operator+(const polynomial_r& p1, const polynomial_r& p2);
polynomial_r operator-(const polynomial_r& p);
polynomial_r operator-(const polynomial_r& p1, const polynomial_r& p2);
polynomial_r operator*(const polynomial_r& p, const double a);
polynomial_r operator*(const double a, const polynomial_r& p);


}
}

#endif
