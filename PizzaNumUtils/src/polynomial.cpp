
#include "polynomial.h"
#include <gsl/gsl_poly.h>
#include <stdexcept>

using namespace std;
namespace Pizza {
namespace NumUtils {

polynomial_r::polynomial_r(const std::vector<double>& c_)
: c(c_)
{
  if (c.empty())
    throw runtime_error("polynomial: initialized with empty vector");
}

//polynomial_r& operator=(const polynomial_r& p2);

double polynomial_r::operator()(const double x) const
{
  return gsl_poly_eval(&(c[0]), c.size(), x);
}

polynomial_r polynomial_r::deriv() const
{
  vector<double> d(size()-1);
  for (size_t k=0; k<d.size(); k++)
    d[k] = (k+1)*c[k+1];
  return polynomial_r(d);
}

polynomial_r operator+(const polynomial_r& p1, const polynomial_r& p2)
{
  if (p1.size()>p2.size())
    return p2+p1;
  polynomial_r res(p2);
  for (size_t k=0; k<p1.size(); k++) res[k] += p1[k];
  return res;
}

polynomial_r operator-(const polynomial_r& p)
{
  polynomial_r res(p);
  for (size_t k=0; k<res.size(); k++) res[k] = -res[k];
  return res;
}

polynomial_r operator-(const polynomial_r& p1, const polynomial_r& p2)
{
  return p1 + (-p2);
}

void polynomial_r::operator*=(const double a)
{
  for (size_t k=0; k<c.size(); k++) c[k] *= a;
}

polynomial_r operator*(const polynomial_r& p, const double a)
{
  polynomial_r res(p);
  res *=a;
  return res;
}

polynomial_r operator*(const double a, const polynomial_r& p)
{
  polynomial_r res(p);
  res *=a;
  return res;
}

}
}

