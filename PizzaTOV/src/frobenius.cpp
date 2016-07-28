#include "frobenius.h"
#include <cmath>

using namespace std;
using namespace Pizza;
using namespace TOV;

double frobenius_method::root0() const
{
  const double k = (1.0 - a[0]) / 2.0;
  return k - sqrt(k*k-b[0]);
}

double frobenius_method::root1() const
{
  const double k = (1.0 - a[0]) / 2.0;
  return k + sqrt(k*k-b[0]);
}


polynomial_r frobenius_method::first_solution(const double root) const
{
  vector<double> u(min(a.size(),b.size()));
  u[0] = 1.0;
  for (size_t k=1; k<u.size(); k++) {
    double s = 0;
    for (size_t l=0; l<k; l++)
      s += u[l]*(a[k-l]*(l+root) + b[k-l]);
    u[k] = -s / (b[0] + (k+root)*(k+root-1+a[0]));
  }
  return polynomial_r(u);
}



