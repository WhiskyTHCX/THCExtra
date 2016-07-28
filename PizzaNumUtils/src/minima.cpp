#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <algorithm>
#include <stdexcept>
#include "minima.h"

using namespace std;

namespace Pizza {
namespace NumUtils {

double func_r2r_gsl(double x, void *par)
{
  const function_r2r* obj= (const function_r2r*) par;
  return (*obj)(x);
}

double find_minimum_from_guess(const function_r2r f, double x0, double x1, double m,
	const double abs_acc, const double rel_acc, int max_iter)
{
  if (x0>x1) swap(x0, x1);
  if ((m<=x0) || (m>=x1))
    throw runtime_error("Minimum finding failed: m not in (x0,x1)");
  double fm = f(m);
  if ((f(x0)<=fm) || (f(x1)<=fm))
    throw runtime_error("Minimum finding failed: f(x0)<=f(m) || f(x1)<=f(m)");
  gsl_function F;
  F.function = &func_r2r_gsl;
  F.params   = const_cast<function_r2r*>(&f);

  const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set(s, &F, m, x0, x1);

  int iter = 0;
  do {
    if (iter++ > max_iter)
      throw runtime_error("Minimum finding failed: maximum iterations");
    if (gsl_min_fminimizer_iterate(s) != GSL_SUCCESS)
      throw runtime_error("Minimum finding failed: gsl error");

    x0 = gsl_min_fminimizer_x_lower(s);
    x1 = gsl_min_fminimizer_x_upper(s);

  } while (gsl_min_test_interval(x0, x1, abs_acc, rel_acc) == GSL_CONTINUE);

  m = gsl_min_fminimizer_x_minimum(s);
  gsl_min_fminimizer_free(s);

  return m;
}

double find_maximum_from_guess(const function_r2r f, double x0, double x1, double m,
	const double abs_acc, const double rel_acc, int max_iter)
{
  return find_minimum_from_guess(-f, x0, x1, m, abs_acc, rel_acc, max_iter);
}

double find_minimum(const function_r2r f, double x0, double x1, const int divisions,
	const double abs_acc, const double rel_acc, int max_iter)
{
  if (divisions<=1)
    throw runtime_error("Number of divisions for minimum search must be greater 1");
  int imin(0);
  double fmin(f(x0));
  double dx = (x1-x0) / divisions;
  for (int i=1; i<=divisions; i++) {
    double xi = x0 + dx*i;
    double fi = f(xi);
    if (fi < fmin) {
      imin = i;
      fmin = fi;
    }
  }
  if ((imin==0) || (imin==divisions))
    throw runtime_error("Minimum finding failed: minimum at boundary.");

  double a = x0+(imin-1)*dx;
  double b = x0+(imin+1)*dx;
  double m = x0+imin*dx;
  return find_minimum_from_guess(f, a, b, m, abs_acc, rel_acc, max_iter);
}

double find_maximum(const function_r2r f, double x0, double x1, const int divisions,
	const double abs_acc, const double rel_acc, int max_iter)
{
  return find_minimum(-f, x0, x1, divisions, abs_acc, rel_acc, max_iter);
}

}
}

