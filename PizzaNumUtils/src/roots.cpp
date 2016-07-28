#include "roots.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <stdexcept>

using namespace std;

namespace Pizza {
namespace NumUtils {

double func_1d_gsl(double x, void *par)
{
  const function_r2r* obj= (const function_r2r*) par;
  return (*obj)(x);
}

double findroot(const function_r2r f, double x0, double x1,
	const double abs_acc, const double rel_acc, int max_iter)
{
  if (f(x0)*f(x1) > 0) {
    throw runtime_error("Root finding failed: bad initial interval");
  }
  gsl_function F;
  F.function                = &func_1d_gsl;
  F.params                  = const_cast<function_r2r*>(&f);
  gsl_root_fsolver *solver  = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(solver, &F, x0, x1);
  while (GSL_SUCCESS != gsl_root_test_interval(x0,x1, abs_acc, rel_acc)) {
    int status  = gsl_root_fsolver_iterate (solver);
    x0          = gsl_root_fsolver_x_lower (solver);
    x1          = gsl_root_fsolver_x_upper (solver);
    if (status != GSL_SUCCESS)
      throw runtime_error("Root finding failed: gsl error");
    if (0 == --max_iter)
      throw runtime_error("Root finding failed: maximum iterations");
  }
  const double r  = gsl_root_fsolver_root(solver);
  gsl_root_fsolver_free(solver);
  return r;
}

}
}
