#include "odes.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <stdexcept>
#include <cmath>

using namespace std;
using namespace Pizza;
using namespace NumUtils;


struct p_gsl_odeiv_step {
  gsl_odeiv_step* p;
  p_gsl_odeiv_step(const gsl_odeiv_step_type* typ, int size) {
    p = gsl_odeiv_step_alloc (typ, size);
    if (p==0) throw runtime_error("ODE solver: no memory for GSL structure");
  }
  ~p_gsl_odeiv_step() { gsl_odeiv_step_free(p);}
};

struct p_gsl_odeiv_control {
  gsl_odeiv_control* p;
  p_gsl_odeiv_control(double acc_abs, double acc_rel, const vector<double>& scales) {
    p = gsl_odeiv_control_scaled_new (acc_abs, acc_rel, 1.0, 0.0, &(scales[0]), scales.size());
    if (p==0) throw runtime_error("ODE solver: no memory for GSL structure");
  }
  ~p_gsl_odeiv_control() {gsl_odeiv_control_free(p);}
};

struct p_gsl_odeiv_evolve {
  gsl_odeiv_evolve* p;
  p_gsl_odeiv_evolve(int size) {
    p = gsl_odeiv_evolve_alloc (size);
    if (p==0) throw runtime_error("ODE solver: no memory for GSL structure");
  }
  ~p_gsl_odeiv_evolve() {gsl_odeiv_evolve_free(p);}
};


int gsl_func(double x, const double y[], double dy[], void *par)
{
  const ode_sys_r* ode= (const ode_sys_r*) par;
  if (ode->stop(x, y))
    return GSL_FAILURE;
  ode->deriv(x, y, dy);
  return GSL_SUCCESS;
}

void ode_sys_r::integrate(double acc, std::vector<double>& res_x,
  std::vector<std::vector<double> >&res_y, size_t max_steps) const
{
  vector<double> y_err(size, 0), y(size, 0), w_abs(size, 0);
  double x(0), w_rel(0), w_xrel(0), w_xabs(0), dx(0);

  err_weights(&(w_abs[0]), w_rel, w_xabs, w_xrel);
  const double acc_abs(acc), acc_rel(acc*w_rel), acc_xabs(acc*w_xabs), acc_xrel(acc*w_xrel);

  init(x, &(y[0]), dx);

  const gsl_odeiv_step_type* styp = gsl_odeiv_step_rkf45;
  p_gsl_odeiv_step  stf(styp, size);
  p_gsl_odeiv_control ctrl(acc_abs, acc_rel, w_abs);
  p_gsl_odeiv_evolve evlv(size);
  gsl_odeiv_system sys            = {&gsl_func, NULL, size , const_cast<ode_sys_r*>(this)};


  res_y.resize(size);
  res_x.push_back(x);
  for (size_t i=0; i<size; i++)
    res_y[i].push_back(y[i]);

  for (; true; max_steps--) {
    if (0 == max_steps)
      throw runtime_error("ODE solver: maximum number of iterations reached");
    int status = gsl_odeiv_evolve_apply(evlv.p, ctrl.p, stf.p, &sys, &x, x+100*dx, &dx, &(y[0]));
    if (status != GSL_SUCCESS) break;

    res_x.push_back(x);
    for (size_t i=0; i<size; i++)
      res_y[i].push_back(y[i]);
  }
  for (; dx > acc_xabs + acc_xrel * abs(x); max_steps--) {
    if (0 == max_steps)
      throw runtime_error("ODE solver: maximum number of iterations reached");
    int status = gsl_odeiv_step_apply(stf.p, x, dx, &(y[0]), &(y_err[0]), NULL, NULL, &sys);
    if (status != GSL_SUCCESS) {
      dx /= 10;
      continue;
    }
    x += dx;
    res_x.push_back(x);
    for (size_t i=0; i<size; i++)
      res_y[i].push_back(y[i]);
  }
  if (res_x.size()<5)
    throw runtime_error("ODE solver: something went wrong.");
}

