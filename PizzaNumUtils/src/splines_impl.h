#include "splines.h"
#include <cstring>
#include <gsl/gsl_spline.h>

namespace Pizza {
namespace NumUtils {

struct wrp_interp_accel {
  gsl_interp_accel *p;
  wrp_interp_accel();
  ~wrp_interp_accel();
};

struct wrp_interp_cspline {
  gsl_interp *p;
  wrp_interp_cspline();
  void init(const std::vector<double>& x, const std::vector<double>& y);
  ~wrp_interp_cspline();
};

///Implementation belonging to interface class \ref cubic_spline
class cubic_spline::impl : public functor_r2r {
  mutable wrp_interp_accel acc;
  wrp_interp_cspline interp;
  std::vector<double> x,y;
  opt_const out_l, out_r;

  public:
  impl(const std::vector<double>& x_, const std::vector<double>& y_,
    opt_const out_l_, opt_const out_r_);
  bool in_range(double t) const {return ((x_min() <= t) && (t <= x_max()));}
  double x_min() const {return x[0];}
  double x_max() const {return x.back();}
  double eval(double t) const;
  double derivative(double t) const;
  double integral(double a, double b) const;
  double operator()(double t) const {return eval(t);}
  size_t size() const {return x.size();}
  double xi(size_t i) const {return x.at(i);}
  double yi(size_t i) const {return y.at(i);}
};

}
}

