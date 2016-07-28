#ifndef SPLINES_H
#define SPLINES_H

#include "functors.h"
#include <vector>
#include <cstring>

namespace Pizza {
namespace NumUtils {

///An optional real-valued constant
class opt_const {
  bool valid;
  double value;
  public:
  enum none_t {NONE};

  opt_const(double value_) : valid(true), value(value_) {}
  opt_const(none_t dummy) : valid(false) {}
  void require() const;
  operator double() const;
};

///A cubic spline.
class cubic_spline {
  class impl;
  boost::shared_ptr<impl> pimpl;
  const impl& s() const;

  public:

  cubic_spline() {}
  cubic_spline(const std::vector<double>& x, const std::vector<double>& y,
    opt_const out_left, opt_const out_right);
  cubic_spline(const cubic_spline& that);

  bool in_range(double t) const;
  double x_min() const;
  double x_max() const;
  double operator()(double x) const;
  double derivative(double x) const;
  double integral(double a, double b) const;
  operator function_r2r () const;
  size_t size() const;
  double xi(size_t i) const;
  double yi(size_t i) const;
};

}
}

#endif
