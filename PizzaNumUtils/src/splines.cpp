#include <stdexcept>
#include "splines_impl.h"

using namespace std;
using namespace Pizza;
using namespace NumUtils;

wrp_interp_accel::wrp_interp_accel() :
p(gsl_interp_accel_alloc())
{
  if (p==0)
    throw runtime_error("Cubic spline: no memory");
}

wrp_interp_accel::~wrp_interp_accel()
{
  gsl_interp_accel_free(p);
}

wrp_interp_cspline::wrp_interp_cspline() : p(0) {}

void wrp_interp_cspline::init(const vector<double>& x, const vector<double>& y)
{
  if (x.size() != y.size())
    throw invalid_argument("Cubic spline: array size mismatch");
  if (x.size() < 5)
    throw invalid_argument("Cubic spline: too few points");
  p = gsl_interp_alloc(gsl_interp_cspline, x.size());
  if (p==0)
    throw runtime_error("Cubic spline: no memory");
  gsl_interp_init(p, &(x[0]), &(y[0]), x.size());
}

wrp_interp_cspline::~wrp_interp_cspline()
{
  if (p) gsl_interp_free(p);
}


void opt_const::require() const
{
  if (!valid)
    throw logic_error("Optional const value not set");
}

opt_const::operator double() const
{
  require();
  return value;
}


cubic_spline::impl::impl(const vector<double>& x_, const vector<double>& y_,
  opt_const out_l_, opt_const out_r_) :
x(x_),
y(y_),
out_l(out_l_),
out_r(out_r_)
{
  interp.init(x,y);
}

double cubic_spline::impl::eval(double t) const
{
  if (t<x_min()) {
      return out_l;
  }
  if (t>x_max()) {
    return out_r;
  }
  return gsl_interp_eval(interp.p, &(x[0]), &(y[0]), t, acc.p);
}

double cubic_spline::impl::derivative(double t) const
{
  if (t<x_min()) {
    out_l.require();
    return 0;
  }
  if (t>x_max()) {
    out_r.require();
    return 0;
  }
  return gsl_interp_eval_deriv(interp.p, &(x[0]), &(y[0]), t, acc.p);
}

double cubic_spline::impl::integral(double a, double b) const
{
  if (b < a) {
    return -integral(b,a);
  }

  double res  = 0;
  double ia   = a;
  double ib   = b;

  if (a<x_min()) {
    ia = x_min();
    if (b<x_min()) res += (b-a) * out_l;
    else res += (x_min()-a) * out_l;
  }
  if (b>x_max()) {
    ib  = x_max();
    if (a>x_max()) res += (b-a) * out_r;
    else res += (b-x_max()) * out_r;
  }
  if (ia<ib) {
    res += gsl_interp_eval_integ(interp.p, &(x[0]), &(y[0]), ia, ib, acc.p);
  }

  return res;
}


cubic_spline::cubic_spline(const std::vector<double>& x, const std::vector<double>& y,
    opt_const out_left, opt_const out_right)
: pimpl(new impl(x, y, out_left, out_right)) {}

cubic_spline::cubic_spline(const cubic_spline& that)
: pimpl(that.pimpl) {}

const cubic_spline::impl& cubic_spline::s() const
{
  if (pimpl) return *pimpl;
  throw logic_error("Cubic spline: uninitialized use.");
}

bool cubic_spline::in_range(double t) const {return s().in_range(t);}
double cubic_spline::x_min() const {return s().x_min();}
double cubic_spline::x_max() const {return s().x_max();}
double cubic_spline::operator()(double x) const {return s().eval(x);}
double cubic_spline::derivative(double x) const {return s().derivative(x);}
double cubic_spline::integral(double a, double b) const {return s().integral(a,b);}
size_t cubic_spline::size() const {return s().size();}
double cubic_spline::xi(size_t i) const {return s().xi(i);}
double cubic_spline::yi(size_t i) const {return s().yi(i);}


cubic_spline::operator function_r2r () const
{
  return function_r2r(pimpl);
}

