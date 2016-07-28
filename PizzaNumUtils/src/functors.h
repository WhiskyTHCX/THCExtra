#ifndef FUNCTORS_H
#define FUNCTORS_H
#include <boost/shared_ptr.hpp>

namespace Pizza {
namespace NumUtils {

///Abstract interface describing a real valued function
class functor_r2r {
  public:
  virtual ~functor_r2r() {}
  virtual double operator()(double x) const=0;
};

///Abstract interface describing a real valued binary function
class functor_rr2r {
  public:
  virtual ~functor_rr2r() {}
  virtual double operator()(double x, double y) const=0;
};


class functor_const : public functor_r2r {
  const double c;
  public:
  explicit functor_const(double c_) : c(c_) {}
  virtual double operator()(double x) const {return c;}
};

double operator_neg(double);
double operator_inv(double);
double operator_sum(double, double);
double operator_difference(double, double);
double operator_product(double, double);
double operator_division(double, double);


///A real-valued function object that can be copied by value.
class function_r2r {
  boost::shared_ptr<functor_r2r> f;
  public:
  function_r2r();
  function_r2r(const function_r2r& other);
  function_r2r(functor_r2r* f_);
  function_r2r(boost::shared_ptr<functor_r2r> f_);
  function_r2r(double (*fp)(double));
  function_r2r& operator=(const function_r2r& other);
  double operator()(double x) const;
  function_r2r operator()(const function_r2r& g) const;
};

class function_rr2r {
  boost::shared_ptr<functor_rr2r> f;
  public:
  function_rr2r();
  function_rr2r(const function_rr2r& other);
  function_rr2r(functor_rr2r* f_);
  function_rr2r(boost::shared_ptr<functor_rr2r> f_);
  function_rr2r(double (*fp)(double, double));
  function_rr2r& operator=(const function_rr2r& other);
  double operator()(double x, double y) const;
  function_r2r bind_1st(double x) const;
  function_r2r bind_2nd(double y) const;
};

function_r2r chain(const function_r2r f, const function_r2r g);
function_r2r apply_binary(const function_rr2r, const function_r2r,
                          const function_r2r);

function_r2r operator+(const function_r2r&, const function_r2r&);
function_r2r operator+(const function_r2r&, double);
function_r2r operator+(double, const function_r2r&);
function_r2r operator-(const function_r2r&, const function_r2r&);
function_r2r operator-(const function_r2r&, double);
function_r2r operator-(double, const function_r2r&);
function_r2r operator-(const function_r2r&);
function_r2r operator*(const function_r2r&, const function_r2r&);
function_r2r operator*(const function_r2r&, double);
function_r2r operator*(double, const function_r2r&);
function_r2r operator/(const function_r2r&, const function_r2r&);
function_r2r operator/(const function_r2r&, double);
function_r2r operator/(double, const function_r2r&);

}
}

#endif

