#include "functors.h"
#include <stdexcept>
#include <cmath>

namespace Pizza {
namespace NumUtils {

class functor_fptr_r2r : public functor_r2r {
  double (*f)(double);
  public:
  functor_fptr_r2r(double (*f_)(double))
  : f(f_) {}
  virtual double operator()(double x) const
  {
    return (*f)(x);
  }
};

class functor_fptr_rr2r : public functor_rr2r {
  double (*f)(double, double);
  public:
  functor_fptr_rr2r(double (*f_)(double, double))
  : f(f_) {}
  virtual double operator()(double x, double y) const
  {
    return (*f)(x, y);
  }
};


class functor_chain : public functor_r2r{
  const function_r2r f,g;
  public:
  functor_chain(const function_r2r& f_, const function_r2r& g_)
  : f(f_), g(g_) {}
  virtual double operator()(double x) const
  {
    return f(g(x));
  }
};

class functor_bind1st : public functor_r2r{
  const function_rr2r f;
  const double a;
  public:
  functor_bind1st(const function_rr2r& f_, double a_)
  : f(f_), a(a_) {}
  virtual double operator()(double x) const
  {
    return f(a, x);
  }
};

class functor_bind2nd : public functor_r2r{
  const function_rr2r f;
  const double b;
  public:
  functor_bind2nd(const function_rr2r& f_, double b_)
  : f(f_), b(b_) {}
  virtual double operator()(double x) const
  {
    return f(x, b);
  }
};

class functor_apply_binary : public functor_r2r{
  const function_rr2r b;
  const function_r2r f,g;
  public:
  functor_apply_binary(const function_rr2r& b_, const function_r2r& f_,
                 const function_r2r& g_)
  : b(b_), f(f_), g(g_) {}
  virtual double operator()(double x) const
  {
    return b(f(x), g(x));
  }
};

double operator_neg(double a) {return -a;}
double operator_inv(double a) {return 1.0 / a;}
double operator_sum(double a, double b) {return a+b;}
double operator_difference(double a, double b) {return a-b;}
double operator_product(double a, double b) {return a*b;}
double operator_division(double a, double b) {return a/b;}


function_r2r chain(const function_r2r f, const function_r2r g)
{
  return function_r2r(new functor_chain(f, g));
}

function_r2r apply_binary(const function_rr2r op, const function_r2r a,
                          const function_r2r b)
{
  return function_r2r(new functor_apply_binary(op, a, b));
}

function_r2r apply_binary(const function_rr2r op, const function_r2r a,
                          double b)
{
  return op.bind_2nd(b)(a);
}

function_r2r apply_binary(const function_rr2r op, double a,
                          const function_r2r b)
{
  return op.bind_1st(a)(b);
}


function_r2r function_r2r::operator()(const function_r2r& g) const
{
  return chain(*this, g);
}

function_r2r operator-(const function_r2r& a)
{
  return chain(operator_neg, a);
}

function_r2r operator+(const function_r2r& a, const function_r2r& b)
{
  return apply_binary(operator_sum, a, b);
}

function_r2r operator+(const function_r2r& a, double b)
{
  return apply_binary(operator_sum, a, b);
}

function_r2r operator+(double a, const function_r2r& b)
{
  return apply_binary(operator_sum, a, b);
}

function_r2r operator-(const function_r2r& a, const function_r2r& b)
{
  return apply_binary(operator_difference, a, b);
}

function_r2r operator-(const function_r2r& a, double b)
{
  return apply_binary(operator_difference, a, b);
}

function_r2r operator-(double a, const function_r2r& b)
{
  return apply_binary(operator_difference, a, b);
}

function_r2r operator*(const function_r2r& a, const function_r2r& b)
{
  return apply_binary(operator_product, a, b);
}

function_r2r operator*(const function_r2r& a, double b)
{
  return apply_binary(operator_product, a, b);
}

function_r2r operator*(double a, const function_r2r& b)
{
  return apply_binary(operator_product, a, b);
}

function_r2r operator/(const function_r2r& a, const function_r2r& b)
{
  return apply_binary(operator_division, a, b);
}

function_r2r operator/(const function_r2r& a, double b)
{
  return apply_binary(operator_division, a, b);
}

function_r2r operator/(double a, const function_r2r& b)
{
  return apply_binary(operator_division, a, b);
}

}
}

using namespace Pizza::NumUtils;

function_r2r::function_r2r() : f() {}

function_r2r::function_r2r(const function_r2r& other) : f(other.f) {}

function_r2r::function_r2r(functor_r2r* f_) : f(f_) {}

function_r2r::function_r2r(boost::shared_ptr<functor_r2r> f_) : f(f_) {}

function_r2r::function_r2r(double (*fp)(double))
{
  f.reset(new functor_fptr_r2r(fp));
}

function_r2r& function_r2r::operator=(const function_r2r& other)
{
  f=other.f;
  return *this;
}

double function_r2r::operator()(double x) const {
  if (f) return (*f)(x);
  throw std::runtime_error("Function_r2r: uninitialized use");
}

function_rr2r::function_rr2r() : f() {}

function_rr2r::function_rr2r(const function_rr2r& other) : f(other.f) {}

function_rr2r::function_rr2r(functor_rr2r* f_) : f(f_) {}

function_rr2r::function_rr2r(boost::shared_ptr<functor_rr2r> f_) : f(f_) {}

function_rr2r::function_rr2r(double (*fp)(double, double))
{
  f.reset(new functor_fptr_rr2r(fp));
}

function_rr2r& function_rr2r::operator=(const function_rr2r& other)
{
  f=other.f;
  return *this;
}

double function_rr2r::operator()(double x, double y) const {
  if (f) return (*f)(x, y);
  throw std::runtime_error("Function_rr2r: uninitialized use");
}

function_r2r function_rr2r::bind_1st(double x) const
{
  return function_r2r(new functor_bind1st(*this, x));
}

function_r2r function_rr2r::bind_2nd(double x) const
{
  return function_r2r(new functor_bind2nd(*this, x));
}


