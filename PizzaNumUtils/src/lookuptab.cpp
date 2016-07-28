#include "lookuptab.h"
#include "splines.h"
#include <cmath>
#include <stdexcept>

using namespace Pizza;
using namespace NumUtils;
using namespace std;

lookup_table_range::lookup_table_range()
: x0(2), x1(1)
{}

lookup_table_range::lookup_table_range(const double xmin_, const double xmax_)
{
  init(xmin_, xmax_);
}

void lookup_table_range::init(const double xmin_, const double xmax_)
{
  if (xmin_ >= xmax_)
    throw invalid_argument("Lookup Table: degenerate range");
  x0 = xmin_;
  x1 = xmax_;
}

double lookup_table_range::x_min() const {return x0;}
double lookup_table_range::x_max() const {return x1;}

bool lookup_table_range::contains(const double x) const
{
  return ((x <= x1) && (x >= x0));
}

lookup_table::lookup_table() {}

lookup_table::lookup_table(const double xmin_, const double xmax_, const unsigned int size_,
                     const function_r2r& y_)
{
  init(xmin_, xmax_, size_, y_);
}

lookup_table::lookup_table(const std::vector<double>& x_,
                                  const std::vector<double>& y_)
{
  init(x_, y_);
}


/**
This samples the given function on the specified range.
*/
void lookup_table::init(const double xmin_, const double xmax_,
    const unsigned int size_, const function_r2r& y_)
{
  lookup_table_range::init(xmin_, xmax_);
  dx = (x_max() - x_min()) / (size_ - 1.0);
  y.resize(size_+1);
  for (unsigned int k=0; k<size_; k++) {
    //min is to prevent overstepping of boundary due to roundoff errors.
    double x = min(x_max(), x_min() + dx*k);
    y[k]  = y_(x);
  }
  //Needed when looking up value exactly on the right boundary
  y[size_] = y[size_-1];
}

/**
Convenience method to initialize lookup table from irregularly spaced input,
using cubic spline interpolation.
The lookup table has the same number of entries as the input vectors.
The input x-values should therefore already be more or less regularly spaced,
or accuracy is lost.
*/
void lookup_table::init(const std::vector<double>& x_,
  const std::vector<double>& y_)
{
  if (y_.size() != x_.size())
    throw invalid_argument("LookupTable: column size mismatch.");
  if (!is_strictly_increasing(x_))
    throw invalid_argument("LookupTable: x values not strictly increasing.");

  cubic_spline spl(x_, y_, opt_const::NONE, opt_const::NONE);
  init(x_[0], x_.back(), x_.size(), spl);
}

/**
The result is computed using fast linear interpolation.
If x is outside the tabulated range, an exception is thrown.
*/
double lookup_table::operator()(const double x) const
{
  if (!contains(x))
    throw runtime_error("LookupTable: out of range");
  const double s  = (x - x_min()) / dx;
  const int i     = floor(s);
  const double w  = s - i;
  const double r  = w * y[i+1] + (1.0-w) * y[i];
  return r;
}

bool lookup_table::is_strictly_increasing(const std::vector<double>& x)
{
  for (size_t k=1; k<x.size(); k++) {
    if (x[k] <= x[k-1]) return false;
  }
  return true;
}

lookup_table_loglog::lookup_table_loglog() {}

lookup_table_loglog::lookup_table_loglog(const std::vector<double>& x_,
                                          const std::vector<double>& y_)
{
  init(x_, y_);
}

/**
Initialize logarithmic lookup table from irregularly spaced input.
The logarithm of both x and y values is used to initialize a standard linear
lookup table.
The lookup table has the same number of entries as the input vectors.
The input x-values should therefore already be more or less geometrically spaced,
or accuracy is lost.
*/
void lookup_table_loglog::init(const std::vector<double>& x_, const std::vector<double>& y_)
{
  if (y_.size() != x_.size())
    throw invalid_argument("Lookup Table LogLog: column size mismatch.");
  const size_t n = x_.size();
  vector<double> lx(n),ly(n);
  for (size_t k=0; k<n; k++) {
    if (x_[k] <= 0)
      throw invalid_argument("Lookup Table LogLog: negative x encountered");
    if (y_[k] <= 0)
      throw invalid_argument("Lookup Table LogLog: negative y encountered");
    lx[k] = log(x_[k]);
    ly[k] = log(y_[k]);
  }
  lookup_table_range::init(x_.front(), x_.back());
  tbl.init(lx,ly);
}

double lookup_table_loglog::operator()(const double x) const
{
  if (x <= 0)
    throw runtime_error("Lookup Table LogLog: evaluated at x <= 0");
  return exp(tbl(log(x)));
}


lookup_table_logx::lookup_table_logx() {}

lookup_table_logx::lookup_table_logx(const std::vector<double>& x_,
                                          const std::vector<double>& y_)
{
  init(x_, y_);
}

/**
Initialize semi-logarithmic lookup table from irregularly spaced input.
A standard linear lookup table is initialized with (log(x), y).
The lookup table has the same number of entries as the input vectors.
The input x-values should therefore already be more or less geometrically spaced,
or accuracy is lost.
*/
void lookup_table_logx::init(const std::vector<double>& x_, const std::vector<double>& y_)
{
  if (y_.size() != x_.size())
    throw invalid_argument("Lookup Table LogLog: column size mismatch.");
  const size_t n = x_.size();
  vector<double> lx(n);
  for (size_t k=0; k<n; k++) {
    if (x_[k] <= 0)
      throw invalid_argument("Lookup Table LogLog: negative x encountered");
    lx[k] = log(x_[k]);
  }
  lookup_table_range::init(x_.front(), x_.back());
  tbl.init(lx,y_);
}

double lookup_table_logx::operator()(const double x) const
{
  if (x <= 0)
    throw runtime_error("Lookup Table LogLog: evaluated at x <= 0");
  return tbl(log(x));
}

