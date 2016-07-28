#ifndef LOOKUPTAB_H
#define LOOKUPTAB_H
#include <vector>
#include "functors.h"

namespace Pizza {
namespace NumUtils {

///Valid range of a lookup table
class lookup_table_range {
  double x0;
  double x1;
  public:
  lookup_table_range();
  lookup_table_range(const double xmin_, const double xmax_);
  void init(const double xmin_, const double xmax_);
  ///Lowest tabulated value
  double x_min() const;
  ///Highest tabulated value
  double x_max() const;
  ///If value is inside range
  bool contains(const double x) const;
};

///Lookup table
class lookup_table : public lookup_table_range {
  std::vector<double> y;
  double dx;
  public:
  ///Default constructor.
  lookup_table();
  ///Construct from function
  lookup_table(const double xmin_, const double xmax_, const unsigned int size_,
                     const function_r2r& y_);
  ///Initialize from irregularly spaced data
  lookup_table(const std::vector<double>& x_, const std::vector<double>& y_);

  ///Initialize from function
  void init(const double xmin_, const double xmax_, const unsigned int size_,
             const function_r2r& y_);
  ///Initialize from irregularly spaced data
  void init(const std::vector<double>& x_, const std::vector<double>& y_);
  ///Look up value
  double operator()(const double x) const;
  static bool is_strictly_increasing(const std::vector<double>& x);
};

///Lookup table covering wide range by using logarithms of the values
class lookup_table_loglog : public lookup_table_range {
  lookup_table tbl;
  public:
  ///Default constructor.
  lookup_table_loglog();
  ///Construct from irregularly spaced data
  lookup_table_loglog(const std::vector<double>& x_, const std::vector<double>& y_);
  ///Initialize from irregularly spaced data
  void init(const std::vector<double>& x_, const std::vector<double>& y_);
  ///Look up value
  double operator()(const double x) const;
};

///Lookup table covering wide range in x by using logarithms of x
class lookup_table_logx : public lookup_table_range {
  lookup_table tbl;
  public:
  ///Default constructor.
  lookup_table_logx();
  ///Construct from irregularly spaced data
  lookup_table_logx(const std::vector<double>& x_, const std::vector<double>& y_);
  ///Initialize from irregularly spaced data
  void init(const std::vector<double>& x_, const std::vector<double>& y_);
  ///Look up value
  double operator()(const double x) const;
};

}
}

#endif

