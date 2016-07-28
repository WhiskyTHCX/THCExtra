#ifndef INTERVAL_H
#define INTERVAL_H


namespace Pizza {
namespace NumUtils {


struct interval {
  double xmin;  ///< Minimum
  double xmax;  ///< Maximum
  ///Default constructor: empty range.
  interval() : xmin(1), xmax(-1) {}
  ///Construct from minimum and maximum
  interval(const double xmin_, const double xmax_) : xmin(xmin_), xmax(xmax_) {}
};

struct interval_closed : public interval {
  ///Default constructor: empty range.
  interval_closed() : interval() {}
  ///Construct from minimum and maximum
  interval_closed(const double xmin_, const double xmax_) : interval(xmin_, xmax_) {}
  bool contains(const double x) const {return (x >= xmin) && (x <= xmax);}
  double limit(const double x) const {return (x >= xmax) ? xmax : ((x <= xmin) ? xmin : x);}
};

struct interval_open : public interval {
  ///Default constructor: empty range.
  interval_open() : interval() {}
  ///Construct from minimum and maximum
  interval_open(const double xmin_, const double xmax_) : interval(xmin_, xmax_) {}
  bool contains(const double x) const {return (x > xmin) && (x < xmax);}
};



}
}

#endif
