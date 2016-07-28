#ifndef MODES_H
#define MODES_H

#include "tovsol.h"
#include <boost/shared_ptr.hpp>
#include <cstring>
#include <complex>
#include <string>

namespace Pizza {
namespace TOV {

///Class describing pressure modes of nonrotating stars in Cowling approximation.
class pmode_cowling
{
  class impl;
  boost::shared_ptr<impl> i;
  const impl& si() const;

  public:

  struct pert_coeff
  {
    double r, scalar, radial, angular;
    pert_coeff() {}
    pert_coeff(double r_, double scalar_, double radial_, double angular_);
  };

  struct pert_local
  {
    std::complex<double> dh, dvr, dvth, dvphi;
    pert_local(const pert_coeff& pc, int l, int m, double th, double phi);
  };


  pmode_cowling(const tovsol& star, int l, int m, int n, double acc);
  pert_coeff operator()(const double r) const;
  pert_local operator()(const double r, const double th, const double phi) const;

  int l() const;
  int m() const;
  int n() const;
  double freq() const;
  double kin_energy() const;
  static double freq(const tovsol& star, int l, int n, double acc);
  void save_astext(std::string fn, size_t res=1000) const;
};

}
}
#endif


