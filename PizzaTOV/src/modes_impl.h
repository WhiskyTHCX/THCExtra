#include <vector>
#include <cstring>
#include "modes.h"
#include "sturm.h"
#include "polynomial.h"
#include "splines.h"

namespace Pizza {
namespace TOV {

using NumUtils::polynomial_r;
using NumUtils::cubic_spline;
using EOS_Barotropic::eos_cold;


///Sturm-Liouville problem describing pressure mode eigenfunctions of TOV stars.
class pruessp : public pruess_problem
{

  const tovsol& star;
  const eos_cold eos;
  const int l;
  polynomial_r left_a, left_bc, left_bv,right_a, right_bc, right_bv;

  void get_eos_lim_surf(double& n, double& dhn) const;
  void get_left_coeff();
  void get_right_coeff();

  public:
  pruessp(const tovsol& star_, const int l_);
  double est_ev(const int n) const;
  virtual ~pruessp() {}
  virtual void coefficients(const double x, double& p, double& kc,
      double& kv) const;
  virtual void left_id(const double ev, const double x, double& u, double& du) const;
  virtual void right_id(const double ev, const double x, double& u, double& du) const;
};

///Eigenfunction spherical harmonics coefficients at given radius.
struct ess_pert
{
  double scalar,radial1,radial2,angular;
  ess_pert(double scalar_, double radial1_, double radial2_, double angular_)
    : scalar(scalar_), radial1(radial1_), radial2(radial2_), angular(angular_) {}
};

///Radial dependence of eigenfunction spherical harmonics coefficients
struct ess_pert_spline
{
  cubic_spline scalar,radial1,radial2,angular;
  ess_pert operator()(double r) const;
  ess_pert_spline() {}
  ess_pert_spline(const std::vector<double>& r_, const std::vector<double>& scalar_,
    const std::vector<double>& radial1_, const std::vector<double>& radial2_,
    const std::vector<double>& angular_);
};

///Implementation belonging to user interface \ref pmode_cowling
class pmode_cowling::impl
{
  typedef EOS_Barotropic::eos_cold eos_cold;

  const tovsol star;
  const eos_cold eos;
  const int md_l;
  const int md_m;
  const int md_n;
  double md_freq;
  double r_surf;
  ess_pert_spline pspl;
  pert_coeff ess2coeff(const ess_pert& ep, const double r) const;
  public:

  impl(const tovsol& star_, int l_, int m_, int n_, const double acc);
  int l() const {return md_l;}
  int m() const {return md_m;}
  int n() const {return md_n;}
  double freq() const {return md_freq;}
  double kin_energy() const;
  void save_astext(std::string fn, size_t res) const;
  pert_coeff operator()(const double r) const;
  pert_local operator()(const double r, const double th, const double phi) const;

};

}
}

