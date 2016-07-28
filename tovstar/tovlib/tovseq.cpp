#include "tovseq.h"
#include "roots.h"
#include "minima.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/format.hpp>

using namespace std;
using namespace EOS_Barotropic;
using namespace Pizza::TOV;
using namespace Pizza::NumUtils;

/*
tovseq::tovseq(double rmd_0_, double rmd_1_, const eos_1p &eos_,
                int n_steps, const double gm1_surf_, const double acc_)
: eos(eos_), rmd_0(rmd_0_), rmd_1(rmd_1_)
{
  if (rmd_0>rmd_1) swap(rmd_0,rmd_1);
  vector<double> vrmd, vgrmass, vrad, vebind;
  for (double lrmd=log(rmd_0); lrmd < log(rmd_1); lrmd += (log(rmd_1)-log(rmd_0))/n_steps) {
    double rmd = exp(lrmd);
    tovsol star(rmd, eos, 0, acc_, gm1_surf_);
    vrmd.push_back(rmd);
    vgrmass.push_back(star.grav_mass());
    vrad.push_back(star.radius());
    vebind.push_back(star.binding_energy());
  }
  grav_mass_rmd = cubic_spline(vrmd, vgrmass, opt_const::NONE, opt_const::NONE);
  radius_rmd    = cubic_spline(vrmd, vrad,    opt_const::NONE, opt_const::NONE);
  ebind_rmd     = cubic_spline(vrmd, vebind,  opt_const::NONE, opt_const::NONE);
}

void tovseq::save(string fn, units u_out) const
{
  units u(units::geom_meter() / u_out);
  ofstream f(fn.c_str());
  f << "# Units: " << u_out << endl
    << "# rmd_c    grav_mass   radius bind_energy/c^2" << endl;
  boost::format fmt("%24.15e %24.15e %24.15e %24.15e\n");
  for (size_t i=0; i<radius_rmd.size(); i++) {
    f << fmt % (radius_rmd.xi(i) *u.density())
             % (grav_mass_rmd.yi(i) *u.mass())
             % (radius_rmd.yi(i) *u.length())
             % (ebind_rmd.yi(i) * u.mass());
  }
}
*/


struct gmass_func : public functor_r2r {
  const tov_sequence& ftov;
  gmass_func(const tov_sequence& ftov_) : ftov(ftov_) {}
  virtual double operator()(double rmdc) const
  {
    return ftov.get_star(rmdc).grav_mass();
  }
};

struct bmass_func : public functor_r2r {
  const tov_sequence& ftov;
  bmass_func(const tov_sequence& ftov_) : ftov(ftov_) {}
  virtual double operator()(double rmdc) const
  {
    return ftov.get_star(rmdc).baryonic_mass();
  }
};

tovsol tov_sequence::get_star(double rhoc) const {
  return tovsol(rhoc, eos, 0, acc_star, gm1_surf);
}

tovsol tov_sequence::operator()(double rhoc) const {
  if (!rmd_range.contains(rhoc))
    throw runtime_error("tov_sequence evaluated out of specified range");
  return get_star(rhoc);
}

tov_sequence::tov_sequence(double rmd_0_, double rmd_1_,
          const EOS_Barotropic::eos_1p &eos_,
          double gm1_surf_, int n_steps_,
          double acc_star_, double acc_max_, double acc_root_)
: eos(eos_), gm1_surf(gm1_surf_), rmd_range(rmd_0_, rmd_1_),
  acc_star(acc_star_), acc_root(acc_root_)
{
  const double lrmd0(log(rmd_range.xmin));
  const double lrmd1(log(rmd_range.xmax));

  const function_r2r fbmass(chain(new bmass_func(*this), exp));
  rmdc_max_mb = exp(find_maximum(fbmass, lrmd0, lrmd1,
                      n_steps_, acc_max_, 0, 10000));
  max_mb      = bmass_from_rmdc(rmdc_max_mb);

  const function_r2r fgmass(chain(new gmass_func(*this), exp));
  rmdc_max_mg = exp(find_maximum(fgmass,  lrmd0, lrmd1,
                      n_steps_, acc_max_, 0, 10000));
  max_mg      = gmass_from_rmdc(rmdc_max_mg);
}

double tov_sequence::bmass_from_rmdc(double rmdc) const
{
  return (*this)(rmdc).baryonic_mass();
}

double tov_sequence::gmass_from_rmdc(double rmdc) const
{
  return (*this)(rmdc).grav_mass();
}

double tov_sequence::rmdc_from_bmass(double mb,
                                double rmd0, double rmd1) const
{
  if (mb > max_mb)
    throw runtime_error("Baryon mass > maximum mass specified.");
  if (mb == max_mb)
    return rmdc_max_mb;
  const function_r2r fbmass(chain(new bmass_func(*this), exp));
  double lrhoc = findroot(fbmass - mb, log(rmd0), log(rmd1),
                          acc_root, 0, 10000);
  return exp(lrhoc);
}

double tov_sequence::rmdc_from_gmass(double mg,
                                double rmd0, double rmd1) const
{
  if (mg > max_mg)
    throw runtime_error("Grav. mass > maximum mass specified.");
  if (mg == max_mg)
    return rmdc_max_mg;
  const function_r2r fgmass(chain(new gmass_func(*this), exp));
  double lrhoc = findroot(fgmass - mg, log(rmd0), log(rmd1),
                          acc_root, 0, 10000);
  return exp(lrhoc);
}

double tov_sequence::rmdc_from_bmass_stable(double mb) const
{
  return rmdc_from_bmass(mb, rmd_range.xmin, rmdc_max_mb);
}

double tov_sequence::rmdc_from_bmass_unstable(double mb) const
{
  return rmdc_from_bmass(mb, rmdc_max_mb, rmd_range.xmax);
}

double tov_sequence::rmdc_from_gmass_stable(double mg) const
{
  return rmdc_from_gmass(mg, rmd_range.xmin, rmdc_max_mb);
}

double tov_sequence::rmdc_from_gmass_unstable(double mg) const
{
  return rmdc_from_gmass(mg, rmdc_max_mb, rmd_range.xmax);
}


std::string tov_sequence::fmt_star(const tovsol& star,
                                   Pizza::units uc) const
{
  boost::format fmt("%24.15e %24.15e %24.15e %24.15e\n");
  fmt % (star.rmd_c() * uc.density())
      % (star.grav_mass() * uc.mass())
      % (star.radius() * uc.length())
      % (star.binding_energy() * uc.mass());
  return fmt.str();
}


void tov_sequence::save(std::string fn, std::string fn_stable,
            std::string fn_unstable, std::string fn_max_mg,
            std::string fn_max_mb, int ndivisions,
            Pizza::units uout, Pizza::units useq) const
{
  units uc(useq / uout);
  ofstream s_all(fn.c_str());
  ofstream s_stbl(fn_stable.c_str());
  ofstream s_unst(fn_unstable.c_str());
  ofstream s_mmb(fn_max_mb.c_str());
  ofstream s_mmg(fn_max_mg.c_str());

  string header = string("# Units: ") + uout.to_str()
          +"\n# rmd_c    grav_mass   radius bind_energy/c^2\n";

  s_all  << header;
  s_stbl << header;
  s_unst << header;

  const double lrmd0(log(rmd_range.xmin));
  const double lrmdm(log(rmdc_max_mg));
  const double lrmd1(log(rmd_range.xmax));
  const double dlrmds((lrmdm - lrmd0) / ndivisions);
  const double dlrmdu((lrmd1 - lrmdm) / ndivisions);

  for (int i=0; i<ndivisions; i++) {
    const double rmdc = exp(lrmd0 + dlrmds*i);
    tovsol star = get_star(rmdc);
    string l    = fmt_star(star, uc);
    s_all  << l;
    s_stbl << l;
  }

  {
    const tovsol star_max_mg = get_star(rmdc_max_mg);
    s_mmg << star_max_mg.to_str(uout, useq) << endl;

    string ml    = fmt_star(star_max_mg, uc);
    s_all  << ml;
    s_stbl << ml;
    s_unst << ml;
  }

  for (int i=1; i<=ndivisions; i++) {
    const double rmdc = exp(lrmdm + dlrmdu*i);
    tovsol star = get_star(rmdc);
    string l    = fmt_star(star, uc);
    s_all  << l;
    s_unst << l;
  }

  const tovsol star_max_mb = get_star(rmdc_max_mb);
  s_mmb << star_max_mb.to_str(uout, useq) << endl;
}


