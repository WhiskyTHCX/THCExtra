#include "eos_hybrid.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace whizza;
using namespace EOS_Barotropic;
using namespace std;


eos_hybrid::eos_hybrid(eos_1p& eos_c_,  pz_real gamma_th_, pz_real eps_max_, pz_real eps_min_,
                       pz_real rho_max_ , const range& rgye_)
: eos_c(eos_c_), gamma_th(gamma_th_), eps_max(eps_max_), eps_min(eps_min_)
{
  if ((gamma_th>2) || (gamma_th<=1)) {
    throw runtime_error("eos_hybrid: gamma_thermal must be in the range (1,2]");
  }
  gm1_th = gamma_th - 1.0;
  set_range_rho(range(0.0, rho_max_));
  set_range_ye(rgye_);
  set_range_temp(range(0,0)); //TODO: Implement temperature;
}


eos_hybrid::~eos_hybrid() {}


pz_real eos_hybrid::eps_cold(pz_real rho) const
{
  return eos_c.sed_from_rmd(rho);
}

pz_real eos_hybrid::p_cold(pz_real rho) const
{
  return eos_c.p_from_rmd(rho);
}

pz_real eos_hybrid::hm1_cold(pz_real rho) const
{
  return eos_c.hm1_from_rmd(rho);
}

pz_real eos_hybrid::cs2_cold(pz_real rho) const
{
  return eos_c.csnd2_from_rmd(rho);
}



pz_real eos_hybrid::press_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  pz_real p_c   = p_cold(rho);
  pz_real eps_c = eps_cold(rho);
  pz_real p_th  = gm1_th * rho * (eps - eps_c);
  return p_c + p_th;
}

pz_real eos_hybrid::csnd_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  pz_real cs2_c   = cs2_cold(rho);
  pz_real eps_c   = eps_cold(rho);
  pz_real h_c     = 1.0 + hm1_cold(rho);
  pz_real eps_th  = eps - eps_c;
  pz_real h_th    = gamma_th * eps_th;
  pz_real w       = h_th / (h_c + h_th);
  pz_real cs2     = (1.0 - w) * cs2_c + w * gm1_th;
  return sqrt(cs2);
}

pz_real eos_hybrid::temp_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  throw runtime_error("eos_hybrid: temperature not implemented");
}

void eos_hybrid::press_derivs_from_valid_rho_eps_ye(pz_real& press, pz_real& dpdrho, pz_real& dpdeps,
           const pz_real rho, const pz_real eps, const pz_real ye) const
{
  pz_real p_c     = p_cold(rho);
  pz_real cs2_c   = cs2_cold(rho);
  pz_real eps_c   = eps_cold(rho);
  pz_real h_c     = 1.0 + hm1_cold(rho);
  pz_real eps_th  = eps - eps_c;
  pz_real p_th    = gm1_th * rho * eps_th;
  press           = p_c + p_th;
  dpdrho          = h_c * cs2_c + gm1_th * (eps_th - p_c / rho);
  dpdeps          = gm1_th * rho;
}

pz_real eos_hybrid::entropy_from_valid_rho_temp_ye(const pz_real rho,
  const pz_real temp, const pz_real ye) const
{
  throw logic_error("EOS: entropy from temperature not implemented for EOS_Thermal_Hybrid.");
}

pz_real eos_hybrid::entropy_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  pz_real eps_th = eps - eps_cold(rho);
  return log(eps_th*pow(rho, -gm1_th));
}

pz_real eos_hybrid::eps_from_valid_rho_temp_ye(const pz_real rho, const pz_real temp, const pz_real ye) const
{
  throw runtime_error("eos_hybrid: temperature not implemented");
}

eos_thermal_impl::range eos_hybrid::range_eps_from_valid_rho_ye(const pz_real rho, const pz_real ye) const
{
  if      ( eps_min >= 0.0 ) {
    return range(eps_min, eps_max);
  }
  else if ( eps_min == -1.0 ) {
    pz_real eps_min_cold = eps_cold(rho);
    return range(eps_min_cold, eps_max);
  }
  else {
    throw runtime_error("eos_hybrid: eps_min is wrong.");
  }
}

eos_thermal whizza::make_eos_hybrid(EOS_Barotropic::eos_1p& eos_c,
             pz_real gamma_th, pz_real eps_max, pz_real eps_min, pz_real rho_max,
             const eos_thermal::range& rgye)
{
  return eos_thermal(new eos_hybrid(eos_c, gamma_th, eps_max, eps_min, rho_max,
                                    rgye));
}



