#include "eos_idealgas.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace whizza;
using namespace std;

eos_idealgas::eos_idealgas(pz_real n_, pz_real umass_, const range& rgeps_, const range& rgrho_, const range& rgye_)
: gamma(1.0 + 1.0/n_), gm1(1.0/n_), rgeps(rgeps_)
{
  if (n_ < 0) {
    throw runtime_error("EOS_Thermal_Idealgas: initialized with gamma < 1");
  }
  if (gamma > 2.0) { // Ensure subluminal Soundspeed and P < E
    rgeps.max = min(rgeps.max, 1.0 / (gamma*(gamma - 2.0)));
  }
  //set_range_h(range(1.0 + gamma * rgeps.min, 1.0 + gamma * rgeps.max));
  set_range_rho(rgrho_);
  set_range_ye(rgye_);
  temp_over_eps =  gm1 * umass_;
  set_range_temp(range(temp_over_eps * rgeps.min, temp_over_eps * rgeps.max));
}

eos_idealgas::~eos_idealgas() {}

pz_real eos_idealgas::press_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  return gm1 * rho * eps;
}

pz_real eos_idealgas::csnd_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  return sqrt(gm1 * eps / (eps + 1.0/gamma));
}

pz_real eos_idealgas::temp_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  return  temp_over_eps * eps;
}

void eos_idealgas::press_derivs_from_valid_rho_eps_ye(pz_real& press, pz_real& dpdrho, pz_real& dpdeps,
           const pz_real rho, const pz_real eps, const pz_real ye) const
{
  press  = gm1 * rho * eps;
  dpdrho = gm1 * eps;
  dpdeps = gm1 * rho;
}

pz_real eos_idealgas::entropy_from_valid_rho_temp_ye(const pz_real rho,
  const pz_real temp, const pz_real ye) const
{
  return log(temp*pow(rho, -gm1)/temp_over_eps);
}

pz_real eos_idealgas::entropy_from_valid_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye) const
{
  return log(eps*pow(rho, -gm1));
}

pz_real eos_idealgas::eps_from_valid_rho_temp_ye(const pz_real rho, const pz_real temp, const pz_real ye) const
{
  return temp / temp_over_eps;
}

eos_thermal_impl::range eos_idealgas::range_eps_from_valid_rho_ye(const pz_real rho, const pz_real ye) const
{
  return rgeps;
}

eos_thermal whizza::make_eos_idealgas(pz_real n, pz_real umass,
    const eos_thermal::range& rgeps,
    const eos_thermal::range& rgrho,
    const eos_thermal::range& rgye)
{
  return eos_thermal(new eos_idealgas(n, umass, rgeps, rgrho, rgye));
}

