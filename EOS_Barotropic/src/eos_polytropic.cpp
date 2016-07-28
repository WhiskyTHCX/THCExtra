#include "eos_polytropic.h"
#include <cmath>

using namespace std;
using namespace EOS_Barotropic;

void eos_polytrope::init(double n_, double rmd_p_, double rmd_max_)
{
  n     = n_;
  rmd_p = rmd_p_;
  np1   = n + 1;
  gamma = 1.0 + 1.0 / n;
  invn  = 1.0 / n;

  set_ranges(range(0, rmd_max_));
}


eos_polytrope::eos_polytrope(double n_, double rmd_p_,
                             double rmd_max_, std::string descr_)
: eos_1p_api(true, false, false, descr_)
{
  init(n_, rmd_p_, rmd_max_);
}

/**
Construct polytrope by specifying \f$ \epsilon \f$ and \f$ P \f$
at a given density \f$ \rho \f$, using the formulas
\f{eqnarray*}{
  n      &=& \frac{\epsilon \rho}{P}  \\
  \rho_p &=& \rho \left(\frac{\rho}{P} \right)^n
\f}
*/
eos_polytrope::eos_polytrope(double rmd_m, double sed_m, double p_m,
                             double rmd_max_, std::string descr_)
: eos_1p_api(true, false, false, descr_)
{
  double n_     = sed_m * rmd_m / p_m;
  double rmd_p_ = rmd_m * pow(rmd_m / p_m, n_);
  init(n_, rmd_p_, rmd_max_);
}

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{\rho}{\rho_p}\right)^{1/n} \f$
*/
double eos_polytrope::gm1_from_valid_rmd(const double rmd) const
{
  return np1 * pow( rmd/rmd_p, invn );
}


/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{P}{\rho_p} \right)^\frac{1}{n+1} \f$
*/
double eos_polytrope::gm1_from_valid_p(const double p) const
{
  return np1 * pow(p / rmd_p, 1.0/np1);
}

/**
\return Specific internal energy \f$ \epsilon = \frac{g-1}{\Gamma} \f$
*/
double eos_polytrope::sed_from_valid_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  return gm1 / gamma;
}

/**
\return Internal energy density
\f$ \rho_I = n \rho_p \left(\frac{g-1}{n+1}\right)^{n+1} \f$
*/
double eos_polytrope::ied_from_valid_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  //return sed_from_gm1(gm1)*rmd_from_gm1(gm1);
  return n * rmd_p * pow( gm1 / np1, np1 );
}


/**
\return Pressure \f$ P = \rho_p \left( \frac{g-1}{1+n} \right)^{1+n} \f$
*/
double eos_polytrope::p_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow( gm1 / np1, np1 );
}

/**
\return Rest mass density \f$ \rho = \rho_p \left( \frac{g-1}{1+n} \right)^n \f$
*/
double eos_polytrope::rmd_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow( gm1 / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$
*/
double eos_polytrope::hm1_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{g-1}{ng} \f$
*/
double eos_polytrope::csnd2_from_valid_gm1(const double gm1) const
{
  return gm1/(n*(gm1+1));
}

