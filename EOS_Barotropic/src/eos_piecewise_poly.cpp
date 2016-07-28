#include "eos_piecewise_poly.h"
#include <stdexcept>
#include <cmath>

using namespace std;
using namespace EOS_Barotropic;


eos_poly_piece::eos_poly_piece(double rmd0_, double sed0_,
                 double gamma_, double rmd_p_)
: rmd0(rmd0_),  gamma(gamma_), rmd_p(rmd_p_)
{
  n    = 1.0 / (gamma - 1.0);
  np1  = n + 1.0;
  invn = 1.0 / n;
  dsed = sed0_ - n*pow(rmd0/rmd_p, invn);
  gm10 = gm1_from_rmd(rmd0);
  p0   = p_from_gm1(gm10);
}

double eos_poly_piece::sed_from_rmd(const double rmd) const
{
  return sed_from_gm1(gm1_from_rmd(rmd));
}

double eos_poly_piece::p_from_rmd(const double rmd) const
{
  return p_from_gm1(gm1_from_rmd(rmd));
}

double eos_poly_piece::csnd2_from_rmd(const double rmd) const
{
  return csnd2_from_gm1(gm1_from_rmd(rmd));
}

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{\rho}{\rho_p}\right)^{1/n}
         + \delta\epsilon\f$
*/
double eos_poly_piece::gm1_from_rmd(const double rmd) const
{
  return np1 * pow( rmd/rmd_p, invn ) + dsed;
}

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{P}{\rho_p}
          \right)^\frac{1}{n+1} + \delta\epsilon\f$
*/
double eos_poly_piece::gm1_from_p(const double p) const
{
  return np1 * pow(p / rmd_p, 1.0/np1) + dsed;
}

/**
\return Specific internal energy \f$ \epsilon =
\frac{g-1-\delta\epsilon}{\Gamma} +\delta\epsilon\f$
*/
double eos_poly_piece::sed_from_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  return (gm1 - dsed) / gamma + dsed;
}

/**
\return Internal energy density
\f$ \rho_I = n \rho_p \left(\frac{g-1}{n+1}\right)^{n+1} \f$
*/
double eos_poly_piece::ied_from_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  return sed_from_gm1(gm1)*rmd_from_gm1(gm1);
}


/**
\return Pressure \f$ P =
\rho_p \left( \frac{g-1-\delta\epsilon}{1+n} \right)^{1+n} \f$
*/
double eos_poly_piece::p_from_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow((gm1 - dsed) / np1, np1 );
}

/**
\return Rest mass density \f$ \rho
= \rho_p \left( \frac{g-1-\delta\epsilon}{1+n} \right)^n \f$
*/
double eos_poly_piece::rmd_from_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return rmd_p * pow( (gm1-dsed) / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$
*/
double eos_poly_piece::hm1_from_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{g-1}{ng} \f$
*/
double eos_poly_piece::csnd2_from_gm1(const double gm1) const
{
  double gm1p  = gm1 - dsed;
  return gm1p / ( n * (gm1 + 1.0));
}



eos_piecewise_poly::eos_piecewise_poly(double rmdp0,
  const vector<double>& segm_bound,
  const  vector<double>& segm_gamma,
  double rmd_max_, std::string descr_)
: eos_1p_api(true, false, false, descr_)
{

  if (segm_bound.size() != segm_gamma.size())
    throw runtime_error("eos_piecewise_poly: vector sizes mismatch.");
  if (segm_bound.size()==0)
    throw runtime_error("eos_piecewise_poly: need at least one segment.");
  if (segm_bound[0] != 0)
    throw runtime_error("eos_piecewise_poly: First segment has to start at zero density.");

  vector<double>::const_iterator ibnd = segm_bound.begin();
  vector<double>::const_iterator iga  = segm_gamma.begin();

  segments.push_back(eos_poly_piece(0, 0, *iga, rmdp0));

  while ((++ibnd != segm_bound.end())
         && (++iga != segm_gamma.end()))
  {
    const eos_poly_piece& lseg = segments.back();
    if (lseg.rmd0 >= *ibnd)
      throw runtime_error("eos_piecewise_poly: segment boundary densities not strictly increasing.");

    double sedc = lseg.sed_from_rmd(*ibnd);
    double np   = 1.0 / (*iga - 1.0);
    double et   = np / lseg.n;
    double rmdp = pow(lseg.rmd_p, et) * pow(*ibnd, 1.0-et);
    segments.push_back(eos_poly_piece(*ibnd, sedc, *iga, rmdp));
  }

  set_ranges(range(0, rmd_max_));
}


const eos_poly_piece&
eos_piecewise_poly::segment_for_rmd(double rmd) const
{
  vector<eos_poly_piece>::const_reverse_iterator
    i = segments.rbegin();
  while (i->rmd0 > rmd) {
    if (++i == segments.rend()) return segments[0];
  }
  return *i;
}

const eos_poly_piece&
eos_piecewise_poly::segment_for_p(double p) const
{
  vector<eos_poly_piece>::const_reverse_iterator
    i = segments.rbegin();
  while (i->p0 > p) {
    if (++i == segments.rend()) return segments[0];
  }
  return *i;
}

const eos_poly_piece&
eos_piecewise_poly::segment_for_gm1(double gm1) const
{
  vector<eos_poly_piece>::const_reverse_iterator
    i = segments.rbegin();
  while (i->gm10 > gm1) {
    if (++i == segments.rend()) return segments[0];
  }
  return *i;
}


/**
\return \f$ g-1 \f$ from polytropic piece at given density.
*/
double eos_piecewise_poly::gm1_from_valid_rmd(const double rmd) const
{
  return segment_for_rmd(rmd).gm1_from_rmd(rmd);
}


/**
\return \f$ g-1 \f$ from polytropic piece at given pressure.
*/
double eos_piecewise_poly::gm1_from_valid_p(const double p) const
{
  return segment_for_p(p).gm1_from_p(p);
}

/**
\return Specific internal energy \f$ \epsilon \f$
from polytropic piece at given \f$ g-1 \f$.
*/
double eos_piecewise_poly::sed_from_valid_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).sed_from_gm1(gm1);
}

/**
\return Internal energy density \f$ \rho_I \f$
from polytropic piece at given \f$ g-1 \f$.
*/
double eos_piecewise_poly::ied_from_valid_gm1(
  const double gm1  ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).ied_from_gm1(gm1);
}


/**
\return Pressure \f$ P \f$
from polytropic piece at given \f$ g-1 \f$.
*/
double eos_piecewise_poly::p_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).p_from_gm1(gm1);
}

/**
\return Rest mass density \f$ \rho \f$
from polytropic piece at given \f$ g-1 \f$.
*/
double eos_piecewise_poly::rmd_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return segment_for_gm1(gm1).rmd_from_gm1(gm1);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$.
*/
double eos_piecewise_poly::hm1_from_valid_gm1(
  const double gm1 ///< \f$ g-1 \f$
) const
{
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2  \f$
from polytropic piece at given \f$ g-1 \f$.
*/
double eos_piecewise_poly::csnd2_from_valid_gm1(const double gm1) const
{
  return segment_for_gm1(gm1).csnd2_from_gm1(gm1);
}

