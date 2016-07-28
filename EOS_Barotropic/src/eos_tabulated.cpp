#include "eos_tabulated.h"
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;
using namespace EOS_Barotropic;
using namespace Pizza::NumUtils;


/**
The input tables will be used to set up logarithmic lookup tables.
All values have to be strictly positive.
The input data can be irregularly spaced, since internally the logarithms
are resampled by cubic splines to set up the tables.
However, the data should be spaced roughly in geometric progression to make
ideal use of a given number of points.

Further, the input tables have to be sorted in ascending order of rest mass density.
Pressure <c>p</c> and <c>gm1</c>  have to be strictly increasing, or an exception
will be thrown. Note it can happen that the cubic interpolation does not conserve
monotonicity if the values have a plateau, in which case an exception is thrown, too.
Soundspeeds above the speed of light are strictly forbidden (an exception will be thrown),
as are values \f$ c_s^2 <=0 \f$.
For densities below the lowest one given, a matching polytrope is constructed.
The units of the input tables are the same as for the resulting EOS object.
By convention, the units should satisfy \f$ G=c=1 \f$.
*/
eos_tabulated::eos_tabulated(
    const std::vector<double>& gm1,
    const std::vector<double>& rmd,
    const std::vector<double>& sed,
    const std::vector<double>& p,
    const std::vector<double>& cs2,
    const std::vector<double>& temp,
    const std::vector<double>& efrac,
    bool isentropic_,
    std::string descr_)
: eos_1p_api(isentropic_, !temp.empty(), !efrac.empty(), descr_),
  poly(rmd[0], sed[0], p[0], 1.1*rmd[0])
{
  const size_t tsize = rmd.size();
  if ((tsize != sed.size()) || (tsize != p.size()) ||
      (tsize != cs2.size()) || (tsize != gm1.size()) ||
      (has_temp() && (temp.size() != tsize)) ||
      (has_efrac() && (efrac.size() != tsize)))
  {
    throw invalid_argument("Tabulated EOS: table column size mismatch");
  }

  const double gsc = (poly.gm1_from_rmd(rmd[0]) - gm1[0]) / (1.0 + gm1[0]);

  vector<double> ngm1(tsize), hm1(tsize);
  for (size_t k=0; k < tsize; k++) {
    ngm1[k] = gm1[k] + gsc * (1.0 + gm1[k]);
    hm1[k]  = sed[k] + p[k] / rmd[k];
    if (cs2[k] >= 1.0)
      throw runtime_error("Tabulated EOS: superluminal sound speed");
    if (cs2[k] <= 0.0)
      throw runtime_error("Tabulated EOS: imaginary sound speed");
    if (sed[k] <= 0.0)
      throw runtime_error("Tabulated EOS: specific internal energy <= 0");
    if (p[k] <= 0.0)
      throw runtime_error("Tabulated EOS: pressure <= 0");
    if (gm1[k] <= 0.0)
      throw runtime_error("Tabulated EOS: g <= 1");
  }

  rmd_gm1.init(ngm1, rmd);
  sed_gm1.init(ngm1, sed);
  p_gm1.init(ngm1, p);
  hm1_gm1.init(ngm1, hm1);
  cs2_gm1.init(ngm1, cs2);
  gm1_rmd.init(rmd, ngm1);
  gm1_p.init(p, ngm1);
  temp0  = 0.0;
  efrac0 = 0.0;
  if (has_temp()) {
    temp_gm1.init(ngm1, temp);
    temp0 = temp[0];
  }
  if (has_efrac()) {
    efrac_gm1.init(ngm1, efrac);
    efrac0 = efrac[0];
  }
  set_ranges(range(0,gm1_rmd.x_max()));
}

double eos_tabulated::gm1_from_valid_rmd(const double rmd) const
{
  return (rmd > gm1_rmd.x_min()) ? gm1_rmd(rmd) : poly.gm1_from_rmd(rmd);
}

double eos_tabulated::gm1_from_valid_p(const double p) const
{
  return (p > gm1_p.x_min()) ? gm1_p(p) : poly.gm1_from_p(p);
}

double eos_tabulated::sed_from_valid_gm1(const double gm1) const
{
  return (gm1 > sed_gm1.x_min()) ? sed_gm1(gm1) : poly.sed_from_gm1(gm1);
}

double eos_tabulated::ied_from_valid_gm1(const double gm1) const
{
  return (gm1 > sed_gm1.x_min()) ? (rmd_gm1(gm1)*sed_gm1(gm1)) : poly.ied_from_gm1(gm1);
}

double eos_tabulated::p_from_valid_gm1(const double gm1) const
{
  return (gm1 > p_gm1.x_min()) ? p_gm1(gm1) : poly.p_from_gm1(gm1);
}

double eos_tabulated::rmd_from_valid_gm1(const double gm1) const
{
  return (gm1 > rmd_gm1.x_min()) ? rmd_gm1(gm1) : poly.rmd_from_gm1(gm1);
}

double eos_tabulated::hm1_from_valid_gm1(const double gm1) const
{
  return (gm1 > hm1_gm1.x_min()) ? hm1_gm1(gm1) : poly.hm1_from_gm1(gm1);
}

double eos_tabulated::csnd2_from_valid_gm1(const double gm1) const
{
  return (gm1 > cs2_gm1.x_min()) ? cs2_gm1(gm1) : poly.csnd2_from_gm1(gm1);
}

double eos_tabulated::temp_from_valid_gm1(const double gm1) const
{
  if (!has_temp())
    throw runtime_error("Tabulated EOS: temperature not available.");
  return (gm1 > temp_gm1.x_min()) ? temp_gm1(gm1) : temp0;
}

double eos_tabulated::efrac_from_valid_gm1(const double gm1) const
{
  if (!has_efrac())
    throw runtime_error("Tabulated EOS: electron fraction not available.");
  return (gm1 > efrac_gm1.x_min()) ? efrac_gm1(gm1) : efrac0;
}


