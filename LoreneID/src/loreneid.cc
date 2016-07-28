#include "loreneid.h"
#include "pizza_central.h"
#include <stdexcept>
#include <gsl/gsl_math.h>

using namespace Pizza;
using namespace Pizza::Base;
using namespace Pizza::IDBase;
using namespace Pizza::LoreneID;

loreneid::loreneid(const eos_1p& eos_, bool set_shift_, bool set_lapse_, bool set_temp_, bool set_efrac_,
    bool do_modify_, bool old_behaviour_, double mass_conv_, double rho_cut_, double scale_vel_,
    double spinup_fac_, double scale_above_density_)
: eos(eos_), set_shift(set_shift_), set_lapse(set_lapse_), set_temp(set_temp_), set_efrac(set_efrac_),
  do_modify(do_modify_), old_behaviour(old_behaviour_), mass_conv(mass_conv_),
  rho_cut(rho_cut_), scale_vel(scale_vel_),
  spinup_fac(spinup_fac_), scale_above_density(scale_above_density_)
{}


/// Set all members of initial data partially initialized from LORENE data.
/** We need to convert the units of rho and Kij. Also, we need
    to consider different definitions of formal rest mass density,
    expressed by the parameter mass_conv. The shift vector is defined
    with opposite sign.
    Pressure and specific internal energy are then recomputed from
    the initial data EOS, which should be the same as used in LORENE.
    Points where LORENE produces NANs/INFs or negative density in the hydro
    part, it is set to vacuum with zero velocity, and a crap counter is
    increased. If density is zero, velocity is set to 0 as well.
    If there are INFs/NANs in the spacetime part, an exception is thrown.
    Finally, the velocity is rescaled by a factor scale_vel (defaults to 1).
**/
void loreneid::complete_vars(idvars& id, int& crap_cnt)
{
  const double km_si  = 1e3;
  const units u       = pizza_base_central::get().internal_units;

  for (int i=0; i<3; i++) {
    for (int j=0;j<=i;j++) {
      if (0==gsl_finite(id.klo(i,j)))
        throw std::runtime_error("LoreneID: K_ij contains NaNs");
      if (0==gsl_finite(id.glo(i,j)))
        throw std::runtime_error("LoreneID: g_ij contains NaNs");
    }
  }

  double lorf_orig = id.florentz_from_vel();

  // Convert rho_cut to Lorene units
  double rho_cut_lorene = rho_cut * u.density() * mass_conv;

  bool id_is_crap =((0==gsl_finite(id.rmd)) ||
                    (id.rmd < 0) ||
                    (0==gsl_finite(id.vel(0))) ||
                    (0==gsl_finite(id.vel(1))) ||
                    (0==gsl_finite(id.vel(2)))  ||
                    (0==gsl_finite(lorf_orig)));

  if (id_is_crap) crap_cnt++;

  // if data from LORENE is crap, set to vacuum
  if (id_is_crap || (id.rmd < rho_cut_lorene)) {
    id.rmd      = 0;
    id.sed      = 0;
    id.press    = 0;
    id.vel      = ZERO;
    id.florentz = 1.0;
  }
  else {
    id.rmd           /= u.density(); // convert units
    id.rmd           /= mass_conv;   // consider different atomic mass units
    //id.vel           *= scale_vel;
    id.florentz       = id.florentz_from_vel();
    const double gm1  = eos.gm1_from_rmd(id.rmd);
    id.sed            = eos.sed_from_gm1(gm1);
    id.press          = eos.p_from_gm1(gm1);
  }


  id.klo             *= (u.length() / km_si);
}

/// Optionally, set temperature and electron fraction
void loreneid::optional_id(const region::iterator& i, IDBase::idvars& id, const double omega)
{
  vec_u shift_id;
  if (set_temp || set_efrac) {
    const double gm1 = eos.gm1_from_rmd(id.rmd);
    if (set_temp)
      temp(i)  = eos.temp_from_gm1(gm1);
    if (set_efrac)
      efrac(i) = eos.efrac_from_gm1(gm1);
  }
  if (set_shift) {
    shift(i) >> shift_id;
    for (int d=0; d<3; d++) if (0==gsl_finite(shift_id(d))) {
        throw std::runtime_error("LoreneID: shift contains NaNs");
    }
    shift_id = -shift_id;
    shift(i) << shift_id;
  }
  if (set_lapse) { //Nothing to convert, just check for NAN.
    if (0==gsl_finite(lapse(i)))
      throw std::runtime_error("LoreneID: lapse contains NaNs");
  }
  if (set_shift && set_lapse && do_modify) {
    modify_id(coord(i), id, shift_id, lapse(i), omega);
  }
}

void loreneid::modify_id(const vec_u& x, IDBase::idvars& id, const vec_u& idshift, const double idlapse,
                 const double omega)
{
  if (id.rmd <= scale_above_density) return;
  if (old_behaviour) {
    id.vel *= scale_vel;
  }
  else {
    vec_u worb;
    worb(0)     = -x(1) * omega;
    worb(1)     = x(0) * omega;
    worb(2)     = 0;
    vec_u w     = idlapse * id.vel - idshift;
    vec_u wnew  = (scale_vel - spinup_fac) * worb + spinup_fac * w;
    id.vel      = (wnew + idshift) / idlapse;
  }
  id.florentz = id.florentz_from_vel();
  if (0==gsl_finite(id.florentz)) {
    throw std::runtime_error("LoreneID: v>c after velocity field modifications.");
  }
}


/// Set up complete initial data from the data already copied from LORENE
void loreneid::add_missing(int& crap_count, const double omega)
{
  crap_count = 0;
  idvars id;
  for (region::iterator i(r_all());i;++i) {
    idata(i)  >> id;
    complete_vars(id, crap_count);
    optional_id(i, id, omega);
    idata(i)  <<  id;
  }
}


