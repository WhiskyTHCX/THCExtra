#include "eos_extable.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>

using namespace whizza;
using namespace std;

const pz_real eos_extable::atomic_mass_unit_mev = 931.494061;

/**Constructor. The new validity ranges have to include the old ones, else
an exception is thrown. We set up the ranges in which to use the original EOS
slightly smaller than its validity range in order to avoid corner case
problems when calling the original EOS.
**/
eos_extable::eos_extable(pz_real temp_max_, pz_real rho_max_,
              eos_thermal::range rgye_, const eos_thermal& eos_)
: temp_max(temp_max_), rho_max(rho_max_), eos(eos_)
{
  rho_max_orig  = eos.range_rho().max * 0.9999999;
  temp_max_orig = eos.range_temp().max;
  rgye_orig = eos_thermal::range(eos.range_ye().min+1e-10,
                                 eos.range_ye().max-1e-10);

  if (rgye_.min >= rgye_.max) {
    throw std::runtime_error("EOS_Thermal_Extable: New electron fraction range degenerate (ye_min>=ye_max).");
  }
  if ((rgye_.min > rgye_orig.min) || (rgye_.max < rgye_orig.max)) {
    throw std::runtime_error("EOS_Thermal_Extable: New electron fraction range does not cover original one.");
  }
  if (rho_max_orig > rho_max) {
    throw std::runtime_error("EOS_Thermal_Extable: New maximum density smaller than original one.");
  }
  if (temp_max_orig > temp_max) {
    throw std::runtime_error("EOS_Thermal_Extable: New maximum temperature smaller than original one.");
  }

  range rgrho(eos.range_rho().min, rho_max);
  range rgtemp(eos.range_temp().min, temp_max);

  set_range_rho(rgrho);
  set_range_ye(rgye_);
  set_range_temp(rgtemp);
}

eos_extable::~eos_extable() {}


pz_real eos_extable::ye_bound(pz_real ye) const
{
  if (ye >= rgye_orig.max) return rgye_orig.max;
  if (ye <= rgye_orig.min) return rgye_orig.min;
  return ye;
}

pz_real eos_extable::press0(pz_real ye) const
{
  pz_real eps0 = eos.range_eps(rho_max_orig, ye).min;
  return orig().press_from_valid_rho_eps_ye(rho_max_orig, eps0, ye);
}

/**This computes the change of specific energy along curves
of zero temperature assuming that the pressure stays constant
above the maximum tabulated density, using
\f[ \delta\epsilon(\rho, Y_e)
    = P(\rho_0, \epsilon_\mathrm{min}(\rho_0, Y_e), Y_e)
      \left(\frac{1}{\rho_0} - \frac{1}{\rho} \right) \f]
**/
pz_real eos_extable::delta_eps_min(pz_real rho, pz_real ye) const
{
  return press0(ye) * (1.0/rho_max_orig - 1.0/rho);
}

/**
Compute \f$ \bar{\epsilon} = \epsilon - \delta\epsilon(\rho, Y_e) \f$
*/
pz_real eos_extable::eps_bar(pz_real rho, pz_real eps, pz_real ye) const
{
  return eps - delta_eps_min(rho, ye);
}



/**The pressure is computed using the original EOS inside its validity
range. In the extended range, we set
\f[ P = \begin{cases}
P(\rho, \epsilon_\mathrm{max}(\rho, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho \le \rho_0, \epsilon
> \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
P(\rho_0, \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
where \f$ \rho_0 \f$ is the maximum valid density of the original EOS,
and \f$ \bar{Y}_e \f$ is the electron fraction limited to the original
range.
\warning This is inconsistent with the Maxwell relations.
**/
pz_real eos_extable::press_from_valid_rho_eps_ye(const pz_real rho,
                           const pz_real eps, const pz_real ye) const
{
  pz_real rho_adj(rho), eps_adj(eps), ye_adj(ye_bound(ye));
  if (rho_adj > rho_max_orig) {
    rho_adj = rho_max_orig;
    eps_adj = eps_bar(rho, eps, ye_adj);
  }
  eps_adj = min(eos.range_eps(rho_adj, ye_adj).max, eps_adj);
  return orig().press_from_valid_rho_eps_ye(rho_adj, eps_adj, ye_adj);
}

/**The soundspeed is computed using the original EOS inside its validity
range. In the extended range, we return an incorrect soundspeed given by
\f[c_s = \begin{cases}
c_s(\rho, \epsilon_\mathrm{max}(\rho, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho \le \rho_0, \epsilon
     > \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
c_s(\rho_0, \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
where \f$ \rho_0 \f$ is the maximum valid density of the original EOS,
and \f$ \bar{Y}_e \f$ is the electron fraction limited to the original
range.
**/
pz_real eos_extable::csnd_from_valid_rho_eps_ye(const pz_real rho,
                              const pz_real eps, const pz_real ye) const
{
  pz_real rho_adj(rho), eps_adj(eps), ye_adj(ye_bound(ye));
  if (rho_adj > rho_max_orig) {
    rho_adj = rho_max_orig;
    eps_adj = eps_bar(rho, eps, ye_adj);
  }
  eps_adj = min(eos.range_eps(rho_adj, ye_adj).max, eps_adj);
  return orig().csnd_from_valid_rho_eps_ye(rho_adj, eps_adj, ye_adj);
}


/**The temperature is computed using the original EOS inside its validity
range. In the extended range, we define a fake temperature by
\f[T = \begin{cases}
T_\mathrm{max}(\rho, \bar{Y}_e) + m_B (\epsilon
                - \epsilon_\mathrm{max}(\rho, \bar{Y}_e))
 & \mbox{if } \rho \le \rho_0, \epsilon
    > \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
T(\rho_0, \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
where \f$ m_B \f$ is the atomic mass unit in MeV, \f$ \rho_0 \f$ is the
maximum valid density of the original EOS, and \f$ \bar{Y}_e \f$ is the
electron fraction limited to the original range.
**/
pz_real eos_extable::temp_from_valid_rho_eps_ye(const pz_real rho,
                              const pz_real eps, const pz_real ye) const
{
  pz_real rho_adj(rho), eps_adj(eps), ye_adj(ye_bound(ye));
  if (rho_adj > rho_max_orig) {
    rho_adj = rho_max_orig;
    eps_adj = eps_bar(rho, eps, ye_adj);
  }
  pz_real eps_max = eos.range_eps(rho_adj, ye_adj).max;
  if (eps_adj < eps_max) {
    return orig().temp_from_valid_rho_eps_ye(rho_adj, eps_adj, ye_adj);
  }
  return temp_max_orig + atomic_mass_unit_mev * (eps_adj - eps_max);
}

/**The pressure and the derivatives are computed using the original EOS
inside its validity range. In the extended range, we set
\f[P = \begin{cases}
P(\rho, \epsilon_\mathrm{max}(\rho, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho \le \rho_0, \epsilon
   > \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
P(\rho_0, \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
\f[\frac{\partial P}{\partial\rho} = \begin{cases}
\frac{\partial P}{\partial\rho}(\rho,
     \epsilon_\mathrm{max}(\rho, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho \le \rho_0, \epsilon
   > \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
-\frac{P_0(\bar{Y}_e)}{\rho^2}
 \frac{\partial P}{\partial\epsilon}(\rho,
   \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
\f[\frac{\partial P}{\partial\epsilon} = \begin{cases}
0
 & \mbox{if } \rho \le \rho_0, \epsilon
          > \epsilon_\mathrm{max}(\rho, \bar{Y}_e) \\
\frac{\partial P}{\partial\epsilon}(\rho,
            \bar{\epsilon}(\rho, \epsilon, \bar{Y}_e), \bar{Y}_e)
 & \mbox{if } \rho > \rho_0
\end{cases}
\f]
where \f$ \rho_0 \f$ is the maximum valid density of the original EOS,
and \f$ \bar{Y}_e \f$ is the electron fraction limited to the original
range.
**/
void eos_extable::press_derivs_from_valid_rho_eps_ye(pz_real& press,
           pz_real& dpdrho, pz_real& dpdeps,
           const pz_real rho, const pz_real eps, const pz_real ye) const
{
  pz_real ye_adj(ye_bound(ye));
  if (rho > rho_max_orig) {
    pz_real eps0     = eps_bar(rho, eps, ye_adj);
    pz_real eps0_max = eos.range_eps(rho_max_orig, ye_adj).max;
    if (eps0 < eps0_max) {
      orig().press_derivs_from_valid_rho_eps_ye(press, dpdrho, dpdeps,
                                            rho_max_orig, eps0, ye_adj);
      dpdrho = -dpdeps * press0(ye_adj) / (rho_max_orig*rho_max_orig);
    }
    else {
      press   = orig().press_from_valid_rho_eps_ye(rho_max_orig,
                                                      eps0_max, ye_adj);
      dpdeps  = 0;
      dpdrho  = 0;
    }
  } else {
    pz_real eps_max = eos.range_eps(rho, ye_adj).max;
    if (eps<eps_max) {
      orig().press_derivs_from_valid_rho_eps_ye(press, dpdrho, dpdeps,
                                              rho, eps, ye_adj);
    } else {
      orig().press_derivs_from_valid_rho_eps_ye(press, dpdrho, dpdeps,
                                              rho, eps_max, ye_adj);
      dpdeps=0;
    }
  }
}

/**The specific entropy in the extended range is computed by the same
extrapolation as for the pressure.
\warning This is completely inconsistent thermodynamically. Do not use
values of the entropy from the extended regions for anything.
**/
pz_real eos_extable::entropy_from_valid_rho_eps_ye(const pz_real rho,
                              const pz_real eps, const pz_real ye) const
{
  pz_real rho_adj(rho), eps_adj(eps), ye_adj(ye_bound(ye));
  if (rho_adj > rho_max_orig) {
    rho_adj = rho_max_orig;
    eps_adj = eps_bar(rho, eps, ye_adj);
  }
  eps_adj = min(eos.range_eps(rho_adj, ye_adj).max, eps_adj);
  return orig().entropy_from_valid_rho_eps_ye(rho_adj, eps_adj, ye_adj);
}

pz_real eos_extable::entropy_from_valid_rho_temp_ye(const pz_real rho,
  const pz_real temp, const pz_real ye) const
{
  const pz_real ye_adj(ye_bound(ye));
  const pz_real rho_adj = (rho < rho_max_orig) ? rho : rho_max_orig;
  const pz_real temp_adj = (temp < temp_max_orig) ? temp : temp_max_orig;

  return orig().entropy_from_valid_rho_temp_ye(rho_adj, temp_adj, ye_adj);
}

pz_real eos_extable::press_from_valid_rho_temp_ye(const pz_real rho,
  const pz_real temp, const pz_real ye) const
{
  const pz_real ye_adj(ye_bound(ye));
  const pz_real rho_adj = (rho < rho_max_orig) ? rho : rho_max_orig;
  const pz_real temp_adj = (temp < temp_max_orig) ? temp : temp_max_orig;

  return orig().press_from_valid_rho_temp_ye(rho_adj, temp_adj, ye_adj);
}

pz_real eos_extable::csnd_from_valid_rho_temp_ye(const pz_real rho,
                              const pz_real temp, const pz_real ye) const
{
  const pz_real ye_adj(ye_bound(ye));
  const pz_real rho_adj = (rho < rho_max_orig) ? rho : rho_max_orig;
  const pz_real temp_adj = (temp < temp_max_orig) ? temp : temp_max_orig;

  return orig().csnd_from_valid_rho_temp_ye(rho_adj, temp_adj, ye_adj);
}

/** The range is extended via
\f[
\epsilon(\rho, T, Y_e) = \begin{cases}
\epsilon_\mathrm{max}(\rho, Y_e) + \left( T - T_0 \right) / m_B
 & \mbox{if } \rho \le \rho_0, T \ge T_0\\
\epsilon_\mathrm{max}(\rho_0, Y_e)
+ \delta\epsilon(\rho, Y_e)
+ \left( T - T_0 \right) / m_B
 & \mbox{if } \rho > \rho_0, T \ge T_0\\
\end{cases}
\f]
where \f$T_0 \f$ is the original maximum temperature and
\f$ m_B \f$ is the atomic mass unit in MeV
*/
pz_real eos_extable::eps_from_valid_rho_temp_ye(const pz_real rho,
                            const pz_real temp, const pz_real ye) const
{
  pz_real rho_adj(rho), ye_adj(ye_bound(ye)), deps(0);

  if (rho_adj > rho_max_orig) {
    rho_adj = rho_max_orig;
    deps    = delta_eps_min(rho, ye_adj);
  }
  if (temp < temp_max_orig) {
    return deps + orig().eps_from_valid_rho_temp_ye(rho_adj, temp, ye_adj);
  }
  const pz_real eps_max_orig = eos.range_eps(rho_adj, ye_adj).max;
  return deps + eps_max_orig +  (temp - temp_max_orig) / atomic_mass_unit_mev;

}

/**The validity range of the specific energy is the following
\f[\epsilon_\mathrm{min}'(\rho, Y_e) = \begin{cases}
\epsilon_\mathrm{min}(\rho, \bar{Y}_e)
 & \mbox{if } \rho \le \rho_0\\
\epsilon_\mathrm{min}(\rho_0, \bar{Y}_e)
+ \delta\epsilon(\rho, Y_e)
 & \mbox{if } \rho > \rho_0\\
\end{cases}
\f]
\f[\epsilon_\mathrm{max}'(\rho, Y_e) = \begin{cases}
\epsilon_\mathrm{max}(\rho, \bar{Y}_e)
+ \left( T'_\mathrm{max} - T_\mathrm{max} \right) / m_B
 & \mbox{if } \rho \le \rho_0 \\
\epsilon_\mathrm{max}(\rho, \bar{Y}_e)
+ \delta\epsilon(\rho, Y_e)
+ \left( T'_\mathrm{max} - T_\mathrm{max} \right) / m_B
 & \mbox{if } \rho > \rho_0\\
\end{cases}
\f]
Here, a prime denotes the extended range and no prime the original range.
\f$ m_B \f$ is the atomic mass unit in MeV.
**/
eos_thermal_impl::range eos_extable::range_eps_from_valid_rho_ye(
                             const pz_real rho, const pz_real ye) const
{
  pz_real ye_adj(ye_bound(ye)), deps(0), rho_adj(rho);

  if (rho > rho_max_orig) {
    rho_adj = rho_max_orig;
    deps    = delta_eps_min(rho, ye_adj);
  }
  const range rgeps   = orig().range_eps_from_valid_rho_ye(rho_adj, ye_adj);
  const pz_real eps0  = rgeps.min + deps;
  const pz_real eps1  = rgeps.max + deps
                        + (temp_max - temp_max_orig) / atomic_mass_unit_mev;

  return range(eps0, eps1);
}

eos_thermal whizza::make_eos_extable(pz_real temp_max,
              pz_real rho_max, eos_thermal::range rgye)
{
  return eos_thermal(new eos_extable(temp_max, rho_max, rgye,
                       new eos_table3d()));
}


