#ifndef EOSEXTABLE_H
#define EOSEXTABLE_H

#include <boost/shared_ptr.hpp>
#include "eos_thermal.h"
#include "eos_table3d.h"

namespace whizza {

///Class representing an EOS extending the tabulated one in temperature
/**
The goal of this class is to allow continuing a simulation in case at some
points considered unimportant the EOS validy range is exceeded. For this,
the EOS is extended analytically in a simplistic unphysical manner.
In the extended validity regime, not even thermodynamic consistency holds.

Internally, this eos stores an EOS of type EOS_Thermal_Table3d. Inside the
validity range of the latter, calls are just forwarded to it. For higher
specific energies, the pressure, soundspeed, and entropy stay constant for a
given density and electron fraction. A fake temperature is defined by
increasing the temperature proportional to the increase in specific energy
similar to the ideal gas law.
The valid range of the electron fraction can be increased as well, using
zeroth order extrapolation.

To extend the density range, we first assume that the pressure stays constant
along the zero temperature curve for densities above the maximum tabulated one.
We can then integrate the adiabaticity conditions analytically to get the
internal energy along
\f$ \epsilon_0(\rho, Y_e) = \epsilon(\rho, T=0, Y_e) \f$.
We then construct the pressure
for larger energy densities by assuming it stays constant along curves
\f$ \epsilon_0(\rho, Y_e) + \Delta\epsilon \f$ parametrized by a constant shift
\f$ \Delta\epsilon \f$ with respect to the zero temperature curve.
This assumption has no physical motivation, the only goal is to smoothly extend
the EOS without violating causality or energy condition.
**/
class eos_extable : public eos_thermal_impl {
  static const pz_real atomic_mass_unit_mev;  ///<Atomic mass unit in MeV
  pz_real temp_max;                           ///<Extended temperature range
  pz_real rho_max_orig;                       ///<Maximum density of original EOS
  pz_real temp_max_orig;                      ///<Max temperature of original EOS
  pz_real rho_max;                            ///<Extend density range to this value
  eos_thermal::range rgye_orig;               ///<Electron fraction range of original EOS
  eos_thermal eos;                            ///<The tabulated EOS

  ///Return electron fraction limited to the original EOSs validity range.
  pz_real ye_bound(pz_real ye) const;
  ///Pressure at maximum tabulated density, zero temperature.
  pz_real press0(pz_real ye) const;
  pz_real delta_eps_min(pz_real rho, pz_real ye) const;
  pz_real eps_bar(pz_real rho, pz_real eps, pz_real ye) const;

  /// Tabulated EOS implementation reference.
  const eos_thermal_impl& orig() const {return eos.implementation();}
  public:

  eos_extable(pz_real temp_max_,        ///< New maximum temperature.
              pz_real rho_max_,         ///< New maximum density.
              eos_thermal::range rgye_, ///<New electron fraction range.
              const eos_thermal& eos_   ///<Original tabulated EOS.
             );

  virtual ~eos_extable();

  protected:

  /// Compute pressure
  virtual pz_real press_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute speed of sound.
  virtual pz_real csnd_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute temperature
  virtual pz_real temp_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute pressure and derivatives.
  virtual void press_derivs_from_valid_rho_eps_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& dpdrho,       ///<Partial derivative \f$ \frac{\partial P}{\partial \rho} \f$
    pz_real& dpdeps,       ///<Partial derivative \f$ \frac{\partial P}{\partial \epsilon} \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute specific entropy.
  virtual pz_real entropy_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute specific energy
  virtual pz_real eps_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$ in MeV
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  /// Compute specific entropy.
  virtual pz_real entropy_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real press_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real csnd_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;


  /// Valid range for specifc energy at given density and electron fraction.
  virtual range range_eps_from_valid_rho_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

};

eos_thermal make_eos_extable(pz_real temp_max, pz_real rho_max,
              eos_thermal::range rgye);

}

#endif


