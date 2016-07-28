#ifndef EOSIDEALGAS_H
#define EOSIDEALGAS_H

#include "eos_thermal.h"

namespace whizza {

class eos_idealgas : public eos_thermal_impl {
  pz_real gamma, gm1, temp_over_eps;
  range rgeps;
  public:

  eos_idealgas(pz_real n_, pz_real umass_, const range& rgeps_, const range& rgrho_, const range& rgye_);
  virtual ~eos_idealgas();

  virtual pz_real press_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real csnd_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real temp_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual void press_derivs_from_valid_rho_eps_ye(
    pz_real& press,        ///<Pressure \f$ P \f$
    pz_real& dpdrho,       ///<Partial derivative \f$ \frac{\partial P}{\partial \rho} \f$
    pz_real& dpdeps,       ///<Partial derivative \f$ \frac{\partial P}{\partial \epsilon} \f$
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real entropy_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real entropy_from_valid_rho_eps_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real eps,     ///<Specific internal energy \f$ \epsilon \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual pz_real eps_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$ in MeV
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual range range_eps_from_valid_rho_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;
};

eos_thermal make_eos_idealgas(pz_real n, pz_real umass,
    const eos_thermal::range& rgeps,
    const eos_thermal::range& rgrho,
    const eos_thermal::range& rgye);

}

#endif


