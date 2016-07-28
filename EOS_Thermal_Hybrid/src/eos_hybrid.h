#ifndef EOSHYBRID_H
#define EOSHYBRID_H
#include "eos_barotropic.h"
#include "eos_thermal.h"

namespace whizza {

class eos_hybrid : public eos_thermal_impl {
  EOS_Barotropic::eos_cold eos_c;
  pz_real gamma_th, gm1_th, eps_max, eps_min;


  pz_real eps_cold(pz_real rho) const;
  pz_real p_cold(pz_real rho) const;
  pz_real hm1_cold(pz_real rho) const;
  pz_real cs2_cold(pz_real rho) const;

  public:

  eos_hybrid(EOS_Barotropic::eos_1p& eos_c_,  pz_real gamma_th_,
             pz_real eps_max_, pz_real eps_min_, pz_real rho_max_, const range& rgye_);
  virtual ~eos_hybrid();

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

eos_thermal make_eos_hybrid(EOS_Barotropic::eos_1p& eos_c,
             pz_real gamma_th, pz_real eps_max, pz_real eps_min, pz_real rho_max,
             const eos_thermal::range& rgye);


}

#endif


