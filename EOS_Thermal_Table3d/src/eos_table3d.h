#ifndef EOSTABLE3D_H
#define EOSTABLE3D_H

#include "eos_thermal.h"
#include <string>

namespace whizza {

class eos_table3d : public eos_thermal_impl {
  static bool global_eos_initialized;
  public:
  static void init_global_eos(std::string name);
  eos_table3d();
  virtual ~eos_table3d();

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

  virtual pz_real eps_from_valid_rho_temp_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
    const pz_real temp,    ///<Temperature \f$ T \f$
    const pz_real ye       ///<Electron fraction \f$ Y_e \f$
  ) const;

  virtual range range_eps_from_valid_rho_ye(
    const pz_real rho,     ///<Rest mass density  \f$ \rho \f$
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

};

eos_thermal make_eos_table3d();

}

#endif


