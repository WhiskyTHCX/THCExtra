#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "global_eos_thermal.h"
#include "eos_idealgas.h"

using namespace whizza;

extern "C" int EOS_Idealgas_Register()
{
  DECLARE_CCTK_PARAMETERS;

  try {
    eos_thermal::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max), rgye(ye_min, ye_max);
    eos_thermal eos(new eos_idealgas(index_n, particle_mass, rgeps, rgrho, rgye));
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


