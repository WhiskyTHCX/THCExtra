#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "global_eos_thermal.h"
#include "pizza_unitconv.h"
#include "eos_hybrid.h"
#include "pizza_eos_barotropic_file.h"

using namespace whizza;
using namespace EOS_Barotropic;
using Pizza::units;

extern "C" int EOS_Hybrid_Register()
{
  DECLARE_CCTK_PARAMETERS;

  try {
    eos_thermal::range rgye(ye_min, ye_max);
    eos_1p eos_cold = load_eos_1p(eos_cold_file, units::geom_solar());
    eos_thermal eos(new eos_hybrid(eos_cold,  gamma_thermal, eps_max, eps_min, rho_max, rgye));
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


