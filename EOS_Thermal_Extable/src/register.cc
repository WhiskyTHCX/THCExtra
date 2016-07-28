#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "global_eos_thermal.h"
#include "the_eos_table3d.h"
#include "eos_extable.h"
#include <stdexcept>

using namespace whizza;

extern "C" int EOS_Extable_Register()
{
  DECLARE_CCTK_PARAMETERS;

  try {
    eos_thermal eos_tab = EOS_Thermal_Table3d::parfile_eos::get();
    eos_thermal::range rgye = eos_tab.range_ye();
    if (extend_ye) {
      if (ye_min >= ye_max) {
        throw std::runtime_error("Extended electron fraction range degenerate (ye_min>=ye_max).");
      }
      if ((ye_min >= rgye.min) || (ye_max <= rgye.max)) {
        throw std::runtime_error("New electron fraction range does not cover original one.");
      }
      rgye = eos_thermal::range(ye_min, ye_max);
    }
    pz_real rgrho_max = eos_tab.range_rho().max;
    if (extend_rho) {
      if (rgrho_max >= rho_max) {
        throw std::runtime_error("New maximum density smaller than original one.");
      }
      rgrho_max = rho_max;
    }

    pz_real rgtemp_max = eos_tab.range_temp().max;
    if (rgtemp_max >= temp_max) {
      throw std::runtime_error("New maximum temperature smaller than original one.");
    }
    rgtemp_max = temp_max;

    eos_thermal eos(new eos_extable(rgtemp_max, rgrho_max, rgye,
                                    eos_tab));
    global_eos_thermal::set_eos(eos);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


