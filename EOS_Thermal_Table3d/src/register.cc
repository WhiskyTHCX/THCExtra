#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "global_eos_thermal.h"
#include "eos_table3d.h"
#include "register.h"
#include <string>

namespace EOS_Thermal_Table3d {

whizza::eos_thermal parfile_eos::eos;
bool parfile_eos::initialized = false;

whizza::eos_thermal parfile_eos::get()
{
  DECLARE_CCTK_PARAMETERS;
  if (!initialized) {
    const std::string eos_fqfn = std::string(eos_db_loc) + "/" + eos_folder
                            + "/" + eos_filename;
    whizza::eos_table3d::init_global_eos(eos_fqfn);
    eos = whizza::eos_thermal(new whizza::eos_table3d());
    initialized = true;
  }
  return eos;
}

}


extern "C" int EOS_Table3D_Register(void)
{
  CCTK_RegisterBanner("EOS_Thermal_Table3d: Tabulated hot EOS");
  try {
    whizza::global_eos_thermal::set_eos(EOS_Thermal_Table3d::parfile_eos::get());
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


