#include "cctk.h"
#include "cctk_Parameters.h"
#include "eos_table3d.h"

namespace EOS_Thermal_Table3d {

/// Singleton storing the EOS specified in the simulations parfile
class parfile_eos {
  static whizza::eos_thermal eos; ///< The EOS.
  static bool initialized;        ///< If it has already been created.
  public:
  ///Get the EOS, creating it if not already existent.
  static whizza::eos_thermal get();
};


}

