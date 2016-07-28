#include "pizzacactus.h"
//#include "pizza_idbase.h"
#include "eos_thermal.h"

namespace Pizza {
namespace ID_Switch_EOS {

class idswitcheos : cactus_grid {
  scalar_pc gf_rmd, gf_sed, gf_temp, gf_efrac, gf_press, gf_entropy;
  const whizza::eos_thermal eos;
  bool sync_eps_temp, temp_from_eps, limit_efrac, set_entropy;

  public:
  idswitcheos(
    const whizza::eos_thermal& eos_,
    bool sync_eps_temp_, bool temp_from_eps_,
    bool limit_efrac_, bool set_entropy_
  );

  void prep(const cGH *cgh);
  void adapt();

};

}
}

