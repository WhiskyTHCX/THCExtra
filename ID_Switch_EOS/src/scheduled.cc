#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <sstream>
#include <iomanip>
#include "id_switch_eos.h"
#include "global_eos_thermal.h"

namespace Pizza {
namespace ID_Switch_EOS {

void idswitcheos::prep(const cGH *cctkGH)
{
  cactus_grid::prep(cctkGH);
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart();

  //size_t s          = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  //pz_real *v[]       = {vel, vel+s, vel+2*s};

  gf_rmd.init(m, rho);
  gf_sed.init(m, eps);
  gf_efrac.init(m, Y_e);
  gf_press.init(m, press);
  if (sync_eps_temp)
    gf_temp.init(m, temperature);
  if (set_entropy)
    gf_entropy.init(m, entropy);
}

cactus_single<idswitcheos> idswitcheos_single;

}
}//namespace Pizza

using namespace Pizza;
using namespace Pizza::ID_Switch_EOS;
using namespace std;

extern "C" int ID_Switch_EOS_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;

  try {
    CCTK_RegisterBanner("ID_Switch_EOS: adapt initial data to evolution EOS");
    const bool set_entropy = (string(initial_entropy) == "IDSwitchEOS");
    const whizza::eos_thermal eos(whizza::global_eos_thermal::get_eos());
    idswitcheos_single = new idswitcheos(eos, sync_eps_temp,
                               temp_from_eps, limit_efrac, set_entropy);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" void ID_Switch_EOS_Adapt (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  try {
    idswitcheos_single.prep(CCTK_PASS_CTOC).adapt();
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

