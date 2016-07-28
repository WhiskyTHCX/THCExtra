#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <sstream>
#include <iomanip>
#include "loreneid.h"

namespace Pizza {
namespace LoreneID {

void loreneid::prep(const cGH *cctkGH)
{
  cactus_grid::prep(cctkGH);
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart();

  size_t s          = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  pz_real *v[]       = {vel, vel+s, vel+2*s};
  pz_real *glo[]    = {gxx, gxy, gyy, gxz, gyz, gzz};
  pz_real *klo[]    = {kxx, kxy, kyy, kxz, kyz, kzz};
  pz_real *beta[]  = {betax, betay, betaz};
  idata.init(m, glo, klo,  rho, eps, press, v, w_lorentz);
  // Optional variables
  if (set_lapse)
    lapse.init(m, alp);
  if (set_shift)
    shift.init(m, beta);
  if (set_temp)
    temp.init(m, temperature);
  if (set_efrac)
    efrac.init(m, Y_e);
}

cactus_single<loreneid> idlorene;

}
}//namespace Pizza

using namespace Pizza;
using namespace Pizza::LoreneID;
using namespace std;

extern "C" int LoreneID_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;

  try {
    const bool set_hydro = (string(initial_hydro)=="LoreneBNS");
    const bool set_temp  = (string(initial_temperature)=="LoreneBNS");
    const bool set_efrac = (string(initial_Y_e)=="LoreneBNS");
    const bool set_lapse = (string(initial_lapse)=="LoreneBNS");
    const bool set_shift = (string(initial_shift)=="LoreneBNS");
    const bool set_metric = (string(initial_data)=="LoreneBNS");

    if (!set_hydro) {
      error::incase(set_lapse || set_shift || set_metric || set_temp || set_efrac,
        "LoreneID does not set any quantities without initializing the hydro part as well");
      return 0;
    }

    if (modify_id) {
      error::unless(set_shift && set_lapse,
        "LoreneID can only modify initial data when shift and lapse is also set.");
    }


    CCTK_RegisterBanner("LoreneID: load Lorene initial data");

    error::unless(set_metric,
      "LoreneID always sets metric variables together with hydro");

    EOS_Barotropic::eos_1p eos  = IDBase::pizza_idbase_central::get().eos;
    idlorene    = new loreneid(eos, set_shift, set_lapse, set_temp,
                               set_efrac, modify_id, old_behaviour,
                               mass_conversion, rho_cut, scale_vel,
                               spinup_fac, scale_above_density);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" void LoreneID_Add_Missing (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  try {
    int crap_count(0);
    idlorene.prep(CCTK_PASS_CTOC).add_missing(crap_count, *lorene_bns_omega);
    if (crap_count > 0) {
      stringstream ss;
      ss << "LoreneID hydro part was crap at " << crap_count << " points";
      CCTK_WARN(1, ss.str().c_str());
    }
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

