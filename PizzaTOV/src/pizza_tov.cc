#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "pizza_central.h"
#include "pizza_tov.h"

namespace Pizza {
namespace TOV {

void ptov::prep(const cGH *cctkGH)
{
  cactus_grid::prep(cctkGH);
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart();

  size_t s          = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  pz_real *v[]       = {vel, vel+s, vel+2*s};
  pz_real *glo[]    = {gxx, gxy, gyy, gxz, gyz, gzz};
  pz_real *klo[]    = {kxx, kxy, kyy, kxz, kyz, kzz};
  pz_real *beta[]  = {betax, betay, betaz};
  idata.init(m, glo, klo, rho, eps, press, v, w_lorentz);
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

cactus_single<ptov> pizzatov;

}//namespace Pizza
}

using namespace Pizza;
using namespace TOV;
using namespace std;

extern "C" int PizzaTOV_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;

  try {
    const bool set_hydro = (string(initial_hydro)=="PizzaTOV");
    const bool set_temp  = (string(initial_temperature)=="PizzaTOV");
    const bool set_efrac = (string(initial_Y_e)=="PizzaTOV");
    const bool set_lapse = (string(initial_lapse)=="PizzaTOV");
    const bool set_shift = (string(initial_shift)=="PizzaTOV");
    const bool set_metric = (string(initial_data)=="PizzaTOV");

    if (!set_hydro) {
      error::incase(set_lapse || set_shift || set_metric || set_temp || set_efrac,
        "PizzaTOV does not set any quantities without initializing the hydro part as well");
      return 0;
    }

    CCTK_RegisterBanner("PizzaTOV: TOV star initial data");

    error::unless(set_metric,
      "PizzaTOV always sets metric variables together with hydro");

    if(CCTK_Equals(idtype, "TOV")) {
      units u     = Pizza::Base::pizza_base_central::get().internal_units;
      eos_1p eos  = Pizza::IDBase::pizza_idbase_central::get().eos;
      pizzatov=new ptov(
        eos,
        star_crmd / u.density(), star_gm1_cut,
        pert_amp, pert_l, pert_m, pert_n, kick_amp,
        set_shift, set_lapse, set_temp, set_efrac, out_dir
      );
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" void PizzaTOV_Initial_Data (CCTK_ARGUMENTS)
{
  try {
    pizzatov.prep(CCTK_PASS_CTOC).initial_data();
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

extern "C" void PizzaTOV_Perturb_Data (CCTK_ARGUMENTS)
{
  try {
    pizzatov.prep(CCTK_PASS_CTOC).perturb();
    pizzatov.prep(CCTK_PASS_CTOC).kick();
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}
