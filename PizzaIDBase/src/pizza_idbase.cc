#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "idbase.h"
#include "pizza_central.h"
#include "pizza_unitconv.h"
#include "pizza_eos_barotropic_file.h"
#include <algorithm>
#include <stdexcept>

using namespace Pizza;
using namespace IDBase;
using namespace Base;

using std::copy;

extern "C" int PizzaIDBase_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_RegisterBanner("PizzaIDBase: pizza initial data infrastructure");
  try {
    units u = pizza_base_central::get().internal_units;
    pizza_idbase_central::init(EOS_Barotropic::load_eos_1p(eos_file, u));
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }

  return 0;
}

extern "C" void PizzaIDBase_Save_Unpert (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_OutputVarAsByMethod(cctkGH,"HydroBase::rho","IOHDF5","rho_unpert");
  CCTK_OutputVarAsByMethod(cctkGH,"HydroBase::rho","IOASCII_1D","rho_unpert");
}

void copy_safe(const double*src, size_t s, double* dst)
{
  if ((0==src) || (0==dst))
    throw std::runtime_error("PizzaIDBase: Cactus passed NULL pointer to active timelevels.");
  copy(src, src+s, dst);
}

void copy_tlvls(const double *l0, double* l1, double* l2, size_t s, int nlvl)
{
  if (nlvl>1) {
    copy_safe(l0, s, l1);
    if (nlvl>2) copy_safe(l0,s,l2);
  }
}

extern "C" void PizzaIDBase_Timelevels (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  try {
    int s=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

    int nl_metric = CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric");
    copy_tlvls(gxx, gxx_p, gxx_p_p, s, nl_metric);
    copy_tlvls(gxy, gxy_p, gxy_p_p, s, nl_metric);
    copy_tlvls(gxz, gxz_p, gxz_p_p, s, nl_metric);
    copy_tlvls(gyy, gyy_p, gyy_p_p, s, nl_metric);
    copy_tlvls(gyz, gyz_p, gyz_p_p, s, nl_metric);
    copy_tlvls(gzz, gzz_p, gzz_p_p, s, nl_metric);

    int nl_curv = CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv");
    copy_tlvls(kxx, kxx_p, kxx_p_p, s, nl_curv);
    copy_tlvls(kxy, kxy_p, kxy_p_p, s, nl_curv);
    copy_tlvls(kxz, kxz_p, kxz_p_p, s, nl_curv);
    copy_tlvls(kyy, kyy_p, kyy_p_p, s, nl_curv);
    copy_tlvls(kyz, kyz_p, kyz_p_p, s, nl_curv);
    copy_tlvls(kzz, kzz_p, kzz_p_p, s, nl_curv);

    copy_tlvls(alp, alp_p, alp_p_p, s, CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse"));

    int nl_beta = CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift");
    copy_tlvls(betax, betax_p, betax_p_p, s, nl_beta);
    copy_tlvls(betay, betay_p, betay_p_p, s, nl_beta);
    copy_tlvls(betaz, betaz_p, betaz_p_p, s, nl_beta);

    copy_tlvls(rho, rho_p, rho_p_p, s, CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho"));

    copy_tlvls(vel,vel_p, vel_p_p, 3*s, CCTK_ActiveTimeLevels(cctkGH, "HydroBase::vel"));

    copy_tlvls(w_lorentz, w_lorentz_p, w_lorentz_p_p, s,
               CCTK_ActiveTimeLevels(cctkGH, "HydroBase::w_lorentz"));

    copy_tlvls(eps, eps_p, eps_p_p, s, CCTK_ActiveTimeLevels(cctkGH, "HydroBase::eps"));

    copy_tlvls(press, press_p, press_p_p, s,
                CCTK_ActiveTimeLevels(cctkGH, "HydroBase::press"));
    copy_tlvls(Y_e, Y_e_p, Y_e_p_p, s, CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e"));
    copy_tlvls(temperature, temperature_p, temperature_p_p, s,
                CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature"));

  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
}

