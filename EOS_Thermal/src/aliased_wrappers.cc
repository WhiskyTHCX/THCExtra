#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "global_eos.h"
#include <sstream>
#include <iomanip>

using namespace whizza;

extern "C" {

CCTK_INT Wrapper_Press_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *press, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    *press                 = eos.press_from_rho_eps_ye(rho, eps, ye, eos_errs);
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

CCTK_INT Wrapper_Press_Derivs_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *press, CCTK_REAL *dpdrho, CCTK_REAL *dpdeps, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    eos.press_derivs_from_rho_eps_ye(*press, *dpdrho, *dpdeps, rho, eps, ye, eos_errs);
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

CCTK_INT Wrapper_Entropy_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *entropy, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    *entropy               = eos.entropy_from_rho_eps_ye(rho, eps, ye, eos_errs);
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

CCTK_INT Wrapper_Entropy_From_Rho_Temp_Ye(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye,
        CCTK_REAL *entropy, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    *entropy               = eos.entropy_from_rho_temp_ye(rho, temp, ye, eos_errs);
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

CCTK_INT Wrapper_Press_Csnd2_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *press, CCTK_REAL *csnd2, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    CCTK_REAL csnd(0);
    const eos_thermal& eos = global_eos_thermal::get_eos();
    eos.press_csnd_from_rho_eps_ye(*press, csnd, rho, eps, ye, eos_errs);
    *csnd2 = csnd*csnd;
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


CCTK_INT Wrapper_Csnd2_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *csnd2, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    CCTK_REAL csnd = eos.csnd_from_rho_eps_ye(rho, eps, ye, eos_errs);
    *csnd2         = csnd*csnd;
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

CCTK_INT Wrapper_Temp_From_Rho_Eps_Ye(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye,
        CCTK_REAL *temp, CCTK_INT warn)
{
  try {
    eos_thermal::status eos_errs;
    const eos_thermal& eos = global_eos_thermal::get_eos();
    *temp                  = eos.temp_from_rho_eps_ye(rho, eps, ye, eos_errs);
    if (eos_errs.failed) {
      if (warn) {
        CCTK_WARN(1, eos_errs.err_msg.c_str());
      }
      return -1;
    }
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

void Wrapper_EOS_Warning(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_INT lvl)
{
  std::stringstream ss;
  ss << "Trouble with EOS at position ("<< x << ", " << y << ", " << z
     << "), on ref. level " << lvl;
  CCTK_WARN(1, ss.str().c_str());
}

CCTK_INT EosValidateRhoEpsYe_wrapper(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye)
{
  return global_eos_thermal::get_eos().is_rho_eps_ye_valid(rho, eps, ye) ? 0 : -1;
}

CCTK_INT EosValidateRhoYe_wrapper(CCTK_REAL rho, CCTK_REAL ye)
{
  const eos_thermal& eos = global_eos_thermal::get_eos();
  return (eos.is_rho_valid(rho) && eos.is_ye_valid(ye))  ? 0 : -1;
}


void EosRangeEps_wrapper(CCTK_REAL rho, CCTK_REAL ye, CCTK_REAL* epsmin, CCTK_REAL* epsmax)
{
  eos_thermal::range rgeps = global_eos_thermal::get_eos().range_eps(rho, ye);
  *epsmin      = rgeps.min;
  *epsmax      = rgeps.max;
}


CCTK_INT equationOfState_wrapper(const CCTK_REAL* independent, CCTK_REAL* dependent, CCTK_INT flag)
{
  CCTK_WARN(0, "equationofstate interface is deprecated");
  return 0;
}
} //extern "C"


