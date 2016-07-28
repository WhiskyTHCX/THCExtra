#include <cctk.h>

int NeutrinoEmissionImpl(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * R_nue,
        CCTK_REAL * R_nua,
        CCTK_REAL * R_nux,
        CCTK_REAL * Q_nue,
        CCTK_REAL * Q_nua,
        CCTK_REAL * Q_nux);

int NeutrinoOpacityImpl(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * kappa_0_nue,
        CCTK_REAL * kappa_0_nua,
        CCTK_REAL * kappa_0_nux,
        CCTK_REAL * kappa_1_nue,
        CCTK_REAL * kappa_1_nua,
        CCTK_REAL * kappa_1_nux);

int NeutrinoAbsorptionRateImpl(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * abs_0_nue,
	CCTK_REAL * abs_0_nua,
        CCTK_REAL * abs_0_nux,
        CCTK_REAL * abs_1_nue,
        CCTK_REAL * abs_1_nua,
        CCTK_REAL * abs_1_nux);
