#include <cctk.h>

int NeutrinoEmissionFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * R_nue,
        CCTK_REAL * R_nua,
        CCTK_REAL * R_nux,
        CCTK_REAL * Q_nue,
        CCTK_REAL * Q_nua,
        CCTK_REAL * Q_nux);

int NeutrinoOpacityFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * kappa_0_nue,
        CCTK_REAL * kappa_0_nua,
        CCTK_REAL * kappa_0_nux,
        CCTK_REAL * kappa_1_nue,
        CCTK_REAL * kappa_1_nua,
        CCTK_REAL * kappa_1_nux);

int NeutrinoAbsorptionRateFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * abs_0_nue,
        CCTK_REAL * abs_0_nua,
        CCTK_REAL * abs_0_nux,
        CCTK_REAL * abs_1_nue,
        CCTK_REAL * abs_1_nua,
        CCTK_REAL * abs_1_nux);

int NeutrinoDensityFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * num_nue,
        CCTK_REAL * num_nua,
        CCTK_REAL * num_nux,
        CCTK_REAL * ene_nue,
        CCTK_REAL * ene_nua,
        CCTK_REAL * ene_nux);

int WeakEquilibriumFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL const num_nue,
        CCTK_REAL const num_nua,
        CCTK_REAL const num_nux,
        CCTK_REAL const ene_nue,
        CCTK_REAL const ene_nua,
        CCTK_REAL const ene_nux,
        CCTK_REAL * temp_eq,
        CCTK_REAL * ye_eq,
        CCTK_REAL * num_nue_eq,
        CCTK_REAL * num_nua_eq,
        CCTK_REAL * num_nux_eq,
        CCTK_REAL * ene_nue_eq,
        CCTK_REAL * ene_nua_eq,
        CCTK_REAL * ene_nux_eq);

int AverageAtomicMassFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * Abar);

CCTK_REAL AverageBaryonMassFake();
