#include <float.h>

#include <cctk_Parameters.h>

#include "rates.h"

int NeutrinoEmissionFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * R_nue,
        CCTK_REAL * R_nua,
        CCTK_REAL * R_nux,
        CCTK_REAL * Q_nue,
        CCTK_REAL * Q_nua,
        CCTK_REAL * Q_nux) {
    DECLARE_CCTK_PARAMETERS

    *R_nue = rho * eta_nue;
    *R_nua = rho * eta_nua;
    *R_nux = rho * eta_nux;
    *Q_nue = rho * eta_nue;
    *Q_nua = rho * eta_nua;
    *Q_nux = rho * eta_nux;

    return 0;
}

int NeutrinoOpacityFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * kappa_0_nue,
        CCTK_REAL * kappa_0_nua,
        CCTK_REAL * kappa_0_nux,
        CCTK_REAL * kappa_1_nue,
        CCTK_REAL * kappa_1_nua,
        CCTK_REAL * kappa_1_nux) {
    DECLARE_CCTK_PARAMETERS

    *kappa_0_nue = rho * (kappa_scat_nue + kappa_abs_nue);
    *kappa_0_nua = rho * (kappa_scat_nua + kappa_abs_nua);
    *kappa_0_nux = rho * (kappa_scat_nux);
    *kappa_1_nue = rho * (kappa_scat_nue + kappa_abs_nue);
    *kappa_1_nua = rho * (kappa_scat_nua + kappa_abs_nua);
    *kappa_1_nux = rho * (kappa_scat_nux);

    return 0;
}

int NeutrinoAbsorptionRateFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * abs_0_nue,
        CCTK_REAL * abs_0_nua,
        CCTK_REAL * abs_0_nux,
        CCTK_REAL * abs_1_nue,
        CCTK_REAL * abs_1_nua,
        CCTK_REAL * abs_1_nux) {
    DECLARE_CCTK_PARAMETERS

    *abs_0_nue = rho * kappa_abs_nue;
    *abs_0_nua = rho * kappa_abs_nua;
    *abs_0_nux = 0.0e0;
    *abs_1_nue = rho * kappa_abs_nue;
    *abs_1_nua = rho * kappa_abs_nua;
    *abs_1_nux = 0.0e0;

    return 0;
}

int NeutrinoDensityFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * num_nue,
        CCTK_REAL * num_nua,
        CCTK_REAL * num_nux,
        CCTK_REAL * ene_nue,
        CCTK_REAL * ene_nua,
        CCTK_REAL * ene_nux) {
    DECLARE_CCTK_PARAMETERS

    if(rho*kappa_abs_nue > FLT_EPSILON*eta_nue) {
        *num_nue = eta_nue/(rho*kappa_abs_nue);
        *ene_nue = eta_nue/(rho*kappa_abs_nue);
    }
    else {
        *num_nue = 1.0;
        *ene_nue = 1.0;
    }

    if(rho*kappa_abs_nua > FLT_EPSILON*eta_nua) {
        *num_nua = eta_nua/(rho*kappa_abs_nua);
        *ene_nua = eta_nua/(rho*kappa_abs_nua);
    }
    else {
        *num_nua = 1.0;
        *ene_nua = 1.0;
    }

    *num_nux = 1.0;
    *ene_nux = 1.0;

    return 0;
}

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
        CCTK_REAL * ene_nux_eq) {
    *temp_eq = temp;
    *ye_eq = Y_e;
    return NeutrinoDensityFake(rho, *temp_eq, *ye_eq,
            num_nue_eq, num_nua_eq, num_nux_eq,
            ene_nue_eq, ene_nua_eq, ene_nux_eq);
}

int AverageAtomicMassFake(
        CCTK_REAL const rho,
        CCTK_REAL const temp,
        CCTK_REAL const Y_e,
        CCTK_REAL * Abar) {
    DECLARE_CCTK_PARAMETERS

    *Abar = avg_atomic_mass;

    return 0;
}

CCTK_REAL AverageBaryonMassFake() {
    return 1.0;
}
