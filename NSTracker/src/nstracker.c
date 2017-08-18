#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

int NSTracker_Startup(void) {
    CCTK_RegisterBanner("NSTracker: single NS tracker.");
    return 0;
}

void NSTracker_TrackNS(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if((cctk_iteration % analyze_every) != 0) {
        return;
    }
    if(cctk_levfac[0] != (1 << analyze_reflevel)) {
        return;
    }
    if(verbose) {
        CCTK_INFO("NSTracker_TrackNS");
    }

    CCTK_REAL rho_sum   = 0;
    CCTK_REAL rho_pos_x = 0;
    CCTK_REAL rho_pos_y = 0;
    CCTK_REAL rho_pos_z = 0;

#pragma omp parallel for reduction(+: rho_sum, rho_pos_x, rho_pos_y, rho_pos_z) collapse(2)
    for(int i = cctk_nghostzones[0]; i < cctk_lsh[0] - cctk_nghostzones[0]; ++i)
    for(int j = cctk_nghostzones[1]; j < cctk_lsh[1] - cctk_nghostzones[1]; ++j)
    for(int k = cctk_nghostzones[2]; k < cctk_lsh[2] - cctk_nghostzones[2]; ++k)
    {
        int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        rho_sum   += rho[ijk];
        rho_pos_x += rho[ijk]*x[ijk];
        rho_pos_y += rho[ijk]*y[ijk];
        rho_pos_z += rho[ijk]*z[ijk];
    }

    CCTK_REAL r_rho_sum;
    CCTK_REAL r_rho_pos_x, r_rho_pos_y, r_rho_pos_z;

    int reduction_handle = CCTK_ReductionArrayHandle("sum");
    if(reduction_handle < 0) {
        CCTK_ERROR("Could not get sum reduction handle");
    }

    int ierr;
    ierr = CCTK_ReduceLocScalar(cctkGH, -1, reduction_handle,
            &rho_sum, &r_rho_sum, CCTK_VARIABLE_REAL);
    if(ierr < 0) {
        CCTK_ERROR("CCTK_ReduceLocalArrays failed!");
    }

    ierr = CCTK_ReduceLocScalar(cctkGH, -1, reduction_handle,
            &rho_pos_x, &r_rho_pos_x, CCTK_VARIABLE_REAL);
    if(ierr < 0) {
        CCTK_ERROR("CCTK_ReduceLocalArrays failed!");
    }
    *ns_pos_x = r_rho_pos_x/r_rho_sum;

    ierr = CCTK_ReduceLocScalar(cctkGH, -1, reduction_handle,
            &rho_pos_y, &r_rho_pos_y, CCTK_VARIABLE_REAL);
    if(ierr < 0) {
        CCTK_ERROR("CCTK_ReduceLocalArrays failed!");
    }
    *ns_pos_y = r_rho_pos_y/r_rho_sum;

    if(reflection_z) {
        *ns_pos_z = 0.0;
    }
    else {
        ierr = CCTK_ReduceLocScalar(cctkGH, -1, reduction_handle,
                &rho_pos_z, &r_rho_pos_z, CCTK_VARIABLE_REAL);
        if(ierr < 0) {
            CCTK_ERROR("CCTK_ReduceLocalArrays failed!");
        }
        *ns_pos_z = r_rho_pos_z/r_rho_sum;
    }

    if(verbose) {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "NS position at t = %g: (%g, %g, %g)",
                cctk_time, *ns_pos_x, *ns_pos_y, *ns_pos_z);
        CCTK_INFO(msg);
    }

    if(surface_index >= 0) {
        int const sn = surface_index;
        if(verbose) {
            char msg[BUFSIZ];
            snprintf(msg, BUFSIZ, "Updating spherical surface index=%d", sn);
            CCTK_INFO(msg);
        }

        sf_centroid_x[sn] = *ns_pos_x;
        sf_centroid_y[sn] = *ns_pos_y;
        sf_centroid_z[sn] = *ns_pos_z;

        sf_active[sn] = 1;
        sf_valid[sn] = 1;
    }
}
