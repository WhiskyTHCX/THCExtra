#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_Table.h"

// ----------------------------------------------------------------------------
// From CarpetRegrid2
// Get indexing information for a vector grid array
static void getvectorindex2(cGH const *const cctkGH, char const *const groupname,
                     int *const lsh) {
  assert(groupname);
  assert(lsh);

  int const gi = CCTK_GroupIndex(groupname);
  assert(gi >= 0);

  {
    int const ierr = CCTK_GrouplshGI(cctkGH, 1, lsh, gi);
    assert(!ierr);
  }

  cGroup groupdata;
  {
    int const ierr = CCTK_GroupData(gi, &groupdata);
    assert(!ierr);
  }
  assert(groupdata.vectorgroup);
  assert(groupdata.vectorlength >= 0);
  lsh[1] = groupdata.vectorlength;
}

static inline int index2(int const *const lsh, int const i, int const j) {
  assert(lsh);
  assert(lsh[0] >= 0);
  assert(lsh[1] >= 0);
  assert(i >= 0 && i < lsh[0]);
  assert(j >= 0 && j < lsh[1]);
  return i + lsh[0] * j;
}
// End CarpetRegrid2
// ----------------------------------------------------------------------------

int BHTracker_Startup(void) {
  CCTK_RegisterBanner("BHTracker: single BH tracker.");
  return 0;
}

void BHTracker_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BHTracker_Init");
  }

  if(adapt_resolution) {
    if(adapt_nlevels_min > adapt_nlevels_max) {
      CCTK_PARAMWARN("adapt_nlevels_min should not be larger than adapt_nlevels_max");
    }
    if(index_surface < 0) {
      CCTK_PARAMWARN("If adapt_resolution is yes, index_surface should also be set");
    }
  }
  if(update_surface && index_surface < 0) {
    CCTK_PARAMWARN("If update_surface is yes, index_surface should also be set");
  }
}

void BHTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BHTracker_Init");
  }

  *bh_pos_x = initial_x;
  *bh_pos_y = initial_y;
  *bh_pos_z = initial_z;
  *bh_vel_x = 0.0;
  *bh_vel_y = 0.0;
  *bh_vel_z = 0.0;
  *bh_time  = cctk_time;
}

void BHTracker_TrackBH(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  // Do not track while setting up initial data;
  // time interpolation may fail
  if(cctk_iteration == 0) {
    return;
  }

  if(verbose) {
    CCTK_INFO("BHTracker_TrackBH");
  }

  *bh_pos_x_old = *bh_pos_x;
  *bh_pos_y_old = *bh_pos_y;
  *bh_pos_z_old = *bh_pos_z;

  // Dimensions
  int const dim = 3;

  // Interpolation operator
  int const operator_handle =
      CCTK_InterpHandle ("Lagrange polynomial interpolation");
  if (operator_handle < 0) {
    CCTK_ERROR("Can't get interpolation handle");
  }

  // Interpolation parameter table
  int const order = 4;
  int const param_table_handle = Util_TableCreateFromString ("order=4");
  if (param_table_handle < 0) {
    CCTK_ERROR("Can't create parameter table");
  }

  // Interpolation coordinate system
  int const coordsys_handle = CCTK_CoordSystemHandle ("cart3d");
  if (coordsys_handle < 0) {
    CCTK_ERROR("Can't get coordinate system handle");
  }

  // Only processor 0 interpolates
  int const num_points = CCTK_MyProc(cctkGH) == 0 ? 1 : 0;

  // Interpolation coordinates
  assert (dim == 3);
  CCTK_POINTER_TO_CONST interp_coords[3];
  interp_coords[0] = bh_pos_x_old;
  interp_coords[1] = bh_pos_y_old;
  interp_coords[2] = bh_pos_z_old;

  // Number of interpolation variables
  int const num_vars = 3;

  // Interpolated variables
  assert (num_vars == 3);
  int input_array_indices[3];
  input_array_indices[0] = CCTK_VarIndex ("ADMBase::betax");
  input_array_indices[1] = CCTK_VarIndex ("ADMBase::betay");
  input_array_indices[2] = CCTK_VarIndex ("ADMBase::betaz");

  // Interpolation result types
  assert (num_vars == 3);
  CCTK_INT output_array_type_codes[3];
  output_array_type_codes[0] = CCTK_VARIABLE_REAL;
  output_array_type_codes[1] = CCTK_VARIABLE_REAL;
  output_array_type_codes[2] = CCTK_VARIABLE_REAL;

  // Interpolation result
  CCTK_REAL bh_betax;
  CCTK_REAL bh_betay;
  CCTK_REAL bh_betaz;

  assert (num_vars == 3);
  CCTK_POINTER output_arrays[3];
  output_arrays[0] = &bh_betax;
  output_arrays[1] = &bh_betay;
  output_arrays[2] = &bh_betaz;

  // Interpolate
  int ierr;
  if (CCTK_IsFunctionAliased ("InterpGridArrays")) {
    // TODO: use correct array types
    // (CCTK_POINTER[] vs. CCTK_REAL[])
    ierr = InterpGridArrays
      (cctkGH, dim,
       order,
       num_points,
       interp_coords,
       num_vars, input_array_indices,
       num_vars, output_arrays);
  } else {
    ierr = CCTK_InterpGridArrays
      (cctkGH, dim,
       operator_handle, param_table_handle, coordsys_handle,
       num_points,
       CCTK_VARIABLE_REAL,
       interp_coords,
       num_vars, input_array_indices,
       num_vars, output_array_type_codes, output_arrays);
  }
  if (ierr < 0) {
    CCTK_ERROR("Interpolation error");
  }

  if (CCTK_MyProc(cctkGH) == 0) {
    // Time evolution
    // First order time integrator
    // Michael Koppitz says this works...
    // if it doesn't, we can make it second order accurate
    CCTK_REAL const dt = cctk_time - *bh_time;
    *bh_pos_x = *bh_pos_x_old + dt * (- bh_betax);
    *bh_pos_y = *bh_pos_y_old + dt * (- bh_betay);
    *bh_pos_z = *bh_pos_z_old + dt * (- bh_betaz);
    *bh_vel_x = - bh_betax;
    *bh_vel_y = - bh_betay;
    *bh_vel_z = - bh_betaz;
  }

  // Broadcast result

  CCTK_REAL loc_local[6]; /* 3 components for location, 3 components for velocity */
  if (CCTK_MyProc(cctkGH) == 0) {
    loc_local[0] = *bh_pos_x;
    loc_local[1] = *bh_pos_y;
    loc_local[2] = *bh_pos_z;
    loc_local[3] = *bh_vel_x;
    loc_local[4] = *bh_vel_y;
    loc_local[5] = *bh_vel_z;
  } else {
    loc_local[0] = 0.0;
    loc_local[1] = 0.0;
    loc_local[2] = 0.0;
    loc_local[3] = 0.0;
    loc_local[4] = 0.0;
    loc_local[5] = 0.0;
  }

  CCTK_REAL loc_global[6]; /* 3 components for location, 3 components for velocity */

  int const handle_sum = CCTK_ReductionHandle ("sum");
  if (handle_sum < 0) {
    CCTK_ERROR("Can't get redunction handle");
  }

  int const ierr2 = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, handle_sum,
     loc_local, loc_global, 6, CCTK_VARIABLE_REAL);
  if (ierr2 < 0) {
    CCTK_ERROR("Reduction error");
  }

  *bh_time  = cctk_time;
  *bh_pos_x = loc_global[0];
  *bh_pos_y = loc_global[1];
  *bh_pos_z = loc_global[2];
  *bh_vel_x = loc_global[3];
  *bh_vel_y = loc_global[4];
  *bh_vel_z = loc_global[5];

  if(verbose) {
    char msg[BUFSIZ];
    snprintf(msg, BUFSIZ, "BH position at t = %g: (%g, %g, %g)",
        cctk_time, *bh_pos_x, *bh_pos_y, *bh_pos_z);
    CCTK_INFO(msg);
  }

  if(update_surface) {
    sf_origin_x[index_surface] = *bh_pos_x;
    sf_origin_y[index_surface] = *bh_pos_y;
    sf_origin_z[index_surface] = *bh_pos_z;
  }
}

void BHTracker_UpdateAMR(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(index_region < 0) {
    return;
  }

  if(verbose) {
    CCTK_INFO("BHTracker_UpdateAMR");
  }

  position_x[index_region] = *bh_pos_x;
  position_y[index_region] = *bh_pos_y;
  position_z[index_region] = *bh_pos_z;

  if(adapt_resolution) {
    CCTK_REAL grid_size[30];
    int lsh[2];
    getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);
    for(int rl = 0; rl < 30; ++rl) {
      int const ind = index2(lsh, rl, index_region);
      grid_size[rl] = radius[ind];
    }

    CCTK_REAL const bh_radius = sf_min_radius[index_surface];
    if(verbose) {
      char msg[BUFSIZ];
      snprintf(msg, BUFSIZ, "BH min radius t = %g: %g",
          cctk_time, bh_radius);
      CCTK_INFO(msg);
    }

    int const nl = num_levels[index_region];
    if(nl < 1 || nl > 30) {
      return;
    }

    if(bh_radius < adapt_factor_refine*grid_size[nl-1] &&
        nl < adapt_nlevels_max - 1) {
      num_levels[index_region]++;
      if(verbose) {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Refine grid: num_levels[%d]: %d -> %d",
            index_region, num_levels[index_region]-1, num_levels[index_region]);
        CCTK_INFO(msg);
      }
    }
    else if(bh_radius > adapt_factor_coarse*grid_size[nl-1]
        && nl > adapt_nlevels_min + 1) {
      num_levels[index_region]--;
      if(verbose) {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Coarsen grid: num_levels[%d]: %d -> %d",
            index_region, num_levels[index_region]+1, num_levels[index_region]);
        CCTK_INFO(msg);
      }
    }
  }

  *bh_num_levels = num_levels[index_region];
}
