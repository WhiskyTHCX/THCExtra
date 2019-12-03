#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_Table.h"

#define NBH (10)

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

int BHTrackerMany_Startup(void) {
  CCTK_RegisterBanner("BHTrackerMany: BH tracker.");
  return 0;
}

void BHTrackerMany_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BHTrackerMany_Init");
  }

  for(int i = 0; i < NBH; ++i) {
    if(adapt_resolution[i]) {
      if(adapt_nlevels_min[i] > adapt_nlevels_max[i]) {
        CCTK_PARAMWARN("adapt_nlevels_min should not be larger than adapt_nlevels_max");
      }
      if(index_surface[i] < 0) {
        CCTK_PARAMWARN("If adapt_resolution is yes, index_surface should also be set");
      }
      if(index_region[i] < 0) {
        CCTK_PARAMWARN("If adapt_resolution is yes, index_region should also be set");
      }
    }
    if(update_surface[i] && index_surface[i] < 0) {
      CCTK_PARAMWARN("If update_surface is yes, index_surface should also be set");
    }
  }
}

void BHTrackerMany_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BHTrackerMany_Init");
  }

  for(int i = 0; i < NBH; ++i) {
    bh_time[i]  = cctk_time;
    bh_pos_x[i] = initial_x[i];
    bh_pos_y[i] = initial_y[i];
    bh_pos_z[i] = initial_z[i];
    bh_vel_x[i] = 0.0;
    bh_vel_y[i] = 0.0;
    bh_vel_z[i] = 0.0;
  }
}

void BHTrackerMany_TrackBH(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  // Do not track while setting up initial data;
  // time interpolation may fail
  if(cctk_iteration == 0) {
    return;
  }

  if(verbose) {
    CCTK_INFO("BHTrackerMany_TrackBH");
  }

  // Cycle timelevels
  for(int i = 0; i < NBH; ++i) {
    bh_pos_x_old[i] = bh_pos_x[i];
    bh_pos_y_old[i] = bh_pos_y[i];
    bh_pos_z_old[i] = bh_pos_z[i];
  }

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
  int const num_points = CCTK_MyProc(cctkGH) == 0 ? NBH : 0;

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
  CCTK_REAL bh_betax[NBH];
  CCTK_REAL bh_betay[NBH];
  CCTK_REAL bh_betaz[NBH];

  assert (num_vars == 3);
  CCTK_POINTER output_arrays[3];
  output_arrays[0] = bh_betax;
  output_arrays[1] = bh_betay;
  output_arrays[2] = bh_betaz;

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
    for(int i = 0; i < NBH; ++i) {
      if(track[i]) {
        // Time evolution
        // First order time integrator
        // Michael Koppitz says this works...
        // if it doesn't, we can make it second order accurate
        CCTK_REAL const dt = cctk_time - bh_time[i];
        bh_pos_x[i] = bh_pos_x_old[i] + dt * (- bh_betax[i]);
        bh_pos_y[i] = bh_pos_y_old[i] + dt * (- bh_betay[i]);
        bh_pos_z[i] = bh_pos_z_old[i] + dt * (- bh_betaz[i]);
        bh_vel_x[i] = - bh_betax[i];
        bh_vel_y[i] = - bh_betay[i];
        bh_vel_z[i] = - bh_betaz[i];
      }
    }
  }

  // Broadcast result

  CCTK_REAL loc_local[6*NBH]; /* 3 components for location, 3 components for velocity */
  if (CCTK_MyProc(cctkGH) == 0) {
    for(int i = 0; i < NBH; ++i) {
      loc_local[0*NBH + i] = bh_pos_x[i];
      loc_local[1*NBH + i] = bh_pos_y[i];
      loc_local[2*NBH + i] = bh_pos_z[i];
      loc_local[3*NBH + i] = bh_vel_x[i];
      loc_local[4*NBH + i] = bh_vel_y[i];
      loc_local[5*NBH + i] = bh_vel_z[i];
    }
  } else {
    for(int i = 0; i < NBH; ++i) {
      loc_local[0*NBH + i] = 0.0;
      loc_local[1*NBH + i] = 0.0;
      loc_local[2*NBH + i] = 0.0;
      loc_local[3*NBH + i] = 0.0;
      loc_local[4*NBH + i] = 0.0;
      loc_local[5*NBH + i] = 0.0;
    }
  }

  CCTK_REAL loc_global[6*NBH]; /* 3 components for location, 3 components for velocity */

  int const handle_sum = CCTK_ReductionHandle ("sum");
  if (handle_sum < 0) {
    CCTK_ERROR("Can't get redunction handle");
  }

  int const ierr2 = CCTK_ReduceLocArrayToArray1D
    (cctkGH, -1, handle_sum,
     loc_local, loc_global, 6*NBH, CCTK_VARIABLE_REAL);
  if (ierr2 < 0) {
    CCTK_ERROR("Reduction error");
  }

  for(int i = 0; i < NBH; ++i) {
    bh_time[i]  = cctk_time;
    bh_pos_x[i] = loc_global[0*NBH + i];
    bh_pos_y[i] = loc_global[1*NBH + i];
    bh_pos_z[i] = loc_global[2*NBH + i];
    bh_vel_x[i] = loc_global[3*NBH + i];
    bh_vel_y[i] = loc_global[4*NBH + i];
    bh_vel_z[i] = loc_global[5*NBH + i];
  }

  if(verbose) {
    char msg[BUFSIZ];
    for(int i = 0; i < NBH; ++i) {
      if(track[i]) {
        snprintf(msg, BUFSIZ, "BH[%d] position at t = %g: (%g, %g, %g)",
            i, cctk_time, bh_pos_x[i], bh_pos_y[i], bh_pos_z[i]);
        CCTK_INFO(msg);
      }
    }
  }

  for(int i = 0; i < NBH; ++i) {
    if(update_surface[i]) {
      sf_origin_x[index_surface[i]] = bh_pos_x[i];
      sf_origin_y[index_surface[i]] = bh_pos_y[i];
      sf_origin_z[index_surface[i]] = bh_pos_z[i];
    }
  }
}

void BHTrackerMany_UpdateAMR(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BHTrackerMany_UpdateAMR");
  }

  for(int i = 0; i < NBH; ++i) {
    if(index_region[i] >= 0) {
      position_x[index_region[i]] = bh_pos_x[i];
      position_y[index_region[i]] = bh_pos_y[i];
      position_z[index_region[i]] = bh_pos_z[i];
    }
  }

  for(int i = 0; i < NBH; ++i) {
    if(adapt_resolution[i]) {
      if(!sf_active[index_surface[i]]) {
        continue;
      }

      CCTK_REAL grid_size[30];
      int lsh[2];
      getvectorindex2(cctkGH, "CarpetRegrid2::radii", lsh);
      for(int rl = 0; rl < 30; ++rl) {
        int const ind = index2(lsh, rl, index_region[i]);
        grid_size[rl] = radius[ind];
      }

      CCTK_REAL const bh_radius = sf_min_radius[index_surface[i]];
      if(verbose) {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "BH min radius t = %g: %g",
            cctk_time, bh_radius);
        CCTK_INFO(msg);
      }

      int const nl = num_levels[index_region[i]];
      if(nl < 1 || nl > 30) {
        return;
      }

      if(bh_radius < adapt_factor_refine[i]*grid_size[nl-1] &&
          nl < adapt_nlevels_max[i] - 1) {
        num_levels[index_region[i]]++;
        if(verbose) {
          char msg[BUFSIZ];
          snprintf(msg, BUFSIZ, "Refine grid: num_levels[%d]: %d -> %d",
              index_region[i], num_levels[index_region[i]]-1,
              num_levels[index_region[i]]);
          CCTK_INFO(msg);
        }
      }
      else if(bh_radius > adapt_factor_coarse[i]*grid_size[nl-1]
          && nl > adapt_nlevels_min[i] + 1) {
        num_levels[index_region[i]]--;
        if(verbose) {
          char msg[BUFSIZ];
          snprintf(msg, BUFSIZ, "Coarsen grid: num_levels[%d]: %d -> %d",
              index_region[i], num_levels[index_region[i]]+1,
              num_levels[index_region[i]]);
          CCTK_INFO(msg);
        }
      }

      bh_num_levels[i] = num_levels[index_region[i]];
    }
  }
}
