#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "pizza_central.h"

#include "mpi.h"

#include <bin_ns.h>
#include <algorithm>
#include <stdexcept>

namespace LoreneID {

/// Copy data array and make sure pointers are nonzero
/** dst Destination
    src Source
    s   Size of array
*/
void check_and_copy(CCTK_REAL* dst, const CCTK_REAL* src, size_t s)
{
  if ((0==src) || (0==dst))
    CCTK_ERROR("LoreneID: Cactus passed NULL pointer to active timelevels.");
  copy(src, src+s, dst);
}

/// Makes a copy of a file
void copy_file(char const * src, char const * dst)
{
  // Source and destination are the same... nothing to do
  if (0 == strcmp(src, dst)) {
    return;
  }

  char msg[BUFSIZ];
  FILE * fsrc = fopen(src, "r");
  if (!fsrc) {
    snprintf(msg, BUFSIZ, "File not found: \"%s\"", src);
    CCTK_ERROR(msg);
  }
  FILE * fdst = fopen(dst, "w");
  if (!fdst) {
    fclose(fsrc);
    snprintf(msg, BUFSIZ, "Could not open file: \"%s\"", src);
    CCTK_ERROR(msg);
  }

  int ch;
  while ((ch = getc(fsrc)) != EOF) {
    putc(ch, fdst);
  }

  fclose(fsrc);
  fclose(fdst);
  sync();
}

} // namespace

using namespace std;
using namespace Pizza;
using namespace LoreneID;

/// Load binary neutron star initial data in LORENE format
/** We just copy data provided by LORENE objects to the CACTUS grid,
    conversions of units and definitions are done in LoreneID_Add_Missing
**/
extern "C" void LoreneID_LoadMeudon(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int local_grid_size = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
  CCTK_REAL* velocity[] = {vel, vel+local_grid_size, vel+2*local_grid_size};

  const bool set_lapse = (string(initial_lapse)=="LoreneBNS");
  const bool set_shift = (string(initial_shift)=="LoreneBNS");

  // Set up an array of coordinates in LORENE units

  const units u       = Base::pizza_base_central::get().internal_units;
  // CACTUS length unit in LORENE length units (km)
  const CCTK_REAL length_ca_in_lo = u.length() / 1e3;
  //was 1.47664;
  for(int i=0; i<cctk_lsh[0]; i++) {
    for(int j=0; j<cctk_lsh[1]; j++) {
      for(int k=0; k<cctk_lsh[2]; k++) {
        int i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
        x_ulorene[i3D] = x[i3D] * length_ca_in_lo;
        y_ulorene[i3D] = y[i3D] * length_ca_in_lo;
        z_ulorene[i3D] = z[i3D] * length_ca_in_lo;
      }
    }
  }

  CCTK_INFO("Reading Meudon BNS initial data.");
  try {
    if (0 != access(lorene_bns_file, R_OK)) {
      char msg[BUFSIZ];
      snprintf(msg, BUFSIZ, "File not found: \"%s\"", lorene_bns_file);
      CCTK_ERROR(msg);
    }
    if (copy_lorene_eos) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (0 == rank) {
        copy_file(lorene_tab_file, lorene_eos_name);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // call the Lorene routine to transform spectral data to cartesian
    Bin_NS binary_system(local_grid_size, x_ulorene, y_ulorene, z_ulorene, lorene_bns_file);

    // save some properties for later use.
    *lorene_bns_omega      = binary_system.omega / u.freq();
    *lorene_bns_separation = binary_system.dist * 1e3 / u.length();

    // copy to CACTUS grid
    check_and_copy(gxx, binary_system.g_xx, local_grid_size);
    check_and_copy(gxy, binary_system.g_xy, local_grid_size);
    check_and_copy(gxz, binary_system.g_xz, local_grid_size);
    check_and_copy(gyy, binary_system.g_yy, local_grid_size);
    check_and_copy(gyz, binary_system.g_yz, local_grid_size);
    check_and_copy(gzz, binary_system.g_zz, local_grid_size);

    if (set_lapse) {
      check_and_copy(alp, binary_system.nnn, local_grid_size);
    }
    if (set_shift) {
      check_and_copy(betax, binary_system.beta_x, local_grid_size);
      check_and_copy(betay, binary_system.beta_y, local_grid_size);
      check_and_copy(betaz, binary_system.beta_z, local_grid_size);
    }

    check_and_copy(kxx, binary_system.k_xx, local_grid_size);
    check_and_copy(kxy, binary_system.k_xy, local_grid_size);
    check_and_copy(kxz, binary_system.k_xz, local_grid_size);
    check_and_copy(kyy, binary_system.k_yy, local_grid_size);
    check_and_copy(kyz, binary_system.k_yz, local_grid_size);
    check_and_copy(kzz, binary_system.k_zz, local_grid_size);

    check_and_copy(rho, binary_system.nbar, local_grid_size);

    check_and_copy(velocity[0], binary_system.u_euler_x, local_grid_size);
    check_and_copy(velocity[1], binary_system.u_euler_y, local_grid_size);
    check_and_copy(velocity[2], binary_system.u_euler_z, local_grid_size);
  }
  catch (exception &e) {
    CCTK_ERROR(e.what());
  }
  catch (...) {
    CCTK_ERROR("Unknown exception in LoreneID_LoadMeudon");
  }

}


