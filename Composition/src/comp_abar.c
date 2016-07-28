#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Comp_Abar(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

#pragma omp parallel
    {
        for(int k = 0; k < cctk_lsh[2]; ++k)
        for(int j = 0; j < cctk_lsh[1]; ++j)
        for(int i = 0; i < cctk_lsh[0]; ++i) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ierr = NucleiAbar(rho[ijk], temperature[ijk], Y_e[ijk],
                    &Abar[ijk]); assert(!ierr);
        }
    }
}
