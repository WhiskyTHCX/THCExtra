//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"
#include "finite_difference.h"

#define SQ(X) ((X)*(X))

void BNSAnalysis_MachNumber(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BNSAnalysis_MachNumber");
  }

  CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

  /*
    Compute the Mach number.
    M := Wv/(W_s * c_s)
    See eq 4.218 of Rezzolla
  */
  for(int k = 0; k < cctk_lsh[2]; ++k)
  for(int j = 0; j < cctk_lsh[1]; ++j)
  for(int i = 0; i < cctk_lsh[0]; ++i) {
    int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    CCTK_REAL const sound_lorentz = 1.0/sqrt(1.0 - SQ(csound[ijk]));

    CCTK_REAL v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] +
      gxz[ijk]*velz[ijk];
    CCTK_REAL v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] +
      gyz[ijk]*velz[ijk];
    CCTK_REAL v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] +
      gzz[ijk]*velz[ijk];
    CCTK_REAL const velnorm = sqrt(v_x*velx[ijk] + v_y*vely[ijk] +
        v_z*velz[ijk]);
    mach_number[ijk] = velnorm/csound[ijk];

    mach_number[ijk] = mach_number[ijk]/sound_lorentz*w_lorentz[ijk];
  }
}

extern "C" void BNSAnalysis_Vorticity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
    CCTK_INFO("BNSAnalysis_Vorticity");
  }

  CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

  CCTK_REAL * specific_momentum_x = &specific_momentum[
    0*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * specific_momentum_y = &specific_momentum[
    1*UTILS_GFSIZE(cctkGH)];
  CCTK_REAL * specific_momentum_z = &specific_momentum[
    2*UTILS_GFSIZE(cctkGH)];

  CCTK_REAL const idx = 1.0 / CCTK_DELTA_SPACE(0);
  CCTK_REAL const idy = 1.0 / CCTK_DELTA_SPACE(1);
  CCTK_REAL const idz = 1.0 / CCTK_DELTA_SPACE(2);

  for(int ijk = 0; ijk < UTILS_GFSIZE(cctkGH); ++ijk) {
    CCTK_REAL const enthalpy = 1.0 + eps[ijk] + press[ijk] / rho[ijk];
    CCTK_REAL const v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] +
      gxz[ijk]*velz[ijk];
    CCTK_REAL const v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] +
      gyz[ijk]*velz[ijk];
    CCTK_REAL const v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] +
      gzz[ijk]*velz[ijk];
    specific_momentum_x[ijk] = enthalpy * w_lorentz[ijk] * v_x;
    specific_momentum_y[ijk] = enthalpy * w_lorentz[ijk] * v_y;
    specific_momentum_z[ijk] = enthalpy * w_lorentz[ijk] * v_z;
  }

  /*
    Compute the vorticity tensor.
    \Omega_{uv} = \partial_{v}(hu_{u}) - \partial_{u}(hu_{v})
    See. 3.133 of Rezzolla and Zanotti
  */
  for(int k = 0; k < cctk_lsh[2]; ++k)
  for(int j = 0; j < cctk_lsh[1]; ++j)
  for(int i = 0; i < cctk_lsh[0]; ++i) {
    int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

    vorticity_xy[ijk] =
      idy*adiff_y(cctkGH, specific_momentum_x, i,j,k, fd_order) -
      idx*adiff_x(cctkGH, specific_momentum_y, i,j,k, fd_order);

    vorticity_xz[ijk] =
      idz*adiff_z(cctkGH, specific_momentum_x, i,j,k, fd_order) -
      idx*adiff_x(cctkGH, specific_momentum_z, i,j,k, fd_order);

    vorticity_yz[ijk] =
      idz*adiff_z(cctkGH, specific_momentum_y, i,j,k, fd_order) -
      idy*adiff_y(cctkGH, specific_momentum_z, i,j,k, fd_order);
  }
}
