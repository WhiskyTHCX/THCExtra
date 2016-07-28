#include <tableindex.h>

INTEGER (C_INT) function eps_range(myRho, myYe, epsmin, epsmax) BIND(C, NAME="eps_range")
  USE ISO_C_BINDING
  use eos3Dmod
  use interpolations
  implicit none

  REAL (C_DOUBLE), intent(in), VALUE :: myRho, myYe
  REAL (C_DOUBLE), intent(out)       :: epsmin, epsmax

  real*8  :: rho0, lrho0, ye0, wlrho, wye
  integer :: ilrho, iye
  real*8  :: lepsmin, lepsmax, epsfsafety

  ! Check if rho is in allowed range
  rho0 = myRho
  if(rho0.lt.eos_rhomin.or.rho0.gt.eos_rhomax) then
     eps_range = -1
     return
  endif

  ! Check if ye is in allowed range
  ye0  = myYe
  if(ye0.lt.eos_yemin.or.ye0.gt.eos_yemax) then
     eps_range = -1
     return
  endif

  lrho0 = log10(rho0)

  call interp_weights(lrho0, eos_lrhomin, dlrho, nrho, ilrho, wlrho)
  call interp_weights(ye0, eos_yemin, dye, nye, iye, wye)

  lepsmin = lookup_bilin(ilrho, 1, iye, wlrho, wye, LENERGY_IDX)
  lepsmax = lookup_bilin(ilrho, ntemp, iye, wlrho, wye, LENERGY_IDX)

  epsfsafety = 1.00000000001d0
  epsmin = (10.0d0**lepsmin) * epsfsafety
  epsmax = (10.0d0**lepsmax) / epsfsafety

  eps_range = 0
  return
end function eps_range

subroutine rho_range(rho_min, rho_max) BIND(C, NAME="rho_range")
  USE ISO_C_BINDING
  use eos3Dmod
  implicit none
  REAL (C_DOUBLE), intent(out) :: rho_min, rho_max
  REAL*8 :: rho_fsafety
  rho_fsafety = 10.0d0**(dlrho*1.0d-2)
  rho_max = (eos_rhomax ) / rho_fsafety
  rho_min = (eos_rhomin ) * rho_fsafety
end subroutine rho_range

subroutine ye_range(ye_min, ye_max) BIND(C, NAME="ye_range")
  USE ISO_C_BINDING
  use eos3Dmod
  implicit none
  REAL (C_DOUBLE), intent(out) :: ye_min, ye_max
  REAL*8 :: ye_safety
  ye_safety = dye * 1.0d-5
  ye_min = eos_yemin + ye_safety
  ye_max = eos_yemax - ye_safety
end subroutine ye_range

subroutine temp_range(temp_min, temp_max) BIND(C, NAME="temp_range")
  USE ISO_C_BINDING
  use eos3Dmod
  implicit none
  REAL (C_DOUBLE), intent(out) :: temp_min, temp_max
  REAL*8 :: temp_fsafety
  temp_fsafety = 10.0d0**(dltemp*1.0d-2)

  temp_min = eos_tempmin * temp_fsafety
  temp_max = eos_tempmax / temp_fsafety
end subroutine temp_range

