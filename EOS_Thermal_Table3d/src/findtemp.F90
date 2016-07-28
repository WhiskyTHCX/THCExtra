MODULE findtemp
#include <tableindex.h>

CONTAINS

integer function bisection(lr,lt,ye,leps)
  use eos3dmod
  use interpolations

  implicit none

  real*8, intent(in)  :: lr,ye,leps
  real*8, intent(out) :: lt

  integer :: ilr, iye, lotemp, midtemp, hitemp, maxstep, ind
  real*8  :: wlr, wye, leps_hi, leps_mid, leps_lo

  bisection = 0
  maxstep   = 15

  call interp_weights(lr, eos_lrhomin, dlrho, nrho, ilr, wlr)
  call interp_weights(ye, eos_yemin, dye, nye, iye, wye)

  hitemp  = ntemp
  lotemp  = 1

  leps_hi = lookup_bilin(ilr, hitemp, iye, wlr, wye, LENERGY_IDX)

  if (leps .gt. leps_hi) then
     lt        = eos_ltempmax
     bisection = -1
     return
  end if


  leps_lo = lookup_bilin(ilr, lotemp, iye, wlr, wye, LENERGY_IDX)

  if (leps .lt. leps_lo) then
     lt        = eos_ltempmin
     bisection = -1
     return
  end if

  ind = 0
  do while ((hitemp-lotemp).ge.2)
    ind = ind + 1
    if (ind.ge.maxstep) then
      write(*,*) 'Severe problem: bisection error, too many iterations'
    end if

    midtemp  = (hitemp + lotemp) / 2
    leps_mid = lookup_bilin(ilr, midtemp, iye, wlr, wye, LENERGY_IDX)

    if (leps_mid.le.leps) then
      lotemp  = midtemp
      leps_lo = leps_mid
    else
      hitemp  = midtemp
      leps_hi = leps_mid
    end if
  end do

  lt = logtemp(lotemp) + dltemp * (leps - leps_lo) / (leps_hi - leps_lo)

end function bisection

integer Function GetTemperature(lrho, ye, leps, ltemp)
  use eos3dmod
  implicit none
  real*8, intent(in) ::  lrho, ye, leps
  real*8, intent(out)::  ltemp

  GetTemperature = bisection(lrho, ltemp, ye, leps)
end Function GetTemperature


END MODULE findtemp

