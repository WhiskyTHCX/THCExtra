MODULE lk_interpolations

contains

real*8 function lk_interp_lin(f0, f1, w)
  implicit none
  real*8, intent(in) :: f0,f1,w
  lk_interp_lin = f0*(1-w) + f1*w
  return
end function lk_interp_lin

real*8 function lk_interp_bilin(fx0y0, fx1y0, fx0y1, fx1y1, wx, wy)
  implicit none
  real*8, intent(in) :: fx0y0, fx1y0, fx0y1, fx1y1, wx, wy
  real*8 :: fy0, fy1
  fy0           = lk_interp_lin(fx0y0, fx1y0, wx)
  fy1           = lk_interp_lin(fx0y1, fx1y1, wx)
  lk_interp_bilin  = lk_interp_lin(fy0, fy1, wy)
  return
end function lk_interp_bilin

real*8 function lk_lookup_bilin(ilrho, iltemp, iye, wlrho, wye, ivar)
  use table3d_mod
  implicit none
  integer, intent(in) :: ilrho, iltemp, iye, ivar
  real*8, intent(in) :: wlrho, wye
  real*8 :: fx0y0, fx1y0, fx0y1, fx1y1

  fx0y0 = allvariables(ilrho,   iltemp, iye,   ivar)
  fx1y0 = allvariables(ilrho+1, iltemp, iye,   ivar)
  fx0y1 = allvariables(ilrho,   iltemp, iye+1, ivar)
  fx1y1 = allvariables(ilrho+1, iltemp, iye+1, ivar)

  lk_lookup_bilin = lk_interp_bilin(fx0y0, fx1y0, fx0y1, fx1y1, wlrho, wye)
  return
end function lk_lookup_bilin

real*8 function lookup_trilin(ilrho, iltemp, iye, wlrho, wltemp, wye, ivar)
  implicit none
  integer, intent(in) :: ilrho, iltemp, iye, ivar
  real*8, intent(in) :: wlrho, wltemp, wye
  real*8 :: ltemp0, ltemp1

  ltemp0 = lk_lookup_bilin(ilrho, iltemp,   iye, wlrho, wye, ivar)
  ltemp1 = lk_lookup_bilin(ilrho, iltemp+1, iye, wlrho, wye, ivar)

  lookup_trilin = lk_interp_lin(ltemp0, ltemp1, wltemp)
  return
end function lookup_trilin

subroutine lk_interp_weights(x, xmin, dx, max_ix, ix, wx)
  implicit none
  real*8, intent(in) :: x, xmin, dx
  integer, intent(in) :: max_ix
  integer, intent(out) :: ix
  real*8, intent(out) :: wx
  real*8 :: s

  s    = (x - xmin) / dx
  ix   = max(1, min(max_ix-1 , 1+floor(s)))
  wx   = s  - (ix-1)
end subroutine lk_interp_weights

real*8 function lkLinearInterpolation3d(lrho, ltemp, ye, var_index)
  use table3d_mod
  implicit none
  real*8, intent(in):: lrho, ltemp, ye
  integer, intent(in) :: var_index
  integer ilrho, iltemp, iye
  real*8 wlrho, wltemp, wye

  call lk_interp_weights(lrho, eos_lrhomin, dlrho, nrho, ilrho, wlrho)
  call lk_interp_weights(ltemp, eos_ltempmin, dltemp, ntemp, iltemp, wltemp)
  call lk_interp_weights(ye, eos_yemin, dye, nye, iye, wye)

  lkLinearInterpolation3d = lookup_trilin(ilrho, iltemp, iye, wlrho, wltemp, wye, var_index)
  return
end function lkLinearInterpolation3d

END MODULE lk_interpolations
