#include <tableindex.h>

REAL (C_DOUBLE) FUNCTION tab3d_press(rho, eps, ye, ierr) BIND(C, NAME="tab3d_press")
    USE ISO_C_BINDING
    use eos3Dmod
    use findtemp
    use interpolations

    implicit none
    REAL (C_DOUBLE), INTENT(IN), VALUE :: rho, eps, ye
    INTEGER (C_INT), INTENT(OUT) :: ierr
    real*8 :: lrho, leps, ltemp, pressure

    lrho = log10(rho)
    leps = log10(eps)
    ierr = GetTemperature(lrho, ye, leps, ltemp)
    if (ierr.eq.0) then
      pressure = linearInterpolation3d(lrho, ltemp, ye, LPRESS_IDX)
      pressure = (10.**(pressure))
    else
      pressure = 0
    end if
    tab3d_press = pressure
    return
end FUNCTION tab3d_press


REAL (C_DOUBLE) FUNCTION tab3d_csnd2(rho, eps, ye, ierr) BIND(C, NAME="tab3d_csnd2")
    USE ISO_C_BINDING
    use eos3Dmod
    use findtemp
    use interpolations

    implicit none
    REAL (C_DOUBLE), INTENT(IN), VALUE :: rho, eps, ye
    INTEGER (C_INT), INTENT(OUT) :: ierr
    real*8 :: lrho, leps, ltemp, csnd2

    lrho = log10(rho)
    leps = log10(eps)
    ierr = GetTemperature(lrho, ye, leps, ltemp)
    if (ierr.eq.0) then
      csnd2 = linearInterpolation3d(lrho, ltemp, ye, CS2_IDX)
    else
      csnd2 = 0
    end if
    tab3d_csnd2 = csnd2
    return
end FUNCTION tab3d_csnd2

REAL (C_DOUBLE) FUNCTION tab3d_temp(rho, eps, ye, ierr) BIND(C, NAME="tab3d_temp")
    USE ISO_C_BINDING
    use eos3Dmod
    use findtemp

    implicit none
    REAL (C_DOUBLE), INTENT(IN), VALUE :: rho, eps, ye
    INTEGER (C_INT), INTENT(OUT) :: ierr
    real*8 :: lrho, leps, ltemp, temp

    lrho = log10(rho)
    leps = log10(eps)
    ierr = GetTemperature(lrho, ye, leps, ltemp)
    if (ierr.eq.0) then
      temp  = 1.0d1**ltemp
    else
      temp = 0
    end if
    tab3d_temp = temp
    return
end FUNCTION tab3d_temp


REAL (C_DOUBLE) FUNCTION tab3d_eps(rho, temp, ye) BIND(c, NAME="tab3d_eps")
    USE ISO_C_BINDING
    USE eos3dmod
    USE interpolations
    implicit none
    REAL (C_DOUBLE), INTENT(IN), VALUE :: rho, temp, ye
    real*8 :: lrho, leps, ltemp

    lrho      = log10(rho)
    ltemp     = log10(temp)
    leps      = linearInterpolation3d(lrho, ltemp, ye, LENERGY_IDX)
    tab3d_eps = (1.0d1**leps)
    return
end FUNCTION tab3d_eps

REAL (C_DOUBLE) FUNCTION tab3d_press_from_temp(rho, temp, ye) BIND(c, NAME="tab3d_press_from_temp")
    USE ISO_C_BINDING
    USE eos3dmod
    USE interpolations
    implicit none
    real (C_DOUBLE), INTENT(IN), value :: rho, temp, ye
    real*8 :: lrho, ltemp, lpressure

    lrho      = log10(rho)
    ltemp     = log10(temp)
    lpressure = linearInterpolation3d(lrho, ltemp, ye, LPRESS_IDX)

    tab3d_press_from_temp = (1.0d1**(lpressure))
    return
end FUNCTION tab3d_press_from_temp

REAL (C_DOUBLE) FUNCTION tab3d_csnd2_from_temp(rho, temp, ye) BIND(c, NAME="tab3d_csnd2_from_temp")
    USE ISO_C_BINDING
    USE eos3dmod
    USE interpolations
    implicit none
    real (C_DOUBLE), INTENT(IN), value :: rho, temp, ye
    real*8 :: lrho, ltemp

    lrho      = log10(rho)
    ltemp     = log10(temp)

    tab3d_csnd2_from_temp = linearInterpolation3d(lrho, ltemp, ye, CS2_IDX)
    return
end FUNCTION tab3d_csnd2_from_temp


REAL (C_DOUBLE) FUNCTION tab3d_entropy_from_temp(rho, temp, ye) BIND(c, NAME="tab3d_entropy_from_temp")
    USE ISO_C_BINDING
    USE eos3dmod
    USE interpolations
    implicit none
    real (C_DOUBLE), INTENT(IN), value :: rho, temp, ye
    real*8 :: lrho, ltemp, lentr

    lrho      = log10(rho)
    ltemp     = log10(temp)
    lentr     = linearInterpolation3d(lrho, ltemp, ye, ENTROPY_IDX)
    tab3d_entropy_from_temp = (1.0d1**lentr)
    return
end FUNCTION tab3d_entropy_from_temp

REAL (C_DOUBLE) FUNCTION tab3d_entropy_from_eps(rho, eps, ye, ierr) BIND(C, NAME="tab3d_entropy_from_eps")
    USE ISO_C_BINDING
    use eos3Dmod
    use findtemp
    use interpolations

    implicit none
    REAL (C_DOUBLE), INTENT(IN), VALUE :: rho, eps, ye
    INTEGER (C_INT), INTENT(OUT) :: ierr
    real*8 :: lentr, entropy, lrho, leps, ltemp

    lrho = log10(rho)
    leps = log10(eps)
    ierr = GetTemperature(lrho, ye, leps, ltemp)
    if (ierr.eq.0) then
      lentr    = linearInterpolation3d(lrho, ltemp, ye, ENTROPY_IDX)
      entropy  = (1.0d1**lentr)
    else
      entropy  = 0
    end if
    tab3d_entropy_from_eps = entropy
    return
end FUNCTION tab3d_entropy_from_eps


