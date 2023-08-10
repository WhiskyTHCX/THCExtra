#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! CGS: cm, g, s, MeV
INTEGER FUNCTION NeutrinoDens_cgs(rho, temp, ye,&
                                   n_nue, n_nua, n_nux, en_nue, en_nua, en_nux)

#ifndef FORTRAN_DISABLE_IEEE_ARITHMETIC
    use ieee_arithmetic
#endif
    use table3d_mod
    use lk_interpolations
    use units

    implicit none

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: n_nue, n_nua, n_nux, en_nue, en_nua, en_nux

    REAL*8 :: mass_fact_cgs

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !All the emission rates are expressed per unit of volume

    REAL*8 :: fermi2, fermi3

    REAL*8 :: nb, lrho,ltemp, &
              mu_n, mu_p, mu_e, &
              abar, zbar, xp, xn, xh

    REAL*8 :: eta_nue, eta_nua, eta_nux, &
              eta_n, eta_p, eta_hat, eta_e, &
              eta_np, eta_pn

    NeutrinoDens_cgs = 0

#define WEAK_RATES_ITS_ME
#include "inc/weak_rates_guts.inc"
#undef WEAK_RATES_ITS_ME

    n_nue = 4.0d0 * pi / hc_mevcm**3 * temp**3 * FERMI2(eta_nue)
    n_nua = 4.0d0 * pi / hc_mevcm**3 * temp**3 * FERMI2(eta_nua)
    n_nux = 16.0d0 * pi / hc_mevcm**3 * temp**3 * FERMI2(eta_nux)

    en_nue = 4.0d0 * pi / hc_mevcm**3 * temp**4 * FERMI3(eta_nue)
    en_nua = 4.0d0 * pi / hc_mevcm**3 * temp**4 * FERMI3(eta_nua)
    en_nux = 16.0d0 * pi / hc_mevcm**3 * temp**4 * FERMI3(eta_nux)

#ifndef FORTRAN_DISABLE_IEEE_ARITHMETIC
    if (.not.ieee_is_finite(n_nue)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in n_nue", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif

    if (.not.ieee_is_finite(n_nua)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in n_nua", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif

    if (.not.ieee_is_finite(n_nux)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in n_nux", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif

    if (.not.ieee_is_finite(en_nue)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in en_nue", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif

    if (.not.ieee_is_finite(en_nua)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in en_nua", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif

    if (.not.ieee_is_finite(en_nux)) then
       write(*,*) "NeutrinoDens_cgs: NaN/Inf in en_nux", rho, temp, ye
       NeutrinoDens_cgs= -1
    endif
#endif

    return
END FUNCTION NeutrinoDens_cgs
