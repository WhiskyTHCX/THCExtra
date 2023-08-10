#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

CCTK_INT FUNCTION NeutrinoEmissionImpl(rho, temp, ye, &
                                   emissionRatesRloc_nue, emissionRatesRloc_nua, &
                                   emissionRatesRloc_nux, emissionRatesQloc_nue, &
                                   emissionRatesQloc_nua,emissionRatesQloc_nux)
    use table3d_mod
    use units

    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye !density, temperature, Ye
    CCTK_REAL, INTENT(OUT) :: emissionRatesRloc_nue,emissionRatesRloc_nua, &
                            emissionRatesRloc_nux
    CCTK_REAL, INTENT(OUT) :: emissionRatesQloc_nue,emissionRatesQloc_nua, &
                            emissionRatesQloc_nux

    CCTK_REAL :: R_nue_cgs, R_nua_cgs, R_nux_cgs,&
                 Q_nue_cgs, Q_nua_cgs, Q_nux_cgs

    CCTK_REAL :: rho_cgs, temp0, ye0
    INTEGER   :: enforceTableBounds
    INTEGER   :: Emissions_cgs
    INTEGER   :: err, boundsErr
    CCTK_REAL :: r_cgs2cactus, q_cgs2cactus

    NeutrinoEmissionImpl = 0

    !Convert the density in cgs
    rho_cgs  = rho * cactus2cgsRho
    temp0 = temp
    ye0   = ye

    if ((rho_cgs.lt.rho_min_cgs).or.(temp0.lt.temp_min_mev)) then
      emissionRatesRloc_nue = 0.
      emissionRatesRloc_nua = 0.
      emissionRatesRloc_nux = 0.
      emissionRatesQloc_nue = 0.
      emissionRatesQloc_nua = 0.
      emissionRatesQloc_nux = 0.
      return
    end if

    boundsErr = enforceTableBounds(rho_cgs,temp0,ye0)

    if (boundsErr.eq.-1) then
      NeutrinoEmissionImpl = -1
      return
    end if

    !the temperature is in MeV
    err = Emissions_cgs(rho_cgs, temp0, ye0,&
                       R_nue_cgs, R_nua_cgs, R_nux_cgs,&
                       Q_nue_cgs, Q_nua_cgs, Q_nux_cgs)

    if (err.eq.-1) then
      NeutrinoEmissionImpl = -1
    end if


    !Unit conversion for the number rates:
    !number of neutrinos * cm^-3 * s^-1 to cactus units
    !1.586234651026439e+10
    r_cgs2cactus = 1. / (cgs2cactusTime*cgs2cactusLength**3)

    emissionRatesRloc_nue = R_nue_cgs * (r_cgs2cactus/normfact)
    emissionRatesRloc_nua = R_nua_cgs * (r_cgs2cactus/normfact)
    emissionRatesRloc_nux = R_nux_cgs * (r_cgs2cactus/normfact)

    !Unit conversion for the energy rates:
    !From MeV * cm^-3 * s^-1 to erg * cm^-3 * s^-1
    !and finally to cactus units
    !1.421737093266046e-50
    q_cgs2cactus = mev_to_erg*cgs2cactusenergy/ &
                   (cgs2cactusTime * cgs2cactusLength**3)

    emissionRatesQloc_nue = Q_nue_cgs * q_cgs2cactus
    emissionRatesQloc_nua = Q_nua_cgs * q_cgs2cactus
    emissionRatesQloc_nux = Q_nux_cgs * q_cgs2cactus

    return
END FUNCTION NeutrinoEmissionImpl

CCTK_INT FUNCTION NeutrinoOpacityImpl(rho, temp, ye, &
                                      kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                                      kappa_1_nue, kappa_1_nua, kappa_1_nux)
    use table3d_mod
    use units
    implicit none
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                              kappa_1_nue, kappa_1_nua, kappa_1_nux

    CCTK_REAL :: kappa_0_nue_cgs, kappa_0_nua_cgs,&
                 kappa_0_nux_cgs,&
                 kappa_1_nue_cgs, kappa_1_nua_cgs,&
                 kappa_1_nux_cgs

    CCTK_REAL :: rho_cgs, temp0, ye0
    INTEGER   :: enforceTableBounds
    INTEGER   :: Opacities_cgs
    INTEGER   :: err, boundsErr
    CCTK_REAL :: kappa_cgs2cactus

    NeutrinoOpacityImpl = 0

    !Convert the density in cgs
    rho_cgs  = rho * cactus2cgsRho
    temp0 = temp
    ye0   = ye

    if ((rho_cgs.lt.rho_min_cgs).or.(temp0.lt.temp_min_mev)) then
      kappa_0_nue = 0.
      kappa_0_nua = 0.
      kappa_0_nux = 0.
      kappa_1_nue = 0.
      kappa_1_nua = 0.
      kappa_1_nux = 0.
      return
    end if

    boundsErr = enforceTableBounds(rho_cgs,temp0,ye0)

    if (boundsErr.eq.-1) then
      NeutrinoOpacityImpl = -1
      return
    end if

    err = Opacities_cgs(rho_cgs, temp0, ye0,&
                        kappa_0_nue_cgs, kappa_0_nua_cgs,&
                        kappa_0_nux_cgs,&
                        kappa_1_nue_cgs, kappa_1_nua_cgs,&
                        kappa_1_nux_cgs)
    if (err.eq.-1) then
      write(*,*) rho_cgs, temp0, ye0
      NeutrinoOpacityImpl = -1
    end if

    !Unit conversion
    kappa_cgs2cactus = 1. / (cgs2cactusLength)
    kappa_0_nue = kappa_cgs2cactus * kappa_0_nue_cgs
    kappa_0_nua = kappa_cgs2cactus * kappa_0_nua_cgs
    kappa_0_nux = kappa_cgs2cactus * kappa_0_nux_cgs

    kappa_1_nue = kappa_cgs2cactus * kappa_1_nue_cgs
    kappa_1_nua = kappa_cgs2cactus * kappa_1_nua_cgs
    kappa_1_nux = kappa_cgs2cactus * kappa_1_nux_cgs

    return
END FUNCTION NeutrinoOpacityImpl


CCTK_INT FUNCTION NeutrinoAbsorptionRateImpl(rho, temp, ye,&
                              kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                              kappa_1_nue, kappa_1_nua, kappa_1_nux)
    use table3d_mod
    use units
    IMPLICIT NONE
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, INTENT(IN) :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: kappa_0_nue, kappa_0_nua, kappa_0_nux, &
                              kappa_1_nue, kappa_1_nua, kappa_1_nux

    CCTK_REAL :: kappa_0_nue_cgs, kappa_0_nua_cgs,&
                 kappa_0_nux_cgs,&
                 kappa_1_nue_cgs, kappa_1_nua_cgs,&
                 kappa_1_nux_cgs

    CCTK_REAL :: rho_cgs, temp0, ye0
    INTEGER   :: enforceTableBounds
    INTEGER   :: err, boundsErr
    INTEGER   :: Absorption_cgs
    CCTK_REAL :: kappa_cgs2cactus

    NeutrinoAbsorptionRateImpl = 0

    !Convert the density in cgs
    rho_cgs  = rho * cactus2cgsRho
    temp0 = temp
    ye0   = ye

    if ((rho_cgs.lt.rho_min_cgs).or.(temp0.lt.temp_min_mev)) then
      kappa_0_nue = 0.
      kappa_0_nua = 0.
      kappa_0_nux = 0.
      kappa_1_nue = 0.
      kappa_1_nua = 0.
      kappa_1_nux = 0.
      return
    end if

    boundsErr = enforceTableBounds(rho_cgs,temp0,ye0)

    if (boundsErr.eq.-1) then
      NeutrinoAbsorptionRateImpl = -1
    end if

    err = Absorption_cgs(rho_cgs, temp0, ye0,&
                        kappa_0_nue_cgs, kappa_0_nua_cgs,&
                        kappa_0_nux_cgs,&
                        kappa_1_nue_cgs, kappa_1_nua_cgs,&
                        kappa_1_nux_cgs)

    if (err.eq.-1) then
      NeutrinoAbsorptionRateImpl = -1
    end if

    !Unit conversion
    kappa_cgs2cactus = 1. / (cgs2cactusLength)
    kappa_0_nue = kappa_cgs2cactus * kappa_0_nue_cgs
    kappa_0_nua = kappa_cgs2cactus * kappa_0_nua_cgs
    kappa_0_nux = kappa_cgs2cactus * kappa_0_nux_cgs

    kappa_1_nue = kappa_cgs2cactus * kappa_1_nue_cgs
    kappa_1_nua = kappa_cgs2cactus * kappa_1_nua_cgs
    kappa_1_nux = kappa_cgs2cactus * kappa_1_nux_cgs

    return
END FUNCTION NeutrinoAbsorptionRateImpl

CCTK_INT FUNCTION NeutrinoDensityImpl(rho, temp, ye,&
                                n_nue, n_nua, n_nux, en_nue, en_nua, en_nux)
    use table3d_mod
    use units

    implicit none
    DECLARE_CCTK_PARAMETERS

    INTEGER :: enforceTableBounds
    INTEGER :: ierr, boundsErr, NeutrinoDens_cgs

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye
    CCTK_REAL, INTENT(OUT) :: n_nue, n_nua, n_nux, en_nue, en_nua, en_nux

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !All the emission rates are expressed per unit of volume
    REAL*8 :: rho_cgs, temp0, ye0

    NeutrinoDensityImpl = 0

    ! Conversion to cgs units
    rho_cgs = rho* cactus2cgsRho
    temp0 = temp
    ye0   = ye

    if ((rho_cgs.lt.rho_min_cgs).or.(temp0.lt.temp_min_mev)) then
      n_nue = 0.
      n_nua = 0.
      n_nux = 0.
      en_nue = 0.
      en_nua = 0.
      en_nux = 0.
      return
    end if

    boundsErr = enforceTableBounds(rho_cgs,temp0,ye0)

    if (boundsErr.eq.-1) then
      NeutrinoDensityImpl = -1
      return
    end if

    ! Call CGS backend
    ierr = NeutrinoDens_cgs(rho_cgs, temp0, ye0, n_nue, n_nua, n_nux, &
           en_nue, en_nua, en_nux)
    if(ierr.ne.0) then
        NeutrinoDensityImpl = -1
    end if

    ! Conversion from CGS units
    n_nue = n_nue / (cgs2cactusLength**3 * normfact)
    n_nua = n_nua / (cgs2cactusLength**3 * normfact)
    n_nux = n_nux / (cgs2cactusLength**3 * normfact)

    en_nue = en_nue * (mev_to_erg * cgs2cactusEnergy / cgs2cactusLength**3)
    en_nua = en_nua * (mev_to_erg * cgs2cactusEnergy / cgs2cactusLength**3)
    en_nux = en_nux * (mev_to_erg * cgs2cactusEnergy / cgs2cactusLength**3)

    return
END FUNCTION NeutrinoDensityImpl

CCTK_INT FUNCTION WeakEquilibriumImpl(rho, temp, ye,&
        n_nue, n_nua, n_nux, en_nue, en_nua, en_nux, &
        temp_eq, ye_eq, &
        n_nue_eq, n_nua_eq, n_nux_eq, en_nue_eq, en_nua_eq, en_nux_eq)
    use table3d_mod
    use units
    use weak_equilibrium_mod

    IMPLICIT NONE

    DECLARE_CCTK_PARAMETERS

    CCTK_REAL, INTENT(IN)  :: rho, temp, ye
    CCTK_REAL, INTENT(IN)  :: n_nue, n_nua, n_nux, en_nue, en_nua, en_nux
    CCTK_REAL, INTENT(OUT) :: temp_eq, ye_eq
    CCTK_REAL, INTENT(OUT) :: n_nue_eq, n_nua_eq, n_nux_eq, en_nue_eq, en_nua_eq, en_nux_eq

    CCTK_REAL :: rho0, temp0, ye0, eps0
    CCTK_REAL :: mb, AtomicMassImpl, nb

    INTEGER :: enforceTableBounds
    INTEGER :: ierr, na, boundsErr

    CCTK_REAL, DIMENSION(4) :: y_in
    CCTK_REAL, DIMENSION(4) :: e_in
    CCTK_REAL, DIMENSION(4) :: y_eq
    CCTK_REAL, DIMENSION(4) :: e_eq

    ! Conversion to cgs units
    rho0  = rho*cactus2cgsRho
    temp0 = temp
    ye0   = ye

    ! Do not do anything outside of this range
    if ((rho0.lt.rho_min_cgs).or.(temp0.lt.temp_min_mev)) then
      n_nue_eq  = 0.
      n_nua_eq  = 0.
      n_nux_eq  = 0.
      en_nue_eq = 0.
      en_nua_eq = 0.
      en_nux_eq = 0.
      WeakEquilibriumImpl = 0
      return
    end if

    ! Enforce table bounds
    boundsErr = enforceTableBounds(rho0, temp0, ye0)

    if (boundsErr.eq.-1) then
      WeakEquilibriumImpl = -1
      return
    end if

    ! Compute baryon number density. These are both in Cactus units
    mb = AtomicMassImpl()
    nb = rho0/(cactus2cgsRho*mb)

    ! Compute fractions
    y_in(1) = ye0
    y_in(2) = n_nue/nb
    y_in(3) = n_nua/nb
    y_in(4) = 0.25*n_nux/nb

    ! Compute energy (note that tab3d_eps works in Cactus units)
    eps0 = tab3d_eps(rho0/cactus2cgsRho, temp, ye)*cactus2cgsEps
    e_in(1) = rho0*(clight*clight + eps0)
    e_in(2) = en_nue*(cgs2cactusLength**3/cgs2cactusEnergy)
    e_in(3) = en_nua*(cgs2cactusLength**3/cgs2cactusEnergy)
    e_in(4) = en_nux*(cgs2cactusLength**3/cgs2cactusEnergy)

    ! Compute weak equilibrium
    call weak_equil_wnu(rho0, temp0, y_in, e_in, temp_eq, &
        y_eq, e_eq, na, ierr)
    ye_eq = y_eq(1)

    if (ierr.ne.0) then
        WeakEquilibriumImpl = -1
    else
        WeakEquilibriumImpl = 0
    end if

    ! Convert results to Cactus units
    n_nue_eq  = nb*y_eq(2)
    n_nua_eq  = nb*y_eq(3)
    n_nux_eq  = 4.0*nb*y_eq(4)
    en_nue_eq = e_eq(2)*(cgs2cactusEnergy/cgs2cactusLength**3)
    en_nua_eq = e_eq(3)*(cgs2cactusEnergy/cgs2cactusLength**3)
    en_nux_eq = e_eq(4)*(cgs2cactusEnergy/cgs2cactusLength**3)

    return
END FUNCTION WeakEquilibriumImpl
