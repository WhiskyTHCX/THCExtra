!***********************************************************************
!
!     module: weak_equilibrium_mod
!
!***********************************************************************

#include <weaktableindex.h>
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      module weak_equilibrium_mod

#ifndef FORTRAN_DISABLE_IEEE_ARITHMETIC
      use ieee_arithmetic
#endif
      use table3d_mod
      use lk_interpolations
      use units

      implicit none

!.....EOS interface.....................................................
      interface
          real(c_double) function tab3d_eps(rho, temp, ye) BIND(C, name="tab3d_eps")
              use iso_c_binding
              real (c_double), intent(in), value :: rho, temp, ye
          end function
      end interface

!.....some parameters later used in the calculations....................
      CCTK_REAL, parameter :: eps_lim       = 1.e-7     ! standard tollerance in 2D NR
      integer  , parameter :: n_cut_max     = 8         ! number of bisections of dx
      integer  , parameter :: n_max         = 100       ! Newton-Raphson max number of iterations
      integer  , parameter :: n_at          = 16        ! number of independent initial guesses

!.....deltas to compute numerical derivatives in the EOS tables.........
      CCTK_REAL, parameter :: delta_ye = 0.005
      CCTK_REAL, parameter :: delta_t  = 0.01

!.....some constants....................................................
      CCTK_REAL, parameter :: pi2   = pi*pi                  ! pi**2 [-]
      CCTK_REAL, parameter :: pref1 = 4.e0/3.d0*pi/(hc_mevcm)**3        ! 4/3 *pi/(hc)**3 [MeV^3/cm^3]
      CCTK_REAL, parameter :: pref2 = 4.e0*pi*mev_to_erg/(hc_mevcm)**3  ! 4*pi/(hc)**3 [erg/MeV^4/cm^3]
      CCTK_REAL, parameter :: cnst1 = 7.e0*pi**4/2.e1        ! 7*pi**4/20 [-]
      CCTK_REAL, parameter :: cnst5 = 7.e0*pi**4/6.e1        ! 7*pi**4/60 [-]
      CCTK_REAL, parameter :: cnst6 = 7.e0*pi**4/3.e1        ! 7*pi**4/30 [-]
      CCTK_REAL, parameter :: cnst2 = 7.e0*pi**4/5.e0        ! 7*pi**4/5 [-]
      CCTK_REAL, parameter :: cnst3 = 7.e0*pi**4/1.5e1       ! 7*pi**4/15 [-]
      CCTK_REAL, parameter :: cnst4 = 1.4e1*pi**4/1.5e1      ! 14*pi**4/15 [-]

!.....variable to switch between analytic and numerical solutions of Fermi integrals
      logical, parameter :: fermi_analytics = .true.

!.....number of points in Gauss-Legendre integration....................
      integer  , parameter :: ngl = 64
      CCTK_REAL, parameter :: gl_eps=3.d-14

      ! if I remember correctly, gaulag needs to be in real*8. Thus, I
      ! have hardcoded it here (to be checked if compiles). Since we are
      ! probably not using numerical Fermi integrals, it could be that
      ! it doesn't matter
      real*8, dimension(ngl) :: xgl
      real*8, dimension(ngl) :: wgl
      logical  , save        :: gl_init=.false.

      contains

!=======================================================================
!
!     subroutine: weak_equil_wnu
!
!     This subroutine ...
!
!=======================================================================

      subroutine weak_equil_wnu(rho,T,y_in,e_in,T_eq,y_eq,e_eq,na,ierr)

      implicit none

      CCTK_REAL              , intent(in)  :: rho
      CCTK_REAL              , intent(in)  :: T
      CCTK_REAL, dimension(4), intent(in)  :: y_in
      CCTK_REAL, dimension(4), intent(in)  :: e_in
      CCTK_REAL              , intent(out) :: T_eq
      CCTK_REAL, dimension(4), intent(out) :: y_eq
      CCTK_REAL, dimension(4), intent(out) :: e_eq
      integer                , intent(out) :: na
      integer                , intent(out) :: ierr

!=======================================================================
!
!     input:
!
!     rho  ... fluid density [g/cm^3]
!     T    ... fluid temperature [MeV]
!     y_in ... incoming abundances
!              y_in(1) ... initial electron fraction               [#/baryon]
!              y_in(2) ... initial electron neutrino fraction      [#/baryon]
!              y_in(3) ... initial electron antineutrino fraction  [#/baryon]
!              y_in(4) ... initial heavy flavor neutrino fraction  [#/baryon]
!                          The total one would be 0, so we are
!                          assuming this to be each of the single ones.
!                          Anyway, this value is useless for our
!                          calculations. We could also assume it to be
!                          the total and set it to 0
!     e_eq ... incoming energies
!              e_in(1) ... initial fluid energy, incl rest mass    [erg/cm^3]
!              e_in(2) ... initial electron neutrino energy        [erg/cm^3]
!              e_in(3) ... initial electron antineutrino energy    [erg/cm^3]
!              e_in(4) ... total initial heavy flavor neutrino energy    [erg/cm^3]
!                          This is assumed to be 4 times the energy of
!                          each single heavy flavor neutrino species
!
!     output:
!
!     T_eq ... equilibrium temperature   [MeV]
!     y_eq ... equilibrium abundances    [#/baryons]
!              y_eq(1) ... equilibrium electron fraction              [#/baryon]
!              y_eq(2) ... equilibrium electron neutrino fraction     [#/baryon]
!              y_eq(3) ... equilibrium electron antineutrino fraction [#/baryon]
!              y_eq(4) ... equilibrium heavy flavor neutrino fraction [#/baryon]
!                          see explanation above and change if necessary
!     e_eq ... equilibrium energies
!              e_eq(1) ... equilibrium fluid energy                   [erg/cm^3]
!              e_eq(2) ... equilibrium electron neutrino energy       [erg/cm^3]
!              e_eq(3) ... equilibrium electron antineutrino energy   [erg/cm^3]
!              e_eq(4) ... total equilibrium heavy flavor neutrino energy   [erg/cm^3]
!                          see explanation above and change if necessary
!     na   ... number of attempts in 2D Newton-Raphson
!     ierr ... 0 success in Newton-Raphson
!              1 failure in Newton-Raphson
!
!=======================================================================


!.....guesses for the 2D Newton-Raphson.................................
      CCTK_REAL, dimension(2)      :: x0,x1
      CCTK_REAL, dimension(n_at,2) :: vec_guess

      CCTK_REAL, dimension(3) :: mus
      CCTK_REAL, dimension(3) :: eta
      CCTK_REAL, dimension(3) :: nu_dens

      CCTK_REAL :: lrho
      CCTK_REAL :: ltemp
      CCTK_REAL :: mu_n
      CCTK_REAL :: mu_p
      CCTK_REAL :: mu_e
      CCTK_REAL :: nb
      CCTK_REAL :: mass_fact_cgs

      CCTK_REAL :: yl  ! total lepton mumber
      CCTK_REAL :: u   ! total internal energy (fluid + radiation)

!.....compute the total lepton fraction and internal energy
      yl = y_in(1) + y_in(2) - y_in(3)               ![#/baryon]
      u  = e_in(1) + e_in(2) + e_in(3) + e_in(4)     ![erg/cm^3]

!.....vector with the coefficients for the different guesses............
!     at the moment, to solve the 2D NR we assign guesses for the
!     equilibrium ye and T close to the incoming ones. This array
!     quantifies this closeness. Different guesses are used, one after
!     the other, until a solution is found. Hopefully, the first one
!     works already in most of the cases. The other ones are used as
!     backups

      vec_guess(1,:)  = (/ 1.00e0, 1.00e0  /)
      vec_guess(2,:)  = (/ 0.90e0, 1.25e0  /)
      vec_guess(3,:)  = (/ 0.90e0, 1.10e0  /)
      vec_guess(4,:)  = (/ 0.90e0, 1.00e0  /)
      vec_guess(5,:)  = (/ 0.90e0, 0.90e0  /)
      vec_guess(6,:)  = (/ 0.90e0, 0.75e0  /)
      vec_guess(7,:)  = (/ 0.75e0, 1.25e0  /)
      vec_guess(8,:)  = (/ 0.75e0, 1.10e0  /)
      vec_guess(9,:)  = (/ 0.75e0, 1.00e0  /)
      vec_guess(10,:) = (/ 0.75e0, 0.90e0  /)
      vec_guess(11,:) = (/ 0.75e0, 0.75e0  /)
      vec_guess(12,:) = (/ 0.50e0, 1.25e0  /)
      vec_guess(13,:) = (/ 0.50e0, 1.10e0  /)
      vec_guess(14,:) = (/ 0.50e0, 1.00e0  /)
      vec_guess(15,:) = (/ 0.50e0, 0.90e0  /)
      vec_guess(16,:) = (/ 0.50e0, 0.75e0  /)

      na = 0      ! counter for the number of attempts

      ! ierr is the variable that check if equilibrium has been found:
      ! ierr = 0   equilibrium found
      ! ierr = 1   equilibrium not found
      ierr = 1

      ! here we try different guesses, one after the other, until
      ! success is obtained
      do while (ierr.ne.0.and.na.lt.n_at)

        na = na + 1

!.....make an initial guess............................................
        x0(1) = vec_guess(na,1)*T        ! T guess  [MeV]
        x0(2) = vec_guess(na,2)*y_in(1)  ! ye guess [#/baryon]

!.....call the 2d Newton-Raphson........................................
        call new_raph_2dim(rho,u,yl,x0,x1,ierr)

      end do

!.....assign the output.................................................
      if (ierr.eq.0) then
        ! calculations worked
        T_eq = x1(1)
        y_eq(1) = x1(2)
      else
        ! calculations did not work
        ! write(6,*)'2D Newton-Raphson search did not work!'
        ! write(6,*)'Point log10 density [g/cm^3]: ',log10(rho)
        ! write(6,*)'Point temperature [MeV]: ',T
        ! write(6,*)'Point yl [#/baryon]: ',yl
        ! write(6,*)'Point log10 total energy [erg/cm^3]: ',log10(u)

!.....as backup plan, we assign the initial values to all outputs.......
        T_eq = T        ![MeV]
        y_eq = y_in     ![#/baryon]
        e_eq = e_in     ![erg/cm^3]
        return

!25    format(3es14.6)
!        close(6)
      end if

      ! here we want to compute the total energy and fractions in the
      ! equilibrated state

      !Interpolate the chemical potentials (stored in MeV in the table)
      lrho  = log10(rho)
      ltemp = log10(T_eq)
      mu_n = lkLinearInterpolation3d(lrho, ltemp, y_eq(1), MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, y_eq(1), MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, y_eq(1), MU_E)
      mus(1) = mu_e           ! electron chem pot including rest mass [MeV]
      mus(2) = mu_n - mu_p    ! n-p chem pot including rest masses [MeV]

      ! compute the degeneracy parameters
      call nu_deg_param_trap(t_eq,mus,eta)

      ! compute the density of the trapped neutrinos
      call dens_nu_trap(t_eq,eta,nu_dens)

      !Compute the baryon number density (mass_fact is given in MeV)
      mass_fact_cgs = mass_fact * mev_to_erg / (clight*clight)
      nb = rho / mass_fact_cgs   ![#/cm^3]

      y_eq(2) = nu_dens(1)/nb
      y_eq(3) = nu_dens(2)/nb
      y_eq(4) = nu_dens(3)/nb
      y_eq(1) = yl - y_eq(2) + y_eq(3)

      ! compute the energy density of the trapped neutrinos
      call edens_nu_trap(t_eq,eta,nu_dens)

      e_eq(2) = nu_dens(1)*mev_to_erg              ![erg/cm^3]
      e_eq(3) = nu_dens(2)*mev_to_erg              ![erg/cm^3]
      e_eq(4) = 4.*nu_dens(3)*mev_to_erg           ![erg/cm^3]
      e_eq(1) = u - e_eq(2) - e_eq(3) - e_eq(4)    ![erg/cm^3]

      ! check that the energy is positive
      if (e_eq(1).lt.nb*mass_fact_cgs*clight*clight) then
        ierr = 1
        T_eq = T
        y_eq = y_in
        e_eq = e_in
        return
      end if

      ! check that Y_e is within the range
      if (y_eq(1).lt.0.or.y_eq(1).gt.1) then
        ierr = 1
        T_eq = T
        y_eq = y_in
        e_eq = e_in
        return
      end if

      end subroutine weak_equil_wnu

!=======================================================================

!=======================================================================
!
!     subroutine: new_raph_2dim
!
!     This subroutine ...
!
!=======================================================================

      subroutine new_raph_2dim(rho,u,yl,x0,x1,ierr)

      implicit none

      CCTK_REAL              , intent(in)  :: rho
      CCTK_REAL              , intent(in)  :: u
      CCTK_REAL              , intent(in)  :: yl
      CCTK_REAL, dimension(2), intent(in)  :: x0
      CCTK_REAL, dimension(2), intent(out) :: x1
      integer                , intent(out) :: ierr

!=======================================================================
!
!     input:
!     rho ... density               [g/cm^3]
!     u   ... total internal energy [erg/cm^3]
!     yl  ... lepton number         [erg/cm^3]
!     x0  ... T and ye guess
!        x0(1) ... T               [MeV]
!        x0(2) ... ye              [#/baryon]
!
!     output:
!     x1 ... T and ye at equilibrium
!        x1(1) ... T               [MeV]
!        x1(2) ... ye              [#/baryon]
!     ierr  ...
!
!=======================================================================

      integer :: n_iter,n_cut
      CCTK_REAL :: err,err_old
      CCTK_REAL :: det
      CCTK_REAL, dimension(2)   :: y
      CCTK_REAL, dimension(2)   :: dx1
      CCTK_REAL, dimension(2)   :: x1_tmp
      CCTK_REAL, dimension(2,2) :: J
      CCTK_REAL, dimension(2,2) :: invJ

      INTEGER :: enforceTableBounds
      INTEGER :: tabBoundsFlag

      ! If true then we satisfy the Karush-Kuhn-Tucker conditions.
      ! This means that the equilibrium is out of the table and we have the best possible result.
      LOGICAL :: KKT
      ! Normal to the domain
      CCTK_REAL, dimension(2) :: norm
      CCTK_REAL :: scal
      ! Active component of the gradient
      CCTK_REAL, dimension(2) :: dxa

      ! initialize the solution
      x1 = x0
      KKT = .false.

      ! compute the initial residuals
      call func_eq_weak(rho,u,yl,x1,y)

      ! compute the error from the residuals
      call error_func_eq_weak(yl,u,y,err)

      ! initialize the iteration variable
      n_iter = 0

      ! loop until a low enough residual is found or until  a too
      ! large number of steps has been performed
      do while (err.gt.eps_lim.and.n_iter.le.n_max.and..not.KKT)

        ! compute the Jacobian
        call jacobi_eq_weak(rho,u,yl,x1,J,ierr)
        if (ierr.ne.0) then
          return
        end if

        ! compute and check the determinant of the Jacobian
        det = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        if (det.eq.0.) then
          ierr = 1
          return
          !write(6,*)'Singular determinant in weak equilibrium!'
          !stop
        end if

        ! inverte the Jacobian
        call inv_jacobi(det,J,invJ)

        ! compute the next step
        dx1(1) = - (invJ(1,1)*y(1)+invJ(1,2)*y(2))
        dx1(2) = - (invJ(2,1)*y(1)+invJ(2,2)*y(2))

        ! check if we are the boundary of the table
        if (x1(1) .eq. eos_tempmin) then
            norm(1) = -1.0
        else if (x1(1) .eq. eos_tempmax) then
            norm(1) = 1.0
        else
            norm(1) = 0.0
        endif

        if (x1(2) .eq. eos_yemin) then
            norm(2) = -1.0
        else if (x1(2) .eq. eos_yemax) then
            norm(2) = 1.0
        else
            norm(2) = 0.0
        endif

        ! Take the part of the gradient that is active (pointing within the eos domain)
        scal = norm(1)**2 + norm(2)**2
        if (scal.le.0.5) then         ! this can only happen if norm = (0, 0)
            scal = 1.0
        endif
        dxa(1) = dx1(1) - (dx1(1)*norm(1) + dx1(2)*norm(2))*norm(1)/scal
        dxa(2) = dx1(2) - (dx1(1)*norm(1) + dx1(2)*norm(2))*norm(2)/scal

        if ((dxa(1)**2 + dxa(2)**2) .lt. (eps_lim**2 * (dx1(1)**2 + dx1(2)**2))) then
            KKT = .true.
            ierr = 2
            return
        endif

        n_cut = 0
        err_old = err

        do while (n_cut.le.n_cut_max.and.err.ge.err_old)

          ! the variation of x1 is divided by an powers of 2 if the
          ! error is not decreasing along the gradient direction
          x1_tmp(1) = x1(1) + dx1(1)/2**n_cut
          x1_tmp(2) = x1(2) + dx1(2)/2**n_cut

          ! check if the next step calculation had problems
          if (isnan(x1_tmp(1))) then
            ierr = 1
            return
            !write(*,*)'x1_tmp NaN',x1_tmp(1)
            !write(*,*)'x1',x1(1)
            !write(*,*)'dx1',dx1(1)
            !write(*,*)'J',J
            !stop
          end if

          tabBoundsFlag = enforceTableBounds(rho, x1_tmp(1), x1_tmp(2))

          ! assign the new point
          x1 = x1_tmp

          ! compute the residuals for the new point
          call func_eq_weak(rho,u,yl,x1,y)

          ! compute the error
          call error_func_eq_weak(yl,u,y,err)

          ! update the bisection cut along the gradient
          n_cut = n_cut + 1

        end do

        ! update the iteration
        n_iter = n_iter+1

      end do

      ! if equilibrium has been found, set ierr=0 and return
      ! if too many attempts have been performed, set ierr=1
      if (n_iter.le.n_max) then
        ierr = 0
      else
        ierr = 1
      end if
      return

      end subroutine new_raph_2dim

!=======================================================================

!=======================================================================
!
!     function: func_eq_weak
!
!     This function ...
!
!=======================================================================


      subroutine func_eq_weak(rho,u,yl,x,y)

      implicit none

      CCTK_REAL              , intent(in)  :: rho
      CCTK_REAL              , intent(in)  :: u
      CCTK_REAL              , intent(in)  :: yl
      CCTK_REAL, dimension(2), intent(in)  :: x
      CCTK_REAL, dimension(2), intent(out) :: y

!-----------------------------------------------------------------------
!
!     Input:
!
!     rho ... density                                  [g/cm^3]
!     u   ... total (fluid+radiation) internal energy  [erg/cm^3]
!     yl  ... lepton number                            [#/baryon]
!     x   ...  array with the temperature and ye
!        x(1) ... T                                  [MeV]
!        x(2) ... ye                                 [#/baryon]
!
!     Output:
!
!     y ... array with the function whose zeros we are searching for
!
!-----------------------------------------------------------------------

      CCTK_REAL :: lrho
      CCTK_REAL :: ltemp
      CCTK_REAL :: ye
      CCTK_REAL :: mu_n
      CCTK_REAL :: mu_e
      CCTK_REAL :: mu_p
      CCTK_REAL :: rho_cu
      CCTK_REAL :: eps_cu
      CCTK_REAL :: e

      CCTK_REAL               :: nb      ! baryon density             [baryon/cm^3]
      CCTK_REAL, dimension(2) :: mus     ! chemical potential array   [MeV]
      CCTK_REAL, dimension(3) :: eta_vec ! degeneracy parameter array [-]

      CCTK_REAL               :: mass_fact_cgs
      CCTK_REAL               :: eta,eta2

      !Compute the baryon number density (mass_fact is given in MeV)
      mass_fact_cgs = mass_fact * mev_to_erg / (clight*clight)
      nb = rho / mass_fact_cgs   ![#/cm^3]

      !Interpolate the chemical potentials (stored in MeV in the table)
      lrho  = log10(rho)
      ltemp = log10(x(1))
      ye = x(2)
      mu_n = lkLinearInterpolation3d(lrho, ltemp, ye, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, ye, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, ye, MU_E)
      mus(1) = mu_e
      mus(2) = mu_n - mu_p

      !Call the EOS
      rho_cu = rho*cgs2cactusRho
      eps_cu = tab3d_eps(rho_cu, x(1), ye)
      e = rho*(clight**2 + eps_cu*cactus2cgsEps)

!.....compute the neutrino degeneracy paramater at equilibrium..........
      call nu_deg_param_trap(x(1),mus,eta_vec)
      eta = eta_vec(1)                           ![-]
      eta2 = eta*eta                             ![-]

!.....compute the function..............................................
      y(1) = x(2) + pref1*x(1)**3*eta*(pi2 + eta2)/nb - yl
      y(2) = (e+pref2*x(1)**4*((cnst5+0.5e0*eta2*(pi2+0.5e0    &
     &            *eta2))+cnst6))/u - 1.e0

      end subroutine func_eq_weak

!=======================================================================

!=======================================================================
!
!     function: error_func_eq_weak
!
!     This function ...
!
!=======================================================================

      subroutine error_func_eq_weak(yl,u,y,err)

      implicit none

      CCTK_REAL              , intent(in)  :: yl
      CCTK_REAL              , intent(in)  :: u
      CCTK_REAL, dimension(2), intent(in)  :: y
      CCTK_REAL              , intent(out) :: err

!-----------------------------------------------------------------------
!
!     Input:
!
!     yl ... lepton number                            [#/baryon]
!     u  ... total (fluid+radiation) internal energy  [erg/cm^3]
!     y  ... array with residuals                     [-]
!
!     Output:
!
!     err ... error associated with the residuals     [-]
!
!-----------------------------------------------------------------------

!.....since the first equation is has yl as constant, we normalized the error to it.
!     since the second equation was normalized wrt u, we divide it by 1.
!     the modulus of the two contributions are then summed

      err = abs(y(1)/yl) + abs(y(2)/1.)

      return

      end subroutine error_func_eq_weak



!=======================================================================
!
!     subroutine: jacobi_eq_weak
!
!=======================================================================


      subroutine jacobi_eq_weak(rho,u,yl,x,J,ierr)

      implicit none

      CCTK_REAL                , intent(in)  :: rho
      CCTK_REAL                , intent(in)  :: u
      CCTK_REAL                , intent(in)  :: yl
      CCTK_REAL, dimension(2)  , intent(in)  :: x
      CCTK_REAL, dimension(2,2), intent(out) :: J
      integer                  , intent(out) :: ierr

!-----------------------------------------------------------------------
!
!     Input:
!
!     rho ... density                   [g/cm^3]
!     u   ... total energy density      [erg/cm^3]
!     yl  ... lepton fraction           [#/baryon]
!     x   ... array with T and ye
!        x(1) ... T                     [MeV]
!        x(2) ... ye                    [#/baryon]
!
!     Output:
!
!     J ... Jacobian for the 2D Newton-Raphson
!
!-----------------------------------------------------------------------

      CCTK_REAL :: mass_fact_cgs
      CCTK_REAL :: nb,t,ye
      CCTK_REAL :: eta,eta2
      CCTK_REAL :: dedt,dedye
      CCTK_REAL :: detadt,detadye
      CCTK_REAL, dimension(3) :: eta_vec

      CCTK_REAL :: lrho,ltemp
      CCTK_REAL :: mu_e,mu_p,mu_n

      !integer :: ierr
      !CCTK_REAL :: x1,x2
      CCTK_REAL, dimension(2) :: mus

      !Interpolate the chemical potentials (stored in MeV in the table)
      lrho  = log10(rho)
      ltemp = log10(x(1))
      t = x(1)
      ye = x(2)
      mu_n = lkLinearInterpolation3d(lrho, ltemp, ye, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, ye, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, ye, MU_E)
      mus(1) = mu_e               ! electron chemical potential (w rest mass) [MeV]
      mus(2) = mu_n - mu_p        ! n minus p chemical potential (w rest mass) [MeV]
      ! compute the degeneracy parameters
      call nu_deg_param_trap(t,mus,eta_vec)
      eta = eta_vec(1)       ! electron neutrinos degeneracy parameter
      eta2 = eta*eta

!.....compute the gradients of eta and of the internal energy...........
      call eta_e_gradient(rho,t,ye,eta,detadt,detadye,dedt,dedye,ierr)
      if (ierr.ne.0) then
        return
      end if

      !Compute the baryon number density (mass_fact is given in MeV)
      mass_fact_cgs = mass_fact * mev_to_erg / (clight*clight)
      nb = rho / mass_fact_cgs   ![#/cm^3]

!     J(1,1): df1/dt
!     J(1,2): df1/dye
!     J(2,1): df2/dt
!     J(2,2): df2/dye

!.....compute the Jacobian.............................................
      J(1,1) = pref1/nb*t**2*(3.e0*eta*(pi2+eta2)+t*(pi2+3.e0*eta2)*detadt)
      J(1,2) = 1.e0+pref1/nb*t**3*(pi2+3.e0*eta2)*detadye

      J(2,1) = (dedt+pref2*t**3*(cnst3+cnst4+2.e0*&
     &                   eta2*(pi2+0.5*eta2)+eta*t*(pi2+eta2)*detadt))/u
      J(2,2) = (dedye+pref2*t**4*eta*(pi2+eta2)*detadye)/u

!.....check on the degeneracy parameters and temperature................
      if (isnan(eta)) then
        ierr = 1
        return
        !write(*,*)'eta',eta
      end if

      if (isnan(detadt)) then
        ierr = 1
        return
        !write(*,*)'detadt',detadt
      end if

      if (isnan(t)) then
        ierr = 1
        return
        !write(*,*)'t',t
      end if

      ierr = 0

      end subroutine jacobi_eq_weak

!=======================================================================

!=======================================================================
!
!     subroutine: eta_e_gradient
!
!     this subroutine computes the gradient of the degeneracy parameter
!     and of the fluid internal energy with respect to temperature and ye
!
!=======================================================================

      subroutine eta_e_gradient(rho,t,ye,eta,detadt,detadye,dedt,dedye,ierr)

      implicit none

      CCTK_REAL, intent(in)  :: rho
      CCTK_REAL, intent(in)  :: t
      CCTK_REAL, intent(in)  :: ye
      CCTK_REAL, intent(in)  :: eta

      CCTK_REAL, intent(out) :: detadt
      CCTK_REAL, intent(out) :: detadye
      CCTK_REAL, intent(out) :: dedt
      CCTK_REAL, intent(out) :: dedye
      integer  , intent(out) :: ierr

!-----------------------------------------------------------------------
!
!     Input:
!     rho ... density             [g/cm^3]
!     t   ... temperature         [MeV]
!     ye  ... electron fraction   [#/baryon]
!     eta ... electron neutrino degeneracy parameter at equilibrium [-]
!
!     Output:
!     detadt  ... derivative of eta wrt T (for ye and rho fixed)             [1/MeV]
!     detadye ... derivative of eta wrt ye (for T and rho fixed)             [-]
!     dedt  ... derivative of internal energy wrt T (for ye and rho fixed) [erg/cm^3/MeV]
!     dedye ... derivative of internal energy wrt ye (for T and rho fixed) [erg/cm^3]
!
!-----------------------------------------------------------------------

      CCTK_REAL, dimension(2) :: mus1, mus2

      CCTK_REAL :: ye1, ye2  ! values used for the derivatives
      CCTK_REAL :: t1, t2    ! values used for the derivatives
      CCTK_REAL :: dmuedt,dmuedye
      CCTK_REAL :: dmuhatdt,dmuhatdye
      CCTK_REAL :: x1,x2

      CCTK_REAL :: lrho,ltemp,tv,yev
      CCTK_REAL :: rho_cu, eps_cu
      CCTK_REAL :: mu_e,mu_p,mu_n
      CCTK_REAL :: e1,e2

!.....gradients are computed numerically. To do it, we consider small
!     variations in ye and temperature, and we compute the detivative
!     using finite differencing. The real limitation is that this way
!     relies on the EOS table interpolation procedure

!     the goal of this part is to obtain chemical potentials (mus1 and
!     mus2) in two points close to the point we are considering, first
!     varying wrt ye and then wrt T

!.....vary the electron fraction........................................

      ! these are the two calls to the EOS. The goal here is to get
      ! the fluid internal energy and the chemical potential for
      ! electrons and the difference between neutron and proton
      ! chemical potential (usually called mu_hat) for two points with
      ! slightly different ye

      lrho  = log10(rho)
      ltemp = log10(t)

      ! first, for ye slightly smaller
      ye1 = max(ye - delta_ye, eos_yemin)
      yev = ye1
      mu_n = lkLinearInterpolation3d(lrho, ltemp, yev, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, yev, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, yev, MU_E)

      rho_cu = rho*cgs2cactusRho
      eps_cu = tab3d_eps(rho_cu, t, yev)
      e1 = rho*(clight**2 + eps_cu*cactus2cgsEps)

      mus1(1) = mu_e
      mus1(2) = mu_n - mu_p

      ! second, for ye slightly larger
      ye2 = min(ye + delta_ye, eos_yemax)
      yev = ye2
      mu_n = lkLinearInterpolation3d(lrho, ltemp, yev, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, yev, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, yev, MU_E)

      eps_cu = tab3d_eps(rho_cu, t, yev)
      e2 = rho*(clight**2 + eps_cu*cactus2cgsEps)

      mus2(1) = mu_e
      mus2(2) = mu_n - mu_p

!.....compute numerical derivaties......................................
      dmuedye   = (mus2(1)-mus1(1))/(ye2 - ye1)
      dmuhatdye = (mus2(2)-mus1(2))/(ye2 - ye1)
      dedye     = (e2-e1)/(ye2 - ye1)

!.....vary the temperature..............................................
      t1 = max(t - delta_t, eos_tempmin)
      t2 = min(t + delta_t, eos_tempmax)

      ! these are the two other calls to the EOS. The goal here is to get
      ! the fluid internal energy and the chemical potential for
      ! electrons and the difference between neutron and proton
      ! chemical potential (usually called mu_hat) for two points with
      ! slightly different t

      ! ye is the original one

      ! first, for t slightly smaller
      tv = t1
      ltemp = log10(t-delta_t)
      mu_n = lkLinearInterpolation3d(lrho, ltemp, ye, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, ye, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, ye, MU_E)

      eps_cu = tab3d_eps(rho_cu, tv, ye)
      e1 = rho*(clight**2 + eps_cu*cactus2cgsEps)

      mus1(1) = mu_e
      mus1(2) = mu_n - mu_p

      ! second, for t slightly larger
      tv = t2
      ltemp = log10(t+delta_t)
      mu_n = lkLinearInterpolation3d(lrho, ltemp, ye, MU_N)
      mu_p = lkLinearInterpolation3d(lrho, ltemp, ye, MU_P)
      mu_e = lkLinearInterpolation3d(lrho, ltemp, ye, MU_E)

      eps_cu = tab3d_eps(rho_cu, tv, ye)
      e2 = rho*(clight**2 + eps_cu*cactus2cgsEps)

      mus2(1) = mu_e
      mus2(2) = mu_n - mu_p

!.....compute the derivatives wrt temperature...........................
      dmuedt   = (mus2(1) - mus1(1))/(t2 - t1)
      dmuhatdt = (mus2(2) - mus1(2))/(t2 - t1)
      dedt     = (e2   - e1  )/(t2 - t1)

!.....combine the eta derivatives.......................................
      detadt  = (-eta + dmuedt - dmuhatdt)/t    ![1/MeV]
      detadye = (dmuedye - dmuhatdye)/t         ![-]

!.....check if the derivative has a problem.............................
      if (isnan(detadt)) then
        ierr = 1
        return
        !write(*,*)'problem with eta: ',eta
        !write(*,*)'mue1: ',mus1(1)
        !write(*,*)'mue2: ',mus2(1)
        !write(*,*)'dmuedt: ',dmuedt
        !write(*,*)'dmuhatdt: ',dmuhatdt
        !write(*,*)'rho: ',rho
        !write(*,*)'ye: ',ye
        !write(*,*)'t: ',t
      end if

      ierr = 0

      end subroutine eta_e_gradient

!=======================================================================

!=======================================================================
!
!     subroutine: inv_jacobi
!
!     This subroutine inverts the Jacobian matrix, assuming it to be a
!     2x2 matrix
!
!=======================================================================

      subroutine inv_jacobi(det,J,invJ)

      implicit none

      CCTK_REAL,                 intent(in)  :: det
      CCTK_REAL, dimension(2,2), intent(in)  :: J
      CCTK_REAL, dimension(2,2), intent(out) :: invJ

!=======================================================================
!
!     Input:
!     det ... determinant of the Jacobian matrix
!     J   ... Jacobian matrix
!
!     Output:
!     invJ ... inverse of the Jacobian matrix
!
!=======================================================================

      invJ(1,1) =  J(2,2)
      invJ(2,2) =  J(1,1)
      invJ(1,2) = -J(1,2)
      invJ(2,1) = -J(2,1)
      invJ = invJ/det

      end subroutine inv_jacobi

!=======================================================================

!=======================================================================
!
!     subroutine: nu_deg_param_trap
!
!     In this subroutine, we compute the neutrino degeneracy parameters
!     assuming weak and thermal equilibrium, i.e. using as input the
!     local thermodynamical properties
!
!=======================================================================

      subroutine nu_deg_param_trap(temp_m,chem_pot,eta)

      implicit none

      CCTK_REAL              , intent(in)  :: temp_m
      CCTK_REAL, dimension(2), intent(in)  :: chem_pot
      CCTK_REAL, dimension(3), intent(out) :: eta

!-----------------------------------------------------------------------
!
!     Input:
!     temp_m   ----> local matter temperature [MeV]
!     chem_pot ----> matter chemical potential [MeV]
!                    chem_pot(1): electron chemical potential (w rest mass)
!                    chem_pot(2): n minus p chemical potential (w/o rest mass)
!
!     Output:
!     eta      ----> neutrino degeneracy parameters [-]
!                    eta(1): electron neutrino
!                    eta(2): electron antineutrino
!                    eta(3): mu and tau neutrinos
!
!-----------------------------------------------------------------------

      if (temp_m.gt.0.d0) then
        eta(1) = (chem_pot(1) - chem_pot(2))/temp_m            ![-]
        eta(2) = - eta(1)                                      ![-]
        eta(3) = 0.e0                                          ![-]
      else
        !write(*,*)'Problem with the temperature in computing eta_nu'
        !write(*,*)'temp',temp_m
        eta(:) = 0.e0                                          ![-]
      end if

      end subroutine nu_deg_param_trap

!=======================================================================

!=======================================================================
!
!     subroutine: dens_nu_trap
!
!=======================================================================

      subroutine dens_nu_trap(temp_m,eta_nu,nu_dens)

      implicit none

      CCTK_REAL              , intent(in)  :: temp_m
      CCTK_REAL, dimension(3), intent(in)  :: eta_nu
      CCTK_REAL, dimension(3), intent(out) :: nu_dens

!-----------------------------------------------------------------------
!     In this subroutine, we compute the neutrino densities in equilibrium
!     conditions, using as input the local thermodynamical properties
!
!     Input:
!     temp_m   ----> local matter temperature [MeV]
!     eta_nu   ----> neutrino degeneracy parameter [-]
!
!     Output:
!     nu_dens ----> neutrino density [particles/cm^3]
!
!-----------------------------------------------------------------------

      CCTK_REAL, parameter :: pref=4.d0*pi/(hc_mevcm)**3  ![1/MeV^3/cm^3]
      CCTK_REAL ::  f2
      CCTK_REAL ::  temp_m3
      CCTK_REAL :: fermi2

      integer :: it

      temp_m3 = temp_m * temp_m * temp_m        ![MeV^3]
      do it=1,3
        if (fermi_analytics) then
          f2 = fermi2(eta_nu(it))
          !call f2_analytic(eta_nu(it),f2)
        else
          call fermiint(2.d0,eta_nu(it),f2)
        end if
        nu_dens(it) = pref * temp_m3 * f2      ![#/cm^3]
      end do

      end subroutine dens_nu_trap

!=======================================================================

!=======================================================================
!
!     subroutine: dens_nu_trap
!
!=======================================================================

      subroutine edens_nu_trap(temp_m,eta_nu,enu_dens)

      implicit none

      CCTK_REAL              , intent(in)  :: temp_m
      CCTK_REAL, dimension(3), intent(in)  :: eta_nu
      CCTK_REAL, dimension(3), intent(out) :: enu_dens

!-----------------------------------------------------------------------
!     In this subroutine, we compute the neutrino densities in equilibrium
!     conditions, using as input the local thermodynamical properties
!
!     Input:
!     temp_m   ----> local matter temperature [MeV]
!     eta_nu   ----> neutrino degeneracy parameter [-]
!
!     Output:
!     enu_dens ----> neutrino density [MeV/cm^3]
!
!-----------------------------------------------------------------------

      CCTK_REAL, parameter :: pref=4.d0*pi/(hc_mevcm)**3  ![1/MeV^3/cm^3]
      CCTK_REAL :: f3
      CCTK_REAL :: temp_m4
      CCTK_REAL :: fermi3

      integer :: it

      temp_m4 = temp_m * temp_m * temp_m * temp_m     ![MeV^4]
      do it=1,3
        if (fermi_analytics) then
          f3 = fermi3(eta_nu(it))
          !call f3_analytic(eta_nu(it),f3)
        else
          call fermiint(3.d0,eta_nu(it),f3)
        end if
        enu_dens(it) = pref * temp_m4 * f3
      end do

      end subroutine edens_nu_trap

!=======================================================================

!=======================================================================
!
!     Subroutine: f2_analytic
!
!=======================================================================

      subroutine f2_analytic(eta,f2)

      implicit none

      CCTK_REAL, intent(in)  :: eta
      CCTK_REAL, intent(out) :: f2

      if (eta.gt.1.e-3) then
        f2 = (eta**3/3.e0 + 3.2899*eta)/(1.-exp(-1.8246e0*eta))
      else
        f2 = 2.e0*exp(eta)/(1.e0+0.1092*exp(0.8908e0*eta))
      end if

      end subroutine f2_analytic

!=======================================================================

!=======================================================================
!
!     Subroutine: f3_analytic
!
!=======================================================================

      subroutine f3_analytic(eta,f3)

      implicit none

      CCTK_REAL, intent(in)  :: eta
      CCTK_REAL, intent(out) :: f3

      CCTK_REAL :: eta2,eta4

      if (eta.gt.1.e-3) then
        eta2 = eta*eta
        eta4 = eta2*eta2
        f3 = (eta4/4.e0 + 4.9348e0*eta2 + 1.13644e0)/(1.+exp(-1.9039e0*eta))
      else
        f3 = 6.e0*exp(eta)/(1.e0+0.0559e0*exp(0.9069e0*eta))
      end if

      end subroutine f3_analytic

!=======================================================================

!=======================================================================
!
!     Fermi integral calculation
!
!=======================================================================

      subroutine fermiint(k,eta,f)

!=======================================================================
! This subroutine calculates the Fermi integral function, once the order,
! the point and the order of Gauss-Legendre integration have been
! specified
!=======================================================================

      implicit none

      CCTK_REAL, intent(in) :: k
      CCTK_REAL, intent(in) :: eta

      CCTK_REAL, intent(out) :: f

!.......................................................................
!     Input variables:
!     k     ----> order of the Fermi integral
!     eta   ----> point where to evaluate the Fermi integral function
!
!     Output variables:
!     f     ----> F_k (eta)
!.......................................................................

      integer :: i
      CCTK_REAL :: fxi

!.....initialize function to 0
      f = 0.e0
      if (.not.gl_init) then
        call gauleg(0.d0,1.d0)
        gl_init = .true.
      end if
      do i=1,ngl
         call kernel(k,eta,xgl(i),fxi)
         f = f + wgl(i) * fxi
      end do

      end subroutine fermiint

!=======================================================================

!=======================================================================

      subroutine kernel(k,eta,x,fcompx)
!....................................................................
! This subroutine calculate the kernel for the calculation of the Fermi
! integral of order k, shift factor eta, in the point x, for the
! numerical integration.
!
!     Input:
!     x ------> abscissa
!     k ------> order of the Fermi integral
!     eta ----> shift parameter of the Fermi integral
!
!     Output:
!     fcompx -----> function
!....................................................................

       implicit none
       CCTK_REAL, intent(in) :: x
       CCTK_REAL, intent(in) :: k
       CCTK_REAL, intent(in) :: eta
       CCTK_REAL, intent(out) :: fcompx

       CCTK_REAL :: f
       CCTK_REAL :: t
       CCTK_REAL :: s, a

!....................................................................
!      the parameter a describes how the nodes should be sampled
!      (i.e. instead of [0,1] and [1,infty], we use the intervals
!      [0,eta] and [eta,infty] here), except eta<1.0
!....................................................................
       a = max(1.e0,eta)

       f = a * (x*a)**k * fermi(x*a-eta)
       t = a/x
       s = a * t**k * fermi(t-eta)/(x*x)
       fcompx = f + s

       end subroutine kernel

!=======================================================================

!=======================================================================
!
!
!
!=======================================================================

      function fermi(arg)

      implicit none

      CCTK_REAL, intent(in) :: arg
      CCTK_REAL :: fermi

!-----------------------------------------------------------------------

      CCTK_REAL :: tmp

      if (arg.gt.0.e0) then
        tmp = exp(-arg)
        fermi = tmp/(tmp+1.e0)
      else
        tmp = exp(arg)
        fermi = 1./(tmp+1.e0)
      endif

      end function fermi

!=======================================================================

!=======================================================================
!
!     Subroutine: gauleg
!
!=======================================================================

      subroutine gauleg(x1,x2)

      implicit none

!     Note: iof I remember correctly, this subroutine needs to be in
!     double. I am forcing it to be real*8

      real*8, intent(in) :: x1,x2

      integer :: i,j,m
      real*8 :: p1,p2,p3,pp,xl,xm,z,z1

      m = (ngl+1)/2
      xm = 0.5d0*(x2+x1)
      xl = 0.5d0*(x2-x1)
      do i=1,m
         z = dcos(pi*(dble(i)-0.25d0)/(dble(ngl)+0.5d0))
         z1 = 0.0
         do while(abs(z-z1).gt.gl_eps)
            p1 = 1.0d0
            p2 = 0.0d0
            do j=1,ngl
               p3 = p2
               p2 = p1
               p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)   &
     &              / dble(j)
            end do
            pp = dble(ngl)*(z*p1-p2)/(z*z-1.0d0)
            z1 = z
            z = z1 - p1/pp
         end do
         xgl(i) = xm - xl*z
         xgl(ngl+1-i) = xm + xl*z
         wgl(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
         wgl(ngl+1-i) = wgl(i)
      end do

      end subroutine gauleg

!=======================================================================

      end module weak_equilibrium_mod

!***********************************************************************
