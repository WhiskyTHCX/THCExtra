module eos3dmod

  implicit none

  integer,save :: nrho,ntemp,nye

  real(kind=8) :: energy_shift = 0.0d0

  real(kind=8) :: precision = 1.0d-9

  ! min-max values:
  real(kind=8) :: eos_yemin,eos_yemax
  real(kind=8) :: eos_rhomin, eos_rhomax
  real(kind=8) :: eos_tempmin, eos_tempmax
  real(kind=8) :: eos_lrhomin, eos_lrhomax, eos_ltempmin, eos_ltempmax, &
       eos_lepsmax, eos_lepsmin
  ! Spacings
  real(kind=8) :: dlrho, dye, dltemp
  ! basics
  integer :: nvars
  real(kind=8),allocatable :: allvariables(:,:,:,:)

  ! index variable mapping:
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar
  ! 19 -> gamma

  real(kind=8),allocatable,save :: logrho(:)
  real(kind=8),allocatable,save :: logtemp(:)
  real(kind=8),allocatable,save :: yeTable(:)

  ! constants
  real(kind=8),save :: mev_to_erg = 1.60217733d-6
  real(kind=8),save :: amu_cgs = 1.66053873d-24
  real(kind=8),save :: amu_mev = 931.49432d0
  real(kind=8),save :: pi = 3.14159265358979d0
  real(kind=8),save :: ggrav = 6.672d-8
  real(kind=8),save :: temp_mev_to_kelvin = 1.1604447522806d10
  real(kind=8),save :: clight = 2.99792458d10
  real(kind=8),save :: kb_erg = 1.380658d-16
  real(kind=8),save :: kb_mev = 8.61738568d-11

  !? TODO this part need to be pushed inside the TableReader
  !constants to convert from cactus units (G=c=1, M_sun=1) to cgs*/

  real(kind=8):: cactus2cgsRho    = 6.1762691458861632d+17
  real(kind=8):: cactus2cgsEps    = 8.987551787368178d+20

  real(kind=8):: cgs2cactusRho    = 1.619100425158886e-18
  real(kind=8):: cgs2cactusPress  = 1.8014921788094724d-39
  real(kind=8):: cgs2cactusEps    = 1.112650056053618e-21
  real(kind=8):: cgs2cactusMass   = 5.0278543128934301d-34
  real(kind=8):: cgs2cactusEnergy = 5.5942423830703013d-55
  real(kind=8):: cgs2cactusTime   = 203012.91587112966
  real(kind=8):: cgs2cactusLength = 6.7717819596091924d-06

end module eos3dmod
