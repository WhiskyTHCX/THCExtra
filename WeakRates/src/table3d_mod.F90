module table3d_mod
  implicit none

  integer,save :: nrho,ntemp,nye

  ! min-max va,savelues:
  real(kind=8),save :: eos_yemin,eos_yemax
  real(kind=8),save :: eos_rhomin, eos_rhomax
  real(kind=8),save :: eos_tempmin, eos_tempmax
  real(kind=8),save :: eos_lrhomin, eos_lrhomax, eos_ltempmin, eos_ltempmax, &
       eos_lep,savesmax, eos_lepsmin
  ! Spacings
  real(kind=8),save :: dlrho, dye, dltemp
  ! basics
  integer :: nvars
  real(kind=8),allocatable :: allvariables(:,:,:,:)

  real(kind=8),save :: mass_fact

  real(kind=8),allocatable,save :: logrho(:)
  real(kind=8),allocatable,save :: logtemp(:)
  real(kind=8),allocatable,save :: yeTable(:)
end module table3d_mod
