module units
  ! unit conversion
  real(kind=8),parameter:: cactus2cgsRho    = 6.1762691458861632d+17
  real(kind=8),parameter:: cactus2cgsEps    = 8.987551787368178d+20

  real(kind=8),parameter:: cgs2cactusRho    = 1.619100425158886d-18
  real(kind=8),parameter:: cgs2cactusPress  = 1.8014921788094724d-39
  real(kind=8),parameter:: cgs2cactusEps    = 1.112650056053618d-21
  real(kind=8),parameter:: cgs2cactusMass   = 5.0278543128934301d-34
  real(kind=8),parameter:: cgs2cactusEnergy = 5.5942423830703013d-55
  real(kind=8),parameter:: cgs2cactusTime   = 203012.91587112966d0
  real(kind=8),parameter:: cgs2cactusLength = 6.7717819596091924d-06

  ! fundamental constants
  real(kind=8),parameter :: pi = 3.14159265358979d0
  real(kind=8),parameter :: mev_to_erg = 1.60217733d-6
  real(kind=8),parameter :: erg_to_mev = 6.24150636d5
  real(kind=8),parameter :: amu_cgs = 1.66053873d-24 !Atomic mass unit in g
  real(kind=8),parameter :: amu_mev = 931.49432d0 !Atomic mass unit in MeV
  real(kind=8),parameter :: kb_erg = 1.380658d-16 !Boltzmann constant in erg
  real(kind=8),parameter :: kb_mev = 8.61738568d-11 !Boltzmann constant in MeV
  real(kind=8),parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  real(kind=8),parameter :: fermi_cubed_cgs = 1.d-39
  real(kind=8),parameter :: clight = 2.99792458d10
  real(kind=8),parameter :: h_MeVs = 4.1356943d-21 !Planck constant in MeV
  real(kind=8),parameter :: me_mev = 0.510998910d0 !mass of the electron in MeV
  real(kind=8),parameter :: me_erg = 8.187108692567103d-07 !mass of the electron in erg
  real(kind=8),parameter :: sigma_0 = 1.76d-44 !cross section in unit of cm^2
  real(kind=8),parameter :: alpha = 1.23d0 !dimensionless
  real(kind=8),parameter :: multipl_nuea = 1 !multiplicity factor for nu_e and anti nu_e
  real(kind=8),parameter :: multipl_nux = 2 !multiplicity factor for nu_x
  real(kind=8),parameter :: Qnp = 1.293333d0 !neutron-proton mass difference in MeV
  real(kind=8),parameter :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
  real(kind=8),parameter :: hc_ergvcm = 1.986445683269303d-16 !hc in units of erg*cm/s
  real(kind=8),parameter :: Cv = 0.5d0 + 2.0d0*0.23d0 !vector  const. dimensionless
  real(kind=8),parameter :: Ca = 0.5d0 !axial const. dimensionless
  real(kind=8),parameter :: gamma_0 = 5.565d-2 ! dimensionless
  real(kind=8),parameter :: fsc = 1.0d0/137.036d0 !fine structure constant, dimensionless
  real(kind=8),parameter :: planck = 6.626176d-27 !Planck constant erg s
  real(kind=8),parameter :: avo = 6.0221367d23 !Avogadro's number mol^-1
  real(kind=8),parameter :: Ggrav = 6.6742d-8
  real(kind=8),parameter :: solarMass = 1.9891d33
end module units
