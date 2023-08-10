! #include <table3d_mod.F90>

program weakrates_calc_rates
    ! -----------------------------------------------------------------
    ! Variable declarations
    ! -----------------------------------------------------------------
    use table3d_mod
    use units
    use weak_equilibrium_mod
    use iso_c_binding

    implicit none
    character*200 weak_table_file_name
    character*200 eos_table_file_name
    integer ierr
    integer weak_table_reader
    integer TableReader
    integer Emissions_cgs
    integer Opacities_cgs
    integer Absorption_cgs
    integer NeutrinoDens_cgs

    ! EOS
    real*8 rho, temp, ye
    real(c_double) rho_cu, temp_cu, ye_cu
    real*8 csnd, eps, entropy, press
    real*8 mb, normfact, nb

    ! weak rates
    real*8 R_nue, R_nua, R_nux
    real*8 Q_nue, Q_nua, Q_nux
    real*8 kappa_0_nue, kappa_0_nua, kappa_0_nux
    real*8 kappa_1_nue, kappa_1_nua, kappa_1_nux
    real*8 abs_0_nue, abs_0_nua, abs_0_nux
    real*8 abs_1_nue, abs_1_nua, abs_1_nux
    real*8 n_nue, n_nua, n_nux
    real*8 e_nue, e_nua, e_nux

    ! weak equilibrium interface
    real*8 ynue, ynua, enue, enua, etot
    real*8, dimension(4) :: y_in
    real*8, dimension(4) :: e_in
    real*8 temp_eq
    real*8, dimension(4) :: y_eq
    real*8, dimension(4) :: e_eq
    integer na

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    interface
        real(c_double) function tab3d_csnd2_from_temp(rho, temp, ye) BIND(C)
            use iso_c_binding
            real (c_double), intent(in), value :: rho, temp, ye
        end function
    end interface
    interface
        real(c_double) function tab3d_entropy_from_temp(rho, temp, ye) BIND(C)
            use iso_c_binding
            real (c_double), intent(in), value :: rho, temp, ye
        end function
    end interface
    interface
        real(c_double) function tab3d_press_from_temp(rho, temp, ye) BIND(C)
            use iso_c_binding
            real (c_double), intent(in), value :: rho, temp, ye
        end function
    end interface
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real*8 r_cgs2cactus, q_cgs2cactus, edens_cgs2cactus
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------
    open(111, file='weakrates_calc_rates.inp', status='unknown', &
            form='formatted', action='read')
    read(111, 200) eos_table_file_name
    read(111, 200) weak_table_file_name
    read(111, 201) rho
    read(111, 201) temp
    read(111, 201) ye
    read(111, 201) ynue
    read(111, 201) ynua
    read(111, 201) enue
    read(111, 201) enua
    close(111)

    write(*,*) 'EOS table file name  : ', trim(eos_table_file_name)
    write(*,*) 'Weak table file name : ', trim(weak_table_file_name)
    write(*,*) 'rho  [g/cm^3]        = ', rho
    write(*,*) 'temp [MeV]           = ', temp
    write(*,*) 'Y_e                  = ', ye
    write(*,*) 'Y_nue                = ', ynue
    write(*,*) 'Y_nua                = ', ynua
    write(*,*) 'enue                 = ', enue
    write(*,*) 'enua                 = ', enua

200 format(a200)
201 format(e16.9)
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Read table and compute values
    ! -----------------------------------------------------------------
    ierr = weak_table_reader(weak_table_file_name)
    if(ierr.ne.0) then
        write(6,*) 'Failed reading ', trim(weak_table_file_name)
        stop
    end if

    ierr = TableReader(eos_table_file_name)
    if(ierr.ne.0) then
        write(6,*) 'Failed reading ', trim(eos_table_file_name)
        stop
    end if
    write(6,*) ""

    ! The EOS table is in Cactus units
    rho_cu = rho/cactus2cgsRho
    temp_cu = temp
    ye_cu = ye

    ! Call the EOS
    eps = tab3d_eps(rho_cu, temp_cu, ye_cu)/cgs2cactusEps
    csnd = sqrt(tab3d_csnd2_from_temp(rho_cu, temp_cu, ye_cu))*clight
    entropy = tab3d_entropy_from_temp(rho_cu, temp_cu, ye_cu)
    press = tab3d_press_from_temp(rho_cu, temp_cu, ye_cu)/cgs2cactusPress

    ! Reference baryon mass in Cactus units
    normfact = 1d50
    mb = normfact * cgs2cactusMass * mass_fact * mev_to_erg / (clight*clight)

    ! Baryon number density in CGS
    nb = rho/(mass_fact * mev_to_erg / (clight*clight))

    ! Call weak equilibrium module
    y_in(1) = ye
    y_in(2) = ynue
    y_in(3) = ynua
    y_in(4) = 0.
    e_in(1) = rho*(clight**2 + eps)
    e_in(2) = enue
    e_in(3) = enua
    e_in(4) = 0.
    etot    = e_in(1) + e_in(2) + e_in(3) + e_in(4)
    call weak_equil_wnu(rho, temp, y_in, e_in, &
            temp_eq, y_eq, e_eq, na, ierr)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the weak equilibrium'
        write(6,*) ""
    end if
    if(abs(etot - e_eq(1) - e_eq(2) - e_eq(3) - e_eq(4))/etot.gt.1e-15) then
        write(6,*) 'Energy conservation violated!'
        write(6,*) 'Delta E / E = ', abs(etot - e_eq(1) - e_eq(2) - e_eq(3) - e_eq(4))/etot
        write(6,*) ""
    end if

    ierr = Emissions_cgs(rho, temp, ye, &
        R_nue, R_nua, R_nux, Q_nue, Q_nua, Q_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the emission rates'
        write(6,*) ""
    end if

    ierr = Opacities_cgs(rho, temp, ye, &
            kappa_0_nue, kappa_0_nua, kappa_0_nux, &
            kappa_1_nue, kappa_1_nua, kappa_1_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the opacities'
        write(6,*) ""
    end if

    ierr = Absorption_cgs(rho, temp, ye, &
            abs_0_nue, abs_0_nua, abs_0_nux, &
            abs_1_nue, abs_1_nua, abs_1_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the absorption rates'
        write(6,*) ""
    end if

    ierr = NeutrinoDens_cgs(rho, temp, ye, &
           n_nue, n_nua, n_nux, e_nue, e_nua, e_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the neutrino densities'
        write(6,*) ""
    end if
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Output results CGS
    ! -----------------------------------------------------------------
    write(*,*) '------------ CGS -------------------------------------'

    write(*,*) 'Thermodynamical variables'
    write(*,*) '       eps      = ', eps
    write(*,*) '         e      = ', e_in(1)
    write(*,*) '      csnd      = ', csnd
    write(*,*) '   entropy      = ', entropy
    write(*,*) '     press      = ', press

    write(*,*) 'Weak equilibrium'
    write(*,*) '     temp [MeV] = ', temp_eq
    write(*,*) '         ye     = ', y_eq(1)
    write(*,*) '       ymue     = ', y_eq(2)
    write(*,*) '       ynua     = ', y_eq(3)
    write(*,*) '       ynux     = ', 4*y_eq(4)
    write(*,*) '   e [erg/cm^3] = ', e_eq(1)
    write(*,*) '       enue     = ', e_eq(2)
    write(*,*) '       enua     = ', e_eq(3)
    write(*,*) '       enux     = ', e_eq(4)
    write(*,*) '       znue     = ', e_eq(2)/etot
    write(*,*) '       znua     = ', e_eq(3)/etot
    write(*,*) '       znux     = ', e_eq(4)/etot
    write(*,*) '   nnue [cm^-3] = ', y_eq(2)*nb
    write(*,*) '   nnua [cm^-3] = ', y_eq(3)*nb
    write(*,*) '   nnux [cm^-3] = ', 4*y_eq(4)*nb


    write(*,*) 'Number emission rates [1/(cm^3 s)]'
    write(*,*) '        nue     = ', R_nue
    write(*,*) '        nua     = ', R_nua
    write(*,*) '        nux     = ', R_nux

    write(*,*) 'Energy emission rates [erg/(cm^3 s)]'
    write(*,*) '        nue     = ', Q_nue * mev_to_erg
    write(*,*) '        nua     = ', Q_nua * mev_to_erg
    write(*,*) '        nux     = ', Q_nux * mev_to_erg

    write(*,*) 'Number opacities [1/cm]'
    write(*,*) '        nue     = ', kappa_0_nue
    write(*,*) '        nua     = ', kappa_0_nua
    write(*,*) '        nux     = ', kappa_0_nux

    write(*,*) 'Energy opacities [1/cm]'
    write(*,*) '        nue     = ', kappa_1_nue
    write(*,*) '        nua     = ', kappa_1_nua
    write(*,*) '        nux     = ', kappa_1_nux

    write(*,*) 'Number absorption coefficients [1/cm]'
    write(*,*) '        nue     = ', abs_0_nue
    write(*,*) '        nua     = ', abs_0_nua
    write(*,*) '        nux     = ', abs_0_nux

    write(*,*) 'Energy absorption coefficients [1/cm]'
    write(*,*) '        nue     = ', abs_1_nue
    write(*,*) '        nua     = ', abs_1_nua
    write(*,*) '        nux     = ', abs_1_nux

    write(*,*) 'Neutrino number density [1/cm^3]'
    write(*,*) '        nue     = ', n_nue
    write(*,*) '        nua     = ', n_nua
    write(*,*) '        nux     = ', n_nux

    write(*,*) 'Neutrino energy density [erg/cm^3]'
    write(*,*) '        nue     = ', e_nue * mev_to_erg
    write(*,*) '        nua     = ', e_nua * mev_to_erg
    write(*,*) '        nux     = ', e_nux * mev_to_erg

    write(*,*) '------------ Cactus ----------------------------------'

    write(*,*) '        mb      = ', mb

    write(*,*) 'EOS input'
    write(*,*) '       rho      = ', rho_cu
    write(*,*) '      temp      = ', temp_cu
    write(*,*) '        ye      = ', ye_cu

    write(*,*) 'Thermodynamical variables'
    write(*,*) '       eps      = ', eps*cgs2cactusEps
    write(*,*) '      csnd      = ', csnd/clight
    write(*,*) '   entropy      = ', entropy
    write(*,*) '     press      = ', press*cgs2cactusPress
    write(*,*) '         e      = ', e_in(1)*cgs2cactusEnergy

    write(*,*) 'Number emission rates'
    r_cgs2cactus = 1. / (cgs2cactusTime*cgs2cactusLength**3)
    write(*,*) '        nue     = ', R_nue * r_cgs2cactus
    write(*,*) '        nua     = ', R_nua * r_cgs2cactus
    write(*,*) '        nux     = ', R_nux * r_cgs2cactus

    write(*,*) 'Energy emission rates'
    q_cgs2cactus = mev_to_erg*cgs2cactusenergy/ &
                   (cgs2cactusTime * cgs2cactusLength**3)
    write(*,*) '        nue     = ', Q_nue * q_cgs2cactus
    write(*,*) '        nua     = ', Q_nue * q_cgs2cactus
    write(*,*) '        nux     = ', Q_nue * q_cgs2cactus

    write(*,*) 'Number opacities'
    write(*,*) '        nue     = ', kappa_0_nue / cgs2cactusLength
    write(*,*) '        nua     = ', kappa_0_nua / cgs2cactusLength
    write(*,*) '        nux     = ', kappa_0_nux / cgs2cactusLength

    write(*,*) 'Energy opacities'
    write(*,*) '        nue     = ', kappa_1_nue / cgs2cactusLength
    write(*,*) '        nua     = ', kappa_1_nua / cgs2cactusLength
    write(*,*) '        nux     = ', kappa_1_nux / cgs2cactusLength

    write(*,*) 'Number absorption coefficients'
    write(*,*) '        nue     = ', abs_0_nue / cgs2cactusLength
    write(*,*) '        nua     = ', abs_0_nua / cgs2cactusLength
    write(*,*) '        nux     = ', abs_0_nux / cgs2cactusLength

    write(*,*) 'Energy absorption coefficients'
    write(*,*) '        nue     = ', abs_1_nue / cgs2cactusLength
    write(*,*) '        nua     = ', abs_1_nua / cgs2cactusLength
    write(*,*) '        nux     = ', abs_1_nux / cgs2cactusLength

    write(*,*) 'Neutrino number density'
    write(*,*) '        nue     = ', n_nue / (cgs2cactusLength**3)
    write(*,*) '        nua     = ', n_nue / (cgs2cactusLength**3)
    write(*,*) '        nux     = ', n_nue / (cgs2cactusLength**3)

    write(*,*) 'Neutrino energy density'
    edens_cgs2cactus = mev_to_erg * cgs2cactusEnergy / cgs2cactusLength**3
    write(*,*) '        nue     = ', e_nue * edens_cgs2cactus
    write(*,*) '        nua     = ', e_nua * edens_cgs2cactus
    write(*,*) '        nux     = ', e_nux * edens_cgs2cactus
    ! -----------------------------------------------------------------
end program
