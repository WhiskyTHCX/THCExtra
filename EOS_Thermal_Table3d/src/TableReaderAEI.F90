#include <tableindex.h>
! #########################################################
! TABLE UNITS:
!        density              g/cm^3
!        temperature          MeV
!        ye                   number fraction per baryon
!        energy               erg/g
!        pressure             dyn/cm^2
!        cs2                  units of c^2 (relativistic)
! #########################################################

INTEGER function TableReader(myFilename)
  ! This routine reads the table and initializes
  ! all variables in the module.

  use hdf5
  use eos3dmod
  use hdfOpenField_mod

  implicit none
  character*200,  INTENT(IN) :: myFilename
  character*200 ::  myFieldName
  integer :: error

  write(*,*) "Reading 3D EOS Table", myFilename
  print*,'TabulatedEOS3D: ....loading the table with the EOS'

  !DENSITY
  myFieldName = "density"
  error = hdfOpenField(myFilename,myFieldName,logrho)
  !print*,'Here1'
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, rho"
    TableReader = -1
    return
  endif
  logrho = log10(logrho*cgs2cactusRho)


  !TEMPERATURE
  myFieldName = "temperature"
  error = hdfOpenField(myFilename,myFieldName,logtemp)
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, temp"
    TableReader = -1
    return
  endif
  logtemp = log10(logtemp)

  !Ye
  myFieldName = "ye"
  error = hdfOpenField(myFilename,myFieldName,yeTable)
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, ye"
    TableReader = -1
    return
  endif

  nrho = size(logrho,1)
  ntemp = size(logtemp,1)
  nye = size(yeTable,1)
  write(*,*) 'Number of points in the table:: ', nrho, ntemp, nye

  nvars = 5
  allocate(allvariables(nrho,ntemp,nye,nvars))

  !Pressure
  myFieldName = "pressure"
  error = hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,LPRESS_IDX))
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, PRESS"
    TableReader = -1
    return
  endif
  allvariables(:,:,:,LPRESS_IDX) = log10(allvariables(:,:,:,LPRESS_IDX)*cgs2cactusPress)

  !Internal Energy
  myFieldName = "internalEnergy"
  error = hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,LENERGY_IDX))
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, ENERGY"
    TableReader = -1
    return
  endif
  allvariables(:,:,:,LENERGY_IDX)=log10(allvariables(:,:,:,LENERGY_IDX)*cgs2cactusEps)

  !cs2
  myFieldName = "cs2"
  error = hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,CS2_IDX))
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, CS2"
    TableReader = -1
    return
  endif

  !depsdT_rho
  myFieldName = "depsdT_rho"
  error = hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,DEDT_IDX))
  if (error.ne.0) then
    write(*,*) "TableReader :: problem reading the table, DEDT"
    TableReader = -1
    return
  endif
  allvariables(:,:,:,DEDT_IDX)=log10(allvariables(:,:,:,DEDT_IDX) *cgs2cactusEps)

  !entropy per baryon
  myFieldName = "entropy"
  error = hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,ENTROPY_IDX))
  if (error.ne.0) then
     write(*,*) "TableReader :: problem reading the table, ENTROPY"
     TableReader = -1
     return
  endif
  allvariables(:,:,:,ENTROPY_IDX)=log10(allvariables(:,:,:,ENTROPY_IDX))

  ! set min-max values:
  eos_lrhomin = logrho(1)
  eos_lrhomax = logrho(nrho)
  eos_rhomin = 10.**logrho(1)
  eos_rhomax = 10.**logrho(nrho)
  write(*,*) "lrhomin",eos_lrhomin
  write(*,*) "lrhomax",eos_lrhomax

  eos_ltempmin = logtemp(1)
  eos_ltempmax = logtemp(ntemp)
  eos_tempmin = 10.**eos_ltempmin
  eos_tempmax = 10.**eos_ltempmax
  write(*,*) "ltempmin",eos_ltempmin
  write(*,*) "ltempmax",eos_ltempmax

  eos_yemin = yeTable(1)
  eos_yemax = yeTable(nye)
  write(*,*) "yemin",eos_yemin
  write(*,*) "yemax",eos_yemax

  ! Set the spacings
  dlrho = logrho(nrho)-logrho(nrho-1)
  dye = yeTable(nye) - yeTable(nye-1)
  dltemp = logtemp(ntemp) - logtemp(ntemp-1)
  write(*,*) "Spacing dlrho dye dltemp",dlrho,dye,dltemp

  write(6,*) "Done reading eos table"
  call h5open_f(error)
  TableReader = 0

end function TableReader

INTEGER (C_INT) function TableReaderC(myFilenameC) BIND(C, NAME="tablereaderc")
  USE ISO_C_BINDING
  implicit none
  character(kind=C_CHAR),  dimension(1), INTENT(IN) :: myFilenameC
  character*200 :: myFilename
  integer :: TableReader, i, n

  n = 200
  i = 1
  DO WHILE ((i .LE. n-1) .AND. (myFilenameC(i) .NE. c_null_char))
    myFilename(i:i)= myFilenameC(i)
    i = i + 1
  END DO
  ! myFilename(i:i)= c_null_char
  ! i = i + 1
  DO WHILE (i .LE. n)
    myFilename(i:i) = ' '
    i = i + 1
  END DO
  TableReaderC = TableReader(myFilename)
end function TableReaderC
