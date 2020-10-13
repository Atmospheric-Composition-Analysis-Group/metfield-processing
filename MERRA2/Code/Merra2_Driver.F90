!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Merra2_Driver
!
! !DESCRIPTION: Program Merra2_Driver is the top-level driver for the 
!  GEOS-5.7.x regridding programs.  Merra2_Driver will call routines to 
!  extract, regrid, and save the MERRA2 met data to files for 
!  input to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
PROGRAM Merra2_Driver
!
! !USES:
!
  USE Merra2_A1Module
  USE Merra2_A3CldModule
  USE Merra2_A3DynModule
  USE Merra2_A3MstCModule
  USE Merra2_A3MstEModule
  USE Merra2_CnModule
  USE Merra2_I3Module
  USE Merra2_InputsModule
  USE Merra2_RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Merra2_Initialize()

  ! Initialize GEOS-5 regridding code
  CALL Merra2_RegridInit()

  ! Create the constant data file
  IF ( doMakeCn ) CALL Merra2_MakeCn()

  ! Create the 1-hour average data file
  CALL Merra2_MakeA1()

  ! Create the 3-hour average data files
  CALL Merra2_MakeA3Cld()
  CALL Merra2_MakeA3Dyn()
  CALL Merra2_MakeA3MstC()
  CALL Merra2_MakeA3MstE()

  ! Create the 6-hour instantaneous data file
  CALL Merra2_MakeI3()

  ! Cleanup and quit 
  CALL Merra2_Cleanup()

END PROGRAM Merra2_Driver
!EOP
