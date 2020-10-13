!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Merra2_Driver0
!
! !DESCRIPTION: Program Merra2_Driver0 is a top-level driver for the 
!  MERRA2 regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Merra2_Driver0
!
! !USES:
!
  USE Merra2_A1Module
  USE Merra2_CnModule
  USE Merra2_I3Module
  USE Merra2_InputsModule
  USE Merra2_RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Merra2_Driver1 creates the CN (constant), A1 (1hr time average), and
!  I3 (3hr instantaneous) data files for input into GEOS-Chem.
!                                                                             .
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

  ! Create the 3-hour instantaneous files
  CALL Merra2_MakeI3()

  ! Cleanup and quit 
  CALL Merra2_Cleanup()

END PROGRAM Merra2_Driver0
!EOP
