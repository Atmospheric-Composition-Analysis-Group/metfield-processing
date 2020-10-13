!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Merra2_Driver2
!
! !DESCRIPTION: Program Merra2_Driver2 is the top-level driver for the 
!  GEOS-5.7.x regridding programs.
!\\
!\\
! !INTERFACE:
!
PROGRAM Merra2_Driver2
!
! !USES:
!
  USE Merra2_A3MstCModule
  USE Merra2_A3MstEModule
  USE Merra2_InputsModule
  USE Merra2_RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Merra2_Driver1 creates the A3mstC (3hr time-averaged moist parameters, on
!  level centers) and A3MstE (3hr time-averaged moist parameters on level
!  edges) data files for input into GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  30 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Merra2_Initialize()

  ! Initialize GEOS-5 regridding code
  CALL Merra2_RegridInit()

  ! Create the 3-hour average data file
  CALL Merra2_MakeA3MstC()
  CALL Merra2_MakeA3MstE()

  ! Cleanup and quit 
  CALL Merra2_Cleanup()

END PROGRAM Merra2_Driver2
!EOP
