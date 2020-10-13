!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Merra2_Driver1
!
! !DESCRIPTION: Program Merra2_Driver1 is a top-level driver for the 
!  MERRA2 regridding programs. 
!\\
!\\
! !INTERFACE:
!
PROGRAM Merra2_Driver1
!
! !USES:
!
  USE Merra2_A3CldModule
  USE Merra2_A3DynModule
  USE Merra2_InputsModule
  USE Merra2_RegridModule

  IMPLICIT NONE
!
! !REMARKS:
!  Merra2_Driver1 creates the A3cld (3hr time-averaged cloud parameters) and
!  A3dyn  (3hr time-averaged dynamics parameters) data files for input into 
!  GEOS-Chem.
!                                                                             .
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY: 
!  20 Sep 2013 - R. Yantosca - Now renamed Geos57 to Merra2_ in routine names
!------------------------------------------------------------------------------

  ! Read filenames and fields to process from an input file
  CALL Merra2_Initialize()

  ! Initialize GEOS-5 regridding code
  CALL Merra2_RegridInit()

  ! Create the 3-hour average data files
  CALL Merra2_MakeA3Cld()
  CALL Merra2_MakeA3Dyn()

  ! Cleanup and quit 
  CALL Merra2_Cleanup()

END PROGRAM Merra2_Driver1
!EOP
