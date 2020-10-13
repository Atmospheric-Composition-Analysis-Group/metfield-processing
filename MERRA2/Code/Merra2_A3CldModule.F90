!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Merra2_A3CldModule
!
! !DESCRIPTION: Module Merra2_A3CldModule contains routines to create the 
!  GEOS-Chem average 3-hr data files (w/ cloud parameters) from the MERRA2
!  raw data.
!\\
!\\
! !INTERFACE: 

MODULE Merra2_A3CldModule
! 
! !USES:
!
  ! MERRA2 data modules
  USE CharpakModule
  USE Merra2_InputsModule
  USE Merra2_RegridModule
  USE Merra2_UtilityModule

  ! Modules for writing netCDF
  USE m_netcdf_io_create
  USE m_netcdf_io_define
  USE m_netcdf_io_write
  USE m_netcdf_io_close
  
  ! Modules for reading netCDF
  USE m_netcdf_io_open
  USE m_netcdf_io_close     
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read

  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Merra2_MakeA3Cld
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: Process3dCldNv
  PRIVATE :: Process3dOptDep
  PRIVATE :: RegridTau
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.
!
! !REVISION HISTORY:
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!  30 Jan 2016 - J.-W. Xu    - Rename all CH domain to AS (Asia) domain
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
!  LOGICAL :: use_CfAn       ! Save & regrid CFAN field
!  LOGICAL :: use_CfCu       ! Save & regrid CFCU field
!  LOGICAL :: use_CfLs       ! Save & regrid CFLS field
!
!  LOGICAL :: readOnly_CfAn  ! Only use CFAN to create CLOUD
!  LOGICAL :: readOnly_CfLs  ! Only use CFLS to create CLOUD

  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcOutFileDef
!
! !DESCRIPTION: Subroutine NcOutFileDef pre-defines variable names and 
!  attributes that will be added to the netCDF output files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NcOutFileDef( X,        Y,           Z,    T,    &
                           xMid,     yMid,        zMid, time, &
                           gridName, outFileName, fOut )
!
! !INPUT PARAMETERS:
! 
    INTEGER,          INTENT(IN)    :: X             ! Longitude dimension
    INTEGER,          INTENT(IN)    :: Y             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: Z             ! Latitude dimension
    INTEGER,          INTENT(IN)    :: T             ! Time dimension
    REAL*4,           INTENT(IN)    :: xMid(X)       ! Array of lon centers
    REAL*4,           INTENT(IN)    :: yMid(Y)       ! Array of lat centers
    REAL*4,           INTENT(IN)    :: zMid(Z)       ! Array of alt centers
    INTEGER,          INTENT(IN)    :: time(T)       ! Array of times
    CHARACTER(LEN=*), INTENT(IN)    :: gridName      ! Name of the grid
    CHARACTER(LEN=*), INTENT(IN)    :: outFileName   ! Output file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fOut          ! Output netCDF file ID
!
! !REVISION HISTORY:
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!  13 Aug 2015 - R. Yantosca - If the output file name ends in *.nc4 
!                              then save data to disk in netCDF-4 format
!  08 Sep 2015 - M. Sulprizio- Added global 0.5 x 0.625 grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255) :: sysTime
    CHARACTER(LEN=255) :: lName,   units,   gamap,   DI,    DJ
    CHARACTER(LEN=255) :: delta_t, begin_d, begin_t, incr,  msg, cal
    INTEGER            :: idLon,   idLat,   idLev,   idAp
    INTEGER            :: idBp,    idTime,  vId,     omode, C
    LOGICAL            :: is_nc4

    ! Arrays
    INTEGER            :: var1(1), var4(4)

    !=========================================================================
    ! %%% BEGINNING OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Echo info
    WRITE( 6, 100 ) TRIM( gridName )
100 FORMAT ( '%%% Defining netCDF file vars & attrs for ', a' grid' )

    ! If filename ends in ".nc"  then save as netCDF-3
    ! If filename ends in ".nc4" then save as netCDF-4
    C      = LEN_TRIM( outFileName )
    is_nc4 = ( outFileName(C-3:C) == '.nc4' )

    ! Open netCDF file for writing
    CALL NcCr_Wr( fOut, TRIM( outFileName ), WRITE_NC4=is_nc4 )

    ! Turn filling off
    CALL NcSetFill( fOut, NF_NOFILL, omode )

    !-------------------------------------------------------------------------
    ! Define global attributes and filling mode
    !-------------------------------------------------------------------------
  
    ! Title string
    lName = 'MERRA2 time-averaged 3-hour cloud parameters (A3cld), processed for GEOS-Chem input'
    CALL NcDef_Glob_Attributes( fOut, 'Title',                TRIM( lName ) )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',              TRIM( lName ) )

    ! References
    lName = "www.geos-chem.org; wiki.geos-chem.org"
    CALL NcDef_Glob_Attributes( fOut, 'References',           TRIM( lName ) )

    ! Filename
    lName = NotDir( outFileName )
    CALL NcDef_Glob_Attributes( fOut, 'Filename',             TRIM( lName ) )
    
    ! History
    sysTime = SystemTimeStamp()
    lName = 'File generated on: ' // TRIM( sysTime )
    CALL NcDef_Glob_Attributes( fOut, 'History' ,             TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ProductionDateTime',   TRIM( lName ) )
    CALL NcDef_Glob_Attributes( fOut, 'ModificationDateTime', TRIM( lName ) )

    ! Format
    IF ( is_Nc4 ) THEN
       lName = "NetCDF-4"
    ELSE
       lName = 'NetCDF-3'
    ENDIF
    CALL NcDef_Glob_Attributes( fOut, 'Format' ,              TRIM( lName ) )
                                                              
    ! Format                                                  
    lName = "global" ;                                        
    CALL NcDef_Glob_Attributes( fOut, 'SpatialCoverage',      TRIM( lName ) )
                                                              
    ! Conventions                                             
    lName = 'COARDS'                                          
    CALL NcDef_Glob_Attributes( fOut, 'Conventions',          TRIM( lName ) )
                                                              
    ! Version                                                 
    lName = 'MERRA2'                                        
    CALL NcDef_Glob_Attributes( fOut, 'Version',              TRIM( lName ) )
                                                              
    ! Model                                                   
    lName = TRIM( VersionId )                                         
    CALL NcDef_Glob_Attributes( fOut, 'VersionID',            TRIM( lName ) )
                                                              
    ! NLayers                                                 
    lName = '72'                                              
    CALL NcDef_Glob_Attributes( fOut, 'Nlayers',              TRIM( lName ) )
                                                              
    ! Start Date
    lName = yyyymmdd_string                                        
    CALL NcDef_Glob_Attributes( fOut, 'Start_Date',           TRIM( lName ) )
                                                              
    ! Start Time                                              
    lName = '00:00:00.0'                                      
    CALL NcDef_Glob_Attributes( fOut, 'Start_Time',           TRIM( lName ) )
                                                              
    ! End Date
    lName = yyyymmdd_string
    CALL NcDef_Glob_Attributes( fOut, 'End_Date',             TRIM( lName ) )
                                                              
    ! End Time                                                
    lName = '23:59:59.99999'
    CALL NcDef_Glob_Attributes( fOut, 'End_Time',             TRIM( lName ) )
                                                              
    ! Delta-time                                              
    lName = '030000'                                          
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Time',           TRIM( lName ) )

    ! Pick DI and DJ attributes based on the grid
    SELECT CASE ( TRIM( gridName ) )
       CASE( 'native', 'nested AS', 'nested EU', 'nested NA', 'nested SE' )
          DI = '0.3125'
          DJ = '0.25'
!       CASE ( 'nested 0.5 x 0.625' )
! (lzh,06/21/2014)
       CASE( 'nested AS 05', 'nested EU 05', 'nested NA 05', 'nested SE 05' )
          DI = '0.625'
          DJ = '0.5'
       CASE( '0.5 x 0.625 global' )
          DI = '0.625'
          DJ = '0.5'
       CASE( '2 x 2.5 global' )
          DI = '2.5'
          DJ = '2'
       CASE( '4 x 5 global' )
          DI = '5'
          DJ = '4'
    END SELECT

    ! Delta-lon
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lon',            TRIM( DI    ) )

    ! Delta-lat
    CALL NcDef_Glob_Attributes( fOut, 'Delta_Lat',            TRIM( DJ    ) )

    !-------------------------------------------------------------------------
    ! Define dimensions and index arrays.  NOTE: COARDS specifies that index 
    ! arrays will have the same names as the dimensions that define them.
    !-------------------------------------------------------------------------

    ! netCDF dimension variables
    CALL NcDef_Dimension( fOut, 'time', T,   idTime )
    CALL NcDef_Dimension( fOut, 'lev',  Z,   idLev  )
    CALL NcDef_Dimension( fOut, 'lat',  Y,   idLat  )
    CALL NcDef_Dimension( fOut, 'lon',  X,   idLon  )

    ! Time index array
    var1    = (/ idTime /)
    vId     = 0
    cal     = 'gregorian'
    lName   = 'time'
    units   = UnitsForTime( yyyymmdd )
    delta_t = '0000-00-00 03:00:00'
    begin_d = yyyymmdd_string
    begin_t = '000000'
    incr    = '030000'
    CALL NcDef_Variable      ( fOut, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'calendar',       TRIM( cal     )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName   )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fOut, vId, 'time_increment', TRIM( incr    )  )

    ! Level index array
    var1    = (/ idLev /)
    vId     = vId + 1
    lName   = 'levels'
    units   = '1'
    CALL NcDef_Variable      ( fOut, 'lev', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 

    ! Latitude index array
    var1    = (/ idLat /)
    vId     = vId + 1
    lName   = 'latitude'
    units   = 'degrees_north'
    CALL NcDef_Variable      ( fOut, 'lat', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName )    )
    CALL NcDef_Var_attributes( fOut, vId, 'units',          TRIM( units )    ) 
  
    ! Longitude index array
    vId     = vId + 1
    var1    = (/ idLon /)
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! CLOUD
    IF ( StrPos( 'CLOUD', tavg3_3d_cld_Nv_Data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'cloud_fraction_for_radiation'
       units = '1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'CLOUD', NF_FLOAT, 4, var4, vId      )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! OPTDEPTH
    IF ( StrPos( 'OPTDEPTH', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'in_cloud_optical_thickness'
       units = '1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'OPTDEPTH', NF_FLOAT, 4, var4, vId   )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! QI
    IF ( StrPos( 'QI', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'mass_fraction_of_cloud_ice_water'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QI', NF_FLOAT, 4, var4, vId         )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! QL
    IF ( StrPos( 'QL', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'mass_fraction_of_cloud_liquid_water'
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QL', NF_FLOAT, 4, var4, vId         )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! RH
!    IF ( StrPos( 'RH', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
!       var4  = (/ idLon, idLat, idLev, idTime /)    
!       vId   = vId + 1
!       lName = 'relative_humidity_after_moist'
!       units = '1'
!       gamap = 'GMAO-3D$'
!       CALL NcDef_Variable      ( fOut, 'RH', NF_FLOAT, 4, var4, vId     )
!       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
!       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
!       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
!       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
!       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
!       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
!       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
!       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
!    ENDIF

    ! TAUCLI
    IF ( StrPos( 'TAUCLI', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'in_cloud_optical_thickness_for_ice_clouds'
       units = '1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'TAUCLI', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! TAUCLW
    IF ( StrPos( 'TAUCLW', tavg3_3d_cld_Nv_data_c ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'in_cloud_optical_thickness_for_liquid_clouds'
       units = '1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'TAUCLW', NF_FLOAT, 4, var4, vId     )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    !=========================================================================
    ! %%% END OF NETCDF DEFINITION SECTION %%%
    !=========================================================================

    ! End the definition section
    CALL NcEnd_def( fOut )

    ! Write index arrays
    CALL NcWr( xMid, fOut, 'lon',  (/ 1 /), (/ X   /) )
    CALL NcWr( yMid, fOut, 'lat',  (/ 1 /), (/ Y   /) )
    CALL NcWr( zMid, fOut, 'lev',  (/ 1 /), (/ Z   /) )
    CALL NcWr( time, fOut, 'time', (/ 1 /), (/ T   /) )

    ! Echo info    
    msg = '%%%%%% LEAVING ROUTINE NcOutFileDef %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE NcOutFileDef
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_MakeA3Cld
!
! !DESCRIPTION: Routine Merra2_MakeA3Cld is the the driver routine for 
! \begin{enumerate}
! \item Extracting 3-hr time-averaged data fields (cloud parameters) from 
!       the MERRA2 raw data files (netCDF-4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program Merra2_Driver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_MakeA3Cld()
!
! !REVISION HISTORY:
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!  08 Sep 2015 - M. Sulprizio- Added global 0.5 x 0.625 grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: nFields_3dCldNv
    INTEGER                 :: nFields_3dRadNv
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: allFieldsList
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: allFields     (MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dCldNv(MAX_FLDS)
    CHARACTER(LEN=MAX_CHAR) :: fields_3dRadNv(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info    
    msg = '%%%%%%%%%% ENTERING ROUTINE Merra2_MakeA3Cld %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! List of all the A-3 fields combined
    allFieldsList = TRIM( tavg3_3d_cld_Nv_data_c ) // ',' // &
                    TRIM( tavg3_3d_rad_Nv_data   )

    ! Return the list of fields and number of fields to process
    ! from each of the Merra2_ raw met data files
    CALL GetNFields( tavg3_3d_cld_Nv_data_c, nFields_3dCldNv, fields_3dCldNv )
    CALL GetNFields( tavg3_3d_rad_Nv_data,   nFields_3dRadNv, fields_3dRadNv )
    CALL GetNFields( allFieldsList,          nAllFields,      allFields      )
    
    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( tavg3_3d_cld_Nv_file ), nFields_3dCldNv
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open nested 0625 AS output file
    IF ( doNestAs05 ) THEN
       fName = TRIM( tempDirTmplNestAs05 ) // TRIM( dataTmplNestAs05 )
       gName = 'nested As 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3cld '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestAs05,  J_NestAs05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_as05:I1_as05),          &
                          yMid_05x0625(J0_as05:J1_as05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,        fOut05NestAs          )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3cld '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,       fOut05NestEu          )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3cld '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05, L05x0625,  TIMES_A3,  &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestNa          )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe05 ) THEN
       fName = TRIM( tempDirTmplNestSe05 ) // TRIM( dataTmplNestSe05 )
       gName = 'nested SE 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'A3cld '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestSe05,  J_NestSe05, L05x0625,  TIMES_A3,  &
                          xMid_05x0625(I0_se05:I1_se05),          &
                          yMid_05x0625(J0_se05:J1_se05),          &
                          zMid_05x0625,                a3Mins,    &
                          gName,    fName,      fOut05NestSe         )
    ENDIF

    ! Open 0.5 x 0.625 output file
    IF ( doGlobal05 ) THEN
       fName = TRIM( tempDirTmpl05x0625 ) // TRIM( dataTmpl05x0625 )
       gName = '0.5 x 0.625 global'
       CALL ExpandDate  ( fName,        yyyymmdd,        000000      )      
       CALL StrRepl     ( fName,        '%%%%%%',        'A3cld '    )
       CALL StrCompress ( fName,        RemoveAll=.TRUE.             )
       CALL NcOutFileDef( I05x0625,     J05x0625,        &
                          L05x0625,     TIMES_A3,        &
                          xMid_05x0625, nc_yMid_05x0625, &
                          zMid_05x0625, a3Mins,          &
                          gName,        fName,           fOut05x0625 )
    ENDIF

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fName = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                ) 
       CALL StrRepl     ( fName,     '%%%%%%',    'A3cld '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                    )
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25,      TIMES_A3,  &
                          xMid_2x25, nc_yMid_2x25, zMid_2x25,  a3Mins,    &
                          gName,     fName,        fOut2x25              )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fName = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                )      
       CALL StrRepl     ( fName,     '%%%%%%',    'A3cld '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                    )
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5,       TIMES_A3,  &
                          xMid_4x5,  nc_yMid_4x5,  zMid_4x5,   a3Mins,    &
                          gName,     fName,        fOut4x5               )
    ENDIF

    !=======================================================================
    ! Process data
    !=======================================================================
    CALL Process3dCldNv ( nFields_3dCldNv, fields_3dCldNv ) ! tavg3_3d_cld_Nv
    CALL Process3dOptDep(                                 ) ! optical depths
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing A3cld output files'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Close output files
    IF ( do2x25     ) CALL NcCl( fOut2x25     )
    IF ( do4x5      ) CALL NcCl( fOut4x5      )
    IF ( doGlobal05 ) CALL NcCl( fOut05x0625  )
    IF ( doNestAs05 ) CALL NcCl( fOut05NestAs )
    IF ( doNestEu05 ) CALL NcCl( fOut05NestEu )
    IF ( doNestNa05 ) CALL NcCl( fOut05NestNa )
    IF ( doNestSe05 ) CALL NcCl( fOut05NestSe )
    

    ! Echo info
    msg = '%%%%%%%%%% LEAVING ROUTINE Merra2_MakeA3Cld %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'   
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Merra2_MakeA3Cld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dCldNv
!
! !DESCRIPTION: Subroutine Process3dCldNv regrids the MERRA2 met fields 
!  from the "tavg3\_3d\_cld\_Nv" file and saves output to netCDF file format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dCldNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
!
! !REMARKS:
!  The cloud fraction field CLOUD and cloud optical depth fields TAUCLI, 
!  TAUCLW, and OPTDEPTH are processed separately in routine Process3dOptDep.  
!  This is because these fields must all be regridded together using the
!  algorithm developed by Hongyu Liu (in routine RegridTau).
!
! !REVISION HISTORY:
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!  08 Sep 2015 - M. Sulprizio- Added global 0.5 x 0.625 grid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: F,         H,         L,       LR
    INTEGER                 :: hhmmss 

    ! Variables for netCDF I/O
    INTEGER                 :: X,         Y,         Z,         T
    INTEGER                 :: XNestAs05, YNestAs05, ZNestAs05, TNestAs05
    INTEGER                 :: XNestEu05, YNestEu05, ZNestEu05, TNestEu05
    INTEGER                 :: XNestNa05, YNestNa05, ZNestNa05, TNestNa05
    INTEGER                 :: XNestSe05, YNestSe05, ZNestSe05, TNestSe05
    INTEGER                 :: X05x0625,  Y05x0625,  Z05x0625,  T05x0625
    INTEGER                 :: X2x25,     Y2x25,     Z2x25,     T2x25
    INTEGER                 :: X4x5,      Y4x5,      Z4x5,      T4x5
    INTEGER                 :: st4d(4),   ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Q      ( I05x0625, J05x0625, L05x0625 )
    REAL*4,  TARGET         :: QI     ( I05x0625, J05x0625, L05x0625 )
    REAL*4,  TARGET         :: QL     ( I05x0625, J05x0625, L05x0625 )
    REAL*4                  :: Q2x25  ( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: QI_2x25( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: QL_2x25( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: Q4x5   ( I4x5,     J4x5,     L4x5     )
    REAL*4                  :: QI_4x5 ( I4x5,     J4x5,     L4x5     )
    REAL*4                  :: QL_4x5 ( I4x5,     J4x5,     L4x5     )
    
    ! Pointer arrays
    REAL*4,  POINTER        :: Qflip(:,:,:)
    REAL*4,  POINTER        :: Ptr  (:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

     ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! 2 x 2.5 global grid       
    IF ( do2x25 ) THEN
       CALL NcGet_DimLen( fOut2x25,     'lon',  X2x25     )
       CALL NcGet_DimLen( fOut2x25,     'lat',  Y2x25     ) 
       CALL NcGet_DimLen( fOut2x25,     'lev',  Z2x25     ) 
       CALL NcGet_DimLen( fOut2x25,     'time', T2x25     )
    ENDIF

    ! 4x5 global grid
    IF ( do4x5 ) THEN
       CALL NcGet_DimLen( fOut4x5,      'lon',  X4x5      )
       CALL NcGet_DimLen( fOut4x5,      'lat',  Y4x5      )   
       CALL NcGet_DimLen( fOut4x5,      'lev',  Z4x5      )   
       CALL NcGet_DimLen( fOut4x5,      'time', T4x5      )
    ENDIF
    
    ! 0.5 x 0.625 global grid
    IF ( doGlobal05 ) THEN
       CALL NcGet_DimLen( fOut05x0625,  'lon',  X05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'lat',  Y05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'lev',  Z05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'time', T05x0625  )
    ENDIF

    ! Nested AS grid
    IF ( doNestAs05 ) THEN
       CALL NcGet_DimLen( fOut05NestAs, 'lon',  XNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lat',  YNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lev',  ZNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'time', TNestAs05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestEu05 ) THEN
       CALL NcGet_DimLen( fOut05NestEu, 'lon',  XNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lat',  YNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lev',  ZNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'time', TNestEu05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestNa05 ) THEN
       CALL NcGet_DimLen( fOut05NestNa, 'lon',  XNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lat',  YNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lev',  ZNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'time', TNestNa05 )
    ENDIF

    ! Nested SE grid
    IF ( doNestSe05 ) THEN
       CALL NcGet_DimLen( fOut05NestSe, 'lon',  XNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lat',  YNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lev',  ZNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'time', TNestSe05 )
    ENDIF
    
    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Create input filename from the template
    fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_cld_Nv_file )
    CALL expandDate( fNameInput, yyyymmdd, hhmmss )

    ! Echo info
    msg = '%%% Opening ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Open the netCDF4 file for input
    CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
    ! Get the dimensions from the netCDF file
    CALL NcGet_DimLen( fIn, 'lon',  X )
    CALL NcGet_DimLen( fIn, 'lat',  Y ) 
    CALL NcGet_DimLen( fIn, 'lev',  Z ) 
    CALL NcGet_DimLen( fIn, 'time', T )

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! GMT time of day (hh:mm:ss)
       hhmmss = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

       ! Zero QI, QL arrays
       QI      = 0e0
       QI_2x25 = 0e0
       QI_4x5  = 0e0
       QL      = 0e0
       QL_2x25 = 0e0
       QL_4x5  = 0e0
       
       !====================================================================
       ! Process data
       !====================================================================
       
       ! Loop over data fields
       DO F = 1, nFields
          
          ! Save field name into an 8-char variable. 
          ! This will truncate field names longer than 8 chars.
          name = TRIM( fields(F) )
          
          ! Skip certain fieldnames
          SELECT CASE ( name ) 
             CASE( '' )                                ! Null string
                CYCLE
             CASE( 'TAUCLI', 'TAUCLW', 'OPTDEPTH',  &  
                   'CFAN',   'CFCU',   'CFLS',      &
                   'CLOUD'                         )   ! These fields are
                CYCLE                                  !  procesed elsewhere 
             CASE DEFAULT
                ! Nothing
          END SELECT

          ! Zero data arrays
          Q     = 0e0
          Q2x25 = 0e0
          Q4x5  = 0e0

          !-----------------------------------------------------------------
          ! Read data
          !-----------------------------------------------------------------
          msg = '%%% Reading    ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Start and count index arrays for netCDF
          st4d = (/ 1, 1, 1, H /)
          ct4d = (/ X, Y, Z, 1 /)
          
          ! Read data
          CALL NcRd( Q, fIn, TRIM( name ), st4d, ct4d )

          ! Strip fill values
          WHERE( Q == FILL_VALUE ) Q = 0e0    

          ! Flip data in vertical
          Qflip => Q( :, :, Z:1:-1 )

          !-----------------------------------------------------------
          ! Regrid data fields
          !-----------------------------------------------------------
          msg = '%%% Regridding ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )

          ! Loop over the A-3 times and vertical levels
          DO L = 1, Z

             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridMerra2_To2x25( 0, Qflip(:,:,L), Q2x25(:,:,L) )
             ENDIF

             ! Regrid to 4x5 
             IF ( do4x5 ) THEN
                CALL RegridMerra2_To4x5 ( 0, Qflip(:,:,L), Q4x5(:,:,L)  )
             ENDIF
          ENDDO
                
          !-----------------------------------------------------------
          ! Write netCDF output (all except QI, QL)
          !-----------------------------------------------------------
          msg = '%%% Archiving  ' // name
          WRITE( IU_LOG, '(a)' ) TRIM( msg )
          
          ! Write 2 x 2.5 data
          IF ( do2x25 ) THEN
             st4d = (/ 1,     1,     1,     H  /)
             ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)
             CALL NcWr( Q2x25, fOut2x25, TRIM( name ), st4d, ct4d )
          ENDIF
                
          ! Write 4x5 data
          IF ( do4x5 ) THEN
             st4d  = (/ 1,    1,    1,    H /)
             ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
             CALL NcWr( Q4x5, fOut4x5, TRIM( name ), st4d, ct4d )
          ENDIF

          ! Write 0.5 x 0.625 data
          IF ( doGlobal05 ) THEN
             Ptr  => Qflip(:,:,:)
             st4d = (/ 1,         1,         1,         H /)
             ct4d = (/ X05x0625,  Y05x0625,  Z05x0625 , 1 /)
             CALL NcWr( Ptr, fOut05x0625, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested Asia (point to proper slice of global data)
          IF ( doNestAs05 ) THEN
             Ptr  => Qflip( I0_as05:I1_as05, J0_as05:J1_as05, : )
             st4d = (/ 1,         1,         1,         H /)
             ct4d = (/ XNestAs05, YNestAs05, ZNestAs05, 1 /)
             CALL NcWr( Ptr, fOut05NestAs, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested EU (point to proper slice of global data)
          IF ( doNestEu05 ) THEN
             Ptr  => Qflip( I0_eu05:I1_eu05, J0_eu05:J1_eu05, : )
             st4d = (/ 1,         1,         1,         H /)
             ct4d = (/ XNestEu05, YNestEu05, ZNestEu05, 1 /)
             CALL NcWr( Ptr, fOut05NestEu, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested NA (point to proper slice of global data)
          IF ( doNestNa05 ) THEN
             Ptr  => Qflip( I0_na05:I1_na05, J0_na05:J1_na05, : )
             st4d = (/ 1,         1,         1,         H /)
             ct4d = (/ XNestNa05, YNestNa05, ZNestNa05, 1 /)
             CALL NcWr( Ptr, fOut05NestNa, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Nested SE (point to proper slice of global data)
          IF ( doNestSe05 ) THEN
             Ptr  => Qflip( I0_se05:I1_se05, J0_se05:J1_se05, : )
             st4d = (/ 1,         1,         1,         H /)
             ct4d = (/ XNestSe05, YNestSe05, ZNestSe05, 1 /)
             CALL NcWr( Ptr, fOut05NestSe, TRIM( name ), st4d, ct4d )
             NULLIFY( Ptr )
          ENDIF

          ! Free pointer memory
          NULLIFY( Qflip )
       ENDDO
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================

    ! Close input file
    msg = '%%% Closing ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL NcCl( fIn )       

    ! Echo info    
    msg = '%%%%%% LEAVING ROUTINE Process3dCldNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dCldNv
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process3dOptDep
!
! !DESCRIPTION: Subroutine Process3dOptDep regrids the CLOUD, TAUCLI,
!  TAUCLW, and OPTDEPTH fields using Hongyu Liu's regridding algorithm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Process3dOptDep()
!
! !REMARKS:
!  The OPTDEPTH and CLOUD fields are regridded following the algorithm of
!  Hongyu Liu (cf. routine RegridTau).  OPTDEPTH = TAUCLI + TAUCLW.
!
! !REVISION HISTORY:
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!  08 Sep 2015 - M. Sulprizio- Added global 0.5 x 0.625 grid
!  17 Mar 2016 - M. Sulprizio- Fix bug in std4d to account for multiple data
!                              blocks per file
!  18 Mar 2016 - M. Sulprizio- Add fixes for flipping 0.5x0.625 global and
!                              nested fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop and time variables
    INTEGER                 :: F,         H,         L,         LR
    INTEGER                 :: hhmmss 

    ! Variables for netCDF I/O
    INTEGER                 :: X,         Y,         Z,         T
    INTEGER                 :: XNestAs05, YNestAs05, ZNestAs05, TNestAs05
    INTEGER                 :: XNestEu05, YNestEu05, ZNestEu05, TNestEu05
    INTEGER                 :: XNestNa05, YNestNa05, ZNestNa05, TNestNa05
    INTEGER                 :: XNestSe05, YNestSe05, ZNestSe05, TNestSe05
    INTEGER                 :: X05x0625,  Y05x0625,  Z05x0625,  T05x0625
    INTEGER                 :: X2x25,     Y2x25,     Z2x25,     T2x25
    INTEGER                 :: X4x5,      Y4x5,      Z4x5,      T4x5
    INTEGER                 :: st3d(3),   st4d(4)
    INTEGER                 :: ct3d(3),   ct4d(4)

    ! Data arrays
    REAL*4,  TARGET         :: Cld      ( I05x0625, J05x0625, L05x0625 )
    REAL*4,  TARGET         :: OptD     ( I05x0625, J05x0625, L05x0625 )
    REAL*4,  TARGET         :: TauI     ( I05x0625, J05x0625, L05x0625 )
    REAL*4,  TARGET         :: TauW     ( I05x0625, J05x0625, L05x0625 )
    REAL*4                  :: Cld_2x25 ( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: OptD_2x25( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: TauI_2x25( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: TauW_2x25( I2x25,    J2x25,    L2x25    )
    REAL*4                  :: Cld_4x5  ( I4x5,     J4x5,     L4x5     )
    REAL*4                  :: OptD_4x5 ( I4x5,     J4x5,     L4x5     )
    REAL*4                  :: TauI_4x5 ( I4x5,     J4x5,     L4x5     )
    REAL*4                  :: TauW_4x5 ( I4x5,     J4x5,     L4x5     )
    
    ! Pointer arrays
    REAL*4,  POINTER        :: ptr(:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE Process3dOptDep %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! 2 x 2.5 global grid       
    IF ( do2x25 ) THEN
       CALL NcGet_DimLen( fOut2x25,     'lon',  X2x25     )
       CALL NcGet_DimLen( fOut2x25,     'lat',  Y2x25     ) 
       CALL NcGet_DimLen( fOut2x25,     'lev',  Z2x25     ) 
       CALL NcGet_DimLen( fOut2x25,     'time', T2x25     )
    ENDIF

    ! 4x5 global grid
    IF ( do4x5 ) THEN
       CALL NcGet_DimLen( fOut4x5,      'lon',  X4x5      )
       CALL NcGet_DimLen( fOut4x5,      'lat',  Y4x5      )   
       CALL NcGet_DimLen( fOut4x5,      'lev',  Z4x5      )   
       CALL NcGet_DimLen( fOut4x5,      'time', T4x5      )  
    ENDIF

    ! 0.5 x 0.625 global grid
    IF ( doGlobal05 ) THEN
       CALL NcGet_DimLen( fOut05x0625,  'lon',  X05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'lat',  Y05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'lev',  Z05x0625  )
       CALL NcGet_DimLen( fOut05x0625,  'time', T05x0625  )
    ENDIF

    ! Nested AS grid
    IF ( doNestAs05 ) THEN
       CALL NcGet_DimLen( fOut05NestAs, 'lon',  XNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lat',  YNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'lev',  ZNestAs05 )
       CALL NcGet_DimLen( fOut05NestAs, 'time', TNestAs05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestEu05 ) THEN
       CALL NcGet_DimLen( fOut05NestEu, 'lon',  XNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lat',  YNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'lev',  ZNestEu05 )
       CALL NcGet_DimLen( fOut05NestEu, 'time', TNestEu05 )
    ENDIF

    ! Nested NA grid
    IF ( doNestNa05 ) THEN
       CALL NcGet_DimLen( fOut05NestNa, 'lon',  XNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lat',  YNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'lev',  ZNestNa05 )
       CALL NcGet_DimLen( fOut05NestNa, 'time', TNestNa05 )
    ENDIF

    ! Nested SE grid
    IF ( doNestSe05 ) THEN
       CALL NcGet_DimLen( fOut05NestSe, 'lon',  XNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lat',  YNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'lev',  ZNestSe05 )
       CALL NcGet_DimLen( fOut05NestSe, 'time', TNestSe05 )
    ENDIF
    

    !=======================================================================
    ! Open input file
    !=======================================================================

    ! Create input filename from the template
    fNameInput = TRIM( inputDataDir ) // TRIM( tavg3_3d_cld_Nv_file )
    CALL expandDate( fNameInput, yyyymmdd, hhmmss )

    ! Echo info
    msg = '%%% Opening ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Open the netCDF4 file for input
    CALL NcOp_Rd( fIn, TRIM( fNameInput ) )
       
    ! Get the dimensions from the netCDF file
    CALL NcGet_DimLen( fIn, 'lon',  X )
    CALL NcGet_DimLen( fIn, 'lat',  Y ) 
    CALL NcGet_DimLen( fIn, 'lev',  Z ) 
    CALL NcGet_DimLen( fIn, 'time', T )

    ! Zero arrays
    Cld   = 0e0
    OptD  = 0e0
    TauI  = 0e0
    Tauw  = 0e0

    ! Loop over the number of files per day
    DO H = 1, TIMES_A3

       ! Zero cloud arrays for each hour
       Cld       = 0e0
       Cld_2x25  = 0e0
       Cld_4x5   = 0e0
       TauI      = 0e0
       TauI_2x25 = 0e0
       TauI_4x5  = 0e0
       TauW      = 0e0
       TauW_2x25 = 0e0
       TauW_4x5  = 0e0
       OptD      = 0e0
       OptD_2x25 = 0e0
       OptD_4x5  = 0e0

       ! GMT time of day (hh:mm:ss)
       hhmmss    = ( ( a3mins(H) / 60 ) * 10000 ) + 3000

       !==================================================================
       ! Read TAUCLI, TAUCLW, CLOUD from tavg3_3d_cld_Nv
       !==================================================================

       ! Start and count index arrays for netCDF
       st4d = (/ 1, 1, 1, H /)
       ct4d = (/ X, Y, Z, 1 /)

       ! CLOUD
       msg = '%%% Reading    CLOUD' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( Cld, fIn, 'CLOUD', st4d, ct4d )
       WHERE( Cld == FILL_VALUE ) Cld = 0e0    
       
       ! TAUCLI
       msg = '%%% Reading    TAUCLI' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( TauI, fIn, 'TAUCLI', st4d, ct4d )
       WHERE( TauI == FILL_VALUE ) TauI = 0e0    

       ! TAUCLW
       msg = '%%% Reading    TAUCLW' 
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       CALL NcRd( TauW, fIn, 'TAUCLW', st4d, ct4d )
       WHERE( TauW == FILL_VALUE ) TauI = 0e0    

       ! OPTDEPTH (construct from TAUCLI and TAUCLW)
       OptD = TauI + TauW

       !====================================================================
       ! Regrid cloud & optical depth fields w/ Hongyu Liu's algorithm
       ! NOTE: Algorithm also flips data in the vertical
       !====================================================================

       msg = '%%% Regridding cloud fraction & optical depth fields'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
    
       ! Regrid to 2 x 2.5 (reverse levels of output array Q2x25)
       IF ( do2x25 ) THEN
          CALL RegridTau( OptD,      OptD_2x25, TauI,  TauI_2x25,  &
                          TauW,      TauW_2x25, Cld,   Cld_2x25,   &
                          mapTo2x25, I2x25,     J2x25, L2x25      )
       ENDIF
       
       ! Regrid to 4x5 (reverse levels of output array Q4x5)
       IF ( do4x5 ) THEN
          CALL RegridTau( OptD,      OptD_4x5,  TauI,  TauI_4x5,   &
                          TauW,      TauW_4x5,  Cld,   Cld_4x5,    &
                          mapTo4x5,  I4x5,      J4x5,  L4x5       )
       ENDIF
              
       !====================================================================
       ! Write to netCDF
       !====================================================================
 
       msg = '%%% Archiving  CLOUD, TAUCLI, TAUCLW, OPTDEPTH'
       WRITE( IU_LOG, '(a)' ) TRIM( msg )
       
       !----------------------------------------
       ! 2 x 2.5 GLOBAL GRID
       !----------------------------------------
       IF ( do2x25 ) THEN

          ! netCDF indices
          st4d = (/ 1,     1,     1,     H  /)
          ct4d = (/ X2x25, Y2x25, Z2x25, 1  /)

          ! Write to disk
          CALL NcWr( Cld_2x25,  fOut2x25, 'CLOUD',    st4d, ct4d )
          CALL NcWr( TauI_2x25, fOut2x25, 'TAUCLI',   st4d, ct4d )
          CALL NcWr( TauW_2x25, fOut2x25, 'TAUCLW',   st4d, ct4d )
          CALL NcWr( OptD_2x25, fOut2x25, 'OPTDEPTH', st4d, ct4d )
       ENDIF

       !----------------------------------------
       ! 4 x 5 GLOBAL GRID
       !----------------------------------------
       IF ( do4x5 ) THEN

          ! netCDF indices
          st4d = (/ 1,    1,    1,    H  /)
          ct4d = (/ X4x5, Y4x5, Z4x5, 1  /)

          ! Write to disk
          CALL NcWr( Cld_4x5,  fOut4x5, 'CLOUD',    st4d, ct4d )
          CALL NcWr( TauI_4x5, fOut4x5, 'TAUCLI',   st4d, ct4d )
          CALL NcWr( TauW_4x5, fOut4x5, 'TAUCLW',   st4d, ct4d )
          CALL NcWr( OptD_4x5, fOut4x5, 'OPTDEPTH', st4d, ct4d )
       ENDIF

       !----------------------------------------
       ! 0.5 x 0.625 GLOBAL GRID (flip in vertical)
       !----------------------------------------
       IF ( doGlobal05 ) THEN

          ! netCDF indices
          st4d = (/ 1,         1,         1,         H /)
          ct4d = (/ X05x0625,  Y05x0625,  Z05x0625,  1 /)

          ! CLOUD
          Ptr  => Cld(:,:,Z:1:-1 )
          CALL NcWr( Ptr, fOut05x0625, 'CLOUD',    st4d, ct4d )

          ! TAUCLI 
          Ptr  => TauI(:,:,Z:1:-1 )
          CALL NcWr( Ptr, fOut05x0625, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW 
          Ptr  => TauW(:,:,Z:1:-1 )
          CALL NcWr( Ptr, fOut05x0625, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD(:,:,Z:1:-1 )
          CALL NcWr( Ptr, fOut05x0625, 'OPTDEPTH', st4d, ct4d )

          ! Free pointer memory
          NULLIFY( Ptr )
       ENDIF

       !----------------------------------------
       ! NESTED AS GRID  (flip in vertical)
       !----------------------------------------
       IF ( doNestAs05 ) THEN

          ! netCDF indices
          st4d = (/ 1,         1,         1,         H /)
          ct4d = (/ XNestAs05, YNestAs05, ZNestAs05, 1 /)

          ! CLOUD
          Ptr  => Cld( I0_as05:I1_as05, J0_as05:J1_as05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestAs, 'CLOUD',    st4d, ct4d )

          ! TAUCLI 
          Ptr  => TauI( I0_as05:I1_as05, J0_as05:J1_as05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestAs, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW 
          Ptr  => TauW( I0_as05:I1_as05, J0_as05:J1_as05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestAs, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD( I0_as05:I1_as05, J0_as05:J1_as05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestAs, 'OPTDEPTH', st4d, ct4d )

          ! Free pointer memory
          NULLIFY( Ptr )
       ENDIF

       !----------------------------------------
       ! NESTED EU GRID (flip in vertical)
       !----------------------------------------
       IF ( doNestEu05 ) THEN

          ! netCDF indices
          st4d = (/ 1,         1,         1,         H /)
          ct4d = (/ XNestEu05, YNestEu05, ZNestEu05, 1 /)

          ! CLOUD
          Ptr  => Cld( I0_eu05:I1_eu05, J0_eu05:J1_eu05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestEu, 'CLOUD',    st4d, ct4d )

          ! TAUCLI
          Ptr  => TauI( I0_eu05:I1_eu05, J0_eu05:J1_eu05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestEu, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW
          Ptr  => TauW( I0_eu05:I1_eu05, J0_eu05:J1_eu05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestEu, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD( I0_eu05:I1_eu05, J0_eu05:J1_eu05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestEu, 'OPTDEPTH', st4d, ct4d )

          ! Free pointer memory
          NULLIFY( Ptr )
       ENDIF

       !----------------------------------------
       ! NESTED NA GRID (flip in vertical)
       !----------------------------------------
       IF ( doNestNa05 ) THEN

          ! netCDF indices
          st4d = (/ 1,         1,         1,         H /)
          ct4d = (/ XNestNa05, YNestNa05, ZNestNa05, 1 /)

          ! CLOUD
          Ptr  => Cld( I0_na05:I1_na05, J0_na05:J1_na05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestNa, 'CLOUD',    st4d, ct4d )

          ! TAUCLI
          Ptr  => TauI( I0_na05:I1_na05, J0_na05:J1_na05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestNa, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW
          Ptr  => TauW( I0_na05:I1_na05, J0_na05:J1_na05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestNa, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD( I0_na05:I1_na05, J0_na05:J1_na05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestNa, 'OPTDEPTH', st4d, ct4d )

          ! Free pointer memory
          NULLIFY( Ptr )
       ENDIF

       !----------------------------------------
       ! NESTED SE GRID (flip in vertical)
       !----------------------------------------
       IF ( doNestSe05 ) THEN

          ! netCDF indices
          st4d = (/ 1,         1,         1,         H /)
          ct4d = (/ XNestSe05, YNestSe05, ZNestSe05, 1 /)

          ! CLOUD
          Ptr  => Cld( I0_se05:I1_se05, J0_se05:J1_se05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestSe, 'CLOUD',    st4d, ct4d )

          ! TAUCLI
          Ptr  => TauI( I0_se05:I1_se05, J0_se05:J1_se05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestSe, 'TAUCLI',   st4d, ct4d )

          ! TAUCLW
          Ptr  => TauW( I0_se05:I1_se05, J0_se05:J1_se05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestSe, 'TAUCLW',   st4d, ct4d )

          ! OPTDEPTH
          Ptr  => OptD( I0_se05:I1_se05, J0_se05:J1_se05, Z:1:-1 )
          CALL NcWr( Ptr, fOut05NestSe, 'OPTDEPTH', st4d, ct4d )

          ! Free pointer memory
          NULLIFY( Ptr )
       ENDIF
    ENDDO

    !=======================================================================
    ! Quit
    !=======================================================================

    ! Close input netCDF file
    msg = '%%% Closing ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL NcCl( fIn )

    msg = '%%%%%% LEAVING ROUTINE Process3dOptDep %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Process3dOptDep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegridTau
!
! !DESCRIPTION: This routine regrids the GEOS-5 optical depth and cloud 
!  fraction fields from the 0.5 x 0.666 native resolution grid to a GEOS-Chem 
!  "coarse" grid (e.g. 1 x 1.25, 2 x 2.5, 4 x 5).
!
! !REMARKS:
!  Regridding formulae from Hongyu Liu:
!
!  Eq. 1: Optical depth:
!  ---------------------
!                                                                             .
!           TAU                   {     Fn * Wn            TAUn      }
!      -------------  =  A  =  SUM{ --------------- * -------------  }
!      ( TAU + 7.7 )              {  SUM( Fn * Wn )   ( TAUn + 7.7 ) }
!
!                                   |_____________|   |____________|
!     in the code these are:          Fn_Wn_Norm         optRatio
!                                                        tauiRatio
!                                                        tauwRAtio
!                                                                             .
!  and to solve for TAU from this equation, we manipulate the right hand 
!  side of the equation (called A), as follows:
!                                                                             .
!                                A
!                 TAU = 7.7 * --------
!                              1 - A
!                                                                             .
!  Eq. 2. Cloud fraction:
!  ----------------------
!                                                                             .
!                        SUM( Fn * Wn )
!                  F  =  ----------------   
!                           SUM( Wn )
!                                                                             .
!  In both above equations:
!                                                                             .
!     TAU    = In-cloud optical depth regridded to "coarse" grid box
!     TAUn   = In-cloud optical depth in a 0.5 x 0.666 grid box
 !                                                                             .
!     F      = Cloud fraction regriddded to "coarse" grid box
!     Fn     = Cloud fraction in a 0.5 x 0.666 grid box
!                                                                             .
!     Wn     = Fraction of 0.5 x 0.666 box in the "coarse" grid box
!
! !INTERFACE:
!
  SUBROUTINE RegridTau( optIn,  optOut,  tauiIn, tauiOut, &
                        tauwIn, tauwOut, fIn,    fOut,    &
                        map,    IMX,     JMX,    LMX     )
!
! !INPUT PARAMETERS: 
!
    ! Dimensions of coarse grid
    INTEGER,      INTENT(IN)  :: IMX, JMX, LMX    

    ! Mapping weight object
    TYPE(MapObj), POINTER     :: map(:,:)

    ! Input total in-cloud optical depth (= water OD + ice OD)
    REAL*4,       INTENT(IN)  :: optIn  ( I05x0625, J05x0625, L05x0625 ) 

    ! Input in-cloud ice optical depth
    REAL*4,       INTENT(IN)  :: tauiIn ( I05x0625, J05x0625, L05x0625 )

    ! Input in-cloud water optical depth
    REAL*4,       INTENT(IN)  :: tauwIn ( I05x0625, J05x0625, L05x0625 )
 
    ! Input cloud fractions (total, anvil, convective, large-scale)
    REAL*4,       INTENT(IN)  :: fIn    ( I05x0625, J05x0625, L05x0625 )
!
! !OUTPUT PARAMETERS:
!
    ! Output total in-cloud optical depth
    REAL*4,       INTENT(OUT) :: optOut ( IMX, JMX, LMX )

    ! Output in-cloud ice optical depth
    REAL*4,       INTENT(OUT) :: tauiOut( IMX, JMX, LMX )

    ! Output in-cloud water path optical depth
    REAL*4,       INTENT(OUT) :: tauwOut( IMX, JMX, LMX )

    ! Output cloud fractions (total, anvil, convective, large-scale)
    REAL*4,       INTENT(OUT) :: fOut   ( IMX, JMX, LMX )
!
! !REVISION HISTORY: 
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: I,           J,           L,        LR,   nPoints
    INTEGER :: Nx,          Ny,          X,        Y
    REAL*4  :: sum_Fn_Wn,   sum_Wn,      Fn_Wn
    REAL*4  :: optRatio,    tauiRatio,   tauwRatio
    REAL*4  :: optRHS,      tauiRHS,     tauwRHS
    
    ! Loop over coarse grid boxes
    !$OMP PARALLEL DO                                                      &
    !$OMP DEFAULT( SHARED )                                                &
    !$OMP PRIVATE( I,           J,           L,       nPoints, sum_Fn_Wn ) &
    !$OMP PRIVATE( sum_Wn,      optRHS,      tauiRHS, tauwRHS, Nx        ) &
    !$OMP PRIVATE( Ny,          X,           Y,       Fn_Wn,   optRatio  ) &
    !$OMP PRIVATE( tauiRatio,   tauwRatio,   LR                          )
    DO L = 1, LMX

       ! Reverse level index for output arrays
       LR = LMX - L + 1

    DO J = 1, JMX
    DO I = 1, IMX

       ! Number of "fine" grid boxes in each dimension
       ! that comprise a "coarse" grid box
       nPoints = map(I,J)%nPoints

       !---------------------------------
       ! Regrid cloud fraction & OD
       !---------------------------------

       ! Zero summing variables
       sum_Fn_Wn   = 0e0
       sum_Wn      = 0e0
       optRHS      = 0e0
       tauiRHS     = 0e0
       tauwRHS     = 0e0
       
       ! Loop over "fine" grid boxes
       DO Ny = 1, nPoints
       DO Nx = 1, nPoints
          
          ! Avoid useless clock cycles if the mapping weight is zero
          IF ( map(I,J)%weight(Nx,Ny) > 0d0 ) THEN

             ! Indices of each "fine" grid box that makes up the "coarse" box
             X          = map(I,J)%xInd(Nx)
             Y          = map(I,J)%yInd(Ny)

             ! Sum of the mapping weights over all of the "fine" grid
             ! boxes (X,Y) that make up the "coarse" grid box (I,J)
             sum_Wn     = sum_Wn    + map(I,J)%weight(Nx,Ny)

             ! Compute the cloud fraction * mapping weight
             Fn_Wn      = fIn(X,Y,L) * map(I,J)%weight(Nx,Ny)

             ! Sum of the cloud fraction * mapping weights over all of the 
             ! "fine" grid boxes (X,Y) that make up the "coarse" grid box (I,J)
             sum_Fn_Wn  = sum_Fn_Wn + Fn_Wn

             ! Compute the term TAU / ( TAU + 7.7 ) for the total,
             ! ice-path, and water-path opbical depths.  NOTE: We don't have 
             ! to worry about div by zero due to the +7.7 in the denominator.
             optRatio   = optIn(X,Y,L)  / ( optIn(X,Y,L)  + 7.7e0 )
             tauiRatio  = tauiIn(X,Y,L) / ( tauiIn(X,Y,L) + 7.7e0 )
             tauwRatio  = tauwIn(X,Y,L) / ( tauwIn(X,Y,L) + 7.7e0 )
          
             ! Compute the right hand side of the equation for the total OD, 
             ! ice-path OD, and water-path OD.  For computational expediency, 
             ! we'll divide by SUM( Fn * Wn ) outside this DO loop.
             optRHS     = optRHS  + ( Fn_Wn * optRatio  )
             tauiRHS    = tauiRHS + ( Fn_Wn * tauiRatio )
             tauwRHS    = tauwRHS + ( Fn_Wn * tauwRatio )
          ENDIF
       ENDDO
       ENDDO

       ! Output cloud fraction.  NOTE, we don't have to worry about div by
       ! zero since SUM( Wn ) will always be greater than zero (there is 
       ! always at least 1 "fine" small box in the "coarse" box).
       fOut(I,J,LR) = sum_Fn_Wn / sum_Wn

       ! Total optical depth on the coarse grid
       IF ( IsSafeDiv( optRHS, sum_Fn_Wn-optRHS ) ) THEN
          optOut(I,J,LR) = 7.7e0 * optRHS / ( sum_Fn_Wn - optRHS )
       ELSE
          optOut(I,J,LR) = 0e0
       ENDIF

       ! Ice-path optical depth on the coarse grid
       IF ( IsSafeDiv( tauiRHS, sum_Fn_Wn-tauiRHS ) ) THEN
          tauiOut(I,J,LR) = 7.7e0 * tauiRHS / ( sum_Fn_Wn - tauiRHS )
       ELSE
          tauiOut(I,J,LR) = 0e0
       ENDIF

       ! Water-path optical depth on the coarse grid
       IF ( IsSafeDiv( tauwRHS, sum_Fn_Wn-tauwRHS ) ) THEN
          tauwOut(I,J,LR) = 7.7e0 * tauwRHS / ( sum_Fn_Wn - tauwRHS )
       ELSE
          tauwOut(I,J,LR) = 0e0
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE RegridTau
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsSafeDiv
!
! !DESCRIPTION: Function IsSafeDiv returns TRUE if the numerator N and 
!  denominator D may be divided safely (i.e. without resulting in a 
!  division-by-zero, Not-a-Number (NaN), or Infinity), or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsSafeDiv( N, D ) RESULT( isSafe )
!
! !INPUT PARAMETERS: 
!
    REAL*4, INTENT(IN) :: N        ! Numerator
    REAL*4, INTENT(IN) :: D        ! Denominator
!
! !RETURN VALUE:
!    
    LOGICAL            :: isSafe   ! Returns TRUE if it's safe to divide N/D 
!
! !REVISION HISTORY: 
!  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Local variables
    INTEGER :: MaxExp, MinExp

    ! Maxinum 
    MaxExp = MAXEXPONENT( N )
    MinExp = MINEXPONENT( N )
    
    ! Test if it's safe to divide 
    IF ( ( D                         == 0      )  .or. &
         ( EXPONENT(N) - EXPONENT(D) >= MaxExp )  .or. &
         ( EXPONENT(N) - EXPONENT(D) <= MinExp ) ) THEN
       isSafe = .FALSE.
    ELSE
       isSafe = .TRUE.
    ENDIF

  END FUNCTION IsSafeDiv
!EOC
END MODULE Merra2_A3CldModule

