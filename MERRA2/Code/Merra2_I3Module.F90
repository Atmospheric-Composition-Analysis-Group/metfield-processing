!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Merra2_I3Module
!
! !DESCRIPTION: Module Merra2_I3Module contains routines to create the 
!  GEOS-Chem instantaneous 3-hr data files from the MERRA2 raw data.
!\\
!\\
! !INTERFACE: 
!
MODULE Merra2_I3Module
! 
! !USES:
!
  ! GEOS-5.7.x data modules
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
  PUBLIC  :: Merra2_MakeI3
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: NcOutFileDef
  PRIVATE :: ProcessI33dAsmNv
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
    INTEGER            :: var1(1), var3(3), var4(4)

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
    lName = 'MERRA2 instantaneous 3-hour parameters (I3), processed for GEOS-Chem input'
    CALL NcDef_Glob_Attributes( fOut, 'Title',                TRIM( lName ) )

    ! Contact
    lName = "GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)"
    CALL NcDef_Glob_Attributes( fOut, 'Contact',              TRIM( lName ) )
                                                              
    ! References                                              
    lName = "www.geos-chem.org; wiki.geos-chem.org"           
    CALL NcDef_Glob_Attributes( fOut, 'References',           TRIM( lName ) )
                                                              
    ! Filename                                                
    write(IU_LOG,*) 'About to crash with outFileName:',TRIM(outFileName)
    flush(IU_LOG)
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
    lName = TRIM( VersionID )                                         
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
    var1    = (/ idLon /)
    vId     = vId + 1
    lName   = 'longitude'
    units   = 'degrees_east'
    CALL NcDef_Variable      ( fOut, 'lon', NF_FLOAT, 1, var1, vId           )
    CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName )    )
    CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units )    )

    !-------------------------------------------------------------------------
    ! Define data arrays
    !-------------------------------------------------------------------------

    ! PS
    IF ( StrPos( 'PS', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var3  = (/ idLon, idLat, idTime /)    
       vId   = vId + 1
       lName = 'surface_pressure' 
       units = 'Pa'
       gamap = 'GMAO-2D'
       CALL NcDef_Variable      ( fOut, 'PS',  NF_FLOAT, 3, var3, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! PV (aka EPV)
    IF ( StrPos( 'PV', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'ertels_potential_vorticity' 
       units = 'K m+2 kg-1 s-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'PV', NF_FLOAT, 4, var4, vId        )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! QV
    IF ( StrPos( 'QV', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'specific_humidity' 
       units = 'kg kg-1'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'QV', NF_FLOAT, 4, var4, vId         )
       CALL NcDef_Var_Attributes( fOut, vId, 'long_name',      TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'standard_name',  TRIM( lName ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'units',          TRIM( units ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'gamap_category', TRIM( gamap ) )
       CALL NcDef_Var_Attributes( fOut, vId, 'missing_value',  FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, '_FillValue',     FILL_VALUE    )
       CALL NcDef_Var_Attributes( fOut, vId, 'scale_factor',   1e0           )
       CALL NcDef_Var_Attributes( fOut, vId, 'add_offset',     0e0           )
    ENDIF

    ! T
    IF ( StrPos( 'T', inst3_3d_asm_Nv_Data ) >= 0 ) THEN
       var4  = (/ idLon, idLat, idLev, idTime /)    
       vId   = vId + 1
       lName = 'air_temperature' 
       units = 'K'
       gamap = 'GMAO-3D$'
       CALL NcDef_Variable      ( fOut, 'T', NF_FLOAT, 4, var4, vId    )
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
! !IROUTINE: Merra2_MakeI3
!
! !DESCRIPTION: Routine Merra2_MakeI3 is the the driver routine for 
! \begin{enumerate}
! \item Extracting instantaneous 3-hr data fields (surface values) from 
!       the MERRA2 raw data files (netCDF4 format),
! \item Regridding the fields to GEOS-Chem data resolution, and 
! \item Saving the regridded data to netCDF format.
! \end{enumerate}
! This routine is called directly from the main program Merra2_Driver.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_MakeI3()
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
    INTEGER                 :: F
    INTEGER                 :: nFields
    INTEGER                 :: nAllFields
    CHARACTER(LEN=MAX_CHAR) :: msg
    CHARACTER(LEN=MAX_CHAR) :: fName
    CHARACTER(LEN=MAX_CHAR) :: gName

    ! Arrays
    CHARACTER(LEN=MAX_CHAR) :: fields(MAX_FLDS)

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Echo info
    msg = '%%%%%%%%%% ENTERING ROUTINE Merra2_MakeI3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Return the list of fields and number of fields to process
    ! from each of the MERRA raw met data files
    CALL GetNFields( inst3_3d_asm_Nv_data, nFields, fields )

    ! Total number of fields that we will process
    nAllFields = nFields

    ! Echo info
    WRITE( IU_LOG, 100 ) TRIM( inst3_3d_asm_Nv_file ), nFields
    WRITE( IU_LOG, 110 ) nAllFields

    ! Formats
100 FORMAT( '%%% # of fields from ', a, ' : ',          i5 )
110 FORMAT( '%%% TOTAL # OF FIELDS TO BE REGRIDDED : ', i5 )

    !=======================================================================
    ! Open files for output; define variables, attribute, index arrays
    !=======================================================================

    ! Open 2 x 2.5 output file
    IF ( do2x25 ) THEN
       fname = TRIM( tempDirTmpl2x25 ) // TRIM( dataTmpl2x25 )
       gName = '2 x 2.5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I2x25,     J2x25,        L2x25,      TIMES_A3,   &
                          xMid_2x25, nc_yMid_2x25, zMid_2x25,  a3MinsI,    &
                          gName,     fName   ,     fOut2x25               )
    ENDIF

    ! Open 4 x 5 output file 
    IF ( do4x5 ) THEN
       fname = TRIM( tempDirTmpl4x5 ) // TRIM( dataTmpl4x5 )
       gName = '4 x 5 global'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000                 )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '               )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.                     )
       CALL NcOutFileDef( I4x5,      J4x5,         L4x5,       TIMES_A3,   &
                          xMid_4x5,  nc_yMid_4x5,  zMid_4x5,   a3MinsI,    &
                          gName,     fName,        fOut4x5                )
    ENDIF
    
    ! Open 0.5 x 0.625 output file
    IF ( doGlobal05 ) THEN
       fName = TRIM( tempDirTmpl05x0625 ) // TRIM( dataTmpl05x0625 )
       gName = '0.5 x 0.625 global'
       CALL ExpandDate  ( fName,        yyyymmdd,        000000      )      
       CALL StrRepl     ( fName,        '%%%%%%',        'I3    '    )
       CALL StrCompress ( fName,        RemoveAll=.TRUE.             )
       CALL NcOutFileDef( I05x0625,     J05x0625,        &
                          L05x0625,     TIMES_A3,        &
                          xMid_05x0625, nc_yMid_05x0625, &
                          zMid_05x0625, a3MinsI,         &
                          gName,        fName,           fOut05x0625 )
    ENDIF

    ! Open nested 0625 AS output file
    IF ( doNestAs05 ) THEN
       fName = TRIM( tempDirTmplNestAs05 ) // TRIM( dataTmplNestAs05 )
       gName = 'nested AS 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestAs05,  J_NestAs05, L05x0625,  TIMES_A3,  &
                          xMid_05x0625(I0_as05:I1_as05),          &
                          yMid_05x0625(J0_as05:J1_as05),          &
                          zMid_05x0625,                a3MinsI,    &
                          gName,    fName,        fOut05NestAs          )
    ENDIF

    ! Open nested EU output file
    IF ( doNestEu05 ) THEN
       fName = TRIM( tempDirTmplNestEu05 ) // TRIM( dataTmplNestEu05 )
       gName = 'nested EU 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestEu05,  J_NestEu05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_eu05:I1_eu05),          &
                          yMid_05x0625(J0_eu05:J1_eu05),          &
                          zMid_05x0625,                a3MinsI,    &
                          gName,    fName,       fOut05NestEu          )
    ENDIF

    ! Open nested NA output file
    IF ( doNestNa05 ) THEN
       fName = TRIM( tempDirTmplNestNa05 ) // TRIM( dataTmplNestNa05 )
       gName = 'nested NA 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestNa05,  J_NestNa05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_na05:I1_na05),          &
                          yMid_05x0625(J0_na05:J1_na05),          &
                          zMid_05x0625,                a3MinsI,    &
                          gName,    fName,      fOut05NestNa          )
    ENDIF

    ! Open nested SE output file
    IF ( doNestSe05 ) THEN
       fName = TRIM( tempDirTmplNestSe05 ) // TRIM( dataTmplNestSe05 )
       gName = 'nested SE 05'
       CALL ExpandDate  ( fName,     yyyymmdd,     000000      )
       CALL StrRepl     ( fName,     '%%%%%%',     'I3    '    )
       CALL StrCompress ( fName,     RemoveAll=.TRUE.          )
       CALL NcOutFileDef( I_NestSe05,  J_NestSe05, L05x0625, TIMES_A3,  &
                          xMid_05x0625(I0_se05:I1_se05),          &
                          yMid_05x0625(J0_se05:J1_se05),          &
                          zMid_05x0625,                a3MinsI,    &
                          gName,    fName,      fOut05NestSe         )
    ENDIF

    ! Regrid fields from the various raw data files
    CALL ProcessI33dAsmNv( nFields, fields )
    
    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Close binary files
    msg = '%%% Closing I3 output files'
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
    msg = '%%%%%%%%%% LEAVING ROUTINE Merra2_MakeI3 %%%%%%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE Merra2_MakeI3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ProcessI33dAsmNv
!
! !DESCRIPTION: Subroutine  ProcessI33dAsmNv regrids the MERRA2 met fields 
!  from the "inst3t\_3d\_asm\_Nv" file and saves to netCDF format.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ProcessI33dAsmNv( nFields, fields )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: nFields     ! # of fields to process
    CHARACTER(LEN=*), INTENT(IN) :: fields(:)   ! List of field names
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

    ! Data arrays (NOTE: 2d or 3d refers to spatial dimensions)
    REAL*4,  TARGET         :: Q2d     ( I05x0625, J05x0625            )
    REAL*4                  :: Q2d_2x25( I2x25,    J2x25               )
    REAL*4                  :: Q2d_4x5 ( I4x5,     J4x5                )
    REAL*4,  TARGET         :: Q3d     ( I05x0625, J05x0625, L05x0625  )
    REAL*4                  :: Q3d_2x25( I2x25,    J2x25,    L2x25     )
    REAL*4                  :: Q3d_4x5 ( I4x5,     J4x5,     L4x5      )

    ! Pointers
    REAL*4,  POINTER        :: Ptr_2d(:,:)
    REAL*4,  POINTER        :: Ptr_3d(:,:,:)
    REAL*4,  POINTER        :: Qflip (:,:,:)

    ! Character strings and arrays
    CHARACTER(LEN=8       ) :: name8
    CHARACTER(LEN=9       ) :: name
    CHARACTER(LEN=MAX_CHAR) :: fNameInput
    CHARACTER(LEN=MAX_CHAR) :: fNameNested
    CHARACTER(LEN=MAX_CHAR) :: fName2x25
    CHARACTER(LEN=MAX_CHAR) :: fName4x5
    CHARACTER(LEN=MAX_CHAR) :: msg

    !=======================================================================
    ! Get dimensions from output files
    !=======================================================================

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
    ! NOTE: For constant file, hardwire date to 2011/01/01
    !=======================================================================

    ! Echo info    
    msg = '%%%%%% ENTERING ROUTINE ProcessI33dAsmNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) '%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

    ! Create input filename from the template
    fNameInput = TRIM( inputDataDir ) // TRIM( inst3_3d_asm_Nv_file )
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
       hhmmss = ( a3minsI(H) / 60 ) * 10000

       !====================================================================
       ! Process data
       !====================================================================

       ! Loop over data fields
       DO F = 1, nFields
          
          ! Save field name into an 9-char variable. 
          ! This will truncate field names longer than 8 chars.
          name  = TRIM( fields(F) )
          
          ! Skip if fieldname is empty
          IF ( name == '' ) CYCLE

          ! Zero data arrays
          Q2d      = 0e0
          Q2d_2x25 = 0e0
          Q2d_4x5  = 0e0
          Q3d      = 0e0
          Q3d_2x25 = 0e0
          Q3d_4x5  = 0e0
          
          ! Test field name
          IF ( TRIM( name ) == 'PS' ) THEN
          
             !==============================================================
             ! Special handling for surface pressure data
             ! since this field is defined at the surface only
             !==============================================================

             ! Start and count index arrays for netCDF
             st3d = (/ 1, 1, H /)
             ct3d = (/ X, Y, 1 /)

             !--------------------------------------------------------------
             ! Read data
             !--------------------------------------------------------------
             msg = '%%% Reading     ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             CALL NcRd( Q2d, fIn, TRIM( name ), st3d, ct3d )
             
             ! Replace missing values with zeroes
             WHERE( Q2d == FILL_VALUE ) Q2d = 0e0

!-----------------------------------------------------------------------------
! Prior to 7/28/15:
! Now keep pressure quantities in hPa (bmy, 7/28/15)
!             ! Convert from [Pa] to [hPa]
!             Q2d = Q2d / 100e0
!-----------------------------------------------------------------------------

             !--------------------------------------------------------------
             ! Regrid data
             !--------------------------------------------------------------
             msg = '%%% Regridding  ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Regrid to 2 x 2.5
             IF ( do2x25 ) THEN
                CALL RegridMerra2_to2x25( 0, Q2d, Q2d_2x25 )
             ENDIF
             
             ! Regrid to 4x5 
             IF ( do4x5 ) THEN
                CALL RegridMerra2_To4x5 ( 0, Q2d, Q2d_4x5  )
             ENDIF
 
             !--------------------------------------------------------------
             ! Write netCDF output
             !--------------------------------------------------------------
             msg = '%%% Archiving   ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Write 2 x 2.5 data
             IF ( do2x25 ) THEN
                st3d  = (/ 1,     1,     H  /)
                ct3d  = (/ X2x25, Y2x25, 1  /)
                CALL NcWr( Q2d_2x25, fOut2x25, TRIM( name ), st3d, ct3d )
             ENDIF
       
             ! Write 4x5 data
             IF ( do4x5 ) THEN
                st3d  = (/ 1,    1,    H /)
                ct3d  = (/ X4x5, Y4x5, 1 /)
                CALL NcWr( Q2d_4x5, fOut4x5, TRIM( name ), st3d, ct3d )
             ENDIF

             ! Write 0.5 x 0.625 data
             IF ( doGlobal05 ) THEN
                Ptr_2d => Q2d(:,:)
                st3d = (/ 1,        1,        H /)
                ct3d = (/ X05x0625, Y05x0625, 1 /)
                CALL NcWr( Ptr_2d, fOut05x0625, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

             ! Nested AS (point to proper slice of global data)
             IF ( doNestAs05 ) THEN
                Ptr_2d => Q2d( I0_as05:I1_as05, J0_as05:J1_as05 )
                st3d = (/ 1,         1,         H /)
                ct3d = (/ XNestAs05, YNestAs05, 1 /)
                CALL NcWr( Ptr_2d, fOut05NestAs, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

             ! Nested EU (point to proper slice of global data)
             IF ( doNestEu05 ) THEN
                Ptr_2d  => Q2d( I0_eu05:I1_eu05, J0_eu05:J1_eu05 )
                st3d = (/ 1,         1,         H /)
                ct3d = (/ XNestEu05, YNestEu05, 1 /)
                CALL NcWr( Ptr_2d, fOut05NestEu, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

            ! Nested NA (point to proper slice of global data)
             IF ( doNestNa05 ) THEN
                Ptr_2d  => Q2d( I0_na05:I1_na05, J0_na05:J1_na05 )
                st3d = (/ 1,         1,         H /)
                ct3d = (/ XNestNa05, YNestNa05, 1 /)
                CALL NcWr( Ptr_2d, fOut05NestNa, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

            ! Nested SE (point to proper slice of global data)
             IF ( doNestSe05 ) THEN
                Ptr_2d  => Q2d( I0_se05:I1_se05, J0_se05:J1_se05 )
                st3d = (/ 1,       1,       H /)
                ct3d = (/ XNestSe05, YNestSe05, 1 /)
                CALL NcWr( Ptr_2d, fOut05NestSe, TRIM( name ), st3d, ct3d )
                NULLIFY( Ptr_2d )
             ENDIF

          ELSE
          
             !==============================================================
             ! Process all other fields
             !==============================================================

             ! Start and count index arrays for netCDF
             st4d = (/ 1, 1, 1, H /)
             ct4d = (/ X, Y, Z, 1 /)

             ! Special handling: "EPV" is known as "PV" 
             ! for backwards compatibility
             name8 = name
             IF ( TRIM( name8 ) == 'PV' ) name8 = 'EPV'

             !--------------------------------------------------------------
             ! Read data
             !--------------------------------------------------------------
             msg = '%%% Reading     ' // name8
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             CALL NcRd( Q3d, fIn, TRIM( name8 ), st4d, ct4d )
             
             ! Replace missing values with zeroes
             WHERE( Q3d == FILL_VALUE ) Q3d = 0e0

             ! Flip data in the vertical
             Qflip => Q3d( :, :, Z:1:-1 )

             !--------------------------------------------------------------
             ! Regrid data
             !--------------------------------------------------------------
             msg = '%%% Regridding  ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Loop over vertical levels
             DO L = 1, Z 

                ! Regrid to 2 x 2.5
                IF ( do2x25 ) THEN
                   CALL RegridMerra2_to2x25( 0, Qflip(:,:,L), Q3d_2x25(:,:,L) )
                ENDIF
             
                ! Regrid to 4x5 
                IF ( do4x5 ) THEN
                   CALL RegridMerra2_To4x5 ( 0, Qflip(:,:,L), Q3d_4x5(:,:,L)  )
                ENDIF
                            
             ENDDO

             !--------------------------------------------------------------
             ! Write netCDF output
             !--------------------------------------------------------------
             msg = '%%% Archiving   ' // name
             WRITE( IU_LOG, '(a)' ) TRIM( msg )
             
             ! Write 2 x 2.5 data
             IF ( do2x25 ) THEN
                st4d  = (/ 1,     1,     1,     H  /)
                ct4d  = (/ X2x25, Y2x25, Z2x25, 1  /)
                CALL NcWr( Q3d_2x25, fOut2x25, TRIM( name ), st4d, ct4d )
             ENDIF
       
             ! Write 4x5 data
             IF ( do4x5 ) THEN
                st4d  = (/ 1,    1,    1,    H /)
                ct4d  = (/ X4x5, Y4x5, Z4x5, 1 /)
                CALL NcWr( Q3d_4x5, fOut4x5, TRIM( name ), st4d, ct4d )
             ENDIF

             ! Write 0.5 x 0.625 data
             IF ( doGlobal05 ) THEN
                Ptr_3d => Qflip(:,:,:)
                st4d = (/ 1,        1,        1,        H /)
                ct4d = (/ X05x0625, Y05x0625, Z05x0625, 1 /)
                CALL NcWr( Ptr_3d, fOut05x0625, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested Asia (point to proper slice of global data)
             IF ( doNestAs05 ) THEN
                Ptr_3d => Qflip( I0_as05:I1_as05, J0_as05:J1_as05, : )
                st4d   =  (/ 1,         1,         1,         H /)
                ct4d   =  (/ XNestAs05, YNestAs05, ZNestAs05, 1 /)
                CALL NcWr( Ptr_3d, fOut05NestAs, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested EU (point to proper slice of global data)
             IF ( doNestEu05 ) THEN
                Ptr_3d => Qflip( I0_eu05:I1_eu05, J0_eu05:J1_eu05, : )
                st4d   = (/ 1,         1,         1,         H /)
                ct4d   = (/ XNestEu05, YNestEu05, ZNestEu05, 1 /)
                CALL NcWr( Ptr_3d, fOut05NestEu, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested NA (point to proper slice of global data)
             IF ( doNestNa05 ) THEN
                Ptr_3d => Qflip( I0_na05:I1_na05, J0_na05:J1_na05, : )
                st4d   = (/ 1,         1,         1,         H /)
                ct4d   = (/ XNestNa05, YNestNa05, ZNestNa05, 1 /)
                CALL NcWr( Ptr_3d, fOut05NestNa, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Nested SE (point to proper slice of global data)
             IF ( doNestSe05 ) THEN
                Ptr_3d => Qflip( I0_se05:I1_se05, J0_se05:J1_se05, : )
                st4d   = (/ 1,         1,         1,         H /)
                ct4d   = (/ XNestSe05, YNestSe05, ZNestSe05, 1 /)
                CALL NcWr( Ptr_3d, fOut05NestSe, TRIM( name ), st4d, ct4d )
                NULLIFY( Ptr_3d )
             ENDIF

             ! Free pointer memory
             NULLIFY( Qflip )
          ENDIF
       ENDDO
    ENDDO

    !=======================================================================
    ! Cleanup and Quit
    !=======================================================================

    ! Close input file
    msg = '%%% Closing ' // TRIM( fNameInput )
    WRITE( IU_LOG, '(a)' ) TRIM( msg )
    CALL NcCl( fIn )
   
    ! Echo info 
    msg = '%%%%%% EXITING ROUTINE ProcessI33dAsmNv %%%%%%'
    WRITE( IU_LOG, '(a)' ) TRIM( msg )

  END SUBROUTINE ProcessI33dAsmNv
!EOC
END MODULE Merra2_I3Module

