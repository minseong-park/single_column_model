MODULE Mod_global

  USE NETCDF
  IMPLICIT NONE
   
  INTEGER :: it, iz   !! do parameter

    !! for namelist val
  INTEGER            :: nt,                 &
                        nz,                 &
                        dt,                 &
                        ionum,              &
                        output_interval

  INTEGER            :: dyn_option,         &
                        dz_option

  REAL               :: z_top

  REAL               :: gamma_dry

  REAL               :: dzr

  CHARACTER(LEN=256) :: output_path, &
                        output_name

    ! Declare variables 
    INTEGER                           :: varid
  TYPE varinfo
    INTEGER                           :: varid
    REAL, DIMENSION(:),   ALLOCATABLE :: dz, next_dz,      &
                                         dt,               &
                                         sfc_dt, top_dt  
    REAL, DIMENSION(:,:), ALLOCATABLE :: dout        
    CHARACTER(LEN=256)                :: vname, axis,      &
                                         desc, units
  END TYPE varinfo

  TYPE(varinfo) ::    Temp,      & !! Temperature [K] 
                      q,         & !! Number of water droplets
                      w,         & !! Vertical velocity
                      dz,        & !! Difference z
                      time,      & !! Difference z
                      z,         & !! Height 
                      density,   & !! density               !! by han
                      p,         & !! pressure              !! by han
                      dlnp,      & !! Difference log pressure   !! by han
                      dp,        & !! Difference pressure   !! by han
                      poten_T      !! Potential Temperature [K] 

    ! For nc file 
  INTEGER                       :: ncid,                       &
                                   rec_dimid, lev_dimid,       &
                                   lat_dimid, lon_dimid
  

  INTEGER,            PARAMETER :: dim1     = 1,               &
                                   dim4     = 4
  
  INTEGER, DIMENSION(dim1)      :: dimid1
  INTEGER, DIMENSION(dim4)      :: dimid4
  INTEGER, DIMENSION(dim1)      :: dim1_start, dim1_count
  INTEGER, DIMENSION(dim4)      :: dim4_start, dim4_count
  
  CHARACTER(LEN=256), PARAMETER :: des      = "description"
  CHARACTER(LEN=256), PARAMETER :: un       = "units"
  CHARACTER(LEN=256), PARAMETER :: ax       = "axis"

  CONTAINS

  !!-----------------------------!!
  SUBROUTINE Sub_allocate

    IF (.NOT. ALLOCATED(Temp%dz      )) ALLOCATE(Temp%dz        (nz))
    IF (.NOT. ALLOCATED(Temp%sfc_dt  )) ALLOCATE(Temp%sfc_dt    (nt))
    IF (.NOT. ALLOCATED(Temp%top_dt  )) ALLOCATE(Temp%top_dt    (nt))
    IF (.NOT. ALLOCATED(Temp%next_dz )) ALLOCATE(Temp%next_dz   (nz))
    IF (.NOT. ALLOCATED(Temp%dout    )) ALLOCATE(Temp%dout (nz,nt+1))

    IF (.NOT. ALLOCATED(q%dz         )) ALLOCATE(q%dz           (nz))
    IF (.NOT. ALLOCATED(q%sfc_dt     )) ALLOCATE(q%sfc_dt       (nt))
    IF (.NOT. ALLOCATED(q%top_dt     )) ALLOCATE(q%top_dt       (nt))
    IF (.NOT. ALLOCATED(q%next_dz    )) ALLOCATE(q%next_dz      (nz))
    IF (.NOT. ALLOCATED(q%dout       )) ALLOCATE(q%dout    (nz,nt+1))

    IF (.NOT. ALLOCATED(w%dz         )) ALLOCATE(w%dz         (0:nz))

    IF (.NOT. ALLOCATED(dz%dz        )) ALLOCATE(dz%dz          (nz))
    IF (.NOT. ALLOCATED(z%dz         )) ALLOCATE(z%dz           (nz))

    IF (.NOT. ALLOCATED(density%dz   )) ALLOCATE(density%dz     (nz))  !! byhan
    IF (.NOT. ALLOCATED(dp%dz        )) ALLOCATE(dp%dz          (nz))  !! byhan
    IF (.NOT. ALLOCATED(dlnp%dz      )) ALLOCATE(dlnp%dz        (nz))  !! byhan
    IF (.NOT. ALLOCATED(p%sfc_dt     )) ALLOCATE(p%sfc_dt       (nt))  !! byhan
    IF (.NOT. ALLOCATED(p%dz         )) ALLOCATE(p%dz           (nz))  !! byhan
    IF (.NOT. ALLOCATED(poten_T%dz   )) ALLOCATE(poten_T%dz     (nz))  !! byhan

  END SUBROUTINE Sub_allocate
 
  !!-----------------------------!!
  SUBROUTINE Sub_deallocate
  END SUBROUTINE Sub_deallocate

  !!-----------------------------!!
  SUBROUTINE Sub_nc_attri

    ! Set name of variables.
    Temp%vname     = "T"
    q%vname        = "Q"
    w%vname        = "W"
    z%vname        = "Lev"
    time%vname     = "Time"


    ! Set "Description" attributes.
    Temp%desc      = "Temperature"
    q%desc         = "mass of water droplets"
    w%desc         = "Vertical velocity"
    z%desc         = "Height"

    ! Set "units" attributes.
    Temp%units     = "K"
    q%units        = "kg kg-1"
    w%units        = "m s-1"
    z%units        = "m"
    time%units     = "minutes since 2000-01-01 00:00:00"

    ! Set "axis" attributes.
    z%axis         = "Z"
    time%axis      = "T"

  END SUBROUTINE Sub_nc_attri

  !!-----------------------------!!
  SUBROUTINE CHECK(status)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: status

    !Check errors.
    IF (status .ne. nf90_noerr) THEN
     PRINT *, trim(nf90_strerror(status))
     PRINT *, "    ERROR :: CHECK NC CODE       "
     STOP "##### ERROR: PROGRAM ABORTED. #####"
    END IF

  END SUBROUTINE CHECK

  !!-----------------------------!!
  SUBROUTINE FAIL_MSG(f_msg)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: f_msg

    WRITE (*,'("FAIL: ", a)') f_msg
    STOP "##### ERROR: PROGRAM ABORTED. #####"

  END SUBROUTINE FAIL_MSG

  !!-----------------------------!!
  SUBROUTINE SUCCESS_MSG(s_msg)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: s_msg
    
    WRITE (*,'("SUCCESS: ", a)') s_msg

  END SUBROUTINE SUCCESS_MSG 


END MODULE Mod_global
