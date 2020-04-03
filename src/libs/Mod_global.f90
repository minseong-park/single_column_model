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
                      z            !! Height 

    ! For nc file 
  INTEGER                       :: ncid,                       &
                                   rec_dimid, lev_dimid

  INTEGER,            PARAMETER :: dim1     = 1,               &
                                   dim2     = 2
  
  INTEGER, DIMENSION(dim1)      :: dimid1
  INTEGER, DIMENSION(dim2)      :: dimid2
  INTEGER, DIMENSION(dim1)      :: dim1_start, dim1_count
  INTEGER, DIMENSION(dim2)      :: dim2_start, dim2_count
  
  CHARACTER(LEN=256), PARAMETER :: rc_name  = "Time"            !! time
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
    IF (.NOT. ALLOCATED(Temp%dout    )) ALLOCATE(Temp%dout (nt+1,nz))

    IF (.NOT. ALLOCATED(q%dz         )) ALLOCATE(q%dz           (nz))
    IF (.NOT. ALLOCATED(q%sfc_dt     )) ALLOCATE(q%sfc_dt       (nt))
    IF (.NOT. ALLOCATED(q%top_dt     )) ALLOCATE(q%top_dt       (nt))
    IF (.NOT. ALLOCATED(q%next_dz    )) ALLOCATE(q%next_dz      (nz))
    IF (.NOT. ALLOCATED(q%dout       )) ALLOCATE(q%dout    (nt+1,nz))

    IF (.NOT. ALLOCATED(w%dz         )) THEN  
      IF ( dyn_option .eq. 1) THEN
        ALLOCATE(w%dz   (nz))
      ELSE IF ( dyn_option .eq. 2 ) THEN
        ALLOCATE(w%dz   (0:nz))
      ELSE
        CALL Fail_msg(" dyn_option must be integer // Choose either 1 or 2 ")
      ENDIF
    ENDIF
    IF (.NOT. ALLOCATED(dz%dz     )) ALLOCATE(dz%dz   (nz))
    IF (.NOT. ALLOCATED(z%dz      )) ALLOCATE(z%dz    (nz))

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

    ! Set "Description" attributes.
    Temp%desc      = "Temperature"
    q%desc         = "number of water droplets"
    w%desc         = "Vertical velocity"
    z%desc         = "Height"

    ! Set "units" attributes.
    Temp%units     = "K"
    q%units        = "  "
    w%units        = "m s-1"
    z%units        = "m"

    ! Set "axis" attributes.
    z%axis         = "Z"

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
