PROGRAM main_prog

  USE NETCDF
  USE Mod_global
  USE Mod_read
  USE Mod_init
  USE Mod_dyn_driver    
  USE Mod_integration
  !USE Mod_write
  !
  IMPLICIT NONE
  INTEGER :: nnt

  CALL Sub_read_namelist

  CALL Sub_allocate
  
  CALL Sub_set_grid
  CALL Sub_Cal_W
  CALL Sub_Set_T
  p%sfc_dt = 100000. !! unit = Pa (1hPa = 100 Pa)
  CALL Sub_Cal_P(1, nz, p%sfc_dt(1), Temp%dz, dz%dz, dp%dz, dlnp%dz, p%dz, density%dz)

  CALL Sub_Cal_potental_Temp(1, p%sfc_dt(1), p%dz, Temp%dz, poten_T%dz)

  DO iz =1, SIZE(p%dz)
    print*, z%dz(iz),p%dz(iz), Temp%dz(iz), poten_T%dz(iz)
  ENDDO

  q%dz(20:30)=100.
  q%sfc_dt(:)=0.
  q%top_dt(:)=0.

  IF (dyn_option .eq. 1) THEN
    CALL Sub_Integration_FD
  ELSE IF (dyn_option .eq. 2) THEN
    CALL Sub_Integration_FV
  ENDIF

  ! write(*,*)  nz, nt, w%dz
  ! write(*,*)  dz%dz, z%dz
    do it = 1, nt+1
    write(*,*) temp%dout(it, 1)
    enddo
!    write(*,*) temp%dout(5, :)
!   write(*,*) q%dout(30,:)
!   write(*,*) sum(q%dout(30,:))

!  CALL Sub_write_netcdf ( nz, nt, dz%dz, z%dz,      &
!                          temp%dout, q%dout,        &
!                          w%dz(1:nz),               &
!                          output_path, output_name )     

   nnt = nt + 1
   OPEN(10, FILE='./z.bin', STATUS = 'unknown', &
        FORM='unformatted',ACCESS='direct', RECL=nz*4)

   WRITE(10, REC=1) z%dz
   CLOSE(10)


   OPEN(10, FILE='./w.bin', STATUS = 'unknown', &
        FORM='unformatted',ACCESS='direct', RECL=SIZE(w%dz)*4)

   WRITE(10, REC=1) w%dz
   CLOSE(10)


   OPEN(10, FILE='./temp.bin', STATUS = 'unknown', &
        FORM='unformatted',ACCESS='direct', RECL=nnt*nz*4)

   WRITE(10, REC=1) temp%dout
   CLOSE(10)

   OPEN(10, FILE='./q.bin', STATUS = 'unknown', &
        FORM='unformatted',ACCESS='direct', RECL=nnt*nz*4)
   WRITE(10, REC=1) q%dout(1:nnt, 1:nz)
   CLOSE(10)

  !CALL Sub_deallocate
END PROGRAM main_prog
