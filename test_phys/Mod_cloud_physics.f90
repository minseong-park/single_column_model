PROGRAM Mod_cloud_physics

  IMPLICIT NONE
  INTEGER, PARAMETER :: nz=100, nbin=100
  REAL, PARAMETER :: rmin= 0.1, rratio = 1.05 
  INTEGER, PARAMETER :: size_dist_option = 1
 !! 1 : log-nomal distribution (using sampledata)
 !! 2 : log-nomal distribution (using Population)
 !! 3 : gamma distribution 
  REAL :: rmax
  REAL, DIMENSION(nz, nbin) :: conc
  REAL, DIMENSION(nbin) :: r, m, n_r
  REAL, DIMENSION(nbin+1) :: rb, mb
  INTEGER :: ir
  REAL :: pi
  REAL :: N0, r0, sigma

  pi = ACOS(-1.)
!  pi = 4.0*atan(1.0)

  rmax = rmin*(rratio**nbin)

  rb(1) = rmin
  rb(nbin+1) = rmax
  DO ir = 2, nbin
    rb(ir) = rratio*rb(ir-1)
  ENDDO
!print*, rb

  mb = rb**3

!print*, '    '
!print*, mb     
  DO ir = 1, nbin
    m(ir) = (mb(ir) + mb(ir+1))/2.
  END DO

  r = m**(1./3.)   

!  DO ir = 1, nbin
!    r(ir) = rmin*(rratio**ir)
!  ENDDO

  !--------------------------------------!

  SELECT CASE(size_dist_option)
    CASE(1)

      !! Using sample value
      r0 = SUM(r)/REAL(nbin)
      sigma = SQRT(SUM((r-r0)**2)/REAL(nbin))
      N0 = 100.
      DO ir = 1, nbin
        n_r(ir) = (N0/(SQRT(2*pi)*log(sigma)*r(ir))) *  EXP(-1.*((log(r(ir))-log(r0))**2)/(2*(log(sigma)**2)))
      ENDDO

    CASE(2)
      !! Using Population value
      r0 = -1. 
      sigma = 0.4 
      N0 = 1.

      DO ir = 1, nbin
        n_r(ir) = (N0/(SQRT(2*pi)*sigma*r(ir))) * EXP(((log(r(ir))-r0)**2)/(2.*(sigma**2))*-1.)
      ENDDO
    CASE(3)
      !! gamma distribution
!      N0 = ??
!      mu = ??
!      lambda = ??
!
!      DO ir = 1, nbin
!        n_r(ir) = N0*r(ir)**mu*exp(-1.0*lambda*r(ir))
!      ENDDO
     CASE DEFAULT
       PRINT*, 'please select option of size distribution' 
       PRINT*, '1 : log-nomal distribution (using sampledata)'
       PRINT*, '2 : log-nomal distribution (using Population)'
       PRINT*, '3 : gamma distribution                       '

   END SELECT
END PROGRAM Mod_cloud_physics
