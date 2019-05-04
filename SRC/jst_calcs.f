!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine jst_calcs(dt,np,params,da,area1,area0,
     &                    u12,u11,u10,u1i,
     &                    u22,u21,u20,u2i,
     &                    u32,u31,u30,u3i,
     &                    u52,u51,u50,u5i
     & )
!
!
!
!
!
      integer :: np
      double precision, dimension(1,np) :: params
      double precision, dimension(4) :: da
      double precision :: u12,u11,u10,u1i,
     &                    u22,u21,u20,u2i,
     &                    u32,u31,u30,u3i,
     &                    u52,u51,u50,u5i,
     & area1,area0,dt,
     &                     p2, p1, p0, pi,
     & E2,E4,K2,K4,v1,v0,h
!
!
!
!
!
!     define constants and pressure sensors
!
      K2 = params(1,6)
      K4 = params(1,7)
!
      p2 = 0.4d0*(u52-0.5d0*(u22**2.0d0+u32**2.0d0)/u12)
      p1 = 0.4d0*(u51-0.5d0*(u21**2.0d0+u31**2.0d0)/u11)
      p0 = 0.4d0*(u50-0.5d0*(u20**2.0d0+u30**2.0d0)/u10)
      pi = 0.4d0*(u5i-0.5d0*(u2i**2.0d0+u3i**2.0d0)/u1i)
!
      v1 = DABS(p2-2.0d0*p1+p0)/(DABS(p2)+2.0D0*DABS(p1)+DABS(p0))
      v0 = DABS(p1-2.0d0*p0+pi)/(DABS(p1)+2.0D0*DABS(p0)+DABS(pi))
!
      E2 = K2*DMAX1(v1,v0)
      E4 = DMAX1(0.0d0,K4-E2)
!
      h = 0.5d0*(area1+area0)
!
!
!
!
!
!     now calculate the actual dissipation terms
!
      da(1) = (h/dt)*( E2*(u11-u10) - E4*(u12-3.0d0*u11+3.0d0*u10-u1i) )
      da(2) = (h/dt)*( E2*(u21-u20) - E4*(u22-3.0d0*u21+3.0d0*u20-u2i) )
      da(3) = (h/dt)*( E2*(u31-u30) - E4*(u32-3.0d0*u31+3.0d0*u30-u3i) )
      da(4) = (h/dt)*( E2*(u51-u50) - E4*(u52-3.0d0*u51+3.0d0*u50-u5i) )
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
      end subroutine jst_calcs
