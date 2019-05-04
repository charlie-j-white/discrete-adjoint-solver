!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine postprocess(nx,ny,nl,flow,u1,u2,u3,u5,meshX,meshY,
     &  rmsmax,flow0)
      integer :: nx,ny,nl
      integer :: i,j
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,flow0
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres,mach,rms,
     &         u10,u20,u30,u50
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision :: gam,rmsmax
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,flow0,u10,u20,u30,u50)
!
!
!
!
!
!     calculate corner solutions for file average
!
      i = 0
      j = ny+1
      u1(i,j) = u1(i,j-1)
      u2(i,j) = u2(i,j-1)
      u3(i,j) = u3(i,j-1)
      u5(i,j) = u5(i,j-1)
!
      i = 0
      j = 0
      u1(i,j) = u1(i,j+1)
      u2(i,j) = u2(i,j+1)
      u3(i,j) = u3(i,j+1)
      u5(i,j) = u5(i,j+1)
!
      i = nx+1
      j = ny+1
      u1(i,j) = u1(i,j-1)
      u2(i,j) = u2(i,j-1)
      u3(i,j) = u3(i,j-1)
      u5(i,j) = u5(i,j-1)
!
      i = nx+1
      j = 0
      u1(i,j) = u1(i,j+1)
      u2(i,j) = u2(i,j+1)
      u3(i,j) = u3(i,j+1)
      u5(i,j) = u5(i,j+1)
!
!
!
!
!
!     calculate local mach number and log(rms) at each point
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
      gam = 1.4d0
      rmsmax = rmsmax/(DBLE(nx)*DBLE(ny))
!
      do i = 0,nx+1
      do j = 0,ny+1
!
      mach(i,j) = DSQRT(u2(i,j)**2.0d0+u3(i,j)**2.0d0)/
     &            DSQRT(u1(i,j)*gam*pres(i,j))
!
      rms(i,j) = (u1(i,j) - u10(i,j))**2.0d0 + 
     &           (u2(i,j) - u20(i,j))**2.0d0 + 
     &           (u3(i,j) - u30(i,j))**2.0d0 + 
     &           (u5(i,j) - u50(i,j))**2.0d0
      rms(i,j) = DLOG10(DSQRT(rms(i,j))/rmsmax)
!
!
!
      end do
      end do
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
!
!     write info to file: transect.plt
!
      open(100, file='transect.plt')
      write(100,*) 'VARIABLES = "x" "u1" "u2" "pres"'
101   format(4f12.8)
!
!
      j = nint((dfloat(ny)+1.0d0)/2.0d0)

      do i = 0,nx
            write(100,101) 
     &                  meshX(i,j),
     &                  0.5d0*(u1(i,j)+u1(i+1,j)),
     &                  0.5d0*(u2(i,j)+u2(i+1,j)),
     &                  0.5d0*(pres(i,j)+pres(i+1,j))
      end do
!
      close(100)
!
!
!
!
!
!
!
!
!     write info to file: solution.plt, full 2D tecplot formatting
!
      open(300, file='solution.plt')
      write(300,*) 'VARIABLES = "X" "Y" "RHO" "U" "V" "pres" "mach" 
     & "log(rms)"'
      write(300,*) 'Zone I=',ny+1,', J=',nx+1,', F=POINT'
301   format(8f13.8)
!
      do i = 0,nx
      do j = 0,ny
!
      write(300,301)
     &     meshX(i,j),
     &     meshY(i,j),
     &     0.25d0*(u1(i,j)+u1(i+1,j)+u1(i,j+1)+u1(i+1,j+1)),
     &     0.25d0*(u2(i,j)/u1(i,j)+u2(i+1,j)/u1(i+1,j)
     &         +   u2(i,j+1)/u1(i,j+1)+u2(i+1,j+1)/u1(i+1,j+1)),
     &     0.25d0*(u3(i,j)/u1(i,j)+u3(i+1,j)/u1(i+1,j)
     &         +   u3(i,j+1)/u1(i,j+1)+u3(i+1,j+1)/u1(i+1,j+1)),
     &     0.25d0*(pres(i,j)+pres(i+1,j)+pres(i,j+1)+pres(i+1,j+1)),
     &     0.25d0*(mach(i,j)+mach(i+1,j)+mach(i,j+1)+mach(i+1,j+1)),
     &     0.25d0*(rms(i,j)+rms(i+1,j)+rms(i,j+1)+rms(i+1,j+1))
!
      end do
      end do
!
      close(300)
!
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!
      end subroutine postprocess
