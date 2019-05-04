!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine itinfo(nx,ny,nl,flow,flow0,dt,it,rmsmax,rms,runtype)
!
!
      integer :: nx,ny,nl,it,runtype
      integer :: i,j
!
      double precision :: dtmax,dtmin,k,rms,rmsmax
!
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,flow0
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres,
     &         u10,u20,u30,u50,
     &         dt
!
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,flow0,u10,u20,u30,u50)
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
!
!
      rms = 0.0d0
      dtmax = -10.0d0
      dtmin =  10.0d0
      do j = 1,ny
      do i = 1,nx
!
      k = (u1(i,j)-u10(i,j))**2.0d0+(u2(i,j)-u20(i,j))**2.0d0+
     &    (u3(i,j)-u30(i,j))**2.0d0+(u5(i,j)-u50(i,j))**2.0d0 
      rms = rms + k/(dt(i,j)*dt(i,j))
!
!
      end do
      end do
!
      rms = DSQRT(rms/(nx*ny))
!
      if (rms .gt. rmsmax) then
              rmsmax = rms
      end if
!
!
!   
!
!     print information to screen if necessary
      if (runtype .EQ. 0) then
!  
      if (it == 1) then 
              print*, " "
              print*, "  ITERATION         LOG(RMS) "
      end if
      print*, it, DLOG10(rms/rmsmax)
!
      end if
!
!   
!
!
!
!
!
!
!
      end subroutine itinfo
