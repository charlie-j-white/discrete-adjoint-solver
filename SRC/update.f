!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine update(nx,ny,nl,flow,residual,np,params,
     & dt,flow0)
!
!
      integer :: nx,ny,nl,np
      integer :: i,j
!
      double precision :: gam,CFL,dtmax
!
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,residual
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,
     &         u10,u20,u30,u50,
     &         r1,r2,r3,r5,
     &         dt
!
!
!
!
      gam = params(1,3)
      CFL = params(1,4)
      dtmax = 10000.00000d0
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,flow0,u10,u20,u30,u50)
      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!
!
!
!
!
!
!
!
!
!     do the actual update procedure using first order Euler time
!     integration
!
      do j = 1,ny
      do i = 1,nx
!
!
      u1(i,j) = u10(i,j) - dt(i,j)*r1(i,j)
      u2(i,j) = u20(i,j) - dt(i,j)*r2(i,j)
      u3(i,j) = u30(i,j) - dt(i,j)*r3(i,j)
      u5(i,j) = u50(i,j) - dt(i,j)*r5(i,j)
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
      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
!
!
!
!
      end subroutine update
