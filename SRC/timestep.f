!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine timestep(nx,ny,nl,flow,np,params,meshX,meshY,dt)
!
!
      integer :: nx,ny,nl,np
      integer :: i,j
!
      double precision :: gam,CFL,vx,vy,dx,dy,dtmax,SoS,k
!
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres,
     &         dt
!
!
!
!
      gam = params(1,3)
      CFL = params(1,4)
      dtmax = 10000.00000d0
!
      k = 1.0d0
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
!
!
!     compute timestep: SoS, local speed and local time
!
      do j = 1,ny
      do i = 1,nx
!
      SoS = DSQRT(DABS(gam*pres(i,j)/u1(i,j)))
      vx = DABS(u2(i,j)/u1(i,j)) +DABS(SoS)
      vy = DABS(u3(i,j)/u1(i,j)) +DABS(SoS)
!
      dx = DMAX1(DABS(meshX(i,j) - meshX(i-1,j)),
     &           DABS(meshX(i-1,j) - meshX(i-2,j)))
      dx = DMAX1(DABS(meshX(i+1,j) - meshX(i,j)), dx)*k
      dy = DMAX1(DABS(meshY(i,j) - meshY(i,j-1)),
     &           DABS(meshY(i,j-1) - meshY(i,j-2)))
      dy = DMAX1(DABS(meshY(i,j+1) - meshY(i,j)), dy)*k
!
      dt(i,j) = CFL/(vx/dx + vy/dy)
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
      end subroutine timestep
