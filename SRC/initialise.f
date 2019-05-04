!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine initialise(nx,ny,nl,flow,np,params,meshX,meshY)
      integer :: nx,ny,nl,np
      integer :: i,j
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres
      double precision :: Minf,Prat,gam,Pinf,Pout,Ptot,rinf,rout
!
!
!
!
!
!     create required variables
!
      Minf = params(1,1)
      Prat = params(1,2)
      gam  = params(1,3)
!
      Pinf = 1.0d0/gam
      Pout = Pinf*Prat
      rinf = 1.0d0
      rout = rinf*Prat**(1/gam)
      Ptot = Pinf + 0.5d0*rinf*Minf*Minf
!
!
!
!
!     create pressure vector
!
      do j = 1,ny
      do i = 1,nx
      pres(i,j) = Pinf+i*(Pout-Pinf)/(nx+1)
      end do
      end do
!
!
!
!
!     assign body cells: linearly interpolate from inflow to outflow
!
      do j = 1,ny
      do i = 1,nx
      u1(i,j) = rinf + i*(rout-rinf)/(nx+1)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      u2(i,j) = u1(i,j)*((Ptot-pres(i,j))/(0.5d0*u1(i,j)))**0.5d0
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      u3(i,j) = 0.0000010d0
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      u5(i,j) = pres(i,j)/(gam-1.0d0) + 
     & 0.5d0*u1(i,j)*(Ptot-pres(i,j))/(0.5d0*u1(i,j))
      end do
      end do
!
!
!
!
!
!     apply boundary conditions prior to the residual calculation
!
      call boundaries(nx,ny,nl,np,params,u1,u2,u3,u5,meshX,meshY)
!
!
!
!     COMBINE the local functions u1 u2 u3 u5 into a single one, flow
!
      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
!
!
!
!
!
!
!
      end subroutine initialise
