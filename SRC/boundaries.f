!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine boundaries(nx,ny,nl,np,params,u1,u2,u3,u5,meshX,meshY)
      integer :: nx,ny,nl,np
      integer :: i,j
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1,np) :: 
     &         params 
      double precision :: Minf,Prat,gam,ga1,
     &  Pinf,rinf,Pstag,rstag,
     &  Pext,
     &  k,Vw1,Vw2
      double precision, dimension(1:ny) :: 
     &  Pbc,rbc,Mp2,uint,Pext2
!
!
!
!
!
!
!
!     perform calculations
!-----------------------------------------------------------------------
!
!     take the flow parameters from the input info
!
      Minf = params(1,1)
      Prat = params(1,2)
      gam  = params(1,3)
      ga1  = gam - 1.0d0
!
!
!     calculate non-reflecting inflow and outflow conditions
!
      Pinf = 1.0d0/gam
      rinf = 1.0d0
      Pstag=Pinf*(1.0d0+0.5d0*ga1*Minf*Minf)**(gam/ga1)
      rstag=rinf*(1.0d0+0.5d0*ga1*Minf*Minf)**(1.0d0/ga1)
!
      do i = 1,ny
      uint(i) = u2(1,i)/u1(1,i)
      Mp2(i) = (uint(i)**2.0d0)/
     & (ga1*(1.0d0/ga1 + 0.5d0*Minf**2.0d0-0.5d0*uint(i)**2.0d0))
      Pbc(i) = Pstag/(1.0d0+0.5d0*ga1*Mp2(i))**(gam/ga1)
      rbc(i) = rstag/(1.0d0+0.5d0*ga1*Mp2(i))**(1.0d0/ga1)
      Pext2(i) = Pbc(i)*prat
      end do
!
      Pext = Pinf*prat
!
!
!
!
!
!     inflow conditions;
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,ny
      u1(1-j,i) = rbc(i)
      u2(1-j,i) = rbc(i)*uint(i)
      u3(1-j,i) = 0.0d0
!      u3(1-j,i) = u3(1,i)
      u5(1-j,i) = Pbc(i)/ga1 + 0.5d0*rbc(i)*uint(i)**2.0d0
      end do
      end do
!
!
!
!
!
!
!     outflow conditions; transient for all but u5
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,ny
      u1(nx+j,i) = u1(nx,i)
      u2(nx+j,i) = u2(nx,i)
      u3(nx+j,i) = 0.0d0
!      u3(nx+j,i) = u3(nx,i)
      u5(nx+j,i) = Pext/ga1 + 0.5d0*(u2(nx,i)**2.0d0)/u1(nx,i)
      end do
      end do
!
!
!
!
!
!
!
!     solid wall; use transient solutions
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,nx
!     bottom wall
      u1(i,1-j) = u1(i,1)
      u5(i,1-j) = u5(i,1)
!     top wall
      u1(i,ny+j) = u1(i,ny)
      u5(i,ny+j) = u5(i,ny)
      end do
      end do
!
!
!
!
!
!
!     solid wall; set normal velocity to zero
!-----------------------------------------------------------------------
!
!
      do j = 1,nl
      do i = 1,nx
!
!
!     bottom wall
!
      Vw1 = meshX(i,0) - meshX(i-1,0)
      Vw2 = meshY(i,0) - meshY(i-1,0)
      k = 2.0d0*(Vw1*u2(i,1)+Vw2*u3(i,1))/(Vw1*Vw1+Vw2*Vw2)
!
      u2(i,1-j) = k*Vw1 - u2(i,1)
      u3(i,1-j) = k*Vw2 - u3(i,1)
!
!
!
!     top wall
!
      Vw1 = meshX(i,ny) - meshX(i-1,ny)
      Vw2 = meshY(i,ny) - meshY(i-1,ny)
      k = 2.0d0*(Vw1*u2(i,ny)+Vw2*u3(i,ny))/(Vw1*Vw1+Vw2*Vw2)
!
      u2(i,ny+j) = k*Vw1 - u2(i,ny)
      u3(i,ny+j) = k*Vw2 - u3(i,ny)
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
      end subroutine boundaries
