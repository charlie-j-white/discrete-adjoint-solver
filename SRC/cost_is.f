!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine cost_is(nx,ny,nl,na,np,params,flow,alpha,Jcost)
!
!
      integer :: nx,ny,nl,na,np
      integer :: i
!
      double precision :: ds, Jcost
!
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres
!
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!
!
!     cost function - just integrate pressure along walls. assume line
!     of symmetry
!
      do i = 1,nx
!
      ds = (meshX(i,ny)+meshX(i-1,ny))**2.0d0+
     &     (meshY(i,ny)+meshY(i-1,ny))**2.0d0
      ds = DSQRT(ds)
!
      Jcost = Jcost + 2*ds*pres(i,ny)
!      
      end do
!
!
!
!
!
!
!
!
      end subroutine cost_is
