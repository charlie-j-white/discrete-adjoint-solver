!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      integer :: nx,ny,nl,nh,nb
      integer :: i,j,R
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)), intent(in) 
     & :: flow
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5
!      double precision :: 
!
!
!                                               1  2  3  4  5
!       1 2 3 4 5 6 7 ... 13 14 15     ==>      6  7  8  9  10
!                                               11 12 13 14 15
!
!
!
!
      nh = 2*nl*(nx+ny)
      nb = nx*ny
!
!
!
!
!
!
!
!
!     assign first flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 0
      u1(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = nl*nx
      u1(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx
      u1(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx+nl*ny
      u1(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign second flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 1*nh+0
      u2(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 1*nh+nl*nx
      u2(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx
      u2(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx+nl*ny
      u2(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign third flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 2*nh+0
      u3(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 2*nh+nl*nx
      u3(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx
      u3(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx+nl*ny
      u3(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!     assign fifth flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 3*nh+0
      u5(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 3*nh+nl*nx
      u5(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx
      u5(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx+nl*ny
      u5(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign body cells
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+0*nb
      u1(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+1*nb
      u2(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+2*nb
      u3(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+3*nb
      u5(i,j) = flow(1,R+(j-1)*nx+i)
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
      end subroutine split_fwd
