!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine debug_cell(nx,ny,nl,u1)
      integer :: nx,ny,nl
      integer :: j
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: u1
!
!
!
      print*, "  "
!
      do j = 1-nl,ny+nl
      print*, u1(:,j)
      end do
!
      print*, "  "
!
!
!
      end subroutine debug_cell
