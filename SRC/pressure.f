!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine pressure(nx,ny,nl,pres,u1,u2,u3,u5)
      integer :: nx,ny,nl
      integer :: i,j
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         pres
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl), intent(in) :: 
     &         u1,u2,u3,u5
      double precision :: ga1
!
!
      ga1 = 0.4d0
!
!
!
      do j = 1-nl,ny+nl
      do i = 1-nl,nx+nl
!
      if (u1(i,j).NE.0.0d0) then
      pres(i,j) = ga1*(u5(i,j) - 
     & 0.5d0*(u2(i,j)**2.0d0+u3(i,j)**2.0d0)/u1(i,j))
      end if
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
      end subroutine pressure
