!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine fg_vector(f,g,u1,u2,u3,u5)
      double precision :: u1,u2,u3,u5,p
      double precision, dimension(4) :: f,g
!
!
      p = 0.4d0*(u5-0.5d0*(u2*u2+u3*u3)/u1)
!
      f(1) = u2
      f(2) = u2*u2/u1 + p
      f(3) = u2*u3/u1
      f(4) = u2*(u5+p)/u1
!
      g(1) = u3
      g(2) = u2*u3/u1
      g(3) = u3*u3/u1 + p
      g(4) = u3*(u5+p)/u1
!
!
!
      end subroutine fg_vector
