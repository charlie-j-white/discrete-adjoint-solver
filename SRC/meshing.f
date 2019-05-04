!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!    
      integer:: nx,ny,na,nl,np
!    
      integer :: n,i,j,k
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) ::
     &    meshX,meshY
!
      double precision, dimension(0:na+1) :: xc,yc,a,c,l,z
      double precision, dimension(0:na) :: h,mu,b,d
      double precision, dimension(1:na) :: be
      double precision, dimension(1:nx+1) :: xt,yt,phi,phi_ig,ps,xi,y0
!
      double precision :: Lx,Ly,s_pos,s_hgt,s_wdt
!
!
!
!
!
!
!     define numbers
!-----------------------------------------------------------------------
!
!     the design variables are the spline control points for the
!     y-coordinate of the mesh. they are initialised in the shape of a
!     squared cosine wave with a height of 1.0
!
!       *                                                               *
!   fixed       *                                               *    fixed
!          1 of na                                           na of na
!
!                       *                               *
!                  2 of na                             . . .
!
!                              *                *
!                       3 of na        *
!                                   . . .
!
!     the start and end points should NOT be provided, they are
!     accounted for by the meshing subroutine
!
!
      n = na + 1
!      nx = nx + 1
      Lx = 1.0d0
      Ly = params(1,12)
      s_pos = params(1,9)
      s_hgt = params(1,10)
      s_wdt = params(1,11)
!
!
!     initialise control points:
!
      xc(0) = 0.0d0
      xc(n) = 1.0d0
      yc(0) = 1.0d0
      yc(n) = 1.0d0
!
      do i = 1,na
      xc(i) = DBLE(i)/(na+1.0d0)
      yc(i) = alpha(1,i)
      end do
!
!
!
!
!
!
!
!     create spline polynomial
!-----------------------------------------------------------------------
!
!     first arrays
!
      do i = 0,n
      a(i) = yc(i)
      end do
!
      do i = 0,n-1
      h(i) = xc(i+1) - xc(i)
      end do
!
      do i = 1,n-1
      be(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1)
      end do
!
!
!
!     next stage of arrays
!
      l(0) = 1.0d0
      mu(0) = 0.0d0
      z(0) = 0.0d0
!
      do i = 1,n-1
      l(i) = 2*(xc(i+1)-xc(i-1)) - h(i-1)*mu(i-1)
      mu(i) = h(i)/l(i)
      z(i) = (be(i) - h(i-1)*z(i-1))/l(i)
      end do
!
      l(n) = 1.0d0
      z(n) = 0.0d0
!
!
!     final arrays
!
      c(n) = 0.0d0
!
      do i = 1,n
      j = n-i
      c(j) = z(j) - mu(j)*c(j+1)
      b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2.0d0*c(j))/3.0d0
      d(j) = (c(j+1)-c(j))/(3.0d0*h(j))
      end do
!
!
!
!
!
!     form the actual spline using polynomial coefficients
      do i = 1,nx+1
!
!
!     the actual points that form the spline
      xt(i) = DBLE(i-1)*Lx/(nx)
!
!
!     determine what j value to use; check through control points
      do k = 0,n
      if (xc(k) .GT. xt(i)) then
        j = k-1
        goto 33
      else
        j = 0
      end if
      end do
33    continue
!
!     last value of yt seems buggy, add conditional statement
      if ( nx+1 .EQ. i ) then
        j = na
      end if
!
!
!
!     do the actual spline formula
      yt(i) = a(j)+
     &        b(j)*(xt(i)-xc(j))**1.0d0 +
     &        c(j)*(xt(i)-xc(j))**2.0d0 +
     &        d(j)*(xt(i)-xc(j))**3.0d0
!
!
!
!
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
!
!
!
!
!     adaptive meshing IF PARAMS SAYS SO
!-----------------------------------------------------------------------
      if (1 .EQ. 1) then
!
!     create potential function and integrate it
!
      phi_ig(1) = 0.0d0
!
      do i = 2,nx+1
      phi(i) = 1.0d0/yt(i) + s_hgt*3**(-s_wdt*(i-s_pos*(nx+1))**2)
      phi_ig(i) = phi_ig(i-1) + phi(i)
      end do
!
!
!     scale the potential function
!
      do i = 1,nx+1
      phi_ig(i) = (phi_ig(i) - phi_ig(1))
     &           /(phi_ig(nx+1) - phi_ig(1))
      end do
!
!
!
!
!
!
!     seed the CDF and interpolate the values back to the x axis
!
      do i = 1,nx+1
!
!     calculate seed points to use
      ps(i) = (DBLE(i)-1)/DBLE(nx)
      y0(i) = 0.0d0
!
!     compare psi_ig values to find points for linear interpolation
      do k = 1,nx+1
!
!     if a seed point is equal to an integrated point, avoid /0
      if (phi_ig(k) .EQ. ps(i)) then
        xi(i) = xt(k)
        goto 44
!
!     else perform linear interpolation
      else if (phi_ig(k) .GT. ps(i)) then
        xi(i) = xt(k-1) + (ps(i)-phi_ig(k-1))*
     & (xt(k)-xt(k-1))/(phi_ig(k)-phi_ig(k-1))
        goto 44
      end if
!
!     exit checking loop
      end do
44    continue
!
!     continue seed check loop
      end do
!
!
!
!
!
!
!     write to relevant vector for consistency
!
      do i = 1,nx
      xt(i) = xi(i)
      end do
!
!
!
!
!
!
!     have to now re-fit y according to the spline
      do i = 1,nx+1
!
!     determine what j value to use; check through control points
      do k = 0,n
      if (xc(k) .GT. xt(i)) then
        j = k-1
        goto 55
      else
        j = 0
      end if
      end do
55    continue
!
!     sort out end of yt array
      if ( i .EQ. nx+1) then
        j = na
      end if
!
!
!
!
!     do the actual spline formula
      yt(i) = a(j)+
     &        b(j)*(xt(i)-xc(j))**1.0d0 +
     &        c(j)*(xt(i)-xc(j))**2.0d0 +
     &        d(j)*(xt(i)-xc(j))**3.0d0
!
!
!
      end do
!
!
!
!
!
!
!
!     end of the adaptive meshing section ---------
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
!
!
!
!
!     create mesh_ using linear spacing; scale to account for halo cells
!-----------------------------------------------------------------------
!
!
!
!     meshX: body x cells
!
      do i = 0,nx
      do j = 0-nl,ny+nl
      meshX(i,j) = xt(i+1)
      end do
      end do 
!
!
!     halo x cells
!
      do i = 1,nl
      do j = 0-nl,ny+nl
      meshX(0-i,j) = -DBLE(i)*Lx/(nx+1)
      end do
      end do 
!
      do i = 1,nl
      do j = 0-nl,ny+nl
      meshX(nx+i,j) = Lx + DBLE(i)*Lx/(nx+1)
      end do
      end do 
!
!
!
!
!
!     meshY:
!
      do i = 0-nl,nx+nl
      do j = 0-nl,ny+nl
      meshY(i,j) = Ly*(j*1.0d0/DBLE(ny) - 0.5d0)
      end do
      end do 
!
      do i = 0,nx
      do j = 0-nl,ny+nl
      meshY(i,j) = yt(i+1)*Ly*(j*1.0d0/DBLE(ny) - 0.5d0)
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
      end subroutine meshing
