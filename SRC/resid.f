!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine resid(nx,ny,nl,flow,residual,np,params,meshX,meshY,dt)
!
!
      integer :: nx,ny,nl,np
      integer :: i,j
!
      double precision :: xa,xb,xc,xd,ya,yb,yc,yd,area
      double precision :: pi
!
      double precision, dimension(1,np) :: params
      double precision, dimension(4) ::
     & fa,fb,fc,fd,ga,gb,gc,gd,
     & da,db,dc,dd
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,residual
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,
     &         r1,r2,r3,r5,
     &         volume,dt
!
!
!
!
!
!
!
!
!
!     Visualise the cell in three dimensions to see cross product
!
!                                mesh(i-1,j)         mesh(i,j)
!       z|                             .______________.
!        |   /y                       /       B      ^           
!        |  /                        /              /
!        | /                        /C            A/--------->n_A
!        |/______                  /              /
!               x                 /_______D______/
!                         mesh(i-1,j-1)          mesh(i,j-1)
!
!
!     Make sure to form the correct "face triangles" to get the right
!     normal
!      
!                             xA
!              .__________.<____.    
!     y|        \     B    \    ^  
!      |         \          \   |                        xA   0     yA
!      |          \C        A\  |yA     n.dS = A X k  =  yA X 0  = -xA 
!      |_____      \          \ |                         0   1      0
!           x       \____D_____\|
!
!
!     Area of a triangle between vectors P and Q for use in the cell
!     area calculation
!                                           area = 0.5*|PXQ|
!                                                = 0.5*|xP.yQ-xQ.yP|
!        
!
!
!        
!     When caulation the fluxes for each face use the definitions below
!     for example F_A = F(U(i+1/2,j))   
!        
!
!         *-------------*-------------*-------------*
!         |             |             |             |
!         |             |    i,j+1    |             |    A = (i+1/2,j)
!         |             |             |             |
!         *-------------*------B------*-------------*    B = (i,j+1/2)
!         |             |             |             |
!         |    i-1,j    C     i,j     A    i+1,j    |    C = (i-1/2,j)
!         |             |             |             |
!         *-------------*------D------*-------------*    D = (i,j-1/2)
!         |             |             |             |
!         |             |    i,j-1    |             |
!         |             |             |             |
!         *-------------*-------------*-------------*
!
!     A face: i+1/2 , j
!     B face: i     , j+1/2
!     C face: i-1/2 , j
!     D face: i     , j-1/2
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
!     apply boundary conditions and calculate pressure field
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!
      call boundaries(nx,ny,nl,np,params,u1,u2,u3,u5,meshX,meshY)
!
!      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
      do i = 1,4
      da(i) = 0.0d0
      db(i) = 0.0d0
      dc(i) = 0.0d0
      dd(i) = 0.0d0
      end do
!
!
!
!     calculate cell area
!
      pi = 3.1415926535897932d0
!
      do j = 1-nl,ny+nl
      do i = 1-nl,nx+nl
!
!
      xa = meshX(i,j) - meshX(i,j-1)
      xb = meshX(i-1,j) - meshX(i,j)
      xc = meshX(i-1,j-1) - meshX(i-1,j)
      xd = meshX(i,j-1) - meshX(i-1,j-1)
      ya = meshY(i,j) - meshY(i,j-1)
      yb = meshY(i-1,j) - meshY(i,j)
      yc = meshY(i-1,j-1) - meshY(i-1,j)
      yd = meshY(i,j-1) - meshY(i-1,j-1)
!
!
!
      volume(i,j) = 0.5d0*DABS(xa*yb-xb*ya) + 0.5d0*DABS(xc*yd-xd*yc)
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
!     Perform finite volume calculation over body cells
!
      do j = 1,ny
      do i = 1,nx
!
!     calculate geometric information
!
      xa = (meshX(i,j) - meshX(i,j-1))
      xb = (meshX(i-1,j) - meshX(i,j))
      xc = (meshX(i-1,j-1) - meshX(i-1,j))
      xd = (meshX(i,j-1) - meshX(i-1,j-1))
      ya = (meshY(i,j) - meshY(i,j-1))
      yb = (meshY(i-1,j) - meshY(i,j))
      yc = (meshY(i-1,j-1) - meshY(i-1,j))
      yd = (meshY(i,j-1) - meshY(i-1,j-1))
!
      area = volume(i,j) 
!
!
!
!
!
!
!     create F and G vector
!
!     A face: i+1/2 , j
      call fg_vector(fa,ga,
     & 0.5d0*(u1(i+1,j)+u1(i,j)),
     & 0.5d0*(u2(i+1,j)+u2(i,j)),
     & 0.5d0*(u3(i+1,j)+u3(i,j)),
     & 0.5d0*(u5(i+1,j)+u5(i,j)))
!
!     B face: i     , j+1/2
      call fg_vector(fb,gb,
     & 0.5d0*(u1(i,j+1)+u1(i,j)),
     & 0.5d0*(u2(i,j+1)+u2(i,j)),
     & 0.5d0*(u3(i,j+1)+u3(i,j)),
     & 0.5d0*(u5(i,j+1)+u5(i,j)))
!
!     C face: i-1/2 , j
      call fg_vector(fc,gc,
     & 0.5d0*(u1(i-1,j)+u1(i,j)),
     & 0.5d0*(u2(i-1,j)+u2(i,j)),
     & 0.5d0*(u3(i-1,j)+u3(i,j)),
     & 0.5d0*(u5(i-1,j)+u5(i,j)))
!
!     D face: i     , j-1/2
      call fg_vector(fd,gd,
     & 0.5d0*(u1(i,j-1)+u1(i,j)),
     & 0.5d0*(u2(i,j-1)+u2(i,j)),
     & 0.5d0*(u3(i,j-1)+u3(i,j)),
     & 0.5d0*(u5(i,j-1)+u5(i,j)))
!
!
!
!
!
!     do JST stuff
!
      call jst_calcs(dt(i,j),np,params,
     & da,volume(i+1,j),volume(i,j),
     & u1(i+2,j),u1(i+1,j),u1(i,j),u1(i-1,j),
     & u2(i+2,j),u2(i+1,j),u2(i,j),u2(i-1,j),
     & u3(i+2,j),u3(i+1,j),u3(i,j),u3(i-1,j),
     & u5(i+2,j),u5(i+1,j),u5(i,j),u5(i-1,j)
     & )
!
      call jst_calcs(dt(i,j),np,params,
     & db,volume(i,j+1),volume(i,j),
     & u1(i,j+2),u1(i,j+1),u1(i,j),u1(i,j-1),
     & u2(i,j+2),u2(i,j+1),u2(i,j),u2(i,j-1),
     & u3(i,j+2),u3(i,j+1),u3(i,j),u3(i,j-1),
     & u5(i,j+2),u5(i,j+1),u5(i,j),u5(i,j-1)
     & )
!
      call jst_calcs(dt(i,j),np,params,
     & dc,volume(i,j),volume(i-1,j),
     & u1(i+1,j),u1(i,j),u1(i-1,j),u1(i-2,j),
     & u2(i+1,j),u2(i,j),u2(i-1,j),u2(i-2,j),
     & u3(i+1,j),u3(i,j),u3(i-1,j),u3(i-2,j),
     & u5(i+1,j),u5(i,j),u5(i-1,j),u5(i-2,j)
     & )
!
!
      call jst_calcs(dt(i,j),np,params,
     & dd,volume(i,j),volume(i,j-1),
     & u1(i,j+1),u1(i,j),u1(i,j-1),u1(i,j-2),
     & u2(i,j+1),u2(i,j),u2(i,j-1),u2(i,j-2),
     & u3(i,j+1),u3(i,j),u3(i,j-1),u3(i,j-2),
     & u5(i,j+1),u5(i,j),u5(i,j-1),u5(i,j-2)
     & )
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
!     calculate residual for specific cell; subtract dissipation terms
!     from original FV section
!
      r1(i,j) = (-ga(1)*xa + fa(1)*ya
     &           -gb(1)*xb + fb(1)*yb
     &           -gc(1)*xc + fc(1)*yc
     &           -gd(1)*xd + fd(1)*yd
     & - (da(1) + db(1) - dc(1) - dd(1))
     & )/area 
!
      r2(i,j) = (-ga(2)*xa + fa(2)*ya
     &           -gb(2)*xb + fb(2)*yb
     &           -gc(2)*xc + fc(2)*yc
     &           -gd(2)*xd + fd(2)*yd
     & - (da(2) + db(2) - dc(2) - dd(2))
     & )/area 
!
      r3(i,j) = (-ga(3)*xa + fa(3)*ya
     &           -gb(3)*xb + fb(3)*yb
     &           -gc(3)*xc + fc(3)*yc
     &           -gd(3)*xd + fd(3)*yd
     & - (da(3) + db(3) - dc(3) - dd(3))
     & )/area 
!
      r5(i,j) = (-ga(4)*xa + fa(4)*ya
     &           -gb(4)*xb + fb(4)*yb
     &           -gc(4)*xc + fc(4)*yc
     &           -gd(4)*xd + fd(4)*yd
     & - (da(4) + db(4) - dc(4) - dd(4))
     & )/area 
!
!
!
!
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
      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_rev(nx,ny,nl,residual,r1,r2,r3,r5)
!
!
!
      end subroutine resid
