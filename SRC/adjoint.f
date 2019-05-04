!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine adjoint(nx,ny,nl,flow,residual,np,params,dt,na,alpha,
     & Jcost,SENS)
!
      integer, external :: colour
!
      integer :: nx,ny,nl,np,na
      integer :: nh,nb,nt
      integer :: i,j,Rx,Ry,INFO
      integer :: cerr
!
      double precision :: Jcost,ans
      double precision, dimension(1,np) :: params
!
      double precision, dimension(1,na) :: 
     &         alpha,alpha_s
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,flow_s,
     &         residual,
     &         dt
!
      double precision, dimension(4*nx*ny,na) :: dRda
      double precision, dimension(4*nx*ny,4*nx*ny) :: fluxjac
      integer, dimension(4*nx*ny,4*nx*ny) :: IPIV
      double precision, dimension(1,na) :: dJda,SENS
      double precision, dimension(1,4*nx*ny) :: adj, dJdw
!
      double precision, dimension(4*(nx+2*nl)*(ny+2*nl),na) :: 
     &             dRda_h
      double precision, dimension(4*(nx+2*nl)*(ny+2*nl),
     &                            4*(nx+2*nl)*(ny+2*nl)) :: 
     &             fluxjac_h
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &             dJdw_h
!
!
!
!
!
      nb = nx*ny
      nh = 2*nl*(nx+ny)
      nt = 4*(nx+2*nl)*(ny+2*nl)
      print*, "      adjoint called from main . . ."
!
!
!
!
!
!
!
!
!
!-----------------------------------------------------------------------
!
!     1) Create functions that will give the correct partial derivatives
!
!     resid.f          aresid.f          aMresid.f          aMresid_d.f
!     cost.f           acost.f           aMcost.f           aMcost_d.f
!           ----------->      ----------->       ----------->
!          isolate w and       | merge.sh          Tapenade
!            a by hand         |
!                              `--> edit 
!                             subroutine names
!
!
!     2) Take the partial derivatives from the functions in hold arrays
!        and delete rows with boundary conditions.
!
!
!     3) Solve adjoint equation.
!
!
!     4) Form total derivative.
!
!
!
!
!
!
!
!-----------------------------------------------------------------------
!  
!     get derivatives from the differentiated cost and resdifual
!     functions. should be using
!
!             dJda
!
!             dJdw_h
!
!     which are both multiple input single output functions, and
!
!             dRda_h
!
!             fluxjac_h
!
!     which are both multiple input multiple output functions.
!
!
      print*, "        calculate partial derivatives . . ."
!
!
!
!
!
!
!     ------- get dJda -------
!
!
!
      do i = 1,na
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      alpha_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call cost_is_d(nx, ny, nl, na, np, params, flow, flow_s,
     & alpha, alpha_s, jcost, dJda(1,i))
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
!     ------- get dJdw -------
!
!
      do i = 1,nt
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      flow_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call cost_is_d(nx, ny, nl, na, np, params, flow, flow_s,
     & alpha, alpha_s, jcost, dJdw_h(1,i))
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
!     ------- get dRda -------
!
!
      do i = 1,na
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      alpha_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call aresid_d(nx, ny, nl, flow, flow_s, residual, dRda_h(:,i), 
     & np, params, dt, na, alpha, alpha_s)
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
!     ------- get dRdw -------
!
!
      do i = 1,nt
!
!
!     initialise all seeds as zero
      do j = 1,na
      alpha_s(1,j) = 0.0d0
      end do
      do j = 1,nt
      flow_s(1,j) = 0.d0
      end do
!
!
!     set relevant seed to one
      flow_s(1,i) = 1.0d0
!
!
!     call cost function to relevant part of hold array
      call aresid_d(nx, ny, nl, flow, flow_s, residual, fluxjac_h(:,i), 
     & np, params, dt, na, alpha, alpha_s)
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
!-----------------------------------------------------------------------
!
!     now delete all of the rows with halo or corner cells. note that nt
!     is defined as the storage vector including halo and corner cells,
!     while nw is just the body cells = 4*nx*ny
!
!               dJda     --->   dJda     *no action required*
!              ( 1,na)         ( 1,na)
!
!              dJdw_h    --->    adj
!              ( 1,nt)         ( 1,nw)
!
!              dRda_h    --->   dRda
!              (nt,na)         (nw,na)
!
!             fluxjac_h  --->  fluxjac
!              (nt,nt)         (nw,nw)
!
      nw = 4*nb
!
!
      Rx = 4*nh
      do i = 1,nw
      adj(1,i) = dJdw_h(1,i+Rx)
      end do
!
!
      Rx = 4*nh
      do i = 1,nw
      do j = 1,na
      dRda(i,j) = dRda_h(i+Rx,j)
      end do
      end do
!
!
      Rx = 4*nh
      Ry = 4*nh
      do j = 1,nw
      do i = 1,nw
      fluxjac(i,j) = fluxjac_h(i+Rx,j+Ry)
      end do
      end do
!
!
      print*, "        partial derivatives calculated."
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
!-----------------------------------------------------------------------
!
!     call the C function that create a PNG of the flux jacobian if the
!     line is ucommented
!
      cerr = colour(nw,fluxjac)
!
!
!
!
!
!
!
!-----------------------------------------------------------------------
!
!     solve the adjoint equation. for simplicity LAPack will be used
!     but solution of the sparse flux jacobian can be done more
!     efficiently with different software.
!
!                   fluxjac * adjoint = dJdw
!
!     for greater compatibility with the solver software dJdw has
!     already been stored in the vector call adj.
!
!
!
!
!
!
!     assign the correct values in case needed
!
      do i = 1,nw
      dJdw(1,i) = adj(1,i)
      end do
!
!
!
!
      INFO = -5
      print*, "        begin flux Jacobian factorisation . . ."
      call DGETRF(nw,nw,fluxjac,nw,IPIV,INFO)
      print*, "        matrix solved with status", 
     &  INFO, "; success = '0';"
!
!
!
      print*, "        begin system solution . . ."
      call DGETRS('T',nw,1,fluxjac,nw,IPIV,adj,nw,INFO)
      print*, "        system solved with status", 
     &  INFO, "; success = '0';"
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
!-----------------------------------------------------------------------
!
!     finally, form the total derivative by combining the remaining
!     partial derivatives.
!
!                 SENS = dJda + adjoint * dRda
!
!
!
      do i = 1,na
!
      ans = 0.0d0
      do j = 1,nw
      ans = ans + adj(1,j)*dRda(j,i)
      end do
!
      SENS(1,i) = dJda(1,i) - ans
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
!
!
!
!-----------------------------------------------------------------------
!
      print*, "      adjoint completed in main."
!
      end subroutine adjoint
!-----------------------------------------------------------------------
