c ------------------------------------------------------------------------------------
c ! --- This supplementary subroutine is part of the main subroutine 'UEL_IGA_2D_SGE'
c ------------------------------------------------------------------------------------
      subroutine Strains_2D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                      stran,gstran,dR,ddR)
c
c ------------------------------------------------------------------------------------
c !--- Components of strain tensor (vector) and its gradient
c
c !    stran(1) --> E_xx = d(u_x)/dx
c !    stran(2) --> E_yy = d(u_y)/dy
c !    stran(3) --> E_zz = 0
c !    stran(4) --> E_xy = d(u_x)/dy + d(u_y)/dx
c
c !   gstran(1) --> d(E_xx)/dx
c !   gstran(2) --> d(E_xx)/dy
c !   gstran(3) --> d(E_yy)/dx
c !   gstran(4) --> d(E_yy)/dy
c !   gstran(5) --> d(E_zz)/dx = 0
c !   gstran(6) --> d(E_zz)/dy = 0
c !   gstran(7) --> d(E_xy)/dx
c !   gstran(8) --> d(E_xy)/dy
c ------------------------------------------------------------------------------------
c
      include 'ABA_PARAM.INC'
c
      integer i, j, k
c
      double precision stran, gstran, u_n
      dimension stran(ntens),
     *          gstran(ndim*ntens),
     *          u_n(ndim),
     *          U(NDOFEL)
c
      double precision dR, ddR, dRdx, dRdy, d2Rdx2, d2Rdy2, d2Rdxdy
      dimension dR(NNODE,2), ddR(NNODE,3)
c
      stran(:)  = 0.0D+00
      gstran(:) = 0.0D+00
c
      do j = 1, NNODE
c
         dRdx = dR(j,1)
         dRdy = dR(j,2)
c
         d2Rdx2  = ddR(j,1)
         d2Rdy2  = ddR(j,2)
         d2Rdxdy = ddR(j,3)
c
         k = (j - 1)*ndof
c
         do i = 1, ndof
            u_n(i) = U(i + k)
         end do
c
c !--- Strain components
c
       stran(1) = stran(1) + dRdx*u_n(1)
       stran(2) = stran(2) + dRdy*u_n(2)
       stran(4) = stran(4) + dRdy*u_n(1) + dRdx*u_n(2)
c
c !--- Strain gradient components
c
       gstran(1) = gstran(1) + d2Rdx2*u_n(1)
       gstran(2) = gstran(2) + d2Rdxdy*u_n(1)
c
       gstran(3) = gstran(3) + d2Rdxdy*u_n(2)
       gstran(4) = gstran(4) + d2Rdy2*u_n(2)
c
       gstran(7) = gstran(7) + d2Rdxdy*u_n(1) + d2Rdx2*u_n(2)
       gstran(8) = gstran(8) + d2Rdy2*u_n(1) + d2Rdxdy*u_n(2)
c
      end do
c
      return
      end subroutine Strains_2D