c ------------------------------------------------------------------------------------
c ! --- This supplementary subroutine is part of the main subroutine 'UEL_IGA_3D_SGE'
c ------------------------------------------------------------------------------------
      subroutine Strains_3D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                      stran,gstran,dR,ddR)
c
c  -----------------------------------------------------------------------------------
c !--- Components of strain tensor (vector) and its gradient
c  
c !    stran(1)  --> E_xx = d(u_x)/dx
c !    stran(2)  --> E_yy = d(u_y)/dy
c !    stran(3)  --> E_zz = d(u_z)/dz
c !    stran(4)  --> E_yz = d(u_y)/dz + d(u_z)/dy !
c !    stran(5)  --> E_xz = d(u_x)/dz + d(u_z)/dx !
c !    stran(6)  --> E_xy = d(u_x)/dy + d(u_y)/dx !
c  
c !   gstran(1)  --> d(E_xx)/dx
c !   gstran(2)  --> d(E_xx)/dy
c !   gstran(3)  --> d(E_xx)/dz
c !   gstran(4)  --> d(E_yy)/dx
c !   gstran(5)  --> d(E_yy)/dy
c !   gstran(6)  --> d(E_yy)/dz
c !   gstran(7)  --> d(E_zz)/dx
c !   gstran(8)  --> d(E_zz)/dy
c !   gstran(9)  --> d(E_zz)/dz
c !   gstran(10) --> d(E_yz)/dx
c !   gstran(11) --> d(E_yz)/dy
c !   gstran(12) --> d(E_yz)/dz
c !   gstran(13) --> d(E_xz)/dx
c !   gstran(14) --> d(E_xz)/dy
c !   gstran(15) --> d(E_xz)/dz
c !   gstran(16) --> d(E_xy)/dx
c !   gstran(17) --> d(E_xy)/dy
c !   gstran(18) --> d(E_xy)/dz
c -----------------------------------------------------------------------------------
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
      double precision dR, ddR,
     *                 dRdx, dRdy, dRdz,
     *                 d2Rdx2, d2Rdy2, d2Rdz2, d2Rdydz, d2Rdxdz, d2Rdxdy 
      dimension dR(NNODE,3), ddR(NNODE,6)
c
      stran(:)  = 0.0D+00
      gstran(:) = 0.0D+00
c
      do j = 1, NNODE
c
         dRdx = dR(j,1)
         dRdy = dR(j,2)
         dRdz = dR(j,3)
         d2Rdx2  = ddR(j,1)
         d2Rdy2  = ddR(j,2)
         d2Rdz2  = ddR(j,3)
         d2Rdydz = ddR(j,4)
         d2Rdxdz = ddR(j,5)
         d2Rdxdy = ddR(j,6)
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
       stran(3) = stran(3) + dRdz*u_n(3)
       stran(4) = stran(4) + dRdz*u_n(2) + dRdy*u_n(3)
       stran(5) = stran(5) + dRdz*u_n(1) + dRdx*u_n(3)
       stran(6) = stran(6) + dRdy*u_n(1) + dRdx*u_n(2)
c
c !--- Strain gradient components
c
       gstran(1) = gstran(1) + d2Rdx2*u_n(1)
       gstran(2) = gstran(2) + d2Rdxdy*u_n(1)
       gstran(3) = gstran(3) + d2Rdxdz*u_n(1)
c
       gstran(4) = gstran(4) + d2Rdxdy*u_n(2)
       gstran(5) = gstran(5) + d2Rdy2*u_n(2)
       gstran(6) = gstran(6) + d2Rdydz*u_n(2)
c
       gstran(7) = gstran(7) + d2Rdxdz*u_n(3)
       gstran(8) = gstran(8) + d2Rdydz*u_n(3)
       gstran(9) = gstran(9) + d2Rdz2*u_n(3)
c
       gstran(10) = gstran(10) + d2Rdxdz*u_n(2) + d2Rdxdy*u_n(3)
       gstran(11) = gstran(11) + d2Rdydz*u_n(2) + d2Rdy2*u_n(3)
       gstran(12) = gstran(12) + d2Rdz2*u_n(2) + d2Rdydz*u_n(3)
c
       gstran(13) = gstran(13) + d2Rdxdz*u_n(1) + d2Rdx2*u_n(3)
       gstran(14) = gstran(14) + d2Rdydz*u_n(1) + d2Rdxdy*u_n(3)
       gstran(15) = gstran(15) + d2Rdz2*u_n(1) + d2Rdxdz*u_n(3)
c
       gstran(16) = gstran(16) + d2Rdxdy*u_n(1) + d2Rdx2*u_n(2)
       gstran(17) = gstran(17) + d2Rdy2*u_n(1) + d2Rdxdy*u_n(2)
       gstran(18) = gstran(18) + d2Rdydz*u_n(1) + d2Rdxdz*u_n(2)
c
      end do
c
      return
      end subroutine Strains_3D