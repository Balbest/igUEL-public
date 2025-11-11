      subroutine Force_Stiff_3D(NNODE,ntens,ndim,ndof,stress,gstress,dV,
     1                          dR,ddR,C_St_Matr,A_St_Matr,force,stiff)
c
       include 'ABA_PARAM.INC'
c
      integer j, jdof, jcol, inc_col, i, idof, irow, icol, inc_row
c
      double precision stiff, stiff_n, force, force_n, stress, gstress,
     *                 C_St_Matr, A_St_Matr
      dimension stiff(ndof*NNODE,ndof*NNODE),
     *          stiff_n(ndof,ndof),
     *          force(ndof*NNODE),
     *          force_n(ndof),
     *          stress(ntens),
     *          gstress(ntens*ndim),
     *          C_St_Matr(ntens,ntens),
     *          A_St_Matr(ntens*ndim,ntens*ndim)
c
      double precision Bj, Lj, Bit, Lit
      dimension Bj(6,ndim),
     *          Lj(18,ndim),
     *          Bit(ndim,6),
     *          Lit(ndim,18)
c
      double precision dR, ddR, dV,
     *                 dRjdx, dRjdy, dRjdz,
     *                 d2Rjdx2, d2Rjdy2, d2Rjdz2,
     *                 d2Rjdydz, d2Rjdxdz, d2Rjdxdy,
     *                 dRidx, dRidy, dRidz,
     *                 d2Ridx2, d2Ridy2, d2Ridz2,
     *                 d2Ridydz, d2Ridxdz, d2Ridxdy
      dimension dR(NNODE,3), ddR(NNODE,6)
c
      force  = 0.d0 ! --- Contribution to the vector of internal forces
      stiff  = 0.d0 ! --- Contribution to the stiffness matrix
c
      Bj  = 0.d0 ! --- Strain matrix
      Lj  = 0.d0 ! --- Strain gradient matrix
      BiT = 0.d0 ! --- Transposed strain matrix
      LiT = 0.d0 ! --- Transposed strain gradient matrix
c
      do j = 1, NNODE
c
c !--- first derivatives 1-dx, 2-dy, 3-dz
c !--- second derivatives 1-dx^2, 2-dy^2, 3-dz^2, 4-dydz, 5-dxdz, 6-dxdy
c
         dRjdx = dR(j,1)
         dRjdy = dR(j,2)
         dRjdz = dR(j,3)
c
         d2Rjdx2  = ddR(j,1)
         d2Rjdy2  = ddR(j,2)
         d2Rjdz2  = ddR(j,3)
         d2Rjdydz = ddR(j,4)
         d2Rjdxdz = ddR(j,5)
         d2Rjdxdy = ddR(j,6)
c
         Bj(1,1) = dRjdx
         Bj(2,2) = dRjdy
         Bj(3,3) = dRjdz
         Bj(4,2) = dRjdz
         Bj(4,3) = dRjdy
         Bj(5,1) = dRjdz
         Bj(5,3) = dRjdx
         Bj(6,1) = dRjdy
         Bj(6,2) = dRjdx
c
         Lj(1,1)  = d2Rjdx2
         Lj(2,1)  = d2Rjdxdy
         Lj(3,1)  = d2Rjdxdz
         Lj(4,2)  = d2Rjdxdy
         Lj(5,2)  = d2Rjdy2
         Lj(6,2)  = d2Rjdydz
         Lj(7,3)  = d2Rjdxdz
         Lj(8,3)  = d2Rjdydz
         Lj(9,3)  = d2Rjdz2
         Lj(10,2) = d2Rjdxdz
         Lj(11,2) = d2Rjdydz
         Lj(12,2) = d2Rjdz2
         Lj(10,3) = d2Rjdxdy
         Lj(11,3) = d2Rjdy2
         Lj(12,3) = d2Rjdydz
         Lj(13,1) = d2Rjdxdz
         Lj(14,1) = d2Rjdydz
         Lj(15,1) = d2Rjdz2
         Lj(13,3) = d2Rjdx2
         Lj(14,3) = d2Rjdxdy
         Lj(15,3) = d2Rjdxdz
         Lj(16,1) = d2Rjdxdy
         Lj(17,1) = d2Rjdy2
         Lj(18,1) = d2Rjdydz
         Lj(16,2) = d2Rjdx2
         Lj(17,2) = d2Rjdxdy
         Lj(18,2) = d2Rjdxdz
c
         force_n(1:3) = matmul(transpose(Bj),stress)
     *                + matmul(transpose(Lj),gstress)
c
         inc_col = (j - 1)*ndof
c
         do jdof = 1, ndof
            jcol = jdof + inc_col
            force(jcol) = force(jcol) + force_n(jdof)*dV
         end do
c
         do i = 1, NNODE
c
C !--- first derivatives 1-dx, 2-dy, 3-dz
C !--- second derivatives 1-dx^2, 2-dy^2, 3-dz^2, 4-dydz, 5-dxdz, 6-dxdy
c
            dRidx = dR(i,1)
            dRidy = dR(i,2)
            dRidz = dR(i,3)
c
            d2Ridx2  = ddR(i,1)
            d2Ridy2  = ddR(i,2)
            d2Ridz2  = ddR(i,3)
            d2Ridydz = ddR(i,4)
            d2Ridxdz = ddR(i,5)
            d2Ridxdy = ddR(i,6)
c
            Bit(1,1) = dRidx
            Bit(2,2) = dRidy
            Bit(3,3) = dRidz
            Bit(2,4) = dRidz
            Bit(3,4) = dRidy
            Bit(1,5) = dRidz
            Bit(3,5) = dRidx
            Bit(1,6) = dRidy
            Bit(2,6) = dRidx
c
            Lit(1,1)  = d2Ridx2
            Lit(1,2)  = d2Ridxdy
            Lit(1,3)  = d2Ridxdz
            Lit(2,4)  = d2Ridxdy
            Lit(2,5)  = d2Ridy2
            Lit(2,6)  = d2Ridydz
            Lit(3,7)  = d2Ridxdz
            Lit(3,8)  = d2Ridydz
            Lit(3,9)  = d2Ridz2
            Lit(2,10) = d2Ridxdz
            Lit(2,11) = d2Ridydz
            Lit(2,12) = d2Ridz2
            Lit(3,10) = d2Ridxdy
            Lit(3,11) = d2Ridy2
            Lit(3,12) = d2Ridydz
            Lit(1,13) = d2Ridxdz
            Lit(1,14) = d2Ridydz
            Lit(1,15) = d2Ridz2
            Lit(3,13) = d2Ridx2
            Lit(3,14) = d2Ridxdy
            Lit(3,15) = d2Ridxdz
            Lit(1,16) = d2Ridxdy
            Lit(1,17) = d2Ridy2
            Lit(1,18) = d2Ridydz
            Lit(2,16) = d2Ridx2
            Lit(2,17) = d2Ridxdy
            Lit(2,18) = d2Ridxdz
c
            stiff_n(1:3,1:3) = matmul(Bit, matmul(C_St_Matr,Bj))
     *                       + matmul(Lit, matmul(A_St_Matr,Lj))
c
            inc_row = (i -1)*ndof
c
            do jdof = 1, ndof
               icol = jdof + inc_col
               do idof = 1, ndof
                  irow = idof + inc_row
c
                  stiff(irow,icol) = stiff(irow,icol)
     *                             + stiff_n(idof,jdof)*dV
               end do
            end do
         end do
      end do
c
      return
      end subroutine Force_Stiff_3D