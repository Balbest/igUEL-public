      subroutine Force_Stiff_2D(NNODE,ntens,ndim,ndof,stress,gstress,dV,
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
      double precision Bj, Lj, BiT, LiT
      dimension Bj(4,2),
     *          Lj(8,2),
     *          BiT(2,4),
     *          LiT(2,8)
c
      double precision dR, ddR, dV,
     *                 dRjdx, dRjdy, dRidx, dRidy,
     *                 d2Rjdx2, d2Rjdy2, d2Rjdxdy,
     *                 d2Ridx2, d2Ridy2, d2Ridxdy
      dimension dR(NNODE,2), ddR(NNODE,3)
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
         inc_col = (j - 1)*ndof
c
         dRjdx = dR(j,1)
         dRjdy = dR(j,2)
c
         d2Rjdx2  = ddR(j,1)
         d2Rjdy2  = ddR(j,2)
         d2Rjdxdy = ddR(j,3)
c
         Bj(1,1) = dRjdx
         Bj(2,2) = dRjdy
         Bj(4,1) = dRjdy
         Bj(4,2) = dRjdx
c
         Lj(1,1)  = d2Rjdx2
         Lj(2,1)  = d2Rjdxdy
         Lj(3,2)  = d2Rjdxdy
         Lj(4,2)  = d2Rjdy2
         Lj(7,1)  = d2Rjdxdy
         Lj(7,2)  = d2Rjdx2
         Lj(8,1)  = d2Rjdy2
         Lj(8,2)  = d2Rjdxdy
c
         force_n(1:2) = matmul(transpose(Bj),stress)
     *                + matmul(transpose(Lj),gstress)
c
         do jdof = 1, ndof
            jcol = jdof + inc_col
            force(jcol) = force(jcol) + force_n(jdof)*dV
         end do
c
         do i = 1, NNODE
            inc_row = (i -1)*ndof
c
            dRidx = dR(i,1)
            dRidy = dR(i,2)
c
            d2Ridx2  = ddR(i,1)
            d2Ridy2  = ddR(i,2)
            d2Ridxdy = ddR(i,3)
c
            BiT(1,1) = dRidx
            BiT(2,2) = dRidy
            BiT(1,4) = dRidy
            BiT(2,4) = dRidx
c           
            LiT(1,1)  = d2Ridx2
            LiT(1,2)  = d2Ridxdy
            LiT(2,3)  = d2Ridxdy
            LiT(2,4)  = d2Ridy2
            LiT(1,7)  = d2Ridxdy
            LiT(2,7)  = d2Ridx2
            LiT(1,8)  = d2Ridy2
            LiT(2,8)  = d2Ridxdy
c
            stiff_n(1:2,1:2) = matmul(Bit, matmul(C_St_Matr,Bj))
     *                       + matmul(Lit, matmul(A_St_Matr,Lj))
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
      end subroutine Force_Stiff_2D