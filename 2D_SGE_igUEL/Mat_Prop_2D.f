c ------------------------------------------------------------------------------------
c ! --- This supplementary subroutine is part of the main subroutine 'UEL_IGA_2D_SGE'
c ------------------------------------------------------------------------------------
      subroutine Mat_Prop_2D(ntens,ndim,E_mod,nu_mod,g_par,C_St_Matr,A_St_Matr)
c
      include 'ABA_PARAM.INC'
c
      double precision ZERO, ONE, TWO
      parameter (ZERO = 0.D+00, ONE = 1.D+00, TWO = 2.D+00)
c
      integer i, j, k, ii, jj
c
      double precision E_mod, nu_mod, g_par, L_Lam, L_Mu,
     *                 C_St_Matr, A_St_Matr
      dimension C_St_Matr(ntens,ntens), A_St_Matr(ntens*ndim,ntens*ndim)
c
c !--- Plane strain
c
      L_Lam=nu_mod*E_mod/(ONE+nu_mod)/(ONE-TWO*nu_mod)
      L_Mu=E_mod/TWO/(ONE+nu_mod)
c
      C_St_Matr = ZERO
      A_St_Matr = ZERO
c
      C_St_Matr(1,1) = L_Lam + TWO*L_Mu
      C_St_Matr(1,2) = L_Lam
      C_St_Matr(1,3) = L_Lam
      C_St_Matr(2,1) = C_St_Matr(1,2)
      C_St_Matr(2,2) = L_Lam + TWO*L_Mu
      C_St_Matr(2,3) = L_Lam
      C_St_Matr(3,1) = C_St_Matr(1,3)
      C_St_Matr(3,2) = C_St_Matr(2,3)
      C_St_Matr(3,3) = L_Lam + TWO*L_Mu
      C_St_Matr(4,4) = L_Mu
c
      do i = 1,ntens
         ii = 0 + (i-1)*ndim
         do j = 1,ntens
            jj = 0 + (j-1)*ndim
            do k = 1,ndim
               A_St_Matr(ii+k,jj+k) = g_par**TWO*C_St_Matr(i,j)
            end do
         end do
      end do
c
      return
      end subroutine Mat_Prop_2D