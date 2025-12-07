c -------------------------------------------------------------------------------------
c ! --- This supplementary subroutine is part of the main subroutine 'UEL_IGA_3D_SGE'
c -------------------------------------------------------------------------------------
c ! --- It calculates the 1st and 2nd derivatives of NURBS basis functions in uvw-space
c -------------------------------------------------------------------------------------
      module module_inv
      implicit none
      contains
c ! - Auxiliary function performing a direct calculation of the inverse of a 3×3 matrix. 
      function matinv3(A) result(B)
	  
         real(8), intent(in) :: A(3,3)   !! Matrix
         real(8)             :: B(3,3)   !! Inverse matrix
         real(8)             :: detA
c
c ! ---  Calculate the determinant of the matrix
         detA = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     1          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     2          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
	 
c ! ---  Calculate the inverse of the matrix
         B(1,1) = + 1/detA * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
         B(2,1) = - 1/detA * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
         B(3,1) = + 1/detA * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
         B(1,2) = - 1/detA * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
         B(2,2) = + 1/detA * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
         B(3,2) = - 1/detA * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
         B(1,3) = + 1/detA * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
         B(2,3) = - 1/detA * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
         B(3,3) = + 1/detA * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      end function matinv3
c
c ! - Auxiliary subroutine for n×n matrix inverse using Crout's method
      subroutine Crout_Inv(M_inp, M_inv, M_size)
c
      implicit none
c
      double precision eps
      parameter (eps = 1.d-20)
c
      integer M_size, parity, sing_ind, pivot_index
      dimension  pivot_index(M_size)
c
      double precision M_inp, M_inv, work_vector
      dimension M_inp(M_size,M_size), M_inv(M_size,M_size),
     *          work_vector(M_size)
c
      double precision row_max, temp_swap, abs_value, accum,
     *                 pivot_inv, row_scale
      dimension        row_scale(M_size)
c
      integer row_idx, col_idx, in_idx, pivot_row, nonzero_start, perm_row
c
c ! - Initialization
c
      parity = 1
      sing_ind = 0
c
c !-- Find largest element in each row for scaled pivoting
c
      do row_idx = 1,M_size
         row_max = 0.d0
         do col_idx = 1,M_size
            abs_value = dabs(M_inp(row_idx,col_idx))
            if (abs_value .GT. row_max) row_max = abs_value
         end do
         if(row_max .LT. eps) then
            sing_ind = 1
            return
         end if
         row_scale(row_idx) = 1.d0/row_max
      end do
c
c !-- Crout's method: LU decomposition with partial pivoting
c
      do col_idx = 1,M_size
c
c !---  Calculate U (upper triangular) elements above diagonal
         do row_idx = 1,col_idx-1
            accum = M_inp(row_idx,col_idx)
            do in_idx = 1,row_idx-1
               accum = accum
     *               - M_inp(row_idx,in_idx)*M_inp(in_idx,col_idx)
            end do
            M_inp(row_idx,col_idx) = accum
         end do
c
c !---   Compute L (lower triangular) elements and find best pivot
         row_max = 0.d0
         pivot_row = col_idx
         do row_idx = col_idx,M_size
            accum = M_inp(row_idx,col_idx)
            do in_idx = 1,col_idx-1
               accum = accum
     *               - M_inp(row_idx,in_idx)*M_inp(in_idx,col_idx)
            end do
            M_inp(row_idx,col_idx) = accum
c
c !---      Calculate scaled pivot quality for numerical stability
            abs_value = dabs(accum)
            temp_swap = row_scale(row_idx)*abs_value
            if (temp_swap .GT. row_max) then
                pivot_row = row_idx
                row_max = temp_swap
            end if
         end do
c
c !---   Swap rows if better pivot found
         if (col_idx .NE. pivot_row) then
            do in_idx = 1,M_size
               temp_swap = M_inp(pivot_row,in_idx)
               M_inp(pivot_row,in_idx) = M_inp(col_idx,in_idx)
               M_inp(col_idx,in_idx) = temp_swap
            end do
            parity = -parity
            row_scale(pivot_row) = row_scale(col_idx)
         end if
c
         pivot_index(col_idx) = pivot_row
c
c !---   Protect against division by zero
         abs_value = dabs(M_inp(col_idx,col_idx))
         if (abs_value .LT. eps) M_inp(col_idx,col_idx) = eps
c
c !---   Divide by pivot to complete L matrix
         if (col_idx .NE. M_size) then
             pivot_inv = 1.d0/M_inp(col_idx,col_idx)
             do row_idx = col_idx+1,M_size
                M_inp(row_idx,col_idx) = 
     *          M_inp(row_idx,col_idx)*pivot_inv
             end do
         end if 
      end do
c
c !-- Compute M_inv
c !-- Solve M_inp * M_inv = Identity by solving for each column
c
      do col_idx = 1,M_size
c
c !---   Initialize right-hand side as k-th unit vector
         do row_idx = 1,M_size
            work_vector(row_idx) = 0.d0
         end do
         work_vector(col_idx) = 1.d0
c
c !---   Forward substitution: Solve L*Y = e_k
         nonzero_start = 0
         do row_idx = 1,M_size
            perm_row = pivot_index(row_idx)
            accum = work_vector(perm_row)
            work_vector(perm_row) = work_vector(row_idx)
            if (nonzero_start .NE. 0) then
                do in_idx = nonzero_start,row_idx-1
                   accum = accum
     *                   - M_inp(row_idx,in_idx)*work_vector(in_idx)
                end do
            else if (accum .NE. 0.d0) then
                nonzero_start = row_idx
            end if
            work_vector(row_idx) = accum
         end do
c
c !---   Back substitution: Solve U*X = Y
         do row_idx = M_size,1,-1
            accum = work_vector(row_idx)
            if (row_idx .LT. M_size) then
                do in_idx = row_idx+1,M_size
                   accum = accum
     *                   - M_inp(row_idx,in_idx)*work_vector(in_idx)
                end do
            end if
            work_vector(row_idx) = accum / M_inp(row_idx,row_idx)
         end do
c
         M_inv(:,col_idx) = work_vector
      end do
c
      return
      end subroutine Crout_Inv
c
      end module module_inv
c
c!----------------------MAIN SUBROUTINE STARTS HERE----------------------------
      subroutine NURBS_BF_3D(p,i_el,uu,U_kn,Kn_Num_U,
     1                       q,j_el,vv,V_kn,Kn_Num_V, 
     2                       r,k_el,ww,W_kn,Kn_Num_W, 
     3                       w_CPs,NOCPs,NNODE,COORDS,MCRD,
     4                       RR,dR_dx,ddR_ddx,Djac,JELEM)
c
      use module_inv
c
      include 'ABA_PARAM.INC'
c
      integer NOCPs, Kn_Num_U, Kn_Num_V, Kn_Num_W,
     *        p, i_el, q, j_el, r, k_el
c
      double precision uu, vv, ww, COORDS, w_CPs, U_kn, V_kn, W_kn
      dimension COORDS(MCRD,NNODE), w_CPs(NOCPs),
     *          U_kn(Kn_Num_U), V_kn(Kn_Num_V), W_kn(Kn_Num_W)
c
      double precision RR, dR_dx, ddR_ddx
      dimension RR(NNODE), dR_dx(NNODE,3), ddR_ddx(NNODE,6)

      double precision WS,dWS(3),ddWS(6),NN(3,p+1),MM(3,q+1),LL(3,r+1),
     *	               dx_dxi(3,3),dxi_dx(3,3),
     *	               ddx_ddxi(3,6),dx_dxi_2T(6,6),dxi_dx_2T(6,6),
     *	               dR(NNODE,3), ddR(NNODE,6), AuxV(6), Djac
c
      integer k, No_w, b, c, d, a1, a2
c
      call deriv2(i_el,p,uu,U_kn,Kn_Num_U,NN) ! basis functions in uu
      call deriv2(j_el,q,vv,V_kn,Kn_Num_V,MM) ! basis functions in vv
      call deriv2(k_el,r,ww,W_kn,Kn_Num_W,LL) ! basis functions in ww

      RR = 0.0d0
      dR = 0.0d0
      ddR = 0.0d0
      dR_dx = 0.0d0
      ddR_ddx = 0.0d0
c
      dxi_dx = 0.0d0
      dx_dxi = 0.0d0
      ddx_ddxi = 0.0d0
c
      WS = 0.0d0
      dWS = 0.0d0
      ddWS = 0.0d0
c
      k = 0
      do d = 0,r
         do c = 0,q
            do b = 0,p
               k = k+1
			   No_w  = (Kn_Num_U-p-1)*(j_el-q+c-1) + i_el-p+b +
     *			       (k_el-r+d-1)*(Kn_Num_U-p-1)*(Kn_Num_V-q-1)   !works only for vectors w/o inner repeats
               RR(k) = NN(1,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               WS = WS + RR(k)
c
c ! ---        First parametric derivatives 1-du, 2-dv, 3-dw
               dR(k,1) = NN(2,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               dWS(1) = dWS(1) + dR(k,1)
               dR(k,2) = NN(1,b+1)*MM(2,c+1)*LL(1,d+1)*w_CPs(No_w)
               dWS(2) = dWS(2) + dR(k,2)
               dR(k,3) = NN(1,b+1)*MM(1,c+1)*LL(2,d+1)*w_CPs(No_w)
               dWS(3) = dWS(3) + dR(k,3)
c
c ! ---        Second parametric derivatives 1-dudu, 2-dvdv, 3-dwdw, 4-dvdw, 5-dudw, 6-dudv
               ddR(k,1) = NN(3,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddWS(1) = ddWS(1) + ddR(k,1)
               ddR(k,2) = NN(1,b+1)*MM(3,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddWS(2) = ddWS(2) + ddR(k,2)
               ddR(k,3) = NN(1,b+1)*MM(1,c+1)*LL(3,d+1)*w_CPs(No_w)
               ddWS(3) = ddWS(3) + ddR(k,3)
               ddR(k,4) = NN(1,b+1)*MM(2,c+1)*LL(2,d+1)*w_CPs(No_w)
               ddWS(4) = ddWS(4) + ddR(k,4)
               ddR(k,5) = NN(2,b+1)*MM(1,c+1)*LL(2,d+1)*w_CPs(No_w)
               ddWS(5) = ddWS(5) + ddR(k,5)
               ddR(k,6) = NN(2,b+1)*MM(2,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddWS(6) = ddWS(6) + ddR(k,6)
            end do
         end do
      end do
c
c ! --- Divide by weight sum WS
      do k=1,NNODE
c
           ddR(k,1) = ddR(k,1)/WS - 2*dR(k,1)*dWS(1)/WS**2
     *     -RR(k)*ddWS(1)/WS**2 + 2*RR(k)*dWS(1)**2/WS**3
c
           ddR(k,2) = ddR(k,2)/WS - 2*dR(k,2)*dWS(2)/WS**2
     *     -RR(k)*ddWS(2)/WS**2 + 2*RR(k)*dWS(2)**2/WS**3
c
           ddR(k,3) = ddR(k,3)/WS - 2*dR(k,3)*dWS(3)/WS**2
     *     -RR(k)*ddWS(3)/WS**2 + 2*RR(k)*dWS(3)**2/WS**3
c
           ddR(k,4) = ddR(k,4)/WS - dR(k,2)*dWS(3)/WS**2 -
     *     dR(k,3)*dWS(2)/WS**2
     *     -RR(k)*ddWS(4)/WS**2 + 2*RR(k)*dWS(2)*dWS(3)/WS**3
c
           ddR(k,5) = ddR(k,5)/WS - dR(k,1)*dWS(3)/WS**2 -
     *     dR(k,3)*dWS(1)/WS**2
     *     -RR(k)*ddWS(5)/WS**2 + 2*RR(k)*dWS(1)*dWS(3)/WS**3
c
           ddR(k,6) = ddR(k,6)/WS - dR(k,1)*dWS(2)/WS**2 -
     *     dR(k,2)*dWS(1)/WS**2
     *     -RR(k)*ddWS(6)/WS**2 + 2*RR(k)*dWS(1)*dWS(2)/WS**3
c
           dR(k,1) = dR(k,1)/WS - RR(k)*dWS(1)/WS**2
           dR(k,2) = dR(k,2)/WS - RR(k)*dWS(2)/WS**2
           dR(k,3) = dR(k,3)/WS - RR(k)*dWS(3)/WS**2
c
           RR(k) = RR(k)/WS
      end do
c
c ! --- Jacobian and its derivatives from parameter space to physical and inversion
      k = 0
      do d = 0,r
         do c = 0,q
            do b = 0,p
               k = k+1
			   do a1 = 1,3
			      do a2 = 1,3
				     dx_dxi(a1,a2) = dx_dxi(a1,a2) + COORDS(a1,k)*dR(k,a2)
					 ddx_ddxi(a1,a2) = ddx_ddxi(a1,a2) + COORDS(a1,k)*ddR(k,a2)
					 ddx_ddxi(a1,a2+3) = ddx_ddxi(a1,a2+3) + COORDS(a1,k)*ddR(k,a2+3)
			      end do
			   end do
            end do
         end do
      end do
c ! --- Calculate the determinant of the matrix
      Djac = (dx_dxi(1,1)*dx_dxi(2,2)*dx_dxi(3,3)
     *     -  dx_dxi(1,1)*dx_dxi(2,3)*dx_dxi(3,2)
     *     -  dx_dxi(1,2)*dx_dxi(2,1)*dx_dxi(3,3)
     *     +  dx_dxi(1,2)*dx_dxi(2,3)*dx_dxi(3,1)
     *     +  dx_dxi(1,3)*dx_dxi(2,1)*dx_dxi(3,2)
     *     -  dx_dxi(1,3)*dx_dxi(2,2)*dx_dxi(3,1))
c
      dxi_dx = matinv3(dx_dxi)
c
c ! --- Basis functions and derivatives w.r.t. physical coordinates
      do k = 1,NNODE
	     do a1 = 1,3
	        do a2 = 1,3
		       dR_dx(k,a1) = dR_dx(k,a1) + dxi_dx(a2,a1)*dR(k,a2)
            end do
         end do
      end do
c
      dx_dxi_2T(1,1) = dx_dxi(1,1)**2
      dx_dxi_2T(1,2) = dx_dxi(2,1)**2
      dx_dxi_2T(1,3) = dx_dxi(3,1)**2
      dx_dxi_2T(1,4) = 2*dx_dxi(2,1)*dx_dxi(3,1)
      dx_dxi_2T(1,5) = 2*dx_dxi(1,1)*dx_dxi(3,1)
      dx_dxi_2T(1,6) = 2*dx_dxi(1,1)*dx_dxi(2,1)
c
      dx_dxi_2T(2,1) = dx_dxi(1,2)**2
      dx_dxi_2T(2,2) = dx_dxi(2,2)**2
      dx_dxi_2T(2,3) = dx_dxi(3,2)**2
      dx_dxi_2T(2,4) = 2*dx_dxi(2,2)*dx_dxi(3,2)
      dx_dxi_2T(2,5) = 2*dx_dxi(1,2)*dx_dxi(3,2)
      dx_dxi_2T(2,6) = 2*dx_dxi(1,2)*dx_dxi(2,2)
c
      dx_dxi_2T(3,1) = dx_dxi(1,3)**2
      dx_dxi_2T(3,2) = dx_dxi(2,3)**2
      dx_dxi_2T(3,3) = dx_dxi(3,3)**2
      dx_dxi_2T(3,4) = 2*dx_dxi(2,3)*dx_dxi(3,3)
      dx_dxi_2T(3,5) = 2*dx_dxi(1,3)*dx_dxi(3,3)
      dx_dxi_2T(3,6) = 2*dx_dxi(1,3)*dx_dxi(2,3)
c
      dx_dxi_2T(4,1) = dx_dxi(1,2)*dx_dxi(1,3)
      dx_dxi_2T(4,2) = dx_dxi(2,2)*dx_dxi(2,3)
      dx_dxi_2T(4,3) = dx_dxi(3,2)*dx_dxi(3,3)
      dx_dxi_2T(4,4) = dx_dxi(2,2)*dx_dxi(3,3) + dx_dxi(2,3)*dx_dxi(3,2)
      dx_dxi_2T(4,5) = dx_dxi(1,2)*dx_dxi(3,3) + dx_dxi(1,3)*dx_dxi(3,2)
      dx_dxi_2T(4,6) = dx_dxi(1,2)*dx_dxi(2,3) + dx_dxi(1,3)*dx_dxi(2,2)
c
      dx_dxi_2T(5,1) = dx_dxi(1,1)*dx_dxi(1,3)
      dx_dxi_2T(5,2) = dx_dxi(2,1)*dx_dxi(2,3)
      dx_dxi_2T(5,3) = dx_dxi(3,1)*dx_dxi(3,3)
      dx_dxi_2T(5,4) = dx_dxi(2,1)*dx_dxi(3,3) + dx_dxi(2,3)*dx_dxi(3,1)
      dx_dxi_2T(5,5) = dx_dxi(1,1)*dx_dxi(3,3) + dx_dxi(1,3)*dx_dxi(3,1)
      dx_dxi_2T(5,6) = dx_dxi(1,1)*dx_dxi(2,3) + dx_dxi(1,3)*dx_dxi(2,1)
c
      dx_dxi_2T(6,1) = dx_dxi(1,1)*dx_dxi(1,2)
      dx_dxi_2T(6,2) = dx_dxi(2,1)*dx_dxi(2,2)
      dx_dxi_2T(6,3) = dx_dxi(3,1)*dx_dxi(3,2)
      dx_dxi_2T(6,4) = dx_dxi(2,1)*dx_dxi(3,2) + dx_dxi(2,2)*dx_dxi(3,1)
      dx_dxi_2T(6,5) = dx_dxi(1,1)*dx_dxi(3,2) + dx_dxi(1,2)*dx_dxi(3,1)
      dx_dxi_2T(6,6) = dx_dxi(1,1)*dx_dxi(2,2) + dx_dxi(1,2)*dx_dxi(2,1)
c
      call Crout_Inv(dx_dxi_2T,dxi_dx_2T,6)
c
      do k = 1,NNODE
         AuxV = matmul(transpose(ddx_ddxi),dR_dx(k,:)) !ddx_ddxiT * dR_dx 
         do a1 = 1,6
	        AuxV(a1) = ddR(k,a1) - AuxV(a1)
         end do
         ddR_ddx(k,:) = matmul(dxi_dx_2T,AuxV)
      end do
c
      return
      end subroutine NURBS_BF_3D
c
c ------------------------------------------------------------------------------------
c ! --- Basis functions, their first and second derivatives at uu 
c ! --- Algorithm is described in "The NURBS Book" by Piegl and Tiller, Springer-Verlag Berlin 1997
c ------------------------------------------------------------------------------------
      subroutine deriv2(ii, ppol, uu, U_kn, Kn_Num, Nd)
c
      include 'ABA_PARAM.INC'
      integer :: ii, ppol, Kn_Num
      real(8) :: uu
      real(8), dimension(Kn_Num) :: U_kn
      real(8), dimension(3,ppol+1) :: Nd
      real(8), dimension(ppol+1) :: left, right
      real(8), dimension(ppol+1,ppol+1) :: ndu
      real(8), dimension(2,ppol+1) :: aa
      real(8) :: dd, saved, temp
      integer :: jj, j1, j2, pk, rr, rk, s1, s2, kk
c
      left = 0.0d0
      right = 0.0d0
      ndu = 0.0d0
      ndu(1,1) = 1.0d0
      Nd = 0.0d0
      aa = 0.0d0
c
      do jj = 1,ppol
         left(jj+1) = uu - U_kn(ii+1-jj)
         right(jj+1) = U_kn(ii+jj) - uu
         saved = 0.0d0
         do rr = 0,jj-1
		    ndu(jj+1,rr+1) = right(rr+2) + left(jj-rr+1)
		    temp = ndu(rr+1,jj)/ndu(jj+1,rr+1)
		    ndu(rr+1,jj+1) = saved + right(rr+2)*temp
		    saved = left(jj-rr+1)*temp
         end do
	     ndu(jj+1,jj+1) = saved
      end do
c
      do jj = 0,ppol
         Nd(1,jj+1) = ndu(jj+1,ppol+1)
      end do
      do rr = 0,ppol 
         s1 = 0
         s2 = 1 
         aa(1,1) = 1.0d0
         do kk = 1,2
            dd = 0.0d0
            rk = rr-kk
            pk = ppol-kk
            if (rr >= kk) then
                aa(s2+1,1) = aa(s1+1,1)/ndu(pk+2,rk+1)
                dd = aa(s2+1,1)*ndu(rk+1,pk+1)
            end if
            if (rk >= -1) then
                j1 = 1
            else
                j1 = -rk
            end if
            if ((rr-1) <= pk) then
                j2 = kk-1
            else
                j2 = ppol-rr
            end if
            do jj = j1,j2
               aa(s2+1,jj+1) = (aa(s1+1,jj+1) - aa(s1+1,jj))/
     *			                ndu(pk+2,rk+jj+1)
               dd = dd + aa(s2+1,jj+1)*ndu(rk+jj+1,pk+1)
            end do
            if (rr <= pk) then
                aa(s2+1,kk+1) = -aa(s1+1,kk)/ndu(pk+2,rr+1)
                dd = dd + aa(s2+1,kk+1)*ndu(rr+1,pk+1)
            end if
            Nd(kk+1,rr+1) = dd
            jj = s1
            s1 = s2
            s2 = jj
         end do
      end do
c
      rr = ppol
      do kk = 1,2
         do jj = 0,ppol
            Nd(kk+1,jj+1) = Nd(kk+1,jj+1)*rr
         end do
         rr = rr*(ppol-kk)
      end do
c
      return
      end subroutine deriv2