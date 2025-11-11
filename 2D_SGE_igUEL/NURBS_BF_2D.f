c Copyright (c) 2017 Viacheslav Balobanov, Sergei Khakalo
c calculates the 1st and 2nd derivative of NURBS basis functions in uv-space
c
      module module_inv
      implicit none
      contains
c     Auxiliary function performing a direct calculation of the inverse of a 3×3 matrix.
      function matinv3(A) result(B)
c
         real(8), intent(in) :: A(3,3)   !! Matrix
         real(8)             :: B(3,3)   !! Inverse matrix
         real(8)             :: detinv
c
c        Calculate the inverse determinant of the matrix
         detinv = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     *          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     *          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
c
c        Calculate the inverse of the matrix
         B(1,1) = + 1/detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
         B(2,1) = - 1/detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
         B(3,1) = + 1/detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
         B(1,2) = - 1/detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
         B(2,2) = + 1/detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
         B(3,2) = - 1/detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
         B(1,3) = + 1/detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
         B(2,3) = - 1/detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
         B(3,3) = + 1/detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
      end function matinv3
c
c     Auxiliary function performing a direct calculation of the inverse of a 2×2 matrix.
      function matinv2(A) result(B)
c
         real(8), intent(in) :: A(2,2)   !! Matrix
         real(8)             :: B(2,2)   !! Inverse matrix
         real(8)             :: detinv
c
c        Calculate the inverse determinant of the matrix
         detinv = A(1,1)*A(2,2) - A(1,2)*A(2,1)
c
c        Calculate the inverse of the matrix
         B(1,1) = + 1/detinv * A(2,2)
         B(2,1) = - 1/detinv * A(2,1)
         B(1,2) = - 1/detinv * A(1,2)
         B(2,2) = + 1/detinv * A(1,1)
      end function matinv2
c
      end module module_inv
c
c----------------------MAIN SUBROUTINE STARTS HERE----------------------------
      subroutine NURBS_BF_2D(p,i_el,uu,U_kn,Kn_Num_U,
     1                       q,j_el,vv,V_kn,Kn_Num_V,
     2                       w_CPs,NOCPs,NNODE,COORDS,MCRD,
     3                       ndim,RR,dR_dx,ddR_ddx,Djac)
c
      use module_inv
c
      include 'ABA_PARAM.INC'
c
      integer NOCPs, Kn_Num_U, Kn_Num_V, p, i_el, q, j_el
c
      double precision uu, vv, COORDS, w_CPs, U_kn, V_kn
      dimension COORDS(MCRD,NNODE), w_CPs(NOCPs), U_kn(Kn_Num_U), V_kn(Kn_Num_V)
c
      double precision dW(2), ddW(3), NN(3,p+1), MM(3,q+1),
     *	        dR(NNODE,2), ddR(NNODE,3),
     *	        dx_dxi(ndim,ndim), dxi_dx(ndim,ndim),
     *	        ddx_ddxi(ndim,3), AuxV(3),
     *	        dx_dxi_2T(3,3), dxi_dx_2T(3,3),
     *	        RR(NNODE), dR_dx(NNODE,2), ddR_ddx(NNODE,3),
     *	        Djac, WW
c
      integer k, No_w, b, c, a1, a2
c
      call deriv2(i_el,p,uu,U_kn,Kn_Num_U,NN) !basis functions in uu
      call deriv2(j_el,q,vv,V_kn,Kn_Num_V,MM) !basis functions in vv
c
      RR = 0.0d0
      dR = 0.0d0
      ddR = 0.0d0
      dR_dx = 0.0d0
      ddR_ddx = 0.0d0
c
      dxi_dx = 0.0d0
      dx_dxi = 0.0d0
c
      ddx_ddxi = 0.0d0
c
      WW  = 0.0d0
      dW  = 0.0d0
      ddW = 0.0d0
c
      k = 0
      do c = 0,q
         do b = 0,p
            k = k+1
		No_w  = (Kn_Num_U-p-1)*(j_el-q+c-1) + i_el-p+b   !works only for vectors w/o inner repetitions
            RR(k) = NN(1,b+1)*MM(1,c+1)*w_CPs(No_w)
            WW = WW + RR(k)
c
c           First parametric derivatives 1-du, 2-dv
            dR(k,1) = NN(2,b+1)*MM(1,c+1)*w_CPs(No_w)
            dW(1) = dW(1) + dR(k,1)
            dR(k,2) = NN(1,b+1)*MM(2,c+1)*w_CPs(No_w)
            dW(2) = dW(2) + dR(k,2)
c
c           Second parametric derivatives 1-dudu, 2-dvdv, 3-dudv
            ddR(k,1) = NN(3,b+1)*MM(1,c+1)*w_CPs(No_w)
            ddW(1) = ddW(1) + ddR(k,1)
            ddR(k,2) = NN(1,b+1)*MM(3,c+1)*w_CPs(No_w)
            ddW(2) = ddW(2) + ddR(k,2)
            ddR(k,3) = NN(2,b+1)*MM(2,c+1)*w_CPs(No_w)
            ddW(3) = ddW(3) + ddR(k,3)
         end do
      end do
c
c     Divide by weight sum W
      do k=1,NNODE
         ddR(k,1) = ddR(k,1)/WW - 2*dR(k,1)*dW(1)/WW**2
     *   -RR(k)*ddW(1)/WW**2 + 2*RR(k)*dW(1)**2/WW**3
c
         ddR(k,2) = ddR(k,2)/WW - 2*dR(k,2)*dW(2)/WW**2
     *   -RR(k)*ddW(2)/WW**2 + 2*RR(k)*dW(2)**2/WW**3
c
         ddR(k,3) = ddR(k,3)/WW - dR(k,1)*dW(2)/WW**2 -
     *   dR(k,2)*dW(1)/WW**2
     *   -RR(k)*ddW(3)/WW**2 + 2*RR(k)*dW(1)*dW(2)/WW**3
c
         dR(k,1) = dR(k,1)/WW - RR(k)*dW(1)/WW**2
         dR(k,2) = dR(k,2)/WW - RR(k)*dW(2)/WW**2
c
         RR(k) = RR(k)/WW
      end do
c
c     Jacobian and its derivatives from parameter space to physical and inversion
      k = 0
      do c = 0,q
         do b = 0,p
            k = k+1 ! = (p+1)*(q+1)*d + (p+1)*c + b + 1
			do a1 = 1,ndim
			   do a2 = 1,ndim
				  dx_dxi(a1,a2)   = dx_dxi(a1,a2) + COORDS(a1,k)*dR(k,a2)
				  ddx_ddxi(a1,a2) = ddx_ddxi(a1,a2) + COORDS(a1,k)*ddR(k,a2)
			   end do
				  ddx_ddxi(a1,3)  = ddx_ddxi(a1,3) + COORDS(a1,k)*ddR(k,3)
			end do
         end do
      end do
c     Calculate the inverse determinant of the matrix
      Djac = dx_dxi(1,1)*dx_dxi(2,2) - dx_dxi(1,2)*dx_dxi(2,1)
c
      dxi_dx = matinv2(dx_dxi)
c
c     Basis functions and derivatives w.r.t. physical coordinates
      do k = 1,NNODE
	     do a1 = 1,ndim
	        do a2 = 1,ndim
		       dR_dx(k,a1) = dR_dx(k,a1) + dxi_dx(a2,a1)*dR(k,a2)
            end do
         end do
      end do
c
      dx_dxi_2T(1,1) = dx_dxi(1,1)**2
      dx_dxi_2T(1,2) = dx_dxi(2,1)**2
      dx_dxi_2T(1,3) = 2*dx_dxi(1,1)*dx_dxi(2,1)
c
      dx_dxi_2T(2,1) = dx_dxi(1,2)**2
      dx_dxi_2T(2,2) = dx_dxi(2,2)**2
      dx_dxi_2T(2,3) = 2*dx_dxi(1,2)*dx_dxi(2,2)
c
      dx_dxi_2T(3,1) = dx_dxi(1,1)*dx_dxi(1,2)
      dx_dxi_2T(3,2) = dx_dxi(2,1)*dx_dxi(2,2)
      dx_dxi_2T(3,3) = dx_dxi(1,1)*dx_dxi(2,2) + dx_dxi(1,2)*dx_dxi(2,1)
c
      dxi_dx_2T = matinv3(dx_dxi_2T)
c
      do k = 1,NNODE
         AuxV = matmul(transpose(ddx_ddxi),dR_dx(k,:)) !ddx_ddxiT * dR_dx
         do a1 = 1,3
	        AuxV(a1) = ddR(k,a1) - AuxV(a1)
         end do
         ddR_ddx(k,:) = matmul(dxi_dx_2T,AuxV)
      end do
c
      return
      end subroutine NURBS_BF_2D
c
c ------------------------------------------------------------------------------------
c Basis functions, their first and second derivatives at uu 
c Algorithm described in "The NURBS Book" by Piegl and Tiller, Springer-Verlag Berlin 1997
c
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
