c Copyright (c) 2017 Viacheslav Balobanov
c calculates the 1st and 2nd derivative of NURBS basis functions in uv-space
c
      module module_inv
      implicit none
      contains
c     Auxilary function performing a direct calculation of the inverse of a 3×3 matrix.	  
      function matinv3(A) result(B)
	  
         real(8), intent(in) :: A(3,3)   !! Matrix
         real(8)             :: B(3,3)   !! Inverse matrix
         real(8)             :: detinv
c
c        Calculate the inverse determinant of the matrix
         detinv = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     1          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     2          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
	 
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

c     Auxilary function performing a direct calculation of the inverse of a 6×6 matrix.	  
      function matinv6(M) result(Minv)
	  
         real(8), intent(in) :: M(6,6)   !! Matrix
         real(8)             :: A(3,3)   !! Submatrix of M
         real(8)             :: B(3,3)   !! Submatrix of M
         real(8)             :: C(3,3)   !! Submatrix of M
         real(8)             :: D(3,3)   !! Submatrix of M
         real(8)             :: M_D(3,3) !! The Schur complement of the block D of the matrix M
         real(8)             :: M_Dinv(3,3)!! Inverse matrix
         real(8)             :: Dinv(3,3)!! Inverse matrix
         real(8)             :: M_Dinv__BDinv(3,3)!! Inverse matrix
         real(8)             :: Minv(6,6)!! Inverse matrix
c
c        Define submatrices of M
         A = M(1:3,1:3)
         B = M(1:3,4:6)
         B(1,1) = M(1,4)
         B(1,2) = M(1,5)
         B(1,3) = M(1,6)
         B(2,1) = M(2,4)
         B(2,2) = M(2,5)
         B(2,3) = M(2,6)
         B(3,1) = M(3,4)
         B(3,2) = M(3,5)
         B(3,3) = M(3,6)
         C = M(4:6,1:3)
         D = M(4:6,4:6)
	 
c        The Schur complement of the block D of the matrix M
         Dinv = matinv3(D)
         M_D = A - matmul(B,matmul(Dinv,C))
         M_Dinv = matinv3(M_D)
         M_Dinv__BDinv = matmul(M_Dinv,matmul(B,Dinv))

c        Calculate the blocks of inverse matrix
         Minv(1:3,1:3) =  M_Dinv
         Minv(1:3,4:6) = -matmul(M_Dinv,matmul(B,Dinv))
         Minv(4:6,1:3) = -matmul(Dinv,matmul(C,M_Dinv))
         Minv(4:6,4:6) = Dinv + matmul(Dinv,matmul(C,M_Dinv__BDinv))

      end function matinv6
	  
      subroutine inverse(a,c,n)
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !=========================================================== 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k
c
      integer max_row, s
      double precision aug(n,n),max_val, factor, det

      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0
      aug = a
      s = 1
      ! step 1: forward elimination
c
c
c
c      do k=1, n-1
c
c      max_row = k
c      max_val = ABS(aug(k,k))
c      DO i = k+1, n
c        IF (ABS(aug(i,k)) > max_val) THEN
c          max_row = i
c          max_val = ABS(aug(i,k))
c        END IF
c      END DO
c      
c      ! Swap the k-th row with the row with the maximum absolute value
c      IF (max_row /= k) THEN
c        s = -s
c        DO j = k, n+1
c          factor = aug(k,j)
c          aug(k,j) = aug(max_row,j)
c          aug(max_row,j) = factor
c        END DO
c      END IF
c      
c      ! Eliminate the k-th column below the k-th row
c      IF (aug(k,k) == 0.0) THEN
c        det = 0.0
c      ELSE
c        DO i = k+1, n
c          factor = aug(i,k) / aug(k,k)
c          DO j = k+1, n+1
c            aug(i,j) = aug(i,j) - factor * aug(k,j)
c          END DO
c          aug(i,k) = 0.0
c        END DO
c      END IF
c      end do
c      a = aug
c
c
c
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do
c
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
c
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse
c
      Subroutine LUDCMP(A,C,N)
      IMPLICIT NONE
      integer, parameter :: nmax = 100
      real, parameter :: tiny = 1.D-30
c     
      real*8, dimension(N,N) :: A
      real*8, dimension(N,N) :: C
      real*8, dimension(N) :: B
      integer :: N
      integer :: D, CODE
      integer, dimension(N) :: INDX
      !f2py depend(N) A, indx
c      
      REAL*8  :: AMAX, DUM, SUMM, VV(NMAX)
      INTEGER :: i, j, k, imax, II, LL
c      
      D=1; CODE=0
c      
      DO I=1,N
        AMAX=0.d0
        DO J=1,N
          IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
        END DO ! j loop
        IF(AMAX.LT.TINY) THEN
          CODE = 1
          RETURN
        END IF
        VV(I) = 1.d0 / AMAX
      END DO ! i loop
c      
      DO J=1,N
        DO I=1,J-1
          SUMM = A(I,J)
          DO K=1,I-1
            SUMM = SUMM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUMM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
          SUMM = A(I,J)
          DO K=1,J-1
            SUMM = SUMM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUMM
          DUM = VV(I)*DABS(SUMM)
          IF(DUM.GE.AMAX) THEN
            IMAX = I
            AMAX = DUM
          END IF
        END DO ! i loop  
c        
        IF(J.NE.IMAX) THEN
          DO K=1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
        END IF
c      
        INDX(J) = IMAX
        IF(DABS(A(J,J)) < TINY) A(J,J) = TINY
c      
        IF(J.NE.N) THEN
          DUM = 1.d0 / A(J,J)
          DO I=J+1,N
            A(I,J) = A(I,J)*DUM
          END DO ! i loop
        END IF 
      END DO ! j loop
c
c
c
      DO k=1,N
      B = 0.0
      B(k) = 1.0
      II = 0
      DO I=1,N
        LL = INDX(I)
        SUMM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
          DO J=II,I-1
            SUMM = SUMM - A(I,J)*B(J)
          END DO ! j loop
        ELSE IF(SUMM.NE.0.d0) THEN
          II = I
        END IF
        B(I) = SUMM
      END DO ! i loop
c      
      DO I=N,1,-1
        SUMM = B(I)
        IF(I < N) THEN
          DO J=I+1,N
            SUMM = SUMM - A(I,J)*B(J)
          END DO ! j loop
        END IF
        B(I) = SUMM / A(I,I)
      END DO ! i loop
      C(:,k) = B
      END DO ! k loop
c
c
c
      RETURN
      END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
c      Subroutine LUBKSB(A, N, INDX, B)
c      integer, intent(in) :: N 
c      real*8, intent(in), dimension(N,N) :: A
c      integer, intent(in), dimension(N) :: INDX
c      real*8, intent(inout), dimension(N) :: B
c      !f2py depend(N) A, INDX, B
c      
c      REAL*8  SUMM
c      
c      II = 0
c      
c      DO I=1,N
c        LL = INDX(I)
c        SUMM = B(LL)
c        B(LL) = B(I)
c        IF(II.NE.0) THEN
c          DO J=II,I-1
c            SUMM = SUMM - A(I,J)*B(J)
c          END DO ! j loop
c        ELSE IF(SUMM.NE.0.d0) THEN
c          II = I
c        END IF
c        B(I) = SUMM
c      END DO ! i loop
c      
c      DO I=N,1,-1
c        SUMM = B(I)
c        IF(I < N) THEN
c          DO J=I+1,N
c            SUMM = SUMM - A(I,J)*B(J)
c          END DO ! j loop
c        END IF
c        B(I) = SUMM / A(I,I)
c      END DO ! i loop
c      
c      RETURN
c      END subroutine LUBKSB

      end module module_inv

c----------------------MAIN SUBROUTINE STARTS HERE----------------------------
      subroutine NURBS_BF_3D(p,i_el,uu,U_kn,Kn_Num_U,
     1                       q,j_el,vv,V_kn,Kn_Num_V, 
     2                       r,k_el,ww,W_kn,Kn_Num_W, 
     3                       w_CPs,NOCPs,NNODE,COORDS,MCRD,
     4                       RR,dR_dx,ddR_ddx,Djac,JELEM)
      use module_inv
      include 'ABA_PARAM.INC'
      integer NOCPs,Kn_Num_U,Kn_Num_V,Kn_Num_W,
     *        p,i_el,q,j_el,r,k_el
      double precision uu,vv,ww,COORDS(MCRD,NNODE),w_CPs(NOCPs),
     *                 U_kn(Kn_Num_U),V_kn(Kn_Num_V),W_kn(Kn_Num_W)
      double precision RR, dR_dx, ddR_ddx
      dimension RR(NNODE), dR_dx(NNODE,3), ddR_ddx(NNODE,6)

      double precision WW,dW(3),ddW(6),NN(3,p+1),MM(3,q+1),LL(3,r+1),
     *	               dx_dxi(3,3),dxi_dx(3,3),
     *	               ddx_ddxi(3,6),dx_dxi_2T(6,6),dxi_dx_2T(6,6),
     *	               dR(NNODE,3), ddR(NNODE,6), AuxV(6), Djac
      integer k, No_w, b, c, d, a1, a2
      call deriv2(i_el,p,uu,U_kn,Kn_Num_U,NN) ! basisfunc in uu
      call deriv2(j_el,q,vv,V_kn,Kn_Num_V,MM) ! basisfunc in vv
      call deriv2(k_el,r,ww,W_kn,Kn_Num_W,LL) ! basisfunc in ww

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
      WW = 0.0d0
      dW = 0.0d0
      ddW = 0.0d0
c
      k = 0
      do d = 0,r
         do c = 0,q
            do b = 0,p
               k = k+1
			   No_w = (Kn_Num_U-p-1)*(j_el-q+c-1) + i_el-p+b +
     *			      (k_el-r+d-1)*(Kn_Num_U-p-1)*(Kn_Num_V-q-1)   !works only for vectors w/o inner repeats

               RR(k) = NN(1,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               WW = WW + RR(k)
c
c              First parametric derivatives 1-du, 2-dv, 3-dw
               dR(k,1) = NN(2,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               dW(1) = dW(1) + dR(k,1)
               dR(k,2) = NN(1,b+1)*MM(2,c+1)*LL(1,d+1)*w_CPs(No_w)
               dW(2) = dW(2) + dR(k,2)
               dR(k,3) = NN(1,b+1)*MM(1,c+1)*LL(2,d+1)*w_CPs(No_w)
               dW(3) = dW(3) + dR(k,3)
c
c              Second parametric derivatives 1-dudu, 2-dvdv, 3-dwdw, 4-dvdw, 5-dudw, 6-dudv
               ddR(k,1) = NN(3,b+1)*MM(1,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddW(1) = ddW(1) + ddR(k,1)
               ddR(k,2) = NN(1,b+1)*MM(3,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddW(2) = ddW(2) + ddR(k,2)
               ddR(k,3) = NN(1,b+1)*MM(1,c+1)*LL(3,d+1)*w_CPs(No_w)
               ddW(3) = ddW(3) + ddR(k,3)
               ddR(k,4) = NN(1,b+1)*MM(2,c+1)*LL(2,d+1)*w_CPs(No_w)
               ddW(4) = ddW(4) + ddR(k,4)
               ddR(k,5) = NN(2,b+1)*MM(1,c+1)*LL(2,d+1)*w_CPs(No_w)
               ddW(5) = ddW(5) + ddR(k,5)
               ddR(k,6) = NN(2,b+1)*MM(2,c+1)*LL(1,d+1)*w_CPs(No_w)
               ddW(6) = ddW(6) + ddR(k,6)
            end do
         end do
      end do

c     Divide by weight sum W
      do k=1,NNODE
c
           ddR(k,1) = ddR(k,1)/WW - 2*dR(k,1)*dW(1)/WW**2
     *     -RR(k)*ddW(1)/WW**2 + 2*RR(k)*dW(1)**2/WW**3
c
           ddR(k,2) = ddR(k,2)/WW - 2*dR(k,2)*dW(2)/WW**2
     *     -RR(k)*ddW(2)/WW**2 + 2*RR(k)*dW(2)**2/WW**3
c
           ddR(k,3) = ddR(k,3)/WW - 2*dR(k,3)*dW(3)/WW**2
     *     -RR(k)*ddW(3)/WW**2 + 2*RR(k)*dW(3)**2/WW**3
c
           ddR(k,4) = ddR(k,4)/WW - dR(k,2)*dW(3)/WW**2 -
     *     dR(k,3)*dW(2)/WW**2
     *     -RR(k)*ddW(4)/WW**2 + 2*RR(k)*dW(2)*dW(3)/WW**3
c
           ddR(k,5) = ddR(k,5)/WW - dR(k,1)*dW(3)/WW**2 -
     *     dR(k,3)*dW(1)/WW**2
     *     -RR(k)*ddW(5)/WW**2 + 2*RR(k)*dW(1)*dW(3)/WW**3
c
           ddR(k,6) = ddR(k,6)/WW - dR(k,1)*dW(2)/WW**2 -
     *     dR(k,2)*dW(1)/WW**2
     *     -RR(k)*ddW(6)/WW**2 + 2*RR(k)*dW(1)*dW(2)/WW**3
c
           dR(k,1) = dR(k,1)/WW - RR(k)*dW(1)/WW**2
           dR(k,2) = dR(k,2)/WW - RR(k)*dW(2)/WW**2
           dR(k,3) = dR(k,3)/WW - RR(k)*dW(3)/WW**2
c
           RR(k) = RR(k)/WW
      end do
c
c     Jacobian and its derivatives from parameter space to physical and inversion
      k = 0
      do d = 0,r
         do c = 0,q
            do b = 0,p
               k = k+1 ! = (p+1)*(q+1)*d + (p+1)*c + b + 1
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
c        Calculate the inverse determinant of the matrix
         Djac = (dx_dxi(1,1)*dx_dxi(2,2)*dx_dxi(3,3) - dx_dxi(1,1)*dx_dxi(2,3)*dx_dxi(3,2)
     *          - dx_dxi(1,2)*dx_dxi(2,1)*dx_dxi(3,3) + dx_dxi(1,2)*dx_dxi(2,3)*dx_dxi(3,1)
     *          + dx_dxi(1,3)*dx_dxi(2,1)*dx_dxi(3,2) - dx_dxi(1,3)*dx_dxi(2,2)*dx_dxi(3,1))
c
      dxi_dx = matinv3(dx_dxi)
c     Basis functions and derivatives w.r.t. physical coordinates
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
c      if (JELEM .EQ. 240) then
c      write(7,*)'dx_dxi_2T, before =', dx_dxi_2T
c      end if
c
      dxi_dx_2T = matinv6(dx_dxi_2T)
      if (JELEM .EQ. 9) then
          write(7,*)'JELEM =', JELEM
          write(7,*)'dx_dxi_2T =', dx_dxi_2T
          write(7,*)'dxi_dx_2T, Schur =', dxi_dx_2T
      end if
c
c      call inverse(dx_dxi_2T,dxi_dx_2T,6)
      call LUDCMP(dx_dxi_2T,dxi_dx_2T,6)
c      if (JELEM .EQ. 9) then
c          write(7,*)'dxi_dx_2T, LU =', dxi_dx_2T
c      end if
c
      do k = 1,NNODE
         AuxV = matmul(transpose(ddx_ddxi),dR_dx(k,:)) !ddx_ddxiT * dR_dx 
         do a1 = 1,6
	        AuxV(a1) = ddR(k,a1) - AuxV(a1)
         end do
         ddR_ddx(k,:) = matmul(dxi_dx_2T,AuxV)
      end do
c
      if (JELEM .EQ. 9) then
c      write(7,*)'JELEM =', JELEM
c      write(7,*)'dx_dxi_2T =', dx_dxi_2T
      write(7,*)'dxi_dx_2T, LUDCMP =', dxi_dx_2T
c      write(7,*)'dxi_dx_2T, Schur =', dxi_dx_2T
c      write(7,*)'FIFFERENCE LU-Schur =', dxi_dx_2T - matinv6(dx_dxi_2T)
c      write(7,*)'dx_dxi_2T, after =', dx_dxi_2T
c      write(7,*)'dxi_dx_2T =', dxi_dx_2T
      write(7,*)'AuxV =', AuxV
      write(7,*)'ddR_ddx =', ddR_ddx
      write(7,*)'dR_dx =', dR_dx
      write(7,*)'ddx_ddxi =', ddx_ddxi
      write(7,*)'ddR =', ddR
      end if
c
      return
      end subroutine NURBS_BF_3D
c
c ------------------------------------------------------------------------------------
! Basis functions, their first and second derivatives at uu 
! Algorithm described in "The NURBS Book" by Piegl and Tiller, Springer-Verlag Berlin 1997
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
