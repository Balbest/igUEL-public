! Copyright (c) 2025 Sergei Khakalo, Viacheslav Balobanov
! The codes (main and included subroutines) are distributed under the MIT License
c ------------------------------------------------------------------------------------
!
! igUEL-SGE-Solid: Abaqus user element (UEL) implementation for isogeometric analysis
! of 3D strain-gradient-elastic solids
!
c ------------------------------------------------------------------------------------
c
      module indicator
      implicit none
      integer ind_PP(100000)
      integer ind_close(100000)
      save
      end module
c
      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
c
      use indicator
      include 'ABA_PARAM.INC'
c
      parameter (ndi=3, nshr=3, ntens=6, ndim=3, ndof=3)
c
c ------------------------------------------------------------------------------------
c ! List of parameters/internal variables:
c !
c !    ndi       -- number of direct stress components
c !    nshr      -- number of shear stress components
c !    ntens     -- total number of stress tensor components (=ndi+nshr)
c !    ndim      -- number of spatial dimensions
c !    ndof      -- number of degrees of freedom per control point (node)
c !
c !    stress    -- Cauchy (ordinary) stress
c !    gstress   -- double (higher-order) stress
c !    stran     -- elastic strain
c !    gstran    -- gradient of elastic strain
c !    C_St_Matr -- matrix of (classical) elastic moduli
c !    A_St_Matr -- matrix of (higher-order) elastic moduli
c !    RR        -- values of NURBS basis functions
c !    dR        -- values of 1st derivatives of NURBS basis functions
c !    ddR       -- values of 2nd derivatives of NURBS basis functions
c !    KV_U      -- knot vector in U-direction
c !    KV_V      -- knot vector in V-direction
c !    KV_W      -- knot vector in W-direction
c !    w_CPs     -- weights of control points
c !    gauss_U   -- quadrature (Gauss) nodes in 1-direction
c !    gauss_V   -- quadrature (Gauss) nodes in 2-direction
c !    gauss_W   -- quadrature (Gauss) nodes in 3-direction
c !    w_U       -- weights of quadrature (Gauss) nodes in 1-direction
c !    w_V       -- weights of quadrature (Gauss) nodes in 2-direction
c !    w_W       -- weights of quadrature (Gauss) nodes in 3-direction
c ------------------------------------------------------------------------------------
c
      double precision RHS, AMATRX, SVARS, PROPS, COORDS
      dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     *     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     *     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     *     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     *     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     *     JPROPS(*)
c
      double precision RR, dR, ddR, map, KV_U, KV_V, KV_W, u_ip, v_ip, w_ip, w_CPs,
     *                 w_U, w_V, w_W, Jdet, gauss_U, gauss_V, gauss_W, dV
      dimension RR(NNODE),
     *          dR(NNODE,3),
     *          ddR(NNODE,6),
     *          KV_U(PROPS(11)),
     *          KV_V(PROPS(12)),
     *          KV_W(PROPS(13)),
     *          w_CPs(PROPS(7)),
     *          w_U(PROPS(1)),
     *          w_V(PROPS(2)),
     *          w_W(PROPS(3)),
     *          gauss_U(PROPS(1)),
     *          gauss_V(PROPS(2)),
     *          gauss_W(PROPS(3))
c
      double precision stress, gstress, stran, gstran,
     *                 force, stiff, C_St_Matr, A_St_Matr,
     *                 BF_x, BF_y, BF_z, E_mod, nu_mod, g_par
c
      dimension stress(ntens),
     *          gstress(ndim*ntens),
     *          stran(ntens),
     *          gstran(ndim*ntens),
     *          force(ndof*NNODE),
     *          stiff(ndof*NNODE,ndof*NNODE),
     *          C_St_Matr(ntens,ntens),
     *          A_St_Matr(ntens*ndim,ntens*ndim)
c
      character*256 :: cwd
      character*256 :: folder
c
      integer ninpt_U, ninpt_V, ninpt_W, NOE_U, NOE_V, NOE_W, NOE, NOCPs,
     *        Ppol, Qpol, Rpol, Kn_Num_U, Kn_Num_V, Kn_Num_W, i_U, i_V, i_W,
     *        Num_FE_PP_U, Num_FE_PP_V, Num_FE_PP_W, k1, k2, j,
     *        kintk_zeta, kintk_eta, kintk_xi, kintk, i_el, j_el, k_el,
     *        iAux, istat, El_Output, open_ind, ind_sum
c      
      logical :: file_exists
c
c !--- PID-properties from .inp --------------------------------
c
      ninpt_U  = PROPS(1)            !Number of integration points in U-direction
      ninpt_V  = PROPS(2)            !Number of integration points in V-direction
      ninpt_W  = PROPS(3)            !Number of integration points in W-direction
      NOE_U  = PROPS(4)              !Number of elements in U-direction
      NOE_V  = PROPS(5)              !Number of elements in V-direction
      NOE_W  = PROPS(6)              !Number of elements in W-direction
      NOE  = NOE_U*NOE_V*NOE_W       !Total number of elements
      NOCPs  = PROPS(7)              !Number of control points
      Ppol  = PROPS(8)               !Polynomial order in U-direction
      Qpol  = PROPS(9)               !Polynomial order in V-direction
      Rpol  = PROPS(10)              !Polynomial order in W-direction
      Kn_Num_U  = PROPS(11)          !Num of members in knot vector U
      Kn_Num_V  = PROPS(12)          !Num of members in knot vector V
      Kn_Num_W  = PROPS(13)          !Num of members in knot vector W
      Num_FE_PP_U  = PROPS(14)       !Num of FE-elements dividing each NURBS-element in U-dir for PostProc
      Num_FE_PP_V  = PROPS(15)       !Num of FE-elements dividing each NURBS-element in V-dir for PostProc
      Num_FE_PP_W  = PROPS(16)       !Num of FE-elements dividing each NURBS-element in W-dir for PostProc
	  El_Output  = PROPS(17)         !Element output parameter: 0 - U,S,E; 1 - U,S,E,gS,gE
	  BF_x  = PROPS(18)              !x-component of body force vector
	  BF_y  = PROPS(19)              !y-component of body force vector
	  BF_z  = PROPS(20)              !z-component of body force vector
	  E_mod  = PROPS(21)             !Young's modulus
	  nu_mod  = PROPS(22)            !Poisson's ratio
	  g_par  = PROPS(23)             !Strain gradient length scale parameter
c
c ------------------------------------------------------------------------------------
c
c !--- Read Ks_Ws.dat ------------------------------------------------------
c
c !--- Get path to the current working directory (max 256 symbols)
      istat = 256
      call getoutdir(cwd, istat)
      folder = TRIM(ADJUSTL(cwd))//'/results'
c !--- Get knot vectors and weights
      open_ind = 1000000+JELEM
      open(open_ind,file = TRIM(ADJUSTL(cwd))//'/Ks_Ws.dat')
      read(open_ind,*)(KV_U(j),j=1,Kn_Num_U)
      read(open_ind,*)(KV_V(j),j=1,Kn_Num_V)
      read(open_ind,*)(KV_W(j),j=1,Kn_Num_W)
      read(open_ind,*)(w_CPs(j),j=1,NOCPs)
      close(open_ind)
c
      i_W  = 1 + INT((JELEM-0.1)/NOE_U/NOE_V) !current No of element in zeta-direction (W-direction)
      iAux = JELEM - NOE_U*NOE_V*(i_W-1)
      i_V  = 1 + INT((iAux-0.1)/NOE_U) !current No of element in eta-direction (V-direction)
      i_U  = iAux-NOE_U*(i_V-1) !current No of element in xi-direction (U-direction)
c !--- No of knot spans
      i_el = i_U + Ppol !only for vectors w/o inner repeats / knot repetition
      j_el = i_V + Qpol !only for vectors w/o inner repeats / knot repetition
      k_el = i_W + Rpol !only for vectors w/o inner repeats / knot repetition
c
      map = (KV_U(i_el+1)-KV_U(i_el))*(KV_V(j_el+1)-KV_V(j_el))
     *                               *(KV_W(k_el+1)-KV_W(k_el))/8.d0
c
c !--- Gauss points and material elasticity matrices
c
      call GaussLegendre(ninpt_U,gauss_U,w_U)
      call GaussLegendre(ninpt_V,gauss_V,w_V)
      call GaussLegendre(ninpt_W,gauss_W,w_W)
c
      call Mat_Prop_3D(ntens,ndim,E_mod,nu_mod,g_par,C_St_Matr,A_St_Matr)
c
      RHS(:,1)    = 0.d0
      AMATRX(:,:) = 0.d0
c
c !--- Loop over Integration Points (IPs)
c
      do kintk_zeta = 1, ninpt_W       ! start loop over IPs
         do kintk_eta = 1, ninpt_V
            do kintk_xi = 1, ninpt_U
               kintk = kintk_xi + ninpt_U*(kintk_eta-1) + ninpt_U*ninpt_V*(kintk_zeta-1)
               u_ip = (KV_U(i_el+1)+KV_U(i_el) + 
     *                 gauss_U(kintk_xi)*(KV_U(i_el+1)-KV_U(i_el)))/2.d0
               v_ip = (KV_V(j_el+1)+KV_V(j_el) + 
     *                 gauss_V(kintk_eta)*(KV_V(j_el+1)-KV_V(j_el)))/2.d0
               w_ip = (KV_W(k_el+1)+KV_W(k_el) + 
     *                 gauss_W(kintk_zeta)*(KV_W(k_el+1)-KV_W(k_el)))/2.d0
c
c !--- NURBS basis functions and derivatives w.r.t. parametric coordinates---------------------------
C !--- first derivatives 1-dx, 2-dy, 3-dz
C !--- second parametric derivatives 1-dx^2, 2-dy^2, 3-dz^2, 4-dydz, 5-dxdz, 6-dxdy
c
               call NURBS_BF_3D(Ppol,i_el,u_ip,KV_U,Kn_Num_U,
     1                          Qpol,j_el,v_ip,KV_V,Kn_Num_V, 
     2                          Rpol,k_el,w_ip,KV_W,Kn_Num_W, 
     3                          w_CPs,NOCPs,NNODE,COORDS,MCRD,
     4                          RR,dR,ddR,Jdet,JELEM)
c
               dV = w_U(kintk_xi)*w_V(kintk_eta)*w_W(kintk_zeta)*map*Jdet
c
c !--- Strains and gradient of strains
c
      call Strains_3D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                stran,gstran,dR,ddR)
c
c !--- Cauchy (ordinary) stress and double (higher-order) stress
c
       stress  = matmul(C_St_Matr,stran)
       gstress = matmul(A_St_Matr,gstran)
c
c !--- Define internal force vector (force) and element stiffness matrix (stiff) at the current IP
c
      call Force_Stiff_3D(NNODE,ntens,ndim,ndof,stress,gstress,dV,
     1                    dR,ddR,C_St_Matr,A_St_Matr,force,stiff)
c
c !--- Define RHS and AMATRX
c
         do k1 = 1, NNODE
            k2 = (k1-1)*ndof
c
		    RHS(k2+1,1) = RHS(k2+1,1) - force(k2+1) + BF_x*RR(k1)*dV
		    RHS(k2+2,1) = RHS(k2+2,1) - force(k2+2) + BF_y*RR(k1)*dV
		    RHS(k2+3,1) = RHS(k2+3,1) - force(k2+3) + BF_z*RR(k1)*dV
c
         end do
c
		    AMATRX = AMATRX + stiff
c
            end do
         end do
      end do       ! end loop over IPs
c
c
c! --- Output data
c
      if (TIME(1) .EQ. 0) then
          ind_PP(JELEM)    = 0
          ind_close(JELEM) = 0
      end if
c
      if (TIME(1) .EQ. 0 .AND. JELEM .EQ. 1) then
c
          open(3000000,file = trim(folder)//'/Nodes.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(4000000,file = trim(folder)//'/Nodes_Coords.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(5000000,file = trim(folder)//'/Elements.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(6000000,file = trim(folder)//'/Elements_Nodes.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(7000000,file = trim(folder)//'/U_Nodes.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(8000000,file = trim(folder)//'/S_Nodes.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          open(9000000,file = trim(folder)//'/E_Nodes.dat',
     *                      STATUS='REPLACE', ACTION='WRITE')
c
          if (El_Output .EQ. 1) then
              open(10000000,file = trim(folder)//'/gSx_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(11000000,file = trim(folder)//'/gSy_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(12000000,file = trim(folder)//'/gSz_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(13000000,file = trim(folder)//'/gEx_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(14000000,file = trim(folder)//'/gEy_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(15000000,file = trim(folder)//'/gEz_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
          end if
c          
          if (El_Output .EQ. 0) then
              INQUIRE(file = trim(folder)//'/gSx_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(10000000,file = trim(folder)//'/gSx_Nodes.dat',STATUS='OLD')
              close(10000000, status='DELETE')
              end if    
c
              INQUIRE(file = trim(folder)//'/gSy_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(11000000,file = trim(folder)//'/gSy_Nodes.dat',STATUS='OLD')
              close(11000000, status='DELETE')
              end if
c
              INQUIRE(file = trim(folder)//'/gSz_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(12000000,file = trim(folder)//'/gSz_Nodes.dat',STATUS='OLD')
              close(12000000, status='DELETE')
              end if
c              
              INQUIRE(file = trim(folder)//'/gEx_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(13000000,file = trim(folder)//'/gEx_Nodes.dat',STATUS='OLD')
              close(13000000, status='DELETE')
              end if
c              
              INQUIRE(file = trim(folder)//'/gEy_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(14000000,file = trim(folder)//'/gEy_Nodes.dat',STATUS='OLD')
              close(14000000, status='DELETE')
              end if
c              
              INQUIRE(file = trim(folder)//'/gEz_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(15000000,file = trim(folder)//'/gEz_Nodes.dat',STATUS='OLD')
              close(15000000, status='DELETE')
              end if
          end if
      end if
c
c
      ind_PP(JELEM) = ind_PP(JELEM)+1
c
      if (ind_PP(JELEM) .EQ. 3) then
c
      call Output_3D(Kn_Num_U,Kn_Num_V,Kn_Num_W,KV_U,KV_V,KV_W,Ppol,Qpol,Rpol,
     1               JELEM,NNODE,NOE,NOCPs,COORDS,U,NDOFEL,MCRD,TIME,
     2               Num_FE_PP_U,Num_FE_PP_V,Num_FE_PP_W,C_St_Matr,A_St_Matr,
     3               i_el,j_el,k_el,w_CPs,cwd,KINC,ndim,ndof,ntens,El_Output)
c
      ind_close(JELEM) = 1
      ind_sum = 0
      do i = 1, NOE
         ind_sum = ind_sum + ind_close(i)
      end do
c
      if (ind_sum .EQ. NOE) then
          close(3000000)
          close(4000000)
          close(5000000)
          close(6000000)
          close(7000000)
          close(8000000)
          close(9000000)
          if (El_Output .EQ. 1) then
              close(10000000)
              close(11000000)
              close(12000000)
              close(13000000)
              close(14000000)
              close(15000000)
          end if
      end if
c
      end if
c
      return
      end subroutine UEL
c
c! --- Inclusion of additional subroutines
c
      include "./GaussLegendre.f"
      include "./Mat_Prop_3D.f"
      include "./NURBS_BF_3D.f"
      include "./Strains_3D.f"
      include "./Force_Stiff_3D.f"
      include "./Output_3D.f"