! Copyright (c) 2025 Sergei Khakalo, Viacheslav Balobanov
! The codes (main and included subroutines) are distributed under the MIT License
c ------------------------------------------------------------------------------------
!
! igUEL-SGE-Solid: Abaqus user element (UEL) implementation for isogeometric analysis
! of 2D strain-gradient-elastic solids
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
      parameter (ndi=3, nshr=1, ntens=4, ndim=2, ndof=2)
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
c !    w_CPs     -- weights of control points
c !    gauss_U   -- quadrature (Gauss) nodes in 1-direction
c !    gauss_V   -- quadrature (Gauss) nodes in 2-direction
c !    w_U       -- weights of quadrature (Gauss) nodes in 1-direction
c !    w_V       -- weights of quadrature (Gauss) nodes in 2-direction
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
      double precision RR, dR, ddR, map, KV_U, KV_V, u_ip, v_ip, w_CPs,
     *                 w_U, w_V, Jdet, gauss_U, gauss_V, dV
      dimension RR(NNODE),
     *          dR(NNODE,2),
     *          ddR(NNODE,3),
     *          KV_U(PROPS(8)),
     *          KV_V(PROPS(9)),
     *          w_CPs(PROPS(5)),
     *          w_U(PROPS(1)),
     *          w_V(PROPS(2)),
     *          gauss_U(PROPS(1)),
     *          gauss_V(PROPS(2))
c
      double precision stress, gstress, stran, gstran,
     *                 force, stiff, C_St_Matr, A_St_Matr,
     *                 BF_x, BF_y, E_mod, nu_mod, g_par
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
      integer ninpt_U, ninpt_V, NOE_U, NOE_V, NOE, NOCPs,
     *        Ppol, Qpol, Kn_Num_U, Kn_Num_V,
     *        Num_FE_PP_U, Num_FE_PP_V, i_U, i_V, k1, k2, j,
     *        kintk_eta, kintk_xi, kintk, i_el, j_el, istat,
     *        El_Output, open_ind, ind_sum
c      
      logical :: file_exists
c
c !--- PID-properties from .inp --------------------------------
c
      ninpt_U  = PROPS(1)            !Number of integration points in U-direction
      ninpt_V  = PROPS(2)            !Number of integration points in V-direction
      NOE_U  = PROPS(3)              !Number of elements in U-direction
      NOE_V  = PROPS(4)              !Number of elements in V-direction
      NOE  = NOE_U*NOE_V             !Total number of elements
      NOCPs  = PROPS(5)              !Number of control points
      Ppol  = PROPS(6)               !Polynomial order in U-direction
      Qpol  = PROPS(7)               !Polynomial order in V-direction
      Kn_Num_U  = PROPS(8)           !Num of members in knot vector U
      Kn_Num_V  = PROPS(9)           !Num of members in knot vector V
      Num_FE_PP_U  = PROPS(10)       !Num of FE-elements dividing each NURBS-element in U-dir for PostProc
      Num_FE_PP_V  = PROPS(11)       !Num of FE-elements dividing each NURBS-element in V-dir for PostProc
      El_Output  = PROPS(12)         !Element output parameter: 0 - U,S,E; 1 - U,S,E,gS,gE
      BF_x  = PROPS(13)              !x-component of body force vector
      BF_y  = PROPS(14)              !y-component of body force vector
      E_mod  = PROPS(15)             !Young's modulus
      nu_mod  = PROPS(16)            !Poisson's ratio
      g_par  = PROPS(17)             !Strain gradient length scale parameter
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
      read(open_ind,*)(w_CPs(j),j=1,NOCPs)
      close(open_ind)
c
      i_V = 1 + INT((JELEM-0.1)/NOE_U) !current No of element in eta-direction (V-direction)
      i_U = JELEM - NOE_U*(i_V-1) !current No of element in xi-direction (U-direction)
c !--- No of knot spans
      i_el = i_U + Ppol !only for vectors w/o inner repeats / knot repetition
      j_el = i_V + Qpol !only for vectors w/o inner repeats
c
      map = (KV_U(i_el+1)-KV_U(i_el))*(KV_V(j_el+1)-KV_V(j_el))/4.d0
c
c !--- Gauss points and material elasticity matrices
c
      call GaussLegendre(ninpt_U,gauss_U,w_U)
      call GaussLegendre(ninpt_V,gauss_V,w_V)
c
      call Mat_Prop_2D(ntens,ndim,E_mod,nu_mod,g_par,C_St_Matr,A_St_Matr)
c
      RHS(:,1)    = 0.d0
      AMATRX(:,:) = 0.d0
c
c !--- Loop over Integration Points (IPs)
c
      do kintk_eta = 1, ninpt_V       ! start loop over IPs
         do kintk_xi = 1, ninpt_U
            kintk = kintk_xi + ninpt_U*(kintk_eta-1)
            u_ip = (KV_U(i_el+1)+KV_U(i_el) + 
     *              gauss_U(kintk_xi)*(KV_U(i_el+1)-KV_U(i_el)))/2.d0
            v_ip = (KV_V(j_el+1)+KV_V(j_el) + 
     *              gauss_V(kintk_eta)*(KV_V(j_el+1)-KV_V(j_el)))/2.d0
c
c !--- NURBS basis functions and derivatives w.r.t. parametric coordinates
C !--- first derivatives 1-dx, 2-dy
C !--- second parametric derivatives 1-dx^2, 2-dy^2, 3-dxdy
c
            call NURBS_BF_2D(Ppol,i_el,u_ip,KV_U,Kn_Num_U,
     1                       Qpol,j_el,v_ip,KV_V,Kn_Num_V,  
     2                       w_CPs,NOCPs,NNODE,COORDS,MCRD,
     3                       ndim,RR,dR,ddR,Jdet)
c
            dV = w_U(kintk_xi)*w_V(kintk_eta)*map*Jdet
c
c !--- Strains and gradient of strains
c
      call Strains_2D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                stran,gstran,dR,ddR)
c
c !--- Cauchy (ordinary) stress and double (higher-order) stress
c
       stress  = matmul(C_St_Matr,stran)
       gstress = matmul(A_St_Matr,gstran)
c
c !--- Define internal force vector (force) and element stiffness matrix (stiff) at the current IP
c
      call Force_Stiff_2D(NNODE,ntens,ndim,ndof,stress,gstress,dV,
     1                    dR,ddR,C_St_Matr,A_St_Matr,force,stiff)
c
c !--- Define RHS and AMATRX
c
         do k1 = 1, NNODE
            k2 = (k1-1)*ndof
c
		    RHS(k2+1,1) = RHS(k2+1,1) - force(k2+1) + BF_x*RR(k1)*dV
		    RHS(k2+2,1) = RHS(k2+2,1) - force(k2+2) + BF_y*RR(k1)*dV
c
         end do
c
		    AMATRX = AMATRX + stiff
c
         end do
      end do       ! end loop over IPs
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
              open(12000000,file = trim(folder)//'/gEx_Nodes.dat',
     *                           STATUS='REPLACE', ACTION='WRITE')
c
              open(13000000,file = trim(folder)//'/gEy_Nodes.dat',
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
              INQUIRE(file = trim(folder)//'/gEx_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(12000000,file = trim(folder)//'/gEx_Nodes.dat',STATUS='OLD')
		  close(12000000, status='DELETE')
              end if
c              
              INQUIRE(file = trim(folder)//'/gEy_Nodes.dat', EXIST=file_exists)
              if (file_exists) then
                  open(13000000,file = trim(folder)//'/gEy_Nodes.dat',STATUS='OLD')
		  close(13000000, status='DELETE')
              end if
          end if
c
      end if
c
c
      ind_PP(JELEM) = ind_PP(JELEM)+1
c
      if (ind_PP(JELEM) .EQ. 3) then
c
      call Output_2D(Kn_Num_U,Kn_Num_V,KV_U,KV_V,Ppol,Qpol,JELEM,
     1               NNODE,NOE,NOCPs,COORDS,U,NDOFEL,MCRD,TIME,
     2               Num_FE_PP_U,Num_FE_PP_V,C_St_Matr,A_St_Matr,
     3               i_el,j_el,w_CPs,cwd,KINC,ndim,ndof,ntens,El_Output)
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
      include "./Mat_Prop_2D.f"
      include "./NURBS_BF_2D.f"
      include "./Strains_2D.f"
      include "./Force_Stiff_2D.f"
      include "./Output_2D.f"
