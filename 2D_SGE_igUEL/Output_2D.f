      subroutine Output_2D(Kn_Num_U,Kn_Num_V,KV_U,KV_V,Ppol,Qpol,JELEM,
     1                     NNODE,NOE,NOCPs,COORDS,U,NDOFEL,MCRD,TIME,
     2                     Num_FE_PP_U,Num_FE_PP_V,C_St_Matr,A_St_Matr,
     3                     i_el,j_el,w_CPs,cwd,KINC,ndim,ndof,ntens,El_Output)
c
      include 'ABA_PARAM.INC'
c
      integer NOE,NOCPs,Num_FE_PP_U,Num_FE_PP_V,
     *        Ppol,Qpol,Kn_Num_U,Kn_Num_V,i_el,j_el,El_Output
c
      double precision KV_U, KV_V, w_CPs, COORDS
      dimension KV_U(Kn_Num_U), KV_V(Kn_Num_V), w_CPs(NOCPs),
     *          COORDS(MCRD,NNODE), U(NDOFEL), TIME(2)
c
      double precision h_xi, h_eta, NCOORD_U, NCOORD_V
c
      double precision NC_in_IGFE,
     *                 ND_in_IGFE,
     *                 stress_NC,
     *                 strain_NC,
     *                 gstress_NC,
     *                 gstrain_NC
      dimension NC_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ndim),
     *          ND_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ndof),
     *          stress_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ntens),
     *          strain_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ntens),
     *          gstress_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ntens*ndim),
     *          gstrain_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1),ntens*ndim)
c
      double precision R_NC, dR_NC, ddR_NC, Jdet
      dimension R_NC(NNODE), dR_NC(NNODE,2), ddR_NC(NNODE,3)
c
      double precision stran, gstran,
     *                 stress, gstress,
     *                 C_St_Matr, A_St_Matr
      dimension stran(ntens), gstran(ntens*ndim),
     *          stress(ntens), gstress(ntens*ndim),
     *          C_St_Matr(ntens,ntens), A_St_Matr(ntens*ndim,ntens*ndim)
c
      integer i,j,ij,m,mij,k1,i_u,i_v,i_uv,node_shift,
     *        i1,i2,i3,i4,Nodes_PP,FE_PP
c
      integer Nodes_in_IGFE, FE_Nodes
      dimension Nodes_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1)),
     *          FE_Nodes(Num_FE_PP_U*Num_FE_PP_V,5)
c
      character*256 :: cwd
      character*256 :: folder,inc
c
c !-----------------------------------------------------------------------------------
c
      strain_NC = 0.d0 !--- Strains at post-processing (PP) nodes
      stress_NC = 0.d0 !--- Stresses at PP nodes
      gstrain_NC = 0.d0 !--- Gradient of strains at PP nodes
      gstress_NC = 0.d0 !--- Higher-order stresses at PP nodes
c
c !--- PP node coordinates: NCOORD_U, NCOORD_V
c
      h_xi   = 2.d0/Num_FE_PP_U
      h_eta  = 2.d0/Num_FE_PP_V
      Nodes_PP = (Num_FE_PP_U+1)*(Num_FE_PP_V+1)
      FE_PP    =  Num_FE_PP_U*Num_FE_PP_V
c
      i_uv = 0
      do i_v = 1, Num_FE_PP_V+1 ! --- Start loop over PP nodes in V-direction
         NCOORD_V = -1.d0 + h_eta*(i_v-1)
         NCOORD_V = (KV_V(j_el+1)+KV_V(j_el) + 
     *               NCOORD_V*(KV_V(j_el+1)-KV_V(j_el)))/2.d0
c
         do i_u = 1, Num_FE_PP_U+1 ! --- Start loop over PP nodes in U-direction
            NCOORD_U = -1.d0 + h_xi*(i_u-1)
            NCOORD_U = (KV_U(i_el+1)+KV_U(i_el) + 
     *                  NCOORD_U*(KV_U(i_el+1)-KV_U(i_el)))/2.d0
c
            i_uv = i_uv + 1
c
            call NURBS_BF_2D(Ppol,i_el,NCOORD_U,KV_U,Kn_Num_U,
     1                       Qpol,j_el,NCOORD_V,KV_V,Kn_Num_V, 
     2                       w_CPs,NOCPs,NNODE,COORDS,MCRD,
     3                       ndim,R_NC,dR_NC,ddR_NC,Jdet)
c
c !--- Nodes_in_IGFE: (global) numbers of PP nodes in current knot span (JELEM)
c !--- NC_in_IGFE   : (global) coordinates of PP nodes in current knot span (JELEM)
c !--- ND_in_IGFE   : (global) displacements of PP nodes in current knot span (JELEM)
c
            Nodes_in_IGFE(i_uv) = (JELEM-1)*Nodes_PP + i_uv
c
            NC_in_IGFE(i_uv,:) = 0.d0
            ND_in_IGFE(i_uv,:) = 0.d0
c
            do m = 1, ndof
               ij = 0
               do j = 1, Qpol+1
                  do i = 1, Ppol+1
                     ij = ij + 1
                     mij = m + (ij-1)*ndof
                     if (m .le. ndim) then 
                         NC_in_IGFE(i_uv,m) = NC_in_IGFE(i_uv,m) +
     *                                        COORDS(m,ij)*R_NC(ij)
                     end if
                         ND_in_IGFE(i_uv,m) = ND_in_IGFE(i_uv,m) +
     *                                        U(mij)*R_NC(ij)
                  end do
               end do
            end do
c
c !--- stress_NC: stresses at PP nodes in current knot span (JELEM)
c !--- strain_NC: strains at PP nodes in current knot span (JELEM)
c !--- gstress_NC: higher-order stresses at PP nodes in current knot span (JELEM)
c !--- gstrain_NC: gradient of strains at PP nodes in current knot span (JELEM)
c
            call Strains_2D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                      stran,gstran,dR_NC,ddR_NC)
c
            stress  = matmul(C_St_Matr,stran)
            gstress = matmul(A_St_Matr,gstran)
c
            do k1=1,ntens
               stress_NC(i_uv,k1) = stress(k1)
               strain_NC(i_uv,k1) = stran(k1)
c
               gstress_NC(i_uv,k1) = gstress(k1)
               gstrain_NC(i_uv,k1) = gstran(k1)
               gstress_NC(i_uv,k1+ntens) = gstress(k1+ntens)
               gstrain_NC(i_uv,k1+ntens) = gstran(k1+ntens)
            end do
c
         end do ! --- End loop over PP nodes in U-direction
      end do ! --- End loop over PP nodes in V-direction
c
c !-----------------------------------------------------------------------------------
c !--- This part generates .dat files in folder 'results'
c !--- to be used to create ODB file with results
c !-----------------------------------------------------------------------------------
c
      folder = TRIM(ADJUSTL(cwd))//'/results'
      write (inc, "(I0.2)"),KINC
c
c ! --- Frames_description.dat
c ! --- Note: ???
c
      if (JELEM .EQ. 1) then
102   format(A,I7,A,f8.3)
      open(2000000,file = trim(folder)//'/Frames_description.dat')
               write(2000000,102)
     *             'Increment', 0,': Step Time =', 0.000
               write(2000000,102)
     *             'Increment', KINC,': Step Time =', TIME(1)
      close(2000000)
      end if
c
c !-----------------------------------------------------------------------------------
c ! --- Nodes.dat
c ! --- Note: ???
c
 103  format(A,I6,A,<Nodes_PP>(I9,A))
      write(3000000,103)'<',JELEM,'> ',(Nodes_in_IGFE(i),',',i=1,Nodes_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- Nodes_Coords.dat
c ! --- Note: ???
c
 104  format(A,I6,A,<Nodes_PP>(A,I9,A,f25.15,A,f25.15,A,f25.15,A))
      write(4000000,104)'<',JELEM,'> ',
     *                 ('(',Nodes_in_IGFE(i),
     *                   ',',NC_in_IGFE(i,1),
     *                   ',',NC_in_IGFE(i,2),
     *                   ',',0.0,'),',i=1,Nodes_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- Elements.dat
c ! --- Note: ???
c
 105  format(A,I6,A,<FE_PP>(I9,A))
      write(5000000,105)'<',JELEM,'> ',((JELEM-1)*FE_PP + i,',',i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- Elements_Nodes.dat
c ! --- Note: ???
c
      ij = 0
      node_shift = (JELEM-1)*Nodes_PP
      do j = 1,Num_FE_PP_V
         do i = 1,Num_FE_PP_U
            ij = ij + 1
            FE_Nodes(ij,1) = i + (j-1)*(Num_FE_PP_U+1) ! # of local node 1 within JELEM
            FE_Nodes(ij,2) = 1 + i + (j-1)*(Num_FE_PP_U+1) ! # of local node 2 within JELEM
            FE_Nodes(ij,3) = Num_FE_PP_U+2 + i + (j-1)*(Num_FE_PP_U+1) ! # of local node 3 within JELEM
            FE_Nodes(ij,4) = Num_FE_PP_U+1 + i + (j-1)*(Num_FE_PP_U+1) ! # of local node 4 within JELEM
            FE_Nodes(ij,5) = (JELEM-1)*FE_PP + ij ! PP FE global #
         end do
      end do
c
 106  format(A,I6,A,<FE_PP>(A,I9,A,I9,A,I9,A,I9,A,I9,A))
      write(6000000,106)'<',JELEM,'> ',
     *                 ('(',FE_Nodes(i,5),
     *                  ',',FE_Nodes(i,1) + node_shift, ! Global # of local node 1
     *                  ',',FE_Nodes(i,2) + node_shift, ! Global # of local node 2
     *                  ',',FE_Nodes(i,3) + node_shift, ! Global # of local node 3
     *                  ',',FE_Nodes(i,4) + node_shift, ! Global # of local node 4
     *                  '),',i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- U_Nodes.dat
c ! --- Note: ???
c
 107  format(A,I6,A,<Nodes_PP>(A,f25.15,A,f25.15,A,f25.15,A))
      write(7000000,107)'<',JELEM,'> ',
     *                 ('(',ND_in_IGFE(i,1),
     *                  ',',ND_in_IGFE(i,2),
     *                  ',',0.0,'),',i=1,Nodes_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- S_Nodes.dat and E_Nodes.dat
c ! --- Note: ???
c
 108  format(A,I6,A,<4*FE_PP>(A,f25.15,A,f25.15,A,f25.15,A,f25.15,A))
c
c ! --- S11 - s(1), S22 - s(2), S33 - s(3), S12 - s(4)
      write(8000000,108)'<',JELEM,'> ',
     *                 ('(',stress_NC(FE_Nodes(i,1),1),
     *                  ',',stress_NC(FE_Nodes(i,1),2),
     *                  ',',stress_NC(FE_Nodes(i,1),3),
     *                  ',',stress_NC(FE_Nodes(i,1),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,2),1),
     *                  ',',stress_NC(FE_Nodes(i,2),2),
     *                  ',',stress_NC(FE_Nodes(i,2),3),
     *                  ',',stress_NC(FE_Nodes(i,2),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,3),1),
     *                  ',',stress_NC(FE_Nodes(i,3),2),
     *                  ',',stress_NC(FE_Nodes(i,3),3),
     *                  ',',stress_NC(FE_Nodes(i,3),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,4),1),
     *                  ',',stress_NC(FE_Nodes(i,4),2),
     *                  ',',stress_NC(FE_Nodes(i,4),3),
     *                  ',',stress_NC(FE_Nodes(i,4),4),'),',
     *                      i=1,FE_PP)
c
c ! --- E11 - e(1), E22 - e(2), E33 - e(3), E12 - e(4)/2
      write(9000000,108)'<',JELEM,'> ',
     *                 ('(',strain_NC(FE_Nodes(i,1),1),
     *                  ',',strain_NC(FE_Nodes(i,1),2),
     *                  ',',strain_NC(FE_Nodes(i,1),3),
     *                  ',',strain_NC(FE_Nodes(i,1),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,2),1),
     *                  ',',strain_NC(FE_Nodes(i,2),2),
     *                  ',',strain_NC(FE_Nodes(i,2),3),
     *                  ',',strain_NC(FE_Nodes(i,2),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,3),1),
     *                  ',',strain_NC(FE_Nodes(i,3),2),
     *                  ',',strain_NC(FE_Nodes(i,3),3),
     *                  ',',strain_NC(FE_Nodes(i,3),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,4),1),
     *                  ',',strain_NC(FE_Nodes(i,4),2),
     *                  ',',strain_NC(FE_Nodes(i,4),3),
     *                  ',',strain_NC(FE_Nodes(i,4),4)/2.d0,'),',
     *                      i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c
      if (El_Output .EQ. 1) then
c
c !-----------------------------------------------------------------------------------
c ! --- gSx_Nodes.dat and gSy_Nodes.dat
c ! --- gEx_Nodes.dat and gEy_Nodes.dat
c ! --- Note: ???
c
 110  format(A,I6,A,<4*FE_PP>(A,f25.15,A,f25.15,A,f25.15,A,f25.15,A))
c
c ! --- gSx11 - gs(1)
c ! --- gSx22 - gs(3)
c ! --- gSx33 - gs(5)
c ! --- gSx12 - gs(7)
      write(10000000,110)'<',JELEM,'> ',
     *                 ('(',gstress_NC(FE_Nodes(i,1),1),
     *                  ',',gstress_NC(FE_Nodes(i,1),3),
     *                  ',',gstress_NC(FE_Nodes(i,1),5),
     *                  ',',gstress_NC(FE_Nodes(i,1),7),'),',
     *                  '(',gstress_NC(FE_Nodes(i,2),1),
     *                  ',',gstress_NC(FE_Nodes(i,2),3),
     *                  ',',gstress_NC(FE_Nodes(i,2),5),
     *                  ',',gstress_NC(FE_Nodes(i,2),7),'),',
     *                  '(',gstress_NC(FE_Nodes(i,3),1),
     *                  ',',gstress_NC(FE_Nodes(i,3),3),
     *                  ',',gstress_NC(FE_Nodes(i,3),5),
     *                  ',',gstress_NC(FE_Nodes(i,3),7),'),',
     *                  '(',gstress_NC(FE_Nodes(i,4),1),
     *                  ',',gstress_NC(FE_Nodes(i,4),3),
     *                  ',',gstress_NC(FE_Nodes(i,4),5),
     *                  ',',gstress_NC(FE_Nodes(i,4),7),'),',
     *                      i=1,FE_PP)
c
c ! --- gSy11 - gs(2)
c ! --- gSy22 - gs(4)
c ! --- gSy33 - gs(6)
c ! --- gSy12 - gs(8)
      write(11000000,110)'<',JELEM,'> ',
     *                 ('(',gstress_NC(FE_Nodes(i,1),2),
     *                  ',',gstress_NC(FE_Nodes(i,1),4),
     *                  ',',gstress_NC(FE_Nodes(i,1),6),
     *                  ',',gstress_NC(FE_Nodes(i,1),8),'),',
     *                  '(',gstress_NC(FE_Nodes(i,2),2),
     *                  ',',gstress_NC(FE_Nodes(i,2),4),
     *                  ',',gstress_NC(FE_Nodes(i,2),6),
     *                  ',',gstress_NC(FE_Nodes(i,2),8),'),',
     *                  '(',gstress_NC(FE_Nodes(i,3),2),
     *                  ',',gstress_NC(FE_Nodes(i,3),4),
     *                  ',',gstress_NC(FE_Nodes(i,3),6),
     *                  ',',gstress_NC(FE_Nodes(i,3),8),'),',
     *                  '(',gstress_NC(FE_Nodes(i,4),2),
     *                  ',',gstress_NC(FE_Nodes(i,4),4),
     *                  ',',gstress_NC(FE_Nodes(i,4),6),
     *                  ',',gstress_NC(FE_Nodes(i,4),8),'),',
     *                      i=1,FE_PP)
c
c ! --- gEx11 - ge(1)
c ! --- gEx22 - ge(3)
c ! --- gEx33 - ge(5)
c ! --- gEx12 - ge(7)/2
      write(12000000,110)'<',JELEM,'> ',
     *                 ('(',gstrain_NC(FE_Nodes(i,1),1),
     *                  ',',gstrain_NC(FE_Nodes(i,1),3),
     *                  ',',gstrain_NC(FE_Nodes(i,1),5),
     *                  ',',gstrain_NC(FE_Nodes(i,1),7)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,2),1),
     *                  ',',gstrain_NC(FE_Nodes(i,2),3),
     *                  ',',gstrain_NC(FE_Nodes(i,2),5),
     *                  ',',gstrain_NC(FE_Nodes(i,2),7)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,3),1),
     *                  ',',gstrain_NC(FE_Nodes(i,3),3),
     *                  ',',gstrain_NC(FE_Nodes(i,3),5),
     *                  ',',gstrain_NC(FE_Nodes(i,3),7)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,4),1),
     *                  ',',gstrain_NC(FE_Nodes(i,4),3),
     *                  ',',gstrain_NC(FE_Nodes(i,4),5),
     *                  ',',gstrain_NC(FE_Nodes(i,4),7)/2.d0,'),',
     *                      i=1,FE_PP)
c
c ! --- gEy11 - ge(2)
c ! --- gEy22 - ge(4)
c ! --- gEy33 - ge(6)
c ! --- gEy12 - ge(8)/2
      write(13000000,110)'<',JELEM,'> ',
     *                 ('(',gstrain_NC(FE_Nodes(i,1),2),
     *                  ',',gstrain_NC(FE_Nodes(i,1),4),
     *                  ',',gstrain_NC(FE_Nodes(i,1),6),
     *                  ',',gstrain_NC(FE_Nodes(i,1),8)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,2),2),
     *                  ',',gstrain_NC(FE_Nodes(i,2),4),
     *                  ',',gstrain_NC(FE_Nodes(i,2),6),
     *                  ',',gstrain_NC(FE_Nodes(i,2),8)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,3),2),
     *                  ',',gstrain_NC(FE_Nodes(i,3),4),
     *                  ',',gstrain_NC(FE_Nodes(i,3),6),
     *                  ',',gstrain_NC(FE_Nodes(i,3),8)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,4),2),
     *                  ',',gstrain_NC(FE_Nodes(i,4),4),
     *                  ',',gstrain_NC(FE_Nodes(i,4),6),
     *                  ',',gstrain_NC(FE_Nodes(i,4),8)/2.d0,'),',
     *                      i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c
      end if
c
      return
      end subroutine Output_2D