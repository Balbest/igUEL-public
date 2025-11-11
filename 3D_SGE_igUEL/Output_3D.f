      subroutine Output_3D(Kn_Num_U,Kn_Num_V,Kn_Num_W,KV_U,KV_V,KV_W,Ppol,Qpol,Rpol,
     1                     JELEM,NNODE,NOE,NOCPs,COORDS,U,NDOFEL,MCRD,TIME,
     2                     Num_FE_PP_U,Num_FE_PP_V,Num_FE_PP_W,C_St_Matr,A_St_Matr,
     3                     i_el,j_el,k_el,w_CPs,cwd,KINC,ndim,ndof,ntens,El_Output)
c
      include 'ABA_PARAM.INC'
c
      integer NOE, NOCPs, Num_FE_PP_U, Num_FE_PP_V, Num_FE_PP_W,
     *        Ppol, Qpol, Rpol, Kn_Num_U, Kn_Num_V, Kn_Num_W,
     *        i_el, j_el, k_el, El_Output
c
      double precision KV_U, KV_V, KV_W, w_CPs, COORDS
      dimension KV_U(Kn_Num_U), KV_V(Kn_Num_V), KV_W(Kn_Num_W),
     *          w_CPs(NOCPs), COORDS(MCRD,NNODE), U(NDOFEL), TIME(2)
c
      double precision h_xi, h_eta, h_zeta,
     *                 NCOORD_U, NCOORD_V, NCOORD_W
c
      double precision NC_in_IGFE,
     *                 ND_in_IGFE,
     *                 stress_NC,
     *                 strain_NC,
     *                 gstress_NC,
     *                 gstrain_NC
      dimension NC_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ndim),
     *          ND_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ndof),
     *          stress_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ntens),
     *          strain_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ntens),
     *          gstress_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ntens*ndim),
     *          gstrain_NC((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1),ntens*ndim)
c
      double precision R_NC, dR_NC, ddR_NC, Jdet
      dimension R_NC(NNODE), dR_NC(NNODE,3), ddR_NC(NNODE,6)
c
      double precision stran, gstran,
     *                 stress, gstress,
     *                 C_St_Matr, A_St_Matr
      dimension stran(ntens), gstran(ntens*ndim),
     *          stress(ntens), gstress(ntens*ndim),
     *          C_St_Matr(ntens,ntens), A_St_Matr(ntens*ndim,ntens*ndim)
c
      integer i,j,k,ijk,m,mijk,k1,i_u,i_v,i_w,i_uvw,node_shift,
     *        i1,i2,i3,i4,i5,i6,i7,i8,Nodes_PP,FE_PP
c
      integer Nodes_in_IGFE, FE_Nodes
      dimension Nodes_in_IGFE((Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1)),
     *          FE_Nodes(Num_FE_PP_U*Num_FE_PP_V*Num_FE_PP_W,9)
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
c !--- PP node coordinates: NCOORD_U, NCOORD_V, NCOORD_W
c
      h_xi   = 2.d0/Num_FE_PP_U
      h_eta  = 2.d0/Num_FE_PP_V
      h_zeta = 2.d0/Num_FE_PP_W
      Nodes_PP = (Num_FE_PP_U+1)*(Num_FE_PP_V+1)*(Num_FE_PP_W+1)
      FE_PP    =  Num_FE_PP_U*Num_FE_PP_V*Num_FE_PP_W
c
      i_uvw = 0
      do i_w = 1, Num_FE_PP_W+1 ! --- Start loop over PP nodes in W-direction
         NCOORD_W = -1.d0 + h_zeta*(i_w-1)
         NCOORD_W = (KV_W(k_el+1)+KV_W(k_el) + 
     *               NCOORD_W*(KV_W(k_el+1)-KV_W(k_el)))/2.d0
c
         do i_v = 1, Num_FE_PP_V+1 ! --- Start loop over PP nodes in V-direction
            NCOORD_V = -1.d0 + h_eta*(i_v-1)
            NCOORD_V = (KV_V(j_el+1)+KV_V(j_el) + 
     *                  NCOORD_V*(KV_V(j_el+1)-KV_V(j_el)))/2.d0
c
            do i_u = 1, Num_FE_PP_U+1 ! --- Start loop over PP nodes in U-direction
               NCOORD_U = -1.d0 + h_xi*(i_u-1)
               NCOORD_U = (KV_U(i_el+1)+KV_U(i_el) + 
     *                     NCOORD_U*(KV_U(i_el+1)-KV_U(i_el)))/2.d0
c
               i_uvw = i_uvw + 1
c
               call NURBS_BF_3D(Ppol,i_el,NCOORD_U,KV_U,Kn_Num_U,
     1                          Qpol,j_el,NCOORD_V,KV_V,Kn_Num_V, 
     2                          Rpol,k_el,NCOORD_W,KV_W,Kn_Num_W, 
     3                          w_CPs,NOCPs,NNODE,COORDS,MCRD,
     4                          R_NC,dR_NC,ddR_NC,Jdet,JELEM)
c
c
c      if (JELEM .EQ. 9) then
c      write(7,*)'JELEM =', JELEM
c      write(7,*)'ddR_NC =', ddR_NC
c      end if
c
c !--- Nodes_in_IGFE: (global) numbers of PP nodes in current knot span (JELEM)
c !--- NC_in_IGFE   : (global) coordinates of PP nodes in current knot span (JELEM)
c !--- ND_in_IGFE   : (global) displacements of PP nodes in current knot span (JELEM)
c
               Nodes_in_IGFE(i_uvw) = (JELEM-1)*Nodes_PP + i_uvw
c
               NC_in_IGFE(i_uvw,:) = 0.d0
               ND_in_IGFE(i_uvw,:) = 0.d0
c
               do m = 1, ndof
				  ijk = 0
                  do k = 1, Rpol+1
                     do j = 1, Qpol+1
                        do i = 1, Ppol+1
                           ijk = ijk + 1
                           mijk = m + (ijk-1)*ndof
                           if (m .le. ndim) then 
                               NC_in_IGFE(i_uvw,m) = NC_in_IGFE(i_uvw,m) +
     *                                               COORDS(m,ijk)*R_NC(ijk)
                           end if
                               ND_in_IGFE(i_uvw,m) = ND_in_IGFE(i_uvw,m) +
     *                                               U(mijk)*R_NC(ijk)
                        end do
                     end do
                  end do
               end do
c
c !--- stress_NC: stresses at PP nodes in current knot span (JELEM)
c !--- strain_NC: strains at PP nodes in current knot span (JELEM)
c !--- gstress_NC: higher-order stresses at PP nodes in current knot span (JELEM)
c !--- gstrain_NC: gradient of strains at PP nodes in current knot span (JELEM)
c
               call Strains_3D(ntens,ndof,ndim,NNODE,NDOFEL,U,
     1                         stran,gstran,dR_NC,ddR_NC)
c
               stress  = matmul(C_St_Matr,stran)
               gstress = matmul(A_St_Matr,gstran)
c
               do k1 = 1, ntens
                  stress_NC(i_uvw,k1) = stress(k1)
                  strain_NC(i_uvw,k1) = stran(k1)
c
                  gstress_NC(i_uvw,k1) = gstress(k1)
                  gstrain_NC(i_uvw,k1) = gstran(k1)
                  gstress_NC(i_uvw,k1+ntens) = gstress(k1+ntens)
                  gstrain_NC(i_uvw,k1+ntens) = gstran(k1+ntens)
                  gstress_NC(i_uvw,k1+2*ntens) = gstress(k1+2*ntens)
                  gstrain_NC(i_uvw,k1+2*ntens) = gstran(k1+2*ntens)
               end do
c
c
c      if (JELEM .EQ. 9) then
c      write(7,*)'JELEM =', JELEM
c      write(7,*)'dR_NC =', dR_NC
c      write(7,*)'ddR_NC =', ddR_NC
c      write(7,*)'PPPs, gstran =', gstran
c      write(7,*)'PPPs, gstrain_NC =', gstrain_NC
c      end if
c
c
            end do ! --- End loop over PP nodes in U-direction
         end do ! --- End loop over PP nodes in V-direction
      end do ! --- End loop over PP nodes in W-direction
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
     *                   ',',NC_in_IGFE(i,3),'),',i=1,Nodes_PP)
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
      ijk = 0
      node_shift = (JELEM-1)*Nodes_PP
      do k = 1,Num_FE_PP_W
         do j = 1,Num_FE_PP_V
            do i = 1,Num_FE_PP_U
               ijk = ijk + 1
               FE_Nodes(ijk,1) = i + (j-1)*(Num_FE_PP_U+1)
     *                         +(k-1)*(Num_FE_PP_U+1)*(Num_FE_PP_V+1) ! # of local node 1 within JELEM
               FE_Nodes(ijk,2) = 1 + i + (j-1)*(Num_FE_PP_U+1)
     *                         +(k-1)*(Num_FE_PP_U+1)*(Num_FE_PP_V+1) ! # of local node 2 within JELEM
               FE_Nodes(ijk,3) = Num_FE_PP_U + 2 + i + (j-1)*(Num_FE_PP_U+1)
     *                         +(k-1)*(Num_FE_PP_U+1)*(Num_FE_PP_V+1) ! # of local node 3 within JELEM
               FE_Nodes(ijk,4) = Num_FE_PP_U + 1 + i + (j-1)*(Num_FE_PP_U+1)
     *                         +(k-1)*(Num_FE_PP_U+1)*(Num_FE_PP_V+1) ! # of local node 4 within JELEM
               FE_Nodes(ijk,5) = i + (j-1)*(Num_FE_PP_U+1)
     *                         + k*(Num_FE_PP_U+1)*(Num_FE_PP_V+1)    ! # of local node 5 within JELEM
               FE_Nodes(ijk,6) = 1 + i + (j-1)*(Num_FE_PP_U+1)
     *                         + k*(Num_FE_PP_U+1)*(Num_FE_PP_V+1)    ! # of local node 6 within JELEM
               FE_Nodes(ijk,7) = Num_FE_PP_U + 2 + i + (j-1)*(Num_FE_PP_U+1)
     *                         + k*(Num_FE_PP_U+1)*(Num_FE_PP_V+1)    ! # of local node 7 within JELEM
               FE_Nodes(ijk,8) = Num_FE_PP_U + 1 + i + (j-1)*(Num_FE_PP_U+1)
     *                         + k*(Num_FE_PP_U+1)*(Num_FE_PP_V+1)    ! # of local node 8 within JELEM
               FE_Nodes(ijk,9) = (JELEM-1)*FE_PP + ijk                ! PP FE global #
            end do
         end do
      end do
c
 106  format(A,I6,A,<FE_PP>(A,I9,A,I9,A,I9,A,I9,A,I9,A,I9,A,I9,A,I9,A,I9,A))
      write(6000000,106)'<',JELEM,'> ',
     *                 ('(',FE_Nodes(i,9),
     *                  ',',FE_Nodes(i,1) + node_shift, ! Global # of local node 1
     *                  ',',FE_Nodes(i,2) + node_shift, ! Global # of local node 2
     *                  ',',FE_Nodes(i,3) + node_shift, ! Global # of local node 3
     *                  ',',FE_Nodes(i,4) + node_shift, ! Global # of local node 4
     *                  ',',FE_Nodes(i,5) + node_shift, ! Global # of local node 5
     *                  ',',FE_Nodes(i,6) + node_shift, ! Global # of local node 6
     *                  ',',FE_Nodes(i,7) + node_shift, ! Global # of local node 7
     *                  ',',FE_Nodes(i,8) + node_shift, ! Global # of local node 8
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
     *                  ',',ND_in_IGFE(i,3),'),',i=1,Nodes_PP)
c
c !-----------------------------------------------------------------------------------
c ! --- S_Nodes.dat and E_Nodes.dat
c ! --- Note: ???
c
 108  format(A,I6,A,<8*FE_PP>(A,f25.15,A,f25.15,A,f25.15,A,f25.15,A,f25.15,A,f25.15,A))
c
c! --- S11 - s(1), S22 - s(2), S33 - s(3), S12 - s(6), S13 - s(5), S23 - s(4)
      write(8000000,108)'<',JELEM,'> ',
     *                 ('(',stress_NC(FE_Nodes(i,1),1),
     *                  ',',stress_NC(FE_Nodes(i,1),2),
     *                  ',',stress_NC(FE_Nodes(i,1),3),
     *                  ',',stress_NC(FE_Nodes(i,1),6),
     *                  ',',stress_NC(FE_Nodes(i,1),5),
     *                  ',',stress_NC(FE_Nodes(i,1),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,2),1),
     *                  ',',stress_NC(FE_Nodes(i,2),2),
     *                  ',',stress_NC(FE_Nodes(i,2),3),
     *                  ',',stress_NC(FE_Nodes(i,2),6),
     *                  ',',stress_NC(FE_Nodes(i,2),5),
     *                  ',',stress_NC(FE_Nodes(i,2),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,3),1),
     *                  ',',stress_NC(FE_Nodes(i,3),2),
     *                  ',',stress_NC(FE_Nodes(i,3),3),
     *                  ',',stress_NC(FE_Nodes(i,3),6),
     *                  ',',stress_NC(FE_Nodes(i,3),5),
     *                  ',',stress_NC(FE_Nodes(i,3),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,4),1),
     *                  ',',stress_NC(FE_Nodes(i,4),2),
     *                  ',',stress_NC(FE_Nodes(i,4),3),
     *                  ',',stress_NC(FE_Nodes(i,4),6),
     *                  ',',stress_NC(FE_Nodes(i,4),5),
     *                  ',',stress_NC(FE_Nodes(i,4),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,5),1),
     *                  ',',stress_NC(FE_Nodes(i,5),2),
     *                  ',',stress_NC(FE_Nodes(i,5),3),
     *                  ',',stress_NC(FE_Nodes(i,5),6),
     *                  ',',stress_NC(FE_Nodes(i,5),5),
     *                  ',',stress_NC(FE_Nodes(i,5),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,6),1),
     *                  ',',stress_NC(FE_Nodes(i,6),2),
     *                  ',',stress_NC(FE_Nodes(i,6),3),
     *                  ',',stress_NC(FE_Nodes(i,6),6),
     *                  ',',stress_NC(FE_Nodes(i,6),5),
     *                  ',',stress_NC(FE_Nodes(i,6),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,7),1),
     *                  ',',stress_NC(FE_Nodes(i,7),2),
     *                  ',',stress_NC(FE_Nodes(i,7),3),
     *                  ',',stress_NC(FE_Nodes(i,7),6),
     *                  ',',stress_NC(FE_Nodes(i,7),5),
     *                  ',',stress_NC(FE_Nodes(i,7),4),'),',
     *                  '(',stress_NC(FE_Nodes(i,8),1),
     *                  ',',stress_NC(FE_Nodes(i,8),2),
     *                  ',',stress_NC(FE_Nodes(i,8),3),
     *                  ',',stress_NC(FE_Nodes(i,8),6),
     *                  ',',stress_NC(FE_Nodes(i,8),5),
     *                  ',',stress_NC(FE_Nodes(i,8),4),'),',
     *                      i=1,FE_PP)
c
c! --- E11 - e(1), E22 - e(2), E33 - e(3), E12 - e(6)/2, E13 - e(5)/2, E23 - e(4)/2
      write(9000000,108)'<',JELEM,'> ',
     *                 ('(',strain_NC(FE_Nodes(i,1),1),
     *                  ',',strain_NC(FE_Nodes(i,1),2),
     *                  ',',strain_NC(FE_Nodes(i,1),3),
     *                  ',',strain_NC(FE_Nodes(i,1),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,1),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,1),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,2),1),
     *                  ',',strain_NC(FE_Nodes(i,2),2),
     *                  ',',strain_NC(FE_Nodes(i,2),3),
     *                  ',',strain_NC(FE_Nodes(i,2),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,2),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,2),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,3),1),
     *                  ',',strain_NC(FE_Nodes(i,3),2),
     *                  ',',strain_NC(FE_Nodes(i,3),3),
     *                  ',',strain_NC(FE_Nodes(i,3),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,3),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,3),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,4),1),
     *                  ',',strain_NC(FE_Nodes(i,4),2),
     *                  ',',strain_NC(FE_Nodes(i,4),3),
     *                  ',',strain_NC(FE_Nodes(i,4),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,4),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,4),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,5),1),
     *                  ',',strain_NC(FE_Nodes(i,5),2),
     *                  ',',strain_NC(FE_Nodes(i,5),3),
     *                  ',',strain_NC(FE_Nodes(i,5),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,5),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,5),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,6),1),
     *                  ',',strain_NC(FE_Nodes(i,6),2),
     *                  ',',strain_NC(FE_Nodes(i,6),3),
     *                  ',',strain_NC(FE_Nodes(i,6),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,6),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,6),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,7),1),
     *                  ',',strain_NC(FE_Nodes(i,7),2),
     *                  ',',strain_NC(FE_Nodes(i,7),3),
     *                  ',',strain_NC(FE_Nodes(i,7),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,7),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,7),4)/2.d0,'),',
     *                  '(',strain_NC(FE_Nodes(i,8),1),
     *                  ',',strain_NC(FE_Nodes(i,8),2),
     *                  ',',strain_NC(FE_Nodes(i,8),3),
     *                  ',',strain_NC(FE_Nodes(i,8),6)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,8),5)/2.d0,
     *                  ',',strain_NC(FE_Nodes(i,8),4)/2.d0,'),',
     *                      i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c
      if (El_Output .EQ. 1) then
c
c !-----------------------------------------------------------------------------------
c ! --- gSx_Nodes.dat, gSy_Nodes.dat and gSz_Nodes.dat
c ! --- gEx_Nodes.dat, gEy_Nodes.dat and gEz_Nodes.dat
c ! --- Note: ???
c
 110  format(A,I6,A,<8*FE_PP>(A,f25.15,A,f25.15,A,f25.15,A,f25.15,A,f25.15,A,f25.15,A))
c
c ! --- gSx11 - gs(1)
c ! --- gSx22 - gs(4)
c ! --- gSx33 - gs(7)
c ! --- gSx12 - gs(16)
c ! --- gSx13 - gs(13)
c ! --- gSx23 - gs(10)
      write(10000000,110)'<',JELEM,'> ',
     *                 ('(',gstress_NC(FE_Nodes(i,1),1),
     *                  ',',gstress_NC(FE_Nodes(i,1),4),
     *                  ',',gstress_NC(FE_Nodes(i,1),7),
     *                  ',',gstress_NC(FE_Nodes(i,1),16),
     *                  ',',gstress_NC(FE_Nodes(i,1),13),
     *                  ',',gstress_NC(FE_Nodes(i,1),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,2),1),
     *                  ',',gstress_NC(FE_Nodes(i,2),4),
     *                  ',',gstress_NC(FE_Nodes(i,2),7),
     *                  ',',gstress_NC(FE_Nodes(i,2),16),
     *                  ',',gstress_NC(FE_Nodes(i,2),13),
     *                  ',',gstress_NC(FE_Nodes(i,2),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,3),1),
     *                  ',',gstress_NC(FE_Nodes(i,3),4),
     *                  ',',gstress_NC(FE_Nodes(i,3),7),
     *                  ',',gstress_NC(FE_Nodes(i,3),16),
     *                  ',',gstress_NC(FE_Nodes(i,3),13),
     *                  ',',gstress_NC(FE_Nodes(i,3),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,4),1),
     *                  ',',gstress_NC(FE_Nodes(i,4),4),
     *                  ',',gstress_NC(FE_Nodes(i,4),7),
     *                  ',',gstress_NC(FE_Nodes(i,4),16),
     *                  ',',gstress_NC(FE_Nodes(i,4),13),
     *                  ',',gstress_NC(FE_Nodes(i,4),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,5),1),
     *                  ',',gstress_NC(FE_Nodes(i,5),4),
     *                  ',',gstress_NC(FE_Nodes(i,5),7),
     *                  ',',gstress_NC(FE_Nodes(i,5),16),
     *                  ',',gstress_NC(FE_Nodes(i,5),13),
     *                  ',',gstress_NC(FE_Nodes(i,5),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,6),1),
     *                  ',',gstress_NC(FE_Nodes(i,6),4),
     *                  ',',gstress_NC(FE_Nodes(i,6),7),
     *                  ',',gstress_NC(FE_Nodes(i,6),16),
     *                  ',',gstress_NC(FE_Nodes(i,6),13),
     *                  ',',gstress_NC(FE_Nodes(i,6),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,7),1),
     *                  ',',gstress_NC(FE_Nodes(i,7),4),
     *                  ',',gstress_NC(FE_Nodes(i,7),7),
     *                  ',',gstress_NC(FE_Nodes(i,7),16),
     *                  ',',gstress_NC(FE_Nodes(i,7),13),
     *                  ',',gstress_NC(FE_Nodes(i,7),10),'),',
     *                  '(',gstress_NC(FE_Nodes(i,8),1),
     *                  ',',gstress_NC(FE_Nodes(i,8),4),
     *                  ',',gstress_NC(FE_Nodes(i,8),7),
     *                  ',',gstress_NC(FE_Nodes(i,8),16),
     *                  ',',gstress_NC(FE_Nodes(i,8),13),
     *                  ',',gstress_NC(FE_Nodes(i,8),10),'),',
     *                      i=1,FE_PP)
c
c ! --- gSy11 - gs(2)
c ! --- gSy22 - gs(5)
c ! --- gSy33 - gs(8)
c ! --- gSy12 - gs(17)
c ! --- gSy13 - gs(14)
c ! --- gSy23 - gs(11)
      write(11000000,110)'<',JELEM,'> ',
     *                 ('(',gstress_NC(FE_Nodes(i,1),2),
     *                  ',',gstress_NC(FE_Nodes(i,1),5),
     *                  ',',gstress_NC(FE_Nodes(i,1),8),
     *                  ',',gstress_NC(FE_Nodes(i,1),17),
     *                  ',',gstress_NC(FE_Nodes(i,1),14),
     *                  ',',gstress_NC(FE_Nodes(i,1),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,2),2),
     *                  ',',gstress_NC(FE_Nodes(i,2),5),
     *                  ',',gstress_NC(FE_Nodes(i,2),8),
     *                  ',',gstress_NC(FE_Nodes(i,2),17),
     *                  ',',gstress_NC(FE_Nodes(i,2),14),
     *                  ',',gstress_NC(FE_Nodes(i,2),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,3),2),
     *                  ',',gstress_NC(FE_Nodes(i,3),5),
     *                  ',',gstress_NC(FE_Nodes(i,3),8),
     *                  ',',gstress_NC(FE_Nodes(i,3),17),
     *                  ',',gstress_NC(FE_Nodes(i,3),14),
     *                  ',',gstress_NC(FE_Nodes(i,3),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,4),2),
     *                  ',',gstress_NC(FE_Nodes(i,4),5),
     *                  ',',gstress_NC(FE_Nodes(i,4),8),
     *                  ',',gstress_NC(FE_Nodes(i,4),17),
     *                  ',',gstress_NC(FE_Nodes(i,4),14),
     *                  ',',gstress_NC(FE_Nodes(i,4),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,5),2),
     *                  ',',gstress_NC(FE_Nodes(i,5),5),
     *                  ',',gstress_NC(FE_Nodes(i,5),8),
     *                  ',',gstress_NC(FE_Nodes(i,5),17),
     *                  ',',gstress_NC(FE_Nodes(i,5),14),
     *                  ',',gstress_NC(FE_Nodes(i,5),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,6),2),
     *                  ',',gstress_NC(FE_Nodes(i,6),5),
     *                  ',',gstress_NC(FE_Nodes(i,6),8),
     *                  ',',gstress_NC(FE_Nodes(i,6),17),
     *                  ',',gstress_NC(FE_Nodes(i,6),14),
     *                  ',',gstress_NC(FE_Nodes(i,6),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,7),2),
     *                  ',',gstress_NC(FE_Nodes(i,7),5),
     *                  ',',gstress_NC(FE_Nodes(i,7),8),
     *                  ',',gstress_NC(FE_Nodes(i,7),17),
     *                  ',',gstress_NC(FE_Nodes(i,7),14),
     *                  ',',gstress_NC(FE_Nodes(i,7),11),'),',
     *                  '(',gstress_NC(FE_Nodes(i,8),2),
     *                  ',',gstress_NC(FE_Nodes(i,8),5),
     *                  ',',gstress_NC(FE_Nodes(i,8),8),
     *                  ',',gstress_NC(FE_Nodes(i,8),17),
     *                  ',',gstress_NC(FE_Nodes(i,8),14),
     *                  ',',gstress_NC(FE_Nodes(i,8),11),'),',
     *                      i=1,FE_PP)
c
c ! --- gSz11 - gs(3)
c ! --- gSz22 - gs(6)
c ! --- gSz33 - gs(9)
c ! --- gSz12 - gs(18)
c ! --- gSz13 - gs(15)
c ! --- gSz23 - gs(12)
      write(12000000,110)'<',JELEM,'> ',
     *                 ('(',gstress_NC(FE_Nodes(i,1),3),
     *                  ',',gstress_NC(FE_Nodes(i,1),6),
     *                  ',',gstress_NC(FE_Nodes(i,1),9),
     *                  ',',gstress_NC(FE_Nodes(i,1),18),
     *                  ',',gstress_NC(FE_Nodes(i,1),15),
     *                  ',',gstress_NC(FE_Nodes(i,1),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,2),3),
     *                  ',',gstress_NC(FE_Nodes(i,2),6),
     *                  ',',gstress_NC(FE_Nodes(i,2),9),
     *                  ',',gstress_NC(FE_Nodes(i,2),18),
     *                  ',',gstress_NC(FE_Nodes(i,2),15),
     *                  ',',gstress_NC(FE_Nodes(i,2),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,3),3),
     *                  ',',gstress_NC(FE_Nodes(i,3),6),
     *                  ',',gstress_NC(FE_Nodes(i,3),9),
     *                  ',',gstress_NC(FE_Nodes(i,3),18),
     *                  ',',gstress_NC(FE_Nodes(i,3),15),
     *                  ',',gstress_NC(FE_Nodes(i,3),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,4),3),
     *                  ',',gstress_NC(FE_Nodes(i,4),6),
     *                  ',',gstress_NC(FE_Nodes(i,4),9),
     *                  ',',gstress_NC(FE_Nodes(i,4),18),
     *                  ',',gstress_NC(FE_Nodes(i,4),15),
     *                  ',',gstress_NC(FE_Nodes(i,4),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,5),3),
     *                  ',',gstress_NC(FE_Nodes(i,5),6),
     *                  ',',gstress_NC(FE_Nodes(i,5),9),
     *                  ',',gstress_NC(FE_Nodes(i,5),18),
     *                  ',',gstress_NC(FE_Nodes(i,5),15),
     *                  ',',gstress_NC(FE_Nodes(i,5),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,6),3),
     *                  ',',gstress_NC(FE_Nodes(i,6),6),
     *                  ',',gstress_NC(FE_Nodes(i,6),9),
     *                  ',',gstress_NC(FE_Nodes(i,6),18),
     *                  ',',gstress_NC(FE_Nodes(i,6),15),
     *                  ',',gstress_NC(FE_Nodes(i,6),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,7),3),
     *                  ',',gstress_NC(FE_Nodes(i,7),6),
     *                  ',',gstress_NC(FE_Nodes(i,7),9),
     *                  ',',gstress_NC(FE_Nodes(i,7),18),
     *                  ',',gstress_NC(FE_Nodes(i,7),15),
     *                  ',',gstress_NC(FE_Nodes(i,7),12),'),',
     *                  '(',gstress_NC(FE_Nodes(i,8),3),
     *                  ',',gstress_NC(FE_Nodes(i,8),6),
     *                  ',',gstress_NC(FE_Nodes(i,8),9),
     *                  ',',gstress_NC(FE_Nodes(i,8),18),
     *                  ',',gstress_NC(FE_Nodes(i,8),15),
     *                  ',',gstress_NC(FE_Nodes(i,8),12),'),',
     *                      i=1,FE_PP)
c
c ! --- gEx11 - ge(1)
c ! --- gEx22 - ge(4)
c ! --- gEx33 - ge(7)
c ! --- gEx12 - ge(16)/2
c ! --- gEx13 - ge(13)/2
c ! --- gEx23 - ge(10)/2
      write(13000000,110)'<',JELEM,'> ',
     *                 ('(',gstrain_NC(FE_Nodes(i,1),1),
     *                  ',',gstrain_NC(FE_Nodes(i,1),4),
     *                  ',',gstrain_NC(FE_Nodes(i,1),7),
     *                  ',',gstrain_NC(FE_Nodes(i,1),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,2),1),
     *                  ',',gstrain_NC(FE_Nodes(i,2),4),
     *                  ',',gstrain_NC(FE_Nodes(i,2),7),
     *                  ',',gstrain_NC(FE_Nodes(i,2),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,3),1),
     *                  ',',gstrain_NC(FE_Nodes(i,3),4),
     *                  ',',gstrain_NC(FE_Nodes(i,3),7),
     *                  ',',gstrain_NC(FE_Nodes(i,3),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,4),1),
     *                  ',',gstrain_NC(FE_Nodes(i,4),4),
     *                  ',',gstrain_NC(FE_Nodes(i,4),7),
     *                  ',',gstrain_NC(FE_Nodes(i,4),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,5),1),
     *                  ',',gstrain_NC(FE_Nodes(i,5),4),
     *                  ',',gstrain_NC(FE_Nodes(i,5),7),
     *                  ',',gstrain_NC(FE_Nodes(i,5),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,6),1),
     *                  ',',gstrain_NC(FE_Nodes(i,6),4),
     *                  ',',gstrain_NC(FE_Nodes(i,6),7),
     *                  ',',gstrain_NC(FE_Nodes(i,6),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,7),1),
     *                  ',',gstrain_NC(FE_Nodes(i,7),4),
     *                  ',',gstrain_NC(FE_Nodes(i,7),7),
     *                  ',',gstrain_NC(FE_Nodes(i,7),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),10)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,8),1),
     *                  ',',gstrain_NC(FE_Nodes(i,8),4),
     *                  ',',gstrain_NC(FE_Nodes(i,8),7),
     *                  ',',gstrain_NC(FE_Nodes(i,8),16)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),13)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),10)/2.d0,'),',
     *                      i=1,FE_PP)
c
c ! --- gEy11 - ge(2)
c ! --- gEy22 - ge(5)
c ! --- gEy33 - ge(8)
c ! --- gEy12 - ge(17)/2
c ! --- gEy13 - ge(14)/2
c ! --- gEy23 - ge(11)/2
      write(14000000,110)'<',JELEM,'> ',
     *                 ('(',gstrain_NC(FE_Nodes(i,1),2),
     *                  ',',gstrain_NC(FE_Nodes(i,1),5),
     *                  ',',gstrain_NC(FE_Nodes(i,1),8),
     *                  ',',gstrain_NC(FE_Nodes(i,1),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,2),2),
     *                  ',',gstrain_NC(FE_Nodes(i,2),5),
     *                  ',',gstrain_NC(FE_Nodes(i,2),8),
     *                  ',',gstrain_NC(FE_Nodes(i,2),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,3),2),
     *                  ',',gstrain_NC(FE_Nodes(i,3),5),
     *                  ',',gstrain_NC(FE_Nodes(i,3),8),
     *                  ',',gstrain_NC(FE_Nodes(i,3),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,4),2),
     *                  ',',gstrain_NC(FE_Nodes(i,4),5),
     *                  ',',gstrain_NC(FE_Nodes(i,4),8),
     *                  ',',gstrain_NC(FE_Nodes(i,4),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,5),2),
     *                  ',',gstrain_NC(FE_Nodes(i,5),5),
     *                  ',',gstrain_NC(FE_Nodes(i,5),8),
     *                  ',',gstrain_NC(FE_Nodes(i,5),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,6),2),
     *                  ',',gstrain_NC(FE_Nodes(i,6),5),
     *                  ',',gstrain_NC(FE_Nodes(i,6),8),
     *                  ',',gstrain_NC(FE_Nodes(i,6),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,7),2),
     *                  ',',gstrain_NC(FE_Nodes(i,7),5),
     *                  ',',gstrain_NC(FE_Nodes(i,7),8),
     *                  ',',gstrain_NC(FE_Nodes(i,7),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),11)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,8),2),
     *                  ',',gstrain_NC(FE_Nodes(i,8),5),
     *                  ',',gstrain_NC(FE_Nodes(i,8),8),
     *                  ',',gstrain_NC(FE_Nodes(i,8),17)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),14)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),11)/2.d0,'),',
     *                      i=1,FE_PP)
c
c ! --- gEz11 - ge(3)
c ! --- gEz22 - ge(6)
c ! --- gEz33 - ge(9)
c ! --- gEz12 - ge(18)/2
c ! --- gEz13 - ge(15)/2
c ! --- gEz23 - ge(12)/2
      write(15000000,110)'<',JELEM,'> ',
     *                 ('(',gstrain_NC(FE_Nodes(i,1),3),
     *                  ',',gstrain_NC(FE_Nodes(i,1),6),
     *                  ',',gstrain_NC(FE_Nodes(i,1),9),
     *                  ',',gstrain_NC(FE_Nodes(i,1),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,1),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,2),3),
     *                  ',',gstrain_NC(FE_Nodes(i,2),6),
     *                  ',',gstrain_NC(FE_Nodes(i,2),9),
     *                  ',',gstrain_NC(FE_Nodes(i,2),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,2),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,3),3),
     *                  ',',gstrain_NC(FE_Nodes(i,3),6),
     *                  ',',gstrain_NC(FE_Nodes(i,3),9),
     *                  ',',gstrain_NC(FE_Nodes(i,3),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,3),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,4),3),
     *                  ',',gstrain_NC(FE_Nodes(i,4),6),
     *                  ',',gstrain_NC(FE_Nodes(i,4),9),
     *                  ',',gstrain_NC(FE_Nodes(i,4),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,4),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,5),3),
     *                  ',',gstrain_NC(FE_Nodes(i,5),6),
     *                  ',',gstrain_NC(FE_Nodes(i,5),9),
     *                  ',',gstrain_NC(FE_Nodes(i,5),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,5),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,6),3),
     *                  ',',gstrain_NC(FE_Nodes(i,6),6),
     *                  ',',gstrain_NC(FE_Nodes(i,6),9),
     *                  ',',gstrain_NC(FE_Nodes(i,6),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,6),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,7),3),
     *                  ',',gstrain_NC(FE_Nodes(i,7),6),
     *                  ',',gstrain_NC(FE_Nodes(i,7),9),
     *                  ',',gstrain_NC(FE_Nodes(i,7),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,7),12)/2.d0,'),',
     *                  '(',gstrain_NC(FE_Nodes(i,8),3),
     *                  ',',gstrain_NC(FE_Nodes(i,8),6),
     *                  ',',gstrain_NC(FE_Nodes(i,8),9),
     *                  ',',gstrain_NC(FE_Nodes(i,8),18)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),15)/2.d0,
     *                  ',',gstrain_NC(FE_Nodes(i,8),12)/2.d0,'),',
     *                      i=1,FE_PP)
c
c !-----------------------------------------------------------------------------------
c
      end if
c
      return
      end subroutine Output_3D