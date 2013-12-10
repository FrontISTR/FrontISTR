!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

!C
!C***
!C*** hecmw_adapt_NEW_CELL_351
!C***
!C
!C    create new PRISMs
!C
      subroutine hecmw_adapt_NEW_CELL_351 (hecMESH, icouN)

      use hecmw_util

      implicit REAL*8 (A-H,O-Z)
      dimension NDIV(6), NNTYP(0:11)

      integer(kind=kint)        :: PAR_CEL_TYP
      type (hecmwST_local_mesh) :: hecMESH

!C
!C +----------------+
!C | embedding TET. |
!C +----------------+
!C===
      do i= 0, 11
        NNTYP(i)= 0
      enddo

      do icel0= 1, hecMESH%n_adapt_act_elem_351
        icel_par= hecMESH%adapt_act_elem_351(icel0)
            NPAR= icel_par

        if (hecMESH%elem_ID(2*icel_par).eq.hecMESH%my_rank) then
          inc= 1          
         else
          inc= 0
        endif

         iS= hecMESH%elem_node_index(icel_par-1)
        n01= hecMESH%elem_node_item (iS+1)
        n02= hecMESH%elem_node_item (iS+2)
        n03= hecMESH%elem_node_item (iS+3)
        n11= hecMESH%elem_node_item (iS+4)
        n12= hecMESH%elem_node_item (iS+5)
        n13= hecMESH%elem_node_item (iS+6)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n01, n02, ie01, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n02, n03, ie02, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n03, n01, ie03, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n11, n12, ie11, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n12, n13, ie12, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n13, n11, ie13, 1 )

        NDIV(1)= hecMESH%adapt_iemb(ie01)
        NDIV(2)= hecMESH%adapt_iemb(ie02)
        NDIV(3)= hecMESH%adapt_iemb(ie03)
        NDIV(4)= hecMESH%adapt_iemb(ie11)
        NDIV(5)= hecMESH%adapt_iemb(ie12)
        NDIV(6)= hecMESH%adapt_iemb(ie13)

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)
        PAR_CEL_TYP= hecMESH%elem_type(NPAR)

!C
!C-- init. CHILD. cell array
        iS = hecMESH%adapt_children_index(NPAR-1)
        iS1= iS + 1
        iS2= iS + 2
        iS3= iS + 3
        iS4= iS + 4
        iS5= iS + 5
        iS6= iS + 6
        iS7= iS + 7
        iS8= iS + 8

        hecMESH%adapt_children_item(2*iS1-1)= 0
        hecMESH%adapt_children_item(2*iS2-1)= 0
        hecMESH%adapt_children_item(2*iS3-1)= 0
        hecMESH%adapt_children_item(2*iS4-1)= 0
        hecMESH%adapt_children_item(2*iS5-1)= 0
        hecMESH%adapt_children_item(2*iS6-1)= 0
        hecMESH%adapt_children_item(2*iS7-1)= 0
        hecMESH%adapt_children_item(2*iS8-1)= 0

        hecMESH%adapt_children_item(2*iS1)= -1
        hecMESH%adapt_children_item(2*iS2)= -1
        hecMESH%adapt_children_item(2*iS3)= -1
        hecMESH%adapt_children_item(2*iS4)= -1
        hecMESH%adapt_children_item(2*iS5)= -1
        hecMESH%adapt_children_item(2*iS6)= -1
        hecMESH%adapt_children_item(2*iS7)= -1
        hecMESH%adapt_children_item(2*iS8)= -1

!C
!C== embedding TYPE
        if (NDIVSUM.eq.0) NTYP= 0
!C
!C-- TYP-1
        if (NDIVSUM.eq.2 .and. NDIV(1).eq.1 .and. NDIV(4).eq.1) then
          NTYP= 1
          n04 = hecMESH%adapt_IWK(ie01)
          n14 = hecMESH%adapt_IWK(ie11)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_PRISM (n01, n04, n03, n11, n14, n13, 1)
          call hecmw_adapt_CREATE_NEW_PRISM (n04, n02, n03, n14, n12, n13, 2)
        endif
!C
!C-- TYP-2
        if (NDIVSUM.eq.2 .and. NDIV(2).eq.1 .and. NDIV(5).eq.1) then
          NTYP= 2
          n04 = hecMESH%adapt_IWK(ie02)
          n14 = hecMESH%adapt_IWK(ie12)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_PRISM (n01, n04, n03, n11, n14, n13, 1)
          call hecmw_adapt_CREATE_NEW_PRISM (n01, n02, n04, n11, n12, n14, 2)
        endif
!C
!C-- TYP-3
        if (NDIVSUM.eq.2 .and. NDIV(3).eq.1 .and. NDIV(6).eq.1) then
          NTYP= 3
          n04 = hecMESH%adapt_IWK(ie03)
          n14 = hecMESH%adapt_IWK(ie13)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_PRISM (n01, n02, n04, n11, n12, n14, 1)
          call hecmw_adapt_CREATE_NEW_PRISM (n04, n02, n03, n14, n12, n13, 2)
        endif
!C
!C-- TYP-4
        if (NDIVSUM.eq.6) then
          NTYP= 4
          n04 = hecMESH%adapt_IWK(ie01)
          n05 = hecMESH%adapt_IWK(ie02)
          n06 = hecMESH%adapt_IWK(ie03)
          n14 = hecMESH%adapt_IWK(ie11)
          n15 = hecMESH%adapt_IWK(ie12)
          n16 = hecMESH%adapt_IWK(ie13)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_PRISM (n01, n04, n06, n11, n14, n16, 1)
          call hecmw_adapt_CREATE_NEW_PRISM (n04, n02, n05, n14, n12, n15, 2)
          call hecmw_adapt_CREATE_NEW_PRISM (n06, n05, n03, n16, n15, n13, 3)
          call hecmw_adapt_CREATE_NEW_PRISM (n04, n05, n06, n14, n15, n16, 4)
        endif
!C==

!C
!C-- TYPE of EMBEDDING
           NNTYP          (NTYP)= NNTYP(NTYP) + 1
      enddo

!C===
      return

      contains
        subroutine hecmw_adapt_CREATE_NEW_PRISM (in1,in2,in3,in4,in5,in6, IDchi)

        hecMESH%n_adapt_elem_351_cur= hecMESH%n_adapt_elem_351_cur + 1
        hecMESH%n_adapt_elem_cur    = hecMESH%n_adapt_elem_cur     + 1

        icel = hecMESH%n_adapt_elem_cur
        icouN= icouN + inc

        if (icel.gt.hecMESH%ne_array) then
          call hecmw_adapt_ERROR_EXIT (hecMESH, 61)
        endif

        hecMESH%when_i_was_refined_elem(icel)= hecMESH%n_adapt
        hecMESH%elem_node_index(icel)= hecMESH%elem_node_index(icel-1) + 6

        iS= hecMESH%elem_node_index(icel-1)
        hecMESH%elem_node_item(iS+1)= in1        
        hecMESH%elem_node_item(iS+2)= in2        
        hecMESH%elem_node_item(iS+3)= in3        
        hecMESH%elem_node_item(iS+4)= in4        
        hecMESH%elem_node_item(iS+5)= in5        
        hecMESH%elem_node_item(iS+6)= in6        

        hecMESH%adapt_parent(2*icel-1)= hecMESH%elem_ID(2*NPAR-1)
        hecMESH%adapt_parent(2*icel  )= hecMESH%elem_ID(2*NPAR  )

        hecMESH%elem_ID(2*icel-1)= icouN + hecMESH%ne_internal
        hecMESH%elem_ID(2*icel  )= hecMESH%elem_ID(2*NPAR  )

        hecMESH%elem_mat_ID_item(icel)= hecMESH%elem_mat_ID_item(NPAR)
        hecMESH%section_ID      (icel)= hecMESH%section_ID      (NPAR)

        hecMESH%adapt_type(icel)= 0

        if (NDIVSUM.eq.6) then
          hecMESH%adapt_level(icel)= hecMESH%adapt_level(NPAR) + 2
         else
          hecMESH%adapt_level(icel)= hecMESH%adapt_level(NPAR) + 1
        endif

        iS= hecMESH%adapt_children_index(NPAR-1)

        hecMESH%adapt_children_item(2*(iS+IDchi)-1)= icel
        hecMESH%adapt_children_item(2*(iS+IDchi)-1)= icouN + hecMESH%ne_internal
        hecMESH%adapt_children_item(2*(iS+IDchi)  )= hecMESH%my_rank

        hecMESH%adapt_children_local(iS+IDchi)= icel

        iS= hecMESH%adapt_children_index(icel-1)
        hecMESH%adapt_children_index(icel)= iS + 8

        iS1= iS + 1
        iS2= iS + 2
        iS3= iS + 3
        iS4= iS + 4
        iS5= iS + 5
        iS6= iS + 6
        iS7= iS + 7
        iS8= iS + 8

        hecMESH%adapt_children_item(2*iS1)= -1
        hecMESH%adapt_children_item(2*iS2)= -1
        hecMESH%adapt_children_item(2*iS3)= -1
        hecMESH%adapt_children_item(2*iS4)= -1
        hecMESH%adapt_children_item(2*iS5)= -1
        hecMESH%adapt_children_item(2*iS6)= -1
        hecMESH%adapt_children_item(2*iS7)= -1
        hecMESH%adapt_children_item(2*iS8)= -1

        hecMESH%adapt_children_item(2*iS1-1)=  0
        hecMESH%adapt_children_item(2*iS2-1)=  0
        hecMESH%adapt_children_item(2*iS3-1)=  0
        hecMESH%adapt_children_item(2*iS4-1)=  0
        hecMESH%adapt_children_item(2*iS5-1)=  0
        hecMESH%adapt_children_item(2*iS6-1)=  0
        hecMESH%adapt_children_item(2*iS7-1)=  0
        hecMESH%adapt_children_item(2*iS8-1)=  0

        hecMESH%elem_type        (icel)= PAR_CEL_TYP
        hecMESH%adapt_parent_type(icel)= hecMESH%adapt_type(NPAR)

        end subroutine hecmw_adapt_CREATE_NEW_PRISM
      end subroutine hecmw_adapt_NEW_CELL_351



