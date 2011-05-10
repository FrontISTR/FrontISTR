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
!C*** hecmw_adapt_NEW_CELL_341
!C***
!C
!C    create new TET.
!C
      subroutine hecmw_adapt_NEW_CELL_341 (hecMESH, icouN)

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

      do icel0= 1, hecMESH%n_adapt_act_elem_341
        icel_par= hecMESH%adapt_act_elem_341(icel0)
            NPAR= icel_par

        if (hecMESH%elem_ID(2*icel_par).eq.hecMESH%my_rank) then
          inc= 1          
         else
          inc= 0
        endif

        iS= hecMESH%elem_node_index(icel_par-1)
        n1= hecMESH%elem_node_item (iS+1)
        n2= hecMESH%elem_node_item (iS+2)
        n3= hecMESH%elem_node_item (iS+3)
        n4= hecMESH%elem_node_item (iS+4)
        n5= 0
        n6= 0
        n7= 0
        n8= 0
        n9= 0
        n0= 0

      nnn= hecMESH%n_adapt_edge
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n4, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n4, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n4, ie6, 1 )

        NDIV(1)= hecMESH%adapt_iemb(ie1)
        NDIV(2)= hecMESH%adapt_iemb(ie2)
        NDIV(3)= hecMESH%adapt_iemb(ie3)
        NDIV(4)= hecMESH%adapt_iemb(ie4)
        NDIV(5)= hecMESH%adapt_iemb(ie5)
        NDIV(6)= hecMESH%adapt_iemb(ie6)

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)
        PAR_CEL_TYP= hecMESH%elem_type(NPAR)
!C
!C-- init. CHILD. cell array
        iS= hecMESH%adapt_children_index(NPAR-1)
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

        hecMESH%adapt_children_item(2*iS1-1)= 0     
        hecMESH%adapt_children_item(2*iS2-1)= 0     
        hecMESH%adapt_children_item(2*iS3-1)= 0     
        hecMESH%adapt_children_item(2*iS4-1)= 0     
        hecMESH%adapt_children_item(2*iS5-1)= 0     
        hecMESH%adapt_children_item(2*iS6-1)= 0     
        hecMESH%adapt_children_item(2*iS7-1)= 0     
        hecMESH%adapt_children_item(2*iS8-1)= 0     

!C
!C== embedding TYPE
        if (NDIVSUM.eq.0) NTYP= 0
!C
!C-- TYP-1
        if (NDIVSUM.eq.1 .and. NDIV(1).eq.1) then
          NTYP= 1
          n5  = hecMESH%adapt_IWK(ie1)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n3, n4, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n3, n4, 2)
        endif
!C
!C-- TYP-2
        if (NDIVSUM.eq.1 .and. NDIV(2).eq.1) then
          NTYP= 2
          n5  = hecMESH%adapt_IWK(ie2)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n5, n4, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n3, n4, 2)
        endif
!C
!C-- TYP-3
        if (NDIVSUM.eq.1 .and. NDIV(3).eq.1) then
          NTYP= 3
          n5  = hecMESH%adapt_IWK(ie3)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n3, n5, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n3, n4, 2)
        endif
!C
!C-- TYP-4
        if (NDIVSUM.eq.1 .and. NDIV(4).eq.1) then
          NTYP= 4
          n5  = hecMESH%adapt_IWK(ie4)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n5, n4, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n3, n4, 2)
        endif
!C
!C-- TYP-5
        if (NDIVSUM.eq.1 .and. NDIV(5).eq.1) then
          NTYP= 5
          n5  = hecMESH%adapt_IWK(ie5)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n3, n5, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n3, n4, 2)
        endif
!C
!C-- TYP-6
        if (NDIVSUM.eq.1 .and. NDIV(6).eq.1) then
          NTYP= 6
          n5  = hecMESH%adapt_IWK(ie6)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n3, n5, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n5, n4, 2)
        endif
!C
!C-- TYP-7
        if (NDIVSUM.eq.3 .and. NDIV(1).eq.1 .and.                       &
     &      NDIV(3).eq.1 .and. NDIV(5).eq.1) then
          NTYP= 7
          n5  = hecMESH%adapt_IWK(ie1)
          n6  = hecMESH%adapt_IWK(ie3)
          n7  = hecMESH%adapt_IWK(ie5)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n3, n6, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n3, n7, 2)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n7, n3, n4, 3)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n5, n3, n7, 4)
        endif
!C
!C-- TYP-8
        if (NDIVSUM.eq.3 .and. NDIV(2).eq.1 .and.                       &
     &      NDIV(3).eq.1 .and. NDIV(6).eq.1) then
          NTYP= 8
          n5  = hecMESH%adapt_IWK(ie2)
          n6  = hecMESH%adapt_IWK(ie3)
          n7  = hecMESH%adapt_IWK(ie6)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n5, n6, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n3, n7, 2)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n2, n7, n4, 3)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n2, n5, n7, 4)
        endif
!C
!C-- TYP-9
        if (NDIVSUM.eq.3 .and. NDIV(1).eq.1 .and.                       &
     &      NDIV(2).eq.1 .and. NDIV(4).eq.1) then
          NTYP= 9
          n5  = hecMESH%adapt_IWK(ie1)
          n6  = hecMESH%adapt_IWK(ie2)
          n7  = hecMESH%adapt_IWK(ie4)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n6, n4, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n7, n4, 2)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n7, n3, n4, 3)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n7, n6, n4, 4)
        endif
!C
!C-- TYP-10
        if (NDIVSUM.eq.3 .and. NDIV(4).eq.1 .and.                       &
     &      NDIV(5).eq.1 .and. NDIV(6).eq.1) then
          NTYP= 10
          n5  = hecMESH%adapt_IWK(ie4)
          n6  = hecMESH%adapt_IWK(ie5)
          n7  = hecMESH%adapt_IWK(ie6)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n2, n5, n6, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n3, n7, 2)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n6, n7, n4, 3)
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n7, n6, 4)
        endif
!C
!C-- TYP-11
        if (NDIVSUM.eq.6) then
          NTYP= 11
          n5  = hecMESH%adapt_IWK(ie1)
          n6  = hecMESH%adapt_IWK(ie2)
          n7  = hecMESH%adapt_IWK(ie3)
          n8  = hecMESH%adapt_IWK(ie4)
          n9  = hecMESH%adapt_IWK(ie5)
          n0  = hecMESH%adapt_IWK(ie6)

          hecMESH%adapt_type(NPAR)= NTYP
          call hecmw_adapt_CREATE_NEW_TETRA (n1, n5, n6, n7, 1)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n2, n8, n9, 2)
          call hecmw_adapt_CREATE_NEW_TETRA (n6, n8, n3, n0, 3)
          call hecmw_adapt_CREATE_NEW_TETRA (n7, n9, n0, n4, 4)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n8, n6, n7, 5)
          call hecmw_adapt_CREATE_NEW_TETRA (n5, n8, n7, n9, 6)
          call hecmw_adapt_CREATE_NEW_TETRA (n0, n8, n7, n6, 7)
          call hecmw_adapt_CREATE_NEW_TETRA (n0, n8, n9, n7, 8)
        endif
!C==

!C
!C-- TYPE of EMBEDDING
        NNTYP             (NTYP)= NNTYP(NTYP) + 1
      enddo
!C===

      return

      contains
        subroutine hecmw_adapt_CREATE_NEW_TETRA (in1, in2, in3, in4, IDchi)

        hecMESH%n_adapt_elem_341_cur= hecMESH%n_adapt_elem_341_cur + 1
        hecMESH%n_adapt_elem_cur    = hecMESH%n_adapt_elem_cur     + 1

        icel = hecMESH%n_adapt_elem_cur
        icouN= icouN + inc

        if (icel.gt.hecMESH%ne_array) then
          call hecmw_adapt_ERROR_EXIT (hecMESH, 61)
        endif

        hecMESH%when_i_was_refined_elem(icel)= hecMESH%n_adapt
        hecMESH%elem_node_index(icel)= hecMESH%elem_node_index(icel-1) + 4

        iS= hecMESH%elem_node_index(icel-1)
        hecMESH%elem_node_item(iS+1)= in1        
        hecMESH%elem_node_item(iS+2)= in2        
        hecMESH%elem_node_item(iS+3)= in3        
        hecMESH%elem_node_item(iS+4)= in4        

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


        end subroutine hecmw_adapt_CREATE_NEW_TETRA
      end subroutine hecmw_adapt_NEW_CELL_341



