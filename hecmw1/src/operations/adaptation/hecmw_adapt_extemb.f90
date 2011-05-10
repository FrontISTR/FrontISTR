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
!C*** hecmw_adapt_EXTEMB
!C***
!C
!C    extend EMBEDDED region to NEIGHBORING region
!C
      subroutine hecmw_adapt_EXTEMB (hecMESH)
      use hecmw_util
      use hecmw_adapt_INT_SR
      use hecmw_adapt_INT_SR_REV
      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint), dimension(:), allocatable :: WR, WS
      dimension NDIV(6)

      type (hecmwST_local_mesh) :: hecMESH

      NROW= 0
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      do iedg= 1, hecMESH%n_adapt_edge
        if (hecMESH%adapt_iemb(iedg).eq.1) hecMESH%adapt_iemb(iedg)=2
      enddo

!C
!C-- TETRAHEDRA
      do icel0= 1, hecMESH%n_adapt_act_elem_341
        icel= hecMESH%adapt_act_elem_341(icel0)
          iS= hecMESH%elem_node_index(icel-1)
          n1= hecMESH%elem_node_item (iS+1)
          n2= hecMESH%elem_node_item (iS+2)
          n3= hecMESH%elem_node_item (iS+3)
          n4= hecMESH%elem_node_item (iS+4)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n4, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n4, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n4, ie6, 1 )

        NDIV(1)= 0
        NDIV(2)= 0
        NDIV(3)= 0
        NDIV(4)= 0
        NDIV(5)= 0
        NDIV(6)= 0

        if (hecMESH%adapt_iemb(ie1).eq.2) NDIV(1)= 1
        if (hecMESH%adapt_iemb(ie2).eq.2) NDIV(2)= 1
        if (hecMESH%adapt_iemb(ie3).eq.2) NDIV(3)= 1
        if (hecMESH%adapt_iemb(ie4).eq.2) NDIV(4)= 1
        if (hecMESH%adapt_iemb(ie5).eq.2) NDIV(5)= 1
        if (hecMESH%adapt_iemb(ie6).eq.2) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.ge.1) then
          if (hecMESH%adapt_iemb(ie1).eq.0) hecMESH%adapt_iemb(ie1)= 1
          if (hecMESH%adapt_iemb(ie2).eq.0) hecMESH%adapt_iemb(ie2)= 1
          if (hecMESH%adapt_iemb(ie3).eq.0) hecMESH%adapt_iemb(ie3)= 1
          if (hecMESH%adapt_iemb(ie4).eq.0) hecMESH%adapt_iemb(ie4)= 1
          if (hecMESH%adapt_iemb(ie5).eq.0) hecMESH%adapt_iemb(ie5)= 1
          if (hecMESH%adapt_iemb(ie6).eq.0) hecMESH%adapt_iemb(ie6)= 1
        endif
      enddo

!C
!C-- PRISMs
      do icel0= 1, hecMESH%n_adapt_act_elem_351
        icel= hecMESH%adapt_act_elem_351(icel0)
          iS= hecMESH%elem_node_index(icel-1)
          n1= hecMESH%elem_node_item (iS+1)
          n2= hecMESH%elem_node_item (iS+2)
          n3= hecMESH%elem_node_item (iS+3)
          n4= hecMESH%elem_node_item (iS+4)
          n5= hecMESH%elem_node_item (iS+5)
          n6= hecMESH%elem_node_item (iS+6)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n1, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n4, n5, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n5, n6, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n6, n4, ie6, 1 )

        NDIV(1)= 0
        NDIV(2)= 0
        NDIV(3)= 0
        NDIV(4)= 0
        NDIV(5)= 0
        NDIV(6)= 0

        if (hecMESH%adapt_iemb(ie1).eq.2) NDIV(1)= 1
        if (hecMESH%adapt_iemb(ie2).eq.2) NDIV(2)= 1
        if (hecMESH%adapt_iemb(ie3).eq.2) NDIV(3)= 1
        if (hecMESH%adapt_iemb(ie4).eq.2) NDIV(4)= 1
        if (hecMESH%adapt_iemb(ie5).eq.2) NDIV(5)= 1
        if (hecMESH%adapt_iemb(ie6).eq.2) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.ge.1) then
          if (hecMESH%adapt_iemb(ie1).eq.0) hecMESH%adapt_iemb(ie1)= 1
          if (hecMESH%adapt_iemb(ie2).eq.0) hecMESH%adapt_iemb(ie2)= 1
          if (hecMESH%adapt_iemb(ie3).eq.0) hecMESH%adapt_iemb(ie3)= 1
          if (hecMESH%adapt_iemb(ie4).eq.0) hecMESH%adapt_iemb(ie4)= 1
          if (hecMESH%adapt_iemb(ie5).eq.0) hecMESH%adapt_iemb(ie5)= 1
          if (hecMESH%adapt_iemb(ie6).eq.0) hecMESH%adapt_iemb(ie6)= 1
        endif
      enddo
!C===

!C
!C-- exchange IEMB
      N1= hecMESH%adapt_import_edge_index(hecMESH%n_neighbor_pe)
      N2= hecMESH%adapt_export_edge_index(hecMESH%n_neighbor_pe)
      m = max (N1, N2)
      allocate (WS(m), WR(m))

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV_REV                                &
     &   ( hecMESH%n_adapt_edge, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,  &
     &     hecMESH%adapt_import_edge_index,                             &
     &     hecMESH%adapt_import_edge_item ,                             & 
     &     hecMESH%adapt_export_edge_index,                             &
     &     hecMESH%adapt_export_edge_item ,                             & 
     &     WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,                &
     &     hecMESH%my_rank, 1, m)

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( hecMESH%n_adapt_edge, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,  &
     &     hecMESH%adapt_import_edge_index,                             &
     &     hecMESH%adapt_import_edge_item ,                             & 
     &     hecMESH%adapt_export_edge_index,                             &
     &     hecMESH%adapt_export_edge_item ,                             & 
     &     WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,                &
     &     hecMESH%my_rank, 1, m)

      do iedg= 1, hecMESH%n_adapt_edge
        if (hecMESH%adapt_iemb(iedg).eq.1) hecMESH%adapt_iemb(iedg)=2
      enddo

!C
!C +----------------------+
!C | extend embedded zone |
!C +----------------------+
!C   RULE : set all EDGEs IEMB(ie)=1 if at least ONE
!C          EDGE of the TETRAHEDRON is marked as IEMB(ie)=2
!C===
      do irow= 1, NROW
!C
!C== set IEMB(ie) = 1

!C
!C-- TETRAHEDRA
      do icel0= 1, hecMESH%n_adapt_act_elem_341
        icel= hecMESH%adapt_act_elem_341(icel0)
          iS= hecMESH%elem_node_index(icel-1)
          n1= hecMESH%elem_node_item (iS+1)
          n2= hecMESH%elem_node_item (iS+2)
          n3= hecMESH%elem_node_item (iS+3)
          n4= hecMESH%elem_node_item (iS+4)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n4, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n4, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n4, ie6, 1 )

        NDIV(1)= 0
        NDIV(2)= 0
        NDIV(3)= 0
        NDIV(4)= 0
        NDIV(5)= 0
        NDIV(6)= 0

        if (hecMESH%adapt_iemb(ie1).eq.2) NDIV(1)= 1
        if (hecMESH%adapt_iemb(ie2).eq.2) NDIV(2)= 1
        if (hecMESH%adapt_iemb(ie3).eq.2) NDIV(3)= 1
        if (hecMESH%adapt_iemb(ie4).eq.2) NDIV(4)= 1
        if (hecMESH%adapt_iemb(ie5).eq.2) NDIV(5)= 1
        if (hecMESH%adapt_iemb(ie6).eq.2) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.ge.1) then
          if (hecMESH%adapt_iemb(ie1).eq.0) hecMESH%adapt_iemb(ie1)= 1
          if (hecMESH%adapt_iemb(ie2).eq.0) hecMESH%adapt_iemb(ie2)= 1
          if (hecMESH%adapt_iemb(ie3).eq.0) hecMESH%adapt_iemb(ie3)= 1
          if (hecMESH%adapt_iemb(ie4).eq.0) hecMESH%adapt_iemb(ie4)= 1
          if (hecMESH%adapt_iemb(ie5).eq.0) hecMESH%adapt_iemb(ie5)= 1
          if (hecMESH%adapt_iemb(ie6).eq.0) hecMESH%adapt_iemb(ie6)= 1
        endif
      enddo

!C
!C-- PRISMs
      do icel0= 1, hecMESH%n_adapt_act_elem_351
        icel= hecMESH%adapt_act_elem_351(icel0)
          iS= hecMESH%elem_node_index(icel-1)
          n1= hecMESH%elem_node_item (iS+1)
          n2= hecMESH%elem_node_item (iS+2)
          n3= hecMESH%elem_node_item (iS+3)
          n4= hecMESH%elem_node_item (iS+4)
          n5= hecMESH%elem_node_item (iS+5)
          n6= hecMESH%elem_node_item (iS+6)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n1, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n4, n5, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n5, n6, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n6, n4, ie6, 1 )

        NDIV(1)= 0
        NDIV(2)= 0
        NDIV(3)= 0
        NDIV(4)= 0
        NDIV(5)= 0
        NDIV(6)= 0

        if (hecMESH%adapt_iemb(ie1).eq.2) NDIV(1)= 1
        if (hecMESH%adapt_iemb(ie2).eq.2) NDIV(2)= 1
        if (hecMESH%adapt_iemb(ie3).eq.2) NDIV(3)= 1
        if (hecMESH%adapt_iemb(ie4).eq.2) NDIV(4)= 1
        if (hecMESH%adapt_iemb(ie5).eq.2) NDIV(5)= 1
        if (hecMESH%adapt_iemb(ie6).eq.2) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.ge.1) then
          if (hecMESH%adapt_iemb(ie1).eq.0) hecMESH%adapt_iemb(ie1)= 1
          if (hecMESH%adapt_iemb(ie2).eq.0) hecMESH%adapt_iemb(ie2)= 1
          if (hecMESH%adapt_iemb(ie3).eq.0) hecMESH%adapt_iemb(ie3)= 1
          if (hecMESH%adapt_iemb(ie4).eq.0) hecMESH%adapt_iemb(ie4)= 1
          if (hecMESH%adapt_iemb(ie5).eq.0) hecMESH%adapt_iemb(ie5)= 1
          if (hecMESH%adapt_iemb(ie6).eq.0) hecMESH%adapt_iemb(ie6)= 1
        endif
      enddo
!C==

!C
!C-- exchange IEMB
      WS= 0 
      WR= 0
      call hecmw_adapt_INT_SEND_RECV_REV                                &
     &   ( hecMESH%n_adapt_edge,                                        &
     &     hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                  &
     &     hecMESH%adapt_import_edge_index,                             &
     &     hecMESH%adapt_import_edge_item ,                             & 
     &     hecMESH%adapt_export_edge_index,                             &
     &     hecMESH%adapt_export_edge_item ,                             & 
     &     WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,                &
     &     hecMESH%my_rank, 1, m)

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( hecMESH%n_adapt_edge,                                        &
     &     hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                  &
     &     hecMESH%adapt_import_edge_index,                             &
     &     hecMESH%adapt_import_edge_item ,                             & 
     &     hecMESH%adapt_export_edge_index,                             &
     &     hecMESH%adapt_export_edge_item ,                             & 
     &     WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,                &
     &     hecMESH%my_rank, 1, m)

!C
!C-- set IEMB(ie) = 2
      do iedg= 1, hecMESH%n_adapt_edge
        if ( hecMESH%adapt_iemb(iedg).eq.1 ) then
          hecMESH%adapt_iemb(iedg)= 2
        endif
      enddo

      enddo

      deallocate (WS, WR)
!C=== 

!C
!C-- LOOP OVER ALL EDGEs
      do iedg= 1, hecMESH%n_adapt_edge
        if ( hecMESH%adapt_iemb(iedg).eq.2 ) hecMESH%adapt_iemb(iedg)= 1
      enddo

      return
      end
