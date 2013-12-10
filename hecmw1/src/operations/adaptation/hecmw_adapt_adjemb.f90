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
!C*** hecmw_adapt_ADJEMB
!C***
!C
!C    adjust cell EMBEDDING LEVEL around each NODE
!C
!C    BASIC RULE : 
!C      keep MAX. difference of embedding level of
!C      neighboring cells around each node NOT MORE THAN 2
!C
      subroutine hecmw_adapt_ADJEMB ( hecMESH, NFLAG_INFO)

      use hecmw_util
      use hecmw_adapt_INT_SR
      use hecmw_adapt_INT_SR_REV

      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(:), allocatable :: WR, WS
      dimension NDIV(6)

      integer(kind=kint), pointer :: ADAPT_nodLEVmax (:), ADAPT_LEVcur(:)

      type (hecmwST_local_mesh) :: hecMESH

!C
!C-- INIT.
      allocate (ADAPT_nodLEVmax(hecMESH%nn_array))
      allocate (ADAPT_LEVcur   (hecMESH%ne_array))

      ADAPT_nodLEVmax= 0
      ADAPT_LEVcur   = 0

!C
!C +-------------------------------------------+
!C | find MAX.embedding LEVEL around each node |
!C +-------------------------------------------+
!C   ONE-directional embedding - add +1 to ADAPT_LEV
!C   ALL-directional embedding - add +2 to ADAPT_LEV
!C===

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

        if ( hecMESH%adapt_iemb(ie1).gt.0 ) NDIV(1)= 1
        if ( hecMESH%adapt_iemb(ie2).gt.0 ) NDIV(2)= 1
        if ( hecMESH%adapt_iemb(ie3).gt.0 ) NDIV(3)= 1
        if ( hecMESH%adapt_iemb(ie4).gt.0 ) NDIV(4)= 1
        if ( hecMESH%adapt_iemb(ie5).gt.0 ) NDIV(5)= 1
        if ( hecMESH%adapt_iemb(ie6).gt.0 ) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.eq.0)                   NLEV_ADD= 0        
        if (NDIVSUM.eq.1 .or. NDIVSUM.eq.3) NLEV_ADD= 1        
        if (NDIVSUM.eq.6)                   NLEV_ADD= 2        

        NL= hecMESH%adapt_level(icel) + NLEV_ADD

        ADAPT_LEVcur(icel)= NL

        m1= ADAPT_nodLEVmax(n1)
        m2= ADAPT_nodLEVmax(n2)
        m3= ADAPT_nodLEVmax(n3)
        m4= ADAPT_nodLEVmax(n4)

        ADAPT_nodLEVmax(n1)= max (NL, m1)        
        ADAPT_nodLEVmax(n2)= max (NL, m2)        
        ADAPT_nodLEVmax(n3)= max (NL, m3)        
        ADAPT_nodLEVmax(n4)= max (NL, m4)        

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

        if ( hecMESH%adapt_iemb(ie1).gt.0 ) NDIV(1)= 1
        if ( hecMESH%adapt_iemb(ie2).gt.0 ) NDIV(2)= 1
        if ( hecMESH%adapt_iemb(ie3).gt.0 ) NDIV(3)= 1
        if ( hecMESH%adapt_iemb(ie4).gt.0 ) NDIV(4)= 1
        if ( hecMESH%adapt_iemb(ie5).gt.0 ) NDIV(5)= 1
        if ( hecMESH%adapt_iemb(ie6).gt.0 ) NDIV(6)= 1

        NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

        if (NDIVSUM.eq.0) NLEV_ADD= 0        
        if (NDIVSUM.eq.2) NLEV_ADD= 1        
        if (NDIVSUM.eq.6) NLEV_ADD= 2        

        NL= hecMESH%adapt_level(icel) + NLEV_ADD

        ADAPT_LEVcur(icel)= NL

        m1= ADAPT_nodLEVmax(n1)
        m2= ADAPT_nodLEVmax(n2)
        m3= ADAPT_nodLEVmax(n3)
        m4= ADAPT_nodLEVmax(n4)
        m5= ADAPT_nodLEVmax(n5)
        m6= ADAPT_nodLEVmax(n6)

        ADAPT_nodLEVmax(n1)= max (NL, m1)        
        ADAPT_nodLEVmax(n2)= max (NL, m2)        
        ADAPT_nodLEVmax(n3)= max (NL, m3)        
        ADAPT_nodLEVmax(n4)= max (NL, m4)        
        ADAPT_nodLEVmax(n5)= max (NL, m5)        
        ADAPT_nodLEVmax(n6)= max (NL, m6)        

      enddo
!C===

      if (hecMESH%PETOT.ne.1) then
!C
!C-- exchange ADAPT_nodLEVmax
      N = hecMESH%n_node
      N1= hecMESH%import_index(hecMESH%n_neighbor_pe)
      N2= hecMESH%export_index(hecMESH%n_neighbor_pe)

      m = max (N1, N2)
      allocate (WS(m), WR(m))

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( N, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,               &
     &     hecMESH%import_index, hecMESH%import_item,                   &
     &     hecMESH%export_index, hecMESH%export_item,                   &
     &     WS, WR, ADAPT_nodLEVmax, hecMESH%MPI_COMM, hecMESH%my_rank,  &
     &     1, m)
      deallocate (WS, WR)
      endif

!C
!C +------------------------+
!C | adjust embedding level |
!C +------------------------+
!C===

!C
!C-- TETRAHEDRA
      do icel0= 1, hecMESH%n_adapt_act_elem_341
        icel= hecMESH%adapt_act_elem_341(icel0)
          iS= hecMESH%elem_node_index(icel-1)
          n1= hecMESH%elem_node_item (iS+1)
          n2= hecMESH%elem_node_item (iS+2)
          n3= hecMESH%elem_node_item (iS+3)
          n4= hecMESH%elem_node_item (iS+4)

        NL= ADAPT_LEVcur(icel)

        m1= ADAPT_nodLEVmax(n1)
        m2= ADAPT_nodLEVmax(n2)
        m3= ADAPT_nodLEVmax(n3)
        m4= ADAPT_nodLEVmax(n4)

        if (((m1-NL).gt.2).or.((m2-NL).gt.2).or.((m3-NL).gt.2).or.      &
     &      ((m4-NL).gt.2 )) then

          call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n3, ie2, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n4, ie3, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie4, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n4, ie5, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n4, ie6, 1 )

          hecMESH%adapt_iemb(ie1)= 1
          hecMESH%adapt_iemb(ie2)= 1
          hecMESH%adapt_iemb(ie3)= 1
          hecMESH%adapt_iemb(ie4)= 1
          hecMESH%adapt_iemb(ie5)= 1
          hecMESH%adapt_iemb(ie6)= 1

          NFLAG_INFO= 1
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

        NL= ADAPT_LEVcur(icel)

        m1= ADAPT_nodLEVmax(n1)
        m2= ADAPT_nodLEVmax(n2)
        m3= ADAPT_nodLEVmax(n3)
        m4= ADAPT_nodLEVmax(n4)
        m5= ADAPT_nodLEVmax(n5)
        m6= ADAPT_nodLEVmax(n6)

        if (((m1-NL).gt.2).or.((m2-NL).gt.2).or.((m3-NL).gt.2).or.      &
     &      ((m4-NL).gt.2).or.((m5-NL).gt.2).or.((m6-NL).gt.2)) then

          call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie2, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n1, ie3, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n4, n5, ie4, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n5, n6, ie5, 1 )
          call hecmw_adapt_EDGE_INFO ( hecMESH, n6, n4, ie6, 1 )

          hecMESH%adapt_iemb(ie1)= 1
          hecMESH%adapt_iemb(ie2)= 1
          hecMESH%adapt_iemb(ie3)= 1
          hecMESH%adapt_iemb(ie4)= 1
          hecMESH%adapt_iemb(ie5)= 1
          hecMESH%adapt_iemb(ie6)= 1

          NFLAG_INFO= 1
        endif
      enddo
!C===

!C
!C-- exchange hecMESH%iemb
      N1= hecMESH%adapt_import_edge_index(hecMESH%n_neighbor_pe)
      N2= hecMESH%adapt_export_edge_index(hecMESH%n_neighbor_pe)
      m = max (N1, N2)
      allocate (WS(m), WR(m))

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
      deallocate (WS, WR)

      deallocate (ADAPT_nodLEVmax, ADAPT_LEVcur)

      return
      end
