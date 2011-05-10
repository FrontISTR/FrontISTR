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
!C*** hecmw_adapt_NEW_NODE
!C***
!C
!C    create NEW NODEs in TETRAHEDRAL REGION
!C
      subroutine hecmw_adapt_NEW_NODE (hecMESH)

      use hecmw_util
      use hecmw_adapt_INT_SR
      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint), dimension(:), allocatable :: WR, WS, IW1, IW2

      type (hecmwST_local_mesh) :: hecMESH

!C
!C-- INIT.

!C
!C +----------------+
!C | mid-edge POINT |
!C +----------------+
!C===
      hecMESH%n_adapt_node_cur= hecMESH%n_node
      hecMESH%n_adapt_node_old= hecMESH%n_node
      intCOUNT                = hecMESH%nn_internal

      icou= 0
      do ie= 1, hecMESH%n_adapt_edge
        n1= hecMESH%adapt_edge_node(2*ie-1)
        n2= hecMESH%adapt_edge_node(2*ie  )

        if ( hecMESH%adapt_iemb(ie) .eq. 1) then
          hecMESH%n_adapt_node_cur= hecMESH%n_adapt_node_cur + 1
          if (hecMESH%n_adapt_node_cur.gt.hecMESH%nn_array) then
            call hecmw_adapt_ERROR_EXIT (hecMESH, 62)
          endif

          ieh= hecMESH%adapt_edge_home(ie)
          nod= hecMESH%n_adapt_node_cur

          hecMESH%node_ID(2*nod)= ieh
          if (ieh.eq.hecMESH%my_rank) then
            intCOUNT= intCOUNT + 1
            hecMESH%adapt_mid_edge(ie)= intCOUNT
            hecMESH%node_ID(2*nod-1)  = intCOUNT
           else
            hecMESH%adapt_mid_edge(ie)= nod
            hecMESH%node_ID (2*nod-1) = nod
          endif

          hecMESH%when_i_was_refined_node(nod)= hecMESH%n_adapt
          hecMESH%adapt_IWK(ie)               = nod

          if (ieh.eq.hecMESH%my_rank) icou= icou + 1

          hecMESH%node(3*nod-2)= 0.5d0 * ( hecMESH%node(3*n1-2) +       &
     &                                     hecMESH%node(3*n2-2) )
          hecMESH%node(3*nod-1)= 0.5d0 * ( hecMESH%node(3*n1-1) +       &
     &                                     hecMESH%node(3*n2-1) )
          hecMESH%node(3*nod  )= 0.5d0 * ( hecMESH%node(3*n1  ) +       &
     &                                     hecMESH%node(3*n2  ) )
        endif
      enddo

      call MPI_BARRIER  (hecMESH%MPI_COMM,ierr)

      hecMESH%nn_adapt_internal_cur= hecMESH%nn_internal + icou
!C===

!C
!C +---------------+
!C | GLOBAL node # |
!C +---------------+
!C===

!C
!C-- OLD/NEW global NODE #
      NODTOTg_old= 0
      NODTOTg_cur= 0

      call MPI_allREDUCE (hecMESH%nn_internal, NODTOTg_old, 1,          &
     &                    MPI_INTEGER, MPI_SUM, hecMESH%MPI_COMM, ierr)
      call MPI_allREDUCE (hecMESH%nn_adapt_internal_cur, NODTOTg_cur, 1,&
     &                    MPI_INTEGER, MPI_SUM, hecMESH%MPI_COMM, ierr)

      if (hecMESH%my_rank.eq.0) then
        write (*,'("  total node number (before)", i8  )') NODTOTg_old
        write (*,'("  total node number (curr. )", i8,/)') NODTOTg_cur
      endif

      hecMESH%n_node          = hecMESH%n_adapt_node_cur
      hecMESH%nn_internal     = hecMESH%nn_adapt_internal_cur

!C
!C-- exchange MIDEDG
      allocate (IW1(hecMESH%n_adapt_edge))

      IW1= hecMESH%adapt_mid_edge

      N1= hecMESH%adapt_import_edge_index(hecMESH%n_neighbor_pe)
      N2= hecMESH%adapt_export_edge_index(hecMESH%n_neighbor_pe)
      m = max (N1, N2)
      allocate (WS(m), WR(m))

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( hecMESH%n_adapt_edge,                                        &
     &     hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                  &
     &     hecMESH%adapt_import_edge_index,                             &
     &     hecMESH%adapt_import_edge_item ,                             & 
     &     hecMESH%adapt_export_edge_index,                             &
     &     hecMESH%adapt_export_edge_item ,                             & 
     &     WS, WR, hecMESH%adapt_mid_edge, hecMESH%MPI_COMM,            &
     &     hecMESH%my_rank, 1, m)
      deallocate (WS, WR)

      do ie= 1, hecMESH%n_adapt_edge
        if (hecMESH%adapt_iemb(ie).eq.1 .and.                           &
     &      hecMESH%adapt_edge_home(ie).ne.hecMESH%my_rank) then
          inG= hecMESH%adapt_mid_edge(ie)
          inL= IW1   (ie)
          hecMESH%node_ID(2*inL  )= hecMESH%adapt_edge_home(ie)
          hecMESH%node_ID(2*inL-1)= inG
        endif
      enddo

      do ie= 1, hecMESH%n_adapt_edge
        if (hecMESH%adapt_iemb(ie).eq.1) then
          inL= IW1   (ie)
          hecMESH%adapt_mid_edge(ie)= inL
        endif
      enddo
      deallocate (IW1)
!C===

!C
!C +-------------+
!C | RE-ORDERING |
!C +-------------+
!C===
      allocate (IW1(-1:hecMESH%n_neighbor_pe),                          &
     &          IW2(-1:hecMESH%n_neighbor_pe))
      allocate (hecMESH%adapt_NEWtoOLD_node(hecMESH%n_node),            &
     &          hecMESH%adapt_OLDtoNEW_node(hecMESH%n_node))

      IW1= 0
      IW2= 0
      hecMESH%adapt_NEWtoOLD_node= 0
      hecMESH%adapt_OLDtoNEW_node= 0

      do in= 1, hecMESH%n_node
            ih  = hecMESH%node_ID(2*in)
            ihR = hecMESH%rev_neighbor_pe(ih)
        IW1(ihR)= IW1(ihR) + 1
      enddo

      do neib= 0, hecMESH%n_neighbor_pe
        IW2(neib)= IW2(neib-1) + IW1(neib)
      enddo

      IW1= 0
      do in= 1, hecMESH%n_node
            ih  = hecMESH%node_ID(2*in)
            ihR = hecMESH%rev_neighbor_pe(ih)
        IW1(ihR)= IW1(ihR) + 1
            is  = IW2(ihR-1) + IW1(ihR)
        hecMESH%adapt_OLDtoNEW_node(in)= is
        hecMESH%adapt_NEWtoOLD_node(is)= in
      enddo
!C===

      deallocate (IW1, IW2)

      return
      end






