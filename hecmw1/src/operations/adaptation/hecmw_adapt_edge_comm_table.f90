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
!C*** hecmw_adapt_edge_comm_table
!C***
!C
!C    global edge information
!C    edge-based communication table
!C
      subroutine hecmw_adapt_edge_comm_table (hecMESH)
      
      use  hecmw_util
      use  hecmw_adapt_STACK_SR
      use  hecmw_adapt_ITEM_SR

      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint ), pointer :: wSI(:), wSE(:)
      integer(kind=kint ), pointer :: wiIa(:), wiEa(:), wiIb(:), wiEb(:)
      integer(kind=kint ), dimension(:), allocatable :: IW1, IW2
    
      type (hecmwST_local_mesh) :: hecMESH

!C
!C +-------------------------+
!C | EDGE related parameters |
!C +-------------------------+
!C===
      allocate (hecMESH%adapt_iemb    (hecMESH%n_adapt_edge),           &
     &          hecMESH%adapt_mid_edge(hecMESH%n_adapt_edge))
        hecMESH%adapt_iemb    = 0
        hecMESH%adapt_mid_edge= 0

      allocate (hecMESH%rev_neighbor_pe(0:hecMESH%n_neighbor_pe))

      hecMESH%rev_neighbor_pe(hecMESH%my_rank)= 0
      do neib= 1, hecMESH%n_neighbor_pe
        hecMESH%rev_neighbor_pe(hecMESH%neighbor_pe(neib))= neib
      enddo

      allocate (hecMESH%adapt_act_elem_341(hecMESH%n_adapt_elem_341))
      allocate (hecMESH%adapt_act_elem_351(hecMESH%n_adapt_elem_351))

      icouTA= 0
      icouPA= 0
      do icel= 1, hecMESH%n_elem
        ityp= hecMESH%elem_type(icel)
!C
!C-- TETRAHEDRA : active
        if (ityp.eq.341) then
          if (hecMESH%adapt_type(icel).eq.0) then
                                 icouTA = icouTA + 1
            hecMESH%adapt_act_elem_341(icouTA)= icel
          endif
        endif
!C
!C-- PRISMS : active
        if (ityp.eq.351) then
          if (hecMESH%adapt_type(icel).eq.0) then
                                 icouPA = icouPA + 1
            hecMESH%adapt_act_elem_351(icouPA)= icel
          endif
        endif
      enddo

      hecMESH%n_adapt_act_elem_341= icouTA
      hecMESH%n_adapt_act_elem_351= icouPA
!C===
!C
!C +---------------+
!C | Global EDGE # |
!C +---------------+
!C===
      nnn= 0
      do ie= 1, hecMESH%n_adapt_edge
        if (hecMESH%adapt_edge_home(ie).eq.hecMESH%my_rank) nnn= nnn + 1
      enddo        
      hecMESH%n_adapt_edge_global= hecMESH%n_adapt_act_edge

      call hecmw_allreduce_I (hecMESH, hecMESH%n_adapt_edge_global, 1, hecmw_sum)
!C===

!C
!C +-----------------------------+
!C | prepare EXTERNAL edge info. |
!C +-----------------------------+
!C===

!C
!C-- init.
      NEIBPETOT= hecMESH%n_neighbor_pe

      allocate (IW1(hecMESH%n_neighbor_pe))
      allocate (wSE(0:hecMESH%n_neighbor_pe), wSI(0:hecMESH%n_neighbor_pe))
      IW1 = 0
      wSE = 0
      wSI = 0

!C
!C-- search IMPORT items
      do ie= 1, hecMESH%n_adapt_edge
        ih= hecMESH%adapt_edge_home(ie)
        if (ih.ne.hecMESH%my_rank) then
              ihR = hecMESH%rev_neighbor_pe(ih)
          wSI(ihR)= wSI(ihR) + 1
        endif
      enddo
!C
!C-- exchange INFO. on IMPORT item #
      stime0= MPI_WTIME ()
      call hecmw_adapt_STACK_SEND_RECV                                  &
     &   ( hecMESH%n_neighbor_pe, hecMESH%neighbor_pe, wSI, wSE,        &
     &     hecMESH%MPI_COMM, hecMESH%my_rank)
      etime0= MPI_WTIME ()
      commTIME= commTIME + etime0 - stime0

!C
!C-- IMPORT/EXPORT item #
      do neib= 1, hecMESH%n_neighbor_pe
        wSI(neib)= wSI(neib-1) + wSI(neib)
        wSE(neib)= wSE(neib-1) + wSE(neib)
      enddo

!C
!C-- send as IMPORT/recv. as EXPORT
      allocate (wiIa(wSI(hecMESH%n_neighbor_pe)*4))
      allocate (wiEa(wSE(hecMESH%n_neighbor_pe)*4))
      allocate (wiIb(wSI(hecMESH%n_neighbor_pe)  ))
      allocate (wiEb(wSE(hecMESH%n_neighbor_pe)  ))
! JF
      wiIa=0
      wiEa=0
      wiIb=0
      wiEb=0

      do ie= 1, hecMESH%n_adapt_edge
        ih= hecMESH%adapt_edge_home(ie)

        if (ih.ne.hecMESH%my_rank) then
              ihR     = hecMESH%rev_neighbor_pe(ih)
          IW1(ihR  )  = IW1(ihR) + 1
                 is   = wSI(ihR-1)+IW1(ihR)

          in1= hecMESH%adapt_edge_node(2*ie-1)
          in2= hecMESH%adapt_edge_node(2*ie  )

          wiIa(4*is-3)= hecMESH%node_ID(2*in1  )
          wiIa(4*is-2)= hecMESH%node_ID(2*in1-1)
          wiIa(4*is-1)= hecMESH%node_ID(2*in2  )
          wiIa(4*is  )= hecMESH%node_ID(2*in2-1)
        endif
      enddo

      LEN= max(wSI(hecMESH%n_neighbor_pe),wSE(hecMESH%n_neighbor_pe),   &
     &         hecMESH%n_adapt_edge)

      deallocate (IW1)
        allocate (IW1(LEN*4), IW2(LEN*4))
      IW1 = 0
      IW2 = 0
      stime0= MPI_WTIME ()
      call hecmw_adapt_ITEM_SEND_RECV                                   &      
     &    (LEN, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
     &     wSI, wiIa, wSE, wiEa, IW1, IW2,                              &
     &     hecMESH%MPI_COMM, hecMESH%my_rank, 4)
      etime0= MPI_WTIME ()
      commTIME= commTIME + etime0 - stime0

      deallocate (IW1, IW2)
!C===

!C
!C +---------------------+
!C | EXTERNAL edge info. |
!C +---------------------+
!C===


!C
!C-- find LOCAL edge ID at DESTINY
      do neib= 1, hecMESH%n_neighbor_pe
        is= wSE(neib-1)+1
        ie= wSE(neib  )
        do k= is, ie
          ip1= wiEa(4*k-3)
          in1= wiEa(4*k-2)
          ip2= wiEa(4*k-1)
          in2= wiEa(4*k  )
          call hecmw_adapt_LOCAL_NODE_INFO (ip1,in1,inC1)
          call hecmw_adapt_LOCAL_NODE_INFO (ip2,in2,inC2)
          call hecmw_adapt_EDGE_INFO (hecMESH, inC1,inC2,IE0,1)
          wiEb(k)= IE0
        enddo
      enddo

      LEN= max(wSI(hecMESH%n_neighbor_pe),wSE(hecMESH%n_neighbor_pe),IEDGTOT)
      allocate (IW1(LEN), IW2(LEN))
      IW1 = 0
      IW2 = 0

      stime0= MPI_WTIME ()
      call hecmw_adapt_ITEM_SEND_RECV                                   &      
     &    (LEN, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
     &     wSE, wiEb, wSI, wiIb, IW1, IW2,                              &
     &     hecMESH%MPI_COMM, hecMESH%my_rank, 1)
      etime0= MPI_WTIME ()
      commTIME= commTIME + etime0 - stime0

!C
!C-- reconstruct IMPORT table

      deallocate (IW1)
        allocate (IW1(hecMESH%n_neighbor_pe))
      IW1 = 0

      do ie= 1, hecMESH%n_adapt_edge
        ih= hecMESH%adapt_edge_home(ie)

        if (ih.ne.hecMESH%my_rank) then
              ihR   = hecMESH%rev_neighbor_pe(ih)
          IW1(ihR  )= IW1(ihR) + 1

          wiIb(wSI(ihR-1)+IW1(ihR))= ie
        endif
      enddo

!C
!C-- new ARRAY
      allocate (hecMESH%adapt_import_edge_index(0:hecMESH%n_neighbor_pe))
      allocate (hecMESH%adapt_export_edge_index(0:hecMESH%n_neighbor_pe))

      do neib= 0, hecMESH%n_neighbor_pe
        hecMESH%adapt_import_edge_index(neib)= wSI(neib)
        hecMESH%adapt_export_edge_index(neib)= wSE(neib)
      enddo

      MAXimport= wSI(hecMESH%n_neighbor_pe)
      MAXexport= wSE(hecMESH%n_neighbor_pe)

      allocate (hecMESH%adapt_import_edge_item(MAXimport))
      allocate (hecMESH%adapt_export_edge_item(MAXexport))
      do k= 1, MAXimport
        hecMESH%adapt_import_edge_item(k)= wiIb(k)
      enddo
      do k= 1, MAXexport
        hecMESH%adapt_export_edge_item(k)= wiEb(k)
      enddo

      deallocate (IW1,IW2)
      deallocate (wSE,wSI,wiEa,wiEb,wiIa,wiIb)

!C===
      contains
        subroutine hecmw_adapt_LOCAL_NODE_INFO (ip,in,in0)
          do i= 1, hecMESH%n_node
            if (hecMESH%node_ID(2*i)  .eq.ip .and.                      &
     &          hecMESH%node_ID(2*i-1).eq.in) then
              in0= i
              return
            endif
          enddo
        end subroutine hecmw_adapt_LOCAL_NODE_INFO
      end subroutine hecmw_adapt_EDGE_COMM_TABLE
