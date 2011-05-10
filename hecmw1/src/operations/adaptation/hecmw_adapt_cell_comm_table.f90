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
!C*** hecmw_adapt_cell_comm_table
!C***
!C
!C    cell-based communication table
!C
      subroutine hecmw_adapt_cell_comm_table (hecMESH)
      
      use  hecmw_util
      use  hecmw_adapt_STACK_SR
      use  hecmw_adapt_ITEM_SR

      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint ), pointer :: wSI(:), wSE(:), wiI(:), wiE(:)
      integer(kind=kint ), dimension(:), allocatable :: IW1, IW2

      type (hecmwST_local_mesh) :: hecMESH

!C
!C-- init.
      allocate (IW1(hecMESH%n_neighbor_pe))
      allocate (wSE(0:hecMESH%n_neighbor_pe),                           &
     &          wSI(0:hecMESH%n_neighbor_pe))
      IW1 = 0
      wSE = 0
      wSI = 0

!C
!C-- search IMPORT items
      do icel= 1, hecMESH%n_elem
        ih= hecMESH%elem_ID(2*icel)
        if (ih.ne.hecMESH%my_rank) then
              ihR = hecMESH%rev_neighbor_pe(ih)
          wSI(ihR)= wSI(ihR) + 1
        endif
      enddo
!C
!C-- exchange INFO. on IMPORT item #
      stime0= MPI_WTIME ()
      call hecmw_adapt_STACK_SEND_RECV                                  &
     &   ( hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                  &
     &     wSI, wSE, hecMESH%MPI_COMM, hecMESH%my_rank)
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
      allocate (wiI(wSI(hecMESH%n_neighbor_pe)))
      allocate (wiE(wSE(hecMESH%n_neighbor_pe)))
! JF
      wiI = 0
      wiE = 0

      do icel= 1, hecMESH%n_elem
        ih= hecMESH%elem_ID(2*icel)
        if (ih.ne.hecMESH%my_rank) then
                   ihR             = hecMESH%rev_neighbor_pe(ih)
               IW1(ihR  )          = IW1(ihR) + 1
           wiI(wSI(ihR-1)+IW1(ihR))= hecMESH%elem_ID(2*icel-1)
        endif
      enddo

      LEN= max ( wSI(hecMESH%n_neighbor_pe), wSE(hecMESH%n_neighbor_pe),&
     &           hecMESH%n_elem)
      deallocate (IW1)
        allocate (IW1(LEN), IW2(LEN))
! JF
      IW1 = 0
      IW2 = 0

      stime0= MPI_WTIME ()
      call hecmw_adapt_ITEM_SEND_RECV                                   &      
     &    (LEN, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
     &     wSI, wiI, wSE, wiE, IW1, IW2,                                &
     &     hecMESH%MPI_COMM, hecMESH%my_rank, 1)
      etime0= MPI_WTIME ()
      commTIME= commTIME + etime0 - stime0

!C
!C-- find LOCAL CELL ID at DESTINY
      do neib= 1, hecMESH%n_neighbor_pe
        is= wSE(neib-1)+1
        ie= wSE(neib  )
        ip= hecMESH%neighbor_pe(neib)
        do k= is, ie
          icel= wiE(k)
          call hecmw_adapt_LOCAL_CELL_INFO (icel, icel0)
          wiE(k)= icel0
        enddo
      enddo

!C
!C-- reconstruct IMPORT table
      deallocate (IW1)
        allocate (IW1(hecMESH%n_neighbor_pe))

      IW1= 0
      do icel= 1, hecMESH%n_elem
        ih= hecMESH%elem_ID(2*icel)
        if (ih.ne.hecMESH%my_rank) then
                   ihR             = hecMESH%rev_neighbor_pe(ih)
               IW1(ihR  )          = IW1(ihR) + 1
           wiI(wSI(ihR-1)+IW1(ihR))= icel
        endif
      enddo

!C
!C-- new ARRAY
      allocate (hecMESH%adapt_import_elem_index(0:hecMESH%n_neighbor_pe))
      allocate (hecMESH%adapt_export_elem_index(0:hecMESH%n_neighbor_pe))

      do neib= 0, hecMESH%n_neighbor_pe
        hecMESH%adapt_import_elem_index(neib)= wSI(neib)
        hecMESH%adapt_export_elem_index(neib)= wSE(neib)
      enddo

      MAXimport= wSI(hecMESH%n_neighbor_pe)
      MAXexport= wSE(hecMESH%n_neighbor_pe)

      allocate (hecMESH%adapt_import_elem_item(MAXimport))
      allocate (hecMESH%adapt_export_elem_item(MAXexport))
      do k= 1, MAXimport
        hecMESH%adapt_import_elem_item(k)= wiI(k)
      enddo
      do k= 1, MAXexport
        hecMESH%adapt_export_elem_item(k)= wiE(k)
      enddo

      deallocate (IW1,IW2)
      deallocate (wSE,wSI,wiE,wiI)

!C===
      contains
        subroutine hecmw_adapt_LOCAL_CELL_INFO (icel, in0)
          do i0= 1, hecMESH%ne_internal
            i= hecMESH%elem_internal_list(i0)
            if (icel.eq.hecMESH%elem_ID(2*i-1)) then
              in0= i
              return
            endif
          enddo
        end subroutine hecmw_adapt_LOCAL_CELL_INFO
      end subroutine hecmw_adapt_CELL_COMM_TABLE
