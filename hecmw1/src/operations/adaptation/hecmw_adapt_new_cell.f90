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
!C*** hecmw_adapt_NEW_CELL
!C***
!C
!C    control NEW_CELL_TETRA/PRISM
!C
      subroutine hecmw_adapt_NEW_CELL (hecMESH)

      use hecmw_util
      use hecmw_adapt_INT_SR
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(:), allocatable :: WR, WS
      integer(kind=kint), dimension(:), allocatable :: IW1, IW2
      integer(kind=kint), dimension(:), allocatable :: IW3, IW4

      type (hecmwST_local_mesh) :: hecMESH

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      hecMESH%n_adapt_elem_cur= hecMESH%n_elem
      hecMESH%n_adapt_elem_old= hecMESH%n_elem

      hecMESH%n_adapt_elem_341_cur= hecMESH%n_adapt_elem_341
      hecMESH%n_adapt_elem_351_cur= hecMESH%n_adapt_elem_351

      icouN= 0

      allocate (hecMESH%adapt_children_local(8*hecMESH%ne_array))
      hecMESH%adapt_children_local= 0
!C===

!C
!C +------------------+
!C | create NEW CELLs |
!C +------------------+
!C===
      call hecmw_adapt_NEW_CELL_341 (hecMESH, icouN,                    &
     &                               hecMESH%n_adapt_elem_341_cur)
      call hecmw_adapt_NEW_CELL_351 (hecMESH, icouN,                    &
     &                               hecMESH%n_adapt_elem_351_cur)
!C===

!C
!C +---------------+
!C | GLOBAL CELL # |
!C +---------------+
!C===
      hecMESH%n_elem          = hecMESH%n_adapt_elem_cur
      hecMESH%n_adapt_elem_341= hecMESH%n_adapt_elem_341_cur
      hecMESH%n_adapt_elem_351= hecMESH%n_adapt_elem_351_cur
      
!C
!C-- OLD/NEW global CELL #

      ICELTOTg_old= 0
      ICELTOTg_cur= 0

      nnn1= hecMESH%ne_internal + icouN

      call MPI_allREDUCE (hecMESH%ne_internal, ICELTOTg_old, 1,         &
     &                    MPI_INTEGER, MPI_SUM, hecMESH%MPI_COMM, ierr)
      call MPI_allREDUCE (nnn1               , ICELTOTg_cur, 1,         &
     &                    MPI_INTEGER, MPI_SUM, hecMESH%MPI_COMM, ierr)

      if (hecMESH%my_rank.eq.0) then
       write (*,'("  total cell number (before)", i8  )') ICELTOTg_old
       write (*,'("  total cell number (curr. )", i8,/)') ICELTOTg_cur
      endif
!C===

!C
!C +-------------------------+
!C | exchange CHILDREN-info. |
!C +-------------------------+
!C===
      LEN= max (hecMESH%n_adapt_elem_old,                               &
     &          hecMESH%adapt_import_elem_index(hecMESH%n_neighbor_pe), &
     &          hecMESH%adapt_export_elem_index(hecMESH%n_neighbor_pe))

      allocate (IW1(LEN*8))
      IW1= 0
      do icel= 1, hecMESH%n_adapt_elem_old
        ityp = hecMESH%adapt_type          (icel)
          iS = hecMESH%adapt_children_index(icel-1)
          iS1= iS + 1
          iS2= iS + 2
          iS3= iS + 3
          iS4= iS + 4
          iS5= iS + 5
          iS6= iS + 6
          iS7= iS + 7
          iS8= iS + 8
!C
!C-- TETRAHEDRA
        if (hecMESH%elem_type(icel).eq.341) then
        if (ityp.ne.0) then
          IW1(8*icel-7)= hecMESH%adapt_children_item(2*iS1-1)
          IW1(8*icel-6)= hecMESH%adapt_children_item(2*iS2-1)
          if (ityp.ge.7) then
            IW1(8*icel-5)= hecMESH%adapt_children_item(2*iS3-1)
            IW1(8*icel-4)= hecMESH%adapt_children_item(2*iS4-1)
            if (ityp.eq.11) then
              IW1(8*icel-3)= hecMESH%adapt_children_item(2*iS5-1)
              IW1(8*icel-2)= hecMESH%adapt_children_item(2*iS6-1)
              IW1(8*icel-1)= hecMESH%adapt_children_item(2*iS7-1)
              IW1(8*icel  )= hecMESH%adapt_children_item(2*iS8-1)
            endif
          endif
        endif
        endif
!C
!C-- PRISMs
        if (hecMESH%elem_type(icel).eq.351) then
        if (ityp.ne.0) then
          IW1(8*icel-7)= hecMESH%adapt_children_item(2*iS1-1)
          IW1(8*icel-6)= hecMESH%adapt_children_item(2*iS2-1)
          if (ityp.eq.4) then
            IW1(8*icel-5)= hecMESH%adapt_children_item(2*iS3-1)
            IW1(8*icel-4)= hecMESH%adapt_children_item(2*iS4-1)
          endif
        endif
        endif
      enddo

      allocate (WS(8*LEN), WR(8*LEN))
      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( LEN, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
     &     hecMESH%adapt_import_elem_index,                             &
     &     hecMESH%adapt_import_elem_item ,                             &
     &     hecMESH%adapt_export_elem_index,                             &
     &     hecMESH%adapt_export_elem_item ,                             &
     &     WS, WR, IW1, hecMESH%MPI_COMM, hecMESH%my_rank, 8, LEN)
      deallocate (WS, WR)

      do icel= 1, hecMESH%n_adapt_elem_old
        ityp = hecMESH%adapt_type          (icel)
          iS = hecMESH%adapt_children_index(icel-1)
          iS1= iS + 1
          iS2= iS + 2
          iS3= iS + 3
          iS4= iS + 4
          iS5= iS + 5
          iS6= iS + 6
          iS7= iS + 7
          iS8= iS + 8
        if (ityp.ne.0 .and.                                             &
     &      hecMESH%elem_ID(2*icel).ne.hecMESH%my_rank) then
!C
!C-- TETRAHEDRA
        if (hecMESH%elem_type(icel).eq.341) then
          hecMESH%adapt_children_item(2*iS1-1)= IW1(8*icel-7)
          hecMESH%adapt_children_item(2*iS2-1)= IW1(8*icel-6)
          if (ityp.ge.7) then
            hecMESH%adapt_children_item(2*iS3-1)= IW1(8*icel-5)
            hecMESH%adapt_children_item(2*iS4-1)= IW1(8*icel-4)
            if (ityp.eq.11) then
              hecMESH%adapt_children_item(2*iS5-1)= IW1(8*icel-3)
              hecMESH%adapt_children_item(2*iS6-1)= IW1(8*icel-2)
              hecMESH%adapt_children_item(2*iS7-1)= IW1(8*icel-1)
              hecMESH%adapt_children_item(2*iS8-1)= IW1(8*icel  )
            endif
          endif
        endif
!C
!C-- PRISMs
        if (hecMESH%elem_type(icel).eq.351) then
          hecMESH%adapt_children_item(2*iS1-1)= IW1(8*icel-7)
          hecMESH%adapt_children_item(2*iS2-1)= IW1(8*icel-6)
          if (ityp.eq.4) then
            hecMESH%adapt_children_item(2*iS3-1)= IW1(8*icel-5)
            hecMESH%adapt_children_item(2*iS4-1)= IW1(8*icel-4)
          endif
        endif
        endif
      enddo
      deallocate (IW1)

      allocate (IW1(LEN*8))
      IW1= 0
      do icel= 1, hecMESH%n_adapt_elem_old
        ityp= hecMESH%adapt_type(icel)
          iS= hecMESH%adapt_children_index(icel-1)
          iS1= iS + 1
          iS2= iS + 2
          iS3= iS + 3
          iS4= iS + 4
          iS5= iS + 5
          iS6= iS + 6
          iS7= iS + 7
          iS8= iS + 8
!C
!C-- TETRAHEDRA
        if (hecMESH%elem_type(icel).eq.341) then
        if (ityp.ne.0) then
          IW1(8*icel-7)= hecMESH%adapt_children_item(2*iS1)
          IW1(8*icel-6)= hecMESH%adapt_children_item(2*iS2)
          if (ityp.ge.7) then
            IW1(8*icel-5)= hecMESH%adapt_children_item(2*iS3)
            IW1(8*icel-4)= hecMESH%adapt_children_item(2*iS4)
            if (ityp.eq.11) then
              IW1(8*icel-3)= hecMESH%adapt_children_item(2*iS5)
              IW1(8*icel-2)= hecMESH%adapt_children_item(2*iS6)
              IW1(8*icel-1)= hecMESH%adapt_children_item(2*iS7)
              IW1(8*icel  )= hecMESH%adapt_children_item(2*iS8)
            endif
          endif
        endif
        endif
!C
!C-- PRISMs
        if (hecMESH%elem_type(icel).eq.351) then
        if (ityp.ne.0) then
          IW1(8*icel-7)= hecMESH%adapt_children_item(2*iS1)
          IW1(8*icel-6)= hecMESH%adapt_children_item(2*iS2)
          if (ityp.eq.4) then
            IW1(8*icel-5)= hecMESH%adapt_children_item(2*iS3)
            IW1(8*icel-4)= hecMESH%adapt_children_item(2*iS4)
          endif
        endif
        endif
      enddo

      allocate (WS(8*LEN), WR(8*LEN))
      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                    &
     &   ( LEN, hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,             &
     &     hecMESH%adapt_import_elem_index,                             &
     &     hecMESH%adapt_import_elem_item ,                             &
     &     hecMESH%adapt_export_elem_index,                             &
     &     hecMESH%adapt_export_elem_item ,                             &
     &     WS, WR, IW1, hecMESH%MPI_COMM, hecMESH%my_rank, 8, LEN)
      deallocate (WS, WR)

      do icel= 1, hecMESH%n_adapt_elem_old
        ityp= hecMESH%adapt_type(icel)
          iS= hecMESH%adapt_children_index(icel-1)
          iS1= iS + 1
          iS2= iS + 2
          iS3= iS + 3
          iS4= iS + 4
          iS5= iS + 5
          iS6= iS + 6
          iS7= iS + 7
          iS8= iS + 8
        if (ityp.ne.0 .and.                                             &
     &      hecMESH%elem_ID(2*icel) .ne. hecMESH%my_rank) then
!C
!C-- TETRAHEDRA
        if (hecMESH%elem_type(icel).eq.341) then
          hecMESH%adapt_children_item(2*iS1)= IW1(8*icel-7)
          hecMESH%adapt_children_item(2*iS2)= IW1(8*icel-6)
          if (ityp.ge.7) then
            hecMESH%adapt_children_item(2*iS3)= IW1(8*icel-5)
            hecMESH%adapt_children_item(2*iS4)= IW1(8*icel-4)
            if (ityp.eq.11) then
              hecMESH%adapt_children_item(2*iS5)= IW1(8*icel-3)
              hecMESH%adapt_children_item(2*iS6)= IW1(8*icel-2)
              hecMESH%adapt_children_item(2*iS7)= IW1(8*icel-1)
              hecMESH%adapt_children_item(2*iS8)= IW1(8*icel  )
            endif
          endif
        endif

!C
!C-- PRISMs
        if (hecMESH%elem_type(icel).eq.351) then
          hecMESH%adapt_children_item(2*iS1)= IW1(8*icel-7)
          hecMESH%adapt_children_item(2*iS2)= IW1(8*icel-6)
          if (ityp.eq.4) then
            hecMESH%adapt_children_item(2*iS3)= IW1(8*icel-5)
            hecMESH%adapt_children_item(2*iS4)= IW1(8*icel-4)
          endif
        endif
        endif
      enddo
      deallocate (IW1)
!C===

!C
!C +---------------+
!C | INTERNAL CELL |
!C +---------------+
!C===
      deallocate (hecMESH%elem_internal_list)
      icou= 0
      do icel= 1, hecMESH%n_elem
        if (hecMESH%elem_ID(2*icel).eq.hecMESH%my_rank) icou= icou + 1
      enddo

      hecMESH%ne_internal= icou
      allocate  (hecMESH%elem_internal_list(icou))
      icou= 0
      do icel= 1, hecMESH%n_elem
        if (hecMESH%elem_ID(2*icel).eq.hecMESH%my_rank) then
          icou= icou + 1
          hecMESH%elem_internal_list(icou)= icel
        endif
      enddo
!C===

!C
!C +------------+
!C | REORDERING |
!C +------------+
!C===
      allocate (hecMESH%adapt_OLDtoNEW_elem(hecMESH%n_elem))
      allocate (hecMESH%adapt_NEWtoOLD_elem(hecMESH%n_elem))
      allocate (IW2(hecMESH%n_elem))
      allocate (IW3(hecMESH%n_elem))
      hecMESH%adapt_OLDtoNEW_elem= 0
      hecMESH%adapt_NEWtoOLD_elem= 0
      IW2= 0
      IW3= 0

      icou_341= 0
      icou_351= 0

      do icel= 1, hecMESH%n_elem
        ityp= hecMESH%elem_type(icel)
        if (ityp.eq.341) then
          icou_341= icou_341 + 1
          IW3(icou_341)= icel
        endif        
      enddo

       do icel= 1, hecMESH%n_elem
         ityp= hecMESH%elem_type(icel)
        if (ityp.eq.351) then
          icou_351= icou_351 + 1
          IW3(icou_341+icou_351)= icel
        endif        
      enddo

      do ic0= 1, hecMESH%n_elem
        icel = IW3(ic0)
        hecMESH%adapt_OLDtoNEW_elem(icel)= ic0
        hecMESH%adapt_NEWtoOLD_elem(ic0 )= icel
      enddo

      hecMESH%elem_type_index(1)= icou_341
      hecMESH%elem_type_index(2)= icou_341 + icou_351 

      deallocate (IW2, IW3)
!C===
      return
      end subroutine hecmw_adapt_NEW_CELL



