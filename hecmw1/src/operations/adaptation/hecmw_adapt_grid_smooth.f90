!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

!C
!C***
!C*** hecmw_adapt_GRID_SMOOTH
!C***
!C
!C    grid smoothing for embedding
!C
subroutine hecmw_adapt_GRID_SMOOTH (hecMESH)
  use hecmw_util
  use hecmw_adapt_INT_SR
  use hecmw_adapt_INT_SR_REV

  implicit real*8 (A-H,O-Z)
  integer(kind=kint), dimension(:), allocatable :: WR, WS

  dimension NDIV(6)
  integer(kind=kint), dimension(:), allocatable :: NFLAG_INFO

  type (hecmwST_local_mesh) :: hecMESH

  !C
  !C-- INIT.
  NITERADAP= 0
  NEWtet   = 0
  NEWprism = 0
  allocate (NFLAG_INFO(hecMESH%PETOT))

  NITERADAP_MAX= 10000

  90 continue
  call MPI_BARRIER  (hecMESH%MPI_COMM, ierr)

  NITERADAP= NITERADAP + 1
  NF0      = 0
  if ( NITERADAP .gt. NITERADAP_MAX) then
    call hecmw_adapt_ERROR_EXIT(hecMESH, 7)
  endif

  !C
  !C +------------+
  !C | TETRAHEDRA |
  !C +------------+
  !C
  !C   RULEs :
  !C     a. READ the PAPER for the PATTERNs
  !C     b. CONSECUTIVE DIRECTIONAL REFINEMENT IS NOT ALLOWED
  !C***********************************************************************

  if (hecMESH%my_rank.eq.0)                                                 &
    &    write (*,'("  TETRA iteration=", 2i8)') NITERADAP, NEWtet

  NEWtet= 0
  do 100 icel0= 1, hecMESH%n_adapt_act_elem_341
    icel= hecMESH%adapt_act_elem_341(icel0)
    is= hecMESH%elem_node_index(icel-1)
    n1= hecMESH%elem_node_item (is+1)
    n2= hecMESH%elem_node_item (is+2)
    n3= hecMESH%elem_node_item (is+3)
    n4= hecMESH%elem_node_item (is+4)

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

    !C
    !C +--------------------------+
    !C | ADJUST the CELL DIVISION |
    !C +--------------------------+
    !C===
    if (NDIVSUM .eq. 1) goto 95
    !C
    !C== 2 edges
    if (NDIVSUM .eq. 2) then
      if ( ( (NDIV(1).eq.1) .and. (NDIV(6).eq.1) ) .or.             &
          &         ( (NDIV(2).eq.1) .and. (NDIV(5).eq.1) ) .or.             &
          &         ( (NDIV(3).eq.1) .and. (NDIV(4).eq.1) ) ) then
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        NEWtet    = NEWtet + 1
        NF0       = 1
        goto 95
        !C
        !C-- PATTERN 1-3-5
      else if                                                       &
          &         ( (NDIV(1).eq.1) .and. (NDIV(3).eq.1) ) then
        hecMESH%adapt_iemb(ie5)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(1).eq.1) .and. (NDIV(5).eq.1) ) then
        hecMESH%adapt_iemb(ie3)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(3).eq.1) .and. (NDIV(5).eq.1) ) then
        hecMESH%adapt_iemb(ie1)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
        !C
        !C-- PATTERN 2-3-6
      else if                                                       &
          &         ( (NDIV(2).eq.1) .and. (NDIV(3).eq.1) ) then
        hecMESH%adapt_iemb(ie6)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(2).eq.1) .and. (NDIV(6).eq.1) ) then
        hecMESH%adapt_iemb(ie3)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(3).eq.1) .and. (NDIV(6).eq.1) ) then
        hecMESH%adapt_iemb(ie2)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
        !C
        !C-- PATTERN 1-2-4
      else if                                                       &
          &         ( (NDIV(1).eq.1) .and. (NDIV(2).eq.1) ) then
        hecMESH%adapt_iemb(ie4)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(1).eq.1) .and. (NDIV(4).eq.1) ) then
        hecMESH%adapt_iemb(ie2)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(2).eq.1) .and. (NDIV(4).eq.1) ) then
        hecMESH%adapt_iemb(ie1)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
        !C
        !C-- PATTERN 4-5-6
      else if                                                       &
          &         ( (NDIV(4).eq.1) .and. (NDIV(5).eq.1) ) then
        hecMESH%adapt_iemb(ie6)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(4).eq.1) .and. (NDIV(6).eq.1) ) then
        hecMESH%adapt_iemb(ie5)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else if                                                       &
          &         ( (NDIV(5).eq.1) .and. (NDIV(6).eq.1) ) then
        hecMESH%adapt_iemb(ie4)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      else
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      endif
    endif
    !C==

    !C
    !C== 3 edges
    if (NDIVSUM .eq. 3) then
      if (                                                          &
          &       ((NDIV(1).eq.1).and.(NDIV(3).eq.1).and.(NDIV(5).eq.1)) .or.&
          &       ((NDIV(2).eq.1).and.(NDIV(3).eq.1).and.(NDIV(6).eq.1)) .or.&
          &       ((NDIV(1).eq.1).and.(NDIV(2).eq.1).and.(NDIV(4).eq.1)) .or.&
          &       ((NDIV(4).eq.1).and.(NDIV(5).eq.1).and.(NDIV(6).eq.1)) )   &
          &         then
        goto 95
      else
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        NEWtet   = NEWtet + 1
        NF0      = 1
        goto 95
      endif
    endif
    !C==

    !C
    !C== more than 4 edges
    if (NDIVSUM.eq.4 .or. NDIVSUM.eq.5) then
      hecMESH%adapt_iemb(ie1)= 1
      hecMESH%adapt_iemb(ie2)= 1
      hecMESH%adapt_iemb(ie3)= 1
      hecMESH%adapt_iemb(ie4)= 1
      hecMESH%adapt_iemb(ie5)= 1
      hecMESH%adapt_iemb(ie6)= 1
      NEWtet   = NEWtet + 1
      NF0      = 1
      goto 95
    endif
    !C==

    !C
    !C== check the type of PARENT cell
    95   continue

    NTYP= hecMESH%adapt_parent_type(icel)
    if (NTYP.ne.0 .and. NTYP.ne.11 .and.                            &
        &      NDIVSUM.ne.0.and.NDIVSUM.ne.6) then
      hecMESH%adapt_iemb(ie1)= 1
      hecMESH%adapt_iemb(ie2)= 1
      hecMESH%adapt_iemb(ie3)= 1
      hecMESH%adapt_iemb(ie4)= 1
      hecMESH%adapt_iemb(ie5)= 1
      hecMESH%adapt_iemb(ie6)= 1
      NF0      = 1
    endif

    !C==
    100 continue
    !C***********************************************************************

    !C
    !C +--------+
    !C | PRISMs |
    !C +--------+
    !C
    !C   RULEs :
    !C     a. READ the PAPER for the PATTERNs
    !C     b. CONSECUTIVE DIRECTIONAL REFINEMENT IS NOT ALLOWED
    !C***********************************************************************

    if (hecMESH%my_rank.eq.0)                                                 &
      &    write (*,'("  PRISM iteration=", 2i8)') NITERADAP, NEWprism

    !C
    !C-- ADJUST normal-to-surface direction
    N1= hecMESH%adapt_import_edge_index(hecMESH%n_neighbor_pe)
    N2= hecMESH%adapt_export_edge_index(hecMESH%n_neighbor_pe)
    m = max (N1, N2)
    allocate (WS(m), WR(m))

    icouM= 0
    do layer= 1, hecMESH%n_adapt_act_elem_351
      icou = 0
      do icel0= 1, hecMESH%n_adapt_act_elem_351
        icel= hecMESH%adapt_act_elem_351(icel0)
        is= hecMESH%elem_node_index(icel-1)
        n1= hecMESH%elem_node_item (is+1)
        n2= hecMESH%elem_node_item (is+2)
        n3= hecMESH%elem_node_item (is+3)
        n4= hecMESH%elem_node_item (is+4)
        n5= hecMESH%elem_node_item (is+5)
        n6= hecMESH%elem_node_item (is+6)

        call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie2, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n1, ie3, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n4, n5, ie4, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n5, n6, ie5, 1 )
        call hecmw_adapt_EDGE_INFO ( hecMESH, n6, n4, ie6, 1 )

        if (hecMESH%adapt_iemb(ie1).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie4).eq.0) then
          hecMESH%adapt_iemb(ie4)= 1
          icou     = icou + 1
        endif
        if (hecMESH%adapt_iemb(ie2).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie5).eq.0) then
          hecMESH%adapt_iemb(ie5)= 1
          icou     = icou + 1
        endif
        if (hecMESH%adapt_iemb(ie3).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie6).eq.0) then
          hecMESH%adapt_iemb(ie6)= 1
          icou     = icou + 1
        endif
        if (hecMESH%adapt_iemb(ie4).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie1).eq.0) then
          hecMESH%adapt_iemb(ie1)= 1
          icou     = icou + 1
        endif

        if (hecMESH%adapt_iemb(ie5).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie2).eq.0) then
          hecMESH%adapt_iemb(ie2)= 1
          icou     = icou + 1
        endif
        if (hecMESH%adapt_iemb(ie6).eq.1 .and.                        &
            &        hecMESH%adapt_iemb(ie3).eq.0) then
          hecMESH%adapt_iemb(ie3)= 1
          icou     = icou + 1
        endif
      enddo

      call MPI_allREDUCE ( icouM, icouMmin, 1, MPI_INTEGER,           &
        &                       MPI_MIN, hecMESH%MPI_COMM, ierr)

      if (icou.eq.0 .and. icouMmin.eq.1) exit
      if (icou.eq.0) icouM= 1
      if (icou.ne.0) icouM= 0

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV_REV                              &
        &     ( hecMESH%n_adapt_edge,                                      &
        &       hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                &
        &       hecMESH%adapt_import_edge_index,                           &
        &       hecMESH%adapt_import_edge_item ,                           &
        &       hecMESH%adapt_export_edge_index,                           &
        &       hecMESH%adapt_export_edge_item ,                           &
        &       WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,              &
        &       hecMESH%my_rank, 1, m)

      WS= 0
      WR= 0
      call hecmw_adapt_INT_SEND_RECV                                  &
        &     ( hecMESH%n_adapt_edge,                                      &
        &       hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                &
        &       hecMESH%adapt_import_edge_index,                           &
        &       hecMESH%adapt_import_edge_item ,                           &
        &       hecMESH%adapt_export_edge_index,                           &
        &       hecMESH%adapt_export_edge_item ,                           &
        &       WS, WR, hecMESH%adapt_iemb, hecMESH%MPI_COMM,              &
        &       hecMESH%my_rank, 1, m)

    enddo
    deallocate (WS, WR)

    call MPI_BARRIER  (hecMESH%MPI_COMM,ierr)

    !C
    !C-- real operation for prismatic region
    NEWprism= 0
    do 110 icel0= 1, hecMESH%n_adapt_act_elem_351
      icel= hecMESH%adapt_act_elem_351(icel0)
      is= hecMESH%elem_node_index(icel-1)
      n1= hecMESH%elem_node_item (is+1)
      n2= hecMESH%elem_node_item (is+2)
      n3= hecMESH%elem_node_item (is+3)
      n4= hecMESH%elem_node_item (is+4)
      n5= hecMESH%elem_node_item (is+5)
      n6= hecMESH%elem_node_item (is+6)

      call hecmw_adapt_EDGE_INFO ( hecMESH, n1, n2, ie1, 1 )
      call hecmw_adapt_EDGE_INFO ( hecMESH, n2, n3, ie2, 1 )
      call hecmw_adapt_EDGE_INFO ( hecMESH, n3, n1, ie3, 1 )
      call hecmw_adapt_EDGE_INFO ( hecMESH, n4, n5, ie4, 1 )
      call hecmw_adapt_EDGE_INFO ( hecMESH, n5, n6, ie5, 1 )
      call hecmw_adapt_EDGE_INFO ( hecMESH, n6, n4, ie6, 1 )

      NDIV(1)= hecMESH%adapt_iemb(ie1)
      NDIV(2)= hecMESH%adapt_iemb(ie2)
      NDIV(3)= hecMESH%adapt_iemb(ie3)
      NDIV(4)= hecMESH%adapt_iemb(ie4)
      NDIV(5)= hecMESH%adapt_iemb(ie5)
      NDIV(6)= hecMESH%adapt_iemb(ie6)

      NDIVSUM= NDIV(1)+NDIV(2)+NDIV(3)+NDIV(4)+NDIV(5)+NDIV(6)

      !C
      !C +--------------------------+
      !C | ADJUST the CELL DIVISION |
      !C +--------------------------+
      !C===

      !C
      !C-- 1 edge(s)
      if (NDIVSUM .eq. 1) then
        if (NDIV(1).eq.1) hecMESH%adapt_iemb(ie4)= 1
        if (NDIV(2).eq.1) hecMESH%adapt_iemb(ie5)= 1
        if (NDIV(3).eq.1) hecMESH%adapt_iemb(ie6)= 1
        if (NDIV(4).eq.1) hecMESH%adapt_iemb(ie1)= 1
        if (NDIV(5).eq.1) hecMESH%adapt_iemb(ie2)= 1
        if (NDIV(6).eq.1) hecMESH%adapt_iemb(ie3)= 1
        goto 105
      endif
      !C
      !C-- 2 edges
      if (NDIVSUM.eq.2) then
        if (NDIV(1).eq.1 .and. NDIV(4).eq.1) goto 105
        if (NDIV(2).eq.1 .and. NDIV(5).eq.1) goto 105
        if (NDIV(3).eq.1 .and. NDIV(6).eq.1) goto 105
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        goto 105
      endif

      !C
      !C-- >3 edges
      if (NDIVSUM.ge.3) then
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        goto 105
      endif
      !C
      !C-- check the type of PARENT cell
      105  continue

      NTYP= hecMESH%adapt_parent_type(icel)
      if (NTYP.ne.0 .and. NTYP.ne.4 .and.                             &
          &      NDIVSUM.ne.0.and.NDIVSUM.ne.6) then
        hecMESH%adapt_iemb(ie1)= 1
        hecMESH%adapt_iemb(ie2)= 1
        hecMESH%adapt_iemb(ie3)= 1
        hecMESH%adapt_iemb(ie4)= 1
        hecMESH%adapt_iemb(ie5)= 1
        hecMESH%adapt_iemb(ie6)= 1
        NF0      = 1
      endif

      !C==
      110 continue
      !C***********************************************************************

      !C
      !C-- adjust EMBEDDING level
      call MPI_BARRIER  (hecMESH%MPI_COMM, ierr)
      call hecmw_adapt_ADJEMB (hecMESH, NF0)

      call MPI_BARRIER (hecMESH%MPI_COMM, ierr)
      call MPI_GATHER                                                   &
        &      (NF0, 1, MPI_INTEGER, NFLAG_INFO, 1, MPI_INTEGER, 0,        &
        &       hecMESH%MPI_COMM, ierr)

      if (hecMESH%my_rank.eq.0) then
        icou= 0
        do i= 1, hecMESH%PETOT
          icou= icou + NFLAG_INFO(i)
        enddo
        if (icou.ne.0) NF0= 1
      endif

      call MPI_BCAST  (NF0, 1, MPI_INTEGER, 0, hecMESH%MPI_COMM, ierr)
      if (NF0 .eq. 1) goto 90

      return
end

