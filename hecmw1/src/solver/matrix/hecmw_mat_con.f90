!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_matrix_con
  use hecmw_util
  use hecmw_varray_int
  use hecmw_matrix_contact
  use hecmw_etype
  implicit none

  private

  public :: hecmw_mat_con

contains
  !C***
  !C*** MAT_CON for solver
  !C***
  !C
  subroutine hecmw_mat_con ( hecMESH, hecMAT )
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH

    type (hecmwST_varray_int), allocatable :: CLU(:), CLL(:)
    integer(kind=kint)  :: itype, iS, iE, iSS, iSE, ic_type, icel
    integer(kind=kint)  :: i, j, k, kk
    integer(kind=kint)  :: nn, nid(2048), nid1, nid2

    !C
    !C +-------+
    !C | INIT. |
    !C +-------+
    !C===
    hecMAT%NP= hecMESH%n_node
    hecMAT%N = hecMESH%nn_internal

    !C===
    !C
    !C +----------------------------------------+
    !C | CONNECTIVITY according to ELEMENT TYPE |
    !C +----------------------------------------+
    !C===
    call HECMW_varray_int_initialize_all( CLU, hecMAT%NP, 27 )
    call HECMW_varray_int_initialize_all( CLL, hecMAT%NP, 27 )
    do itype= 1, hecMESH%n_elem_type
      is= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if( hecmw_is_etype_patch(ic_type) ) cycle
      !C element loop
      do icel= is, iE
        iSS = hecMESH%elem_node_index(icel-1)
        iSE = hecMESH%elem_node_index(icel  )
        nn = iSE-iSS
        do j=1,nn
          nid(j) = hecMESH%elem_node_item (iSS+j)
        enddo
        do j=1,nn
          nid1 = nid(j)
          do k=1,nn
            nid2 = nid(k)
            if( nid1 < nid2 ) then
              call HECMW_varray_int_insert_if_not_exists( CLU(nid1), nid2 )
            else if( nid1 > nid2 ) then
              call HECMW_varray_int_insert_if_not_exists( CLL(nid1), nid2 )
            end if
          enddo
        enddo
      enddo
    enddo

    !C===
    !C
    !C +---------------------------+
    !C | ALLOCATE hecMAT structure |
    !C +---------------------------+
    !C===
    allocate (hecMAT%indexL(0:hecMAT%NP), hecMAT%indexU(0:hecMAT%NP))

    hecMAT%indexL = 0
    hecMAT%indexU = 0
    do i = 1, hecMAT%NP
      hecMAT%indexL(i) = hecMAT%indexL(i-1) + HECMW_varray_int_get_nitem(CLL(i))
      hecMAT%indexU(i) = hecMAT%indexU(i-1) + HECMW_varray_int_get_nitem(CLU(i))
    enddo

    hecMAT%NPL = hecMAT%indexL(hecMAT%NP)
    hecMAT%NPU = hecMAT%indexU(hecMAT%NP)

    allocate (hecMAT%itemL(hecMAT%NPL), hecMAT%itemU(hecMAT%NPU))

    do i=1,hecMAT%NP
      do k=1,HECMW_varray_int_get_nitem(CLL(i))
        kk = k + hecMAT%indexL(i-1)
        hecMAT%itemL(kk) = HECMW_varray_int_get_item(CLL(i), k)
      enddo
      do k=1,HECMW_varray_int_get_nitem(CLU(i))
        kk = k + hecMAT%indexU(i-1)
        hecMAT%itemU(kk) = HECMW_varray_int_get_item(CLU(i), k)
      enddo
    enddo

    call hecmw_cmat_init (hecMAT%cmat)

    call HECMW_varray_int_finalize_all( CLU )
    call HECMW_varray_int_finalize_all( CLL )

  end subroutine

end module hecmw_matrix_con
