!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_matrix_con
  use hecmw_util
#ifdef _OPENACC
  use hecmw_varray_int_gpu
#else
  use hecmw_varray_int
#endif
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

#ifdef _OPENACC
    type (hecmwST_varray_int_gpu_container) :: varray_container
    type (hecmwST_varray_int_gpu), allocatable :: CLU(:), CLL(:)
#else
    type (hecmwST_varray_int), allocatable :: CLU(:), CLL(:)
#endif
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
#ifdef _OPENACC
    call HECMW_varray_int_gpu_container_initialize(varray_container)
#endif

    !C===
    !C
    !C +----------------------------------------+
    !C | CONNECTIVITY according to ELEMENT TYPE |
    !C +----------------------------------------+
    !C===
#ifdef _OPENACC
    call HECMW_varray_int_gpu_initialize_all( varray_container, CLU, hecMAT%NP, 27 )
    call HECMW_varray_int_gpu_initialize_all( varray_container, CLL, hecMAT%NP, 27 )
#else
    call HECMW_varray_int_initialize_all( CLU, hecMAT%NP, 27 )
    call HECMW_varray_int_initialize_all( CLL, hecMAT%NP, 27 )
#endif
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
#ifdef _OPENACC
              call HECMW_varray_int_gpu_insert_if_not_exists( varray_container, CLU(nid1), nid2 )
#else
              call HECMW_varray_int_insert_if_not_exists( CLU(nid1), nid2 )
#endif
            else if( nid1 > nid2 ) then
#ifdef _OPENACC
              call HECMW_varray_int_gpu_insert_if_not_exists( varray_container, CLL(nid1), nid2 )
#else
              call HECMW_varray_int_insert_if_not_exists( CLL(nid1), nid2 )
#endif
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
#ifdef _OPENACC
      hecMAT%indexL(i) = hecMAT%indexL(i-1) + HECMW_varray_int_gpu_get_nitem(CLL(i))
      hecMAT%indexU(i) = hecMAT%indexU(i-1) + HECMW_varray_int_gpu_get_nitem(CLU(i))
#else
      hecMAT%indexL(i) = hecMAT%indexL(i-1) + HECMW_varray_int_get_nitem(CLL(i))
      hecMAT%indexU(i) = hecMAT%indexU(i-1) + HECMW_varray_int_get_nitem(CLU(i))
#endif
    enddo

    hecMAT%NPL = hecMAT%indexL(hecMAT%NP)
    hecMAT%NPU = hecMAT%indexU(hecMAT%NP)

    allocate (hecMAT%itemL(hecMAT%NPL), hecMAT%itemU(hecMAT%NPU))

    do i=1,hecMAT%NP
#ifdef _OPENACC
      do k=1,HECMW_varray_int_gpu_get_nitem(CLL(i))
        kk = k + hecMAT%indexL(i-1)
        hecMAT%itemL(kk) = HECMW_varray_int_gpu_get_item(CLL(i), k)
      enddo
      do k=1,HECMW_varray_int_gpu_get_nitem(CLU(i))
        kk = k + hecMAT%indexU(i-1)
        hecMAT%itemU(kk) = HECMW_varray_int_gpu_get_item(CLU(i), k)
      enddo
#else
      do k=1,HECMW_varray_int_get_nitem(CLL(i))
        kk = k + hecMAT%indexL(i-1)
        hecMAT%itemL(kk) = HECMW_varray_int_get_item(CLL(i), k)
      enddo
      do k=1,HECMW_varray_int_get_nitem(CLU(i))
        kk = k + hecMAT%indexU(i-1)
        hecMAT%itemU(kk) = HECMW_varray_int_get_item(CLU(i), k)
      enddo
#endif
    enddo

#ifdef _OPENACC
    call HECMW_varray_int_gpu_finalize_all( varray_container, CLU )
    call HECMW_varray_int_gpu_finalize_all( varray_container, CLL )
    call HECMW_varray_int_gpu_container_finalize(varray_container)
#else
    call HECMW_varray_int_finalize_all( CLU )
    call HECMW_varray_int_finalize_all( CLL )
#endif

  end subroutine

end module hecmw_matrix_con
