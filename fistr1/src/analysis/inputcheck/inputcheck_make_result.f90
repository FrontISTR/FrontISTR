!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module constructs result data for ELEMCHECK visualization/result output
module m_inputcheck_make_result

  implicit none

  private
  public :: inputcheck_make_result
  public :: inputcheck_write_result

contains

  !> Build hecmwST_result_data for visualization via hecmw_visualize
  subroutine inputcheck_make_result(hecMESH, fstrRESULT, elem_vol, elem_asp)
    use hecmw
    use m_fstr
    implicit none
    type(hecmwST_local_mesh), intent(in)     :: hecMESH
    type(hecmwST_result_data), intent(inout) :: fstrRESULT
    real(kind=kreal), intent(in)             :: elem_vol(:)
    real(kind=kreal), intent(in)             :: elem_asp(:)

    integer(kind=kint) :: n_ncomp, n_ecomp, n_nitem, n_eitem
    integer(kind=kint) :: i, j, ie, ig, jS, jE, idx, isect, cid
    integer(kind=kint) :: n_ngrp, n_egrp, n_sgrp

    n_ngrp = hecMESH%node_group%n_grp
    n_egrp = hecMESH%elem_group%n_grp
    n_sgrp = hecMESH%surf_group%n_grp

    ! Node components: node groups
    n_ncomp = n_ngrp
    n_nitem = n_ngrp  ! each is scalar (dof=1)

    ! Element components: elem groups + surf groups + MATERIAL_ID + SECTION_ID + VOLUME + ASPECT_RATIO
    n_ecomp = n_egrp + n_sgrp + 4
    n_eitem = n_ecomp  ! each is scalar (dof=1)

    ! Initialize result
    fstrRESULT%ng_component = 0
    fstrRESULT%nn_component = n_ncomp
    fstrRESULT%ne_component = n_ecomp

    ! Allocate global (none)
    allocate(fstrRESULT%ng_dof(0))
    allocate(fstrRESULT%global_label(0))
    allocate(fstrRESULT%global_val_item(0))

    ! Allocate node data
    allocate(fstrRESULT%nn_dof(n_ncomp))
    allocate(fstrRESULT%node_label(n_ncomp))
    allocate(fstrRESULT%node_val_item(n_nitem * hecMESH%n_node))
    fstrRESULT%nn_dof = 1
    fstrRESULT%node_val_item = 0.0d0

    ! Allocate element data
    allocate(fstrRESULT%ne_dof(n_ecomp))
    allocate(fstrRESULT%elem_label(n_ecomp))
    allocate(fstrRESULT%elem_val_item(n_eitem * hecMESH%n_elem))
    fstrRESULT%ne_dof = 1
    fstrRESULT%elem_val_item = 0.0d0

    ! --- Fill node group flags ---
    do ig = 1, n_ngrp
      fstrRESULT%node_label(ig) = 'NGRP_' // trim(hecMESH%node_group%grp_name(ig))
      jS = hecMESH%node_group%grp_index(ig-1)
      jE = hecMESH%node_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%node_group%grp_item(j)
        if(idx >= 1 .and. idx <= hecMESH%n_node) then
          fstrRESULT%node_val_item(n_nitem*(idx-1) + ig) = 1.0d0
        endif
      enddo
    enddo

    ! --- Fill element group flags ---
    do ig = 1, n_egrp
      fstrRESULT%elem_label(ig) = 'EGRP_' // trim(hecMESH%elem_group%grp_name(ig))
      jS = hecMESH%elem_group%grp_index(ig-1)
      jE = hecMESH%elem_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%elem_group%grp_item(j)
        if(idx >= 1 .and. idx <= hecMESH%n_elem) then
          fstrRESULT%elem_val_item(n_eitem*(idx-1) + ig) = 1.0d0
        endif
      enddo
    enddo

    ! --- Fill surface group (element with face number) ---
    do ig = 1, n_sgrp
      fstrRESULT%elem_label(n_egrp + ig) = 'SGRP_' // trim(hecMESH%surf_group%grp_name(ig))
      jS = hecMESH%surf_group%grp_index(ig-1)
      jE = hecMESH%surf_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%surf_group%grp_item(2*j-1)  ! element ID
        if(idx >= 1 .and. idx <= hecMESH%n_elem) then
          fstrRESULT%elem_val_item(n_eitem*(idx-1) + n_egrp + ig) = &
            dble(hecMESH%surf_group%grp_item(2*j))  ! face number
        endif
      enddo
    enddo

    ! --- Fill MATERIAL_ID ---
    i = n_egrp + n_sgrp + 1
    fstrRESULT%elem_label(i) = 'MATERIAL_ID'
    do ie = 1, hecMESH%n_elem
      isect = hecMESH%section_ID(ie)
      if(isect >= 1 .and. isect <= hecMESH%section%n_sect) then
        cid = hecMESH%section%sect_mat_ID_item(isect)
        fstrRESULT%elem_val_item(n_eitem*(ie-1) + i) = dble(cid)
      endif
    enddo

    ! --- Fill SECTION_ID ---
    i = n_egrp + n_sgrp + 2
    fstrRESULT%elem_label(i) = 'SECTION_ID'
    do ie = 1, hecMESH%n_elem
      fstrRESULT%elem_val_item(n_eitem*(ie-1) + i) = dble(hecMESH%section_ID(ie))
    enddo

    ! --- Fill VOLUME ---
    i = n_egrp + n_sgrp + 3
    fstrRESULT%elem_label(i) = 'VOLUME'
    do ie = 1, hecMESH%n_elem
      fstrRESULT%elem_val_item(n_eitem*(ie-1) + i) = elem_vol(ie)
    enddo

    ! --- Fill ASPECT_RATIO ---
    i = n_egrp + n_sgrp + 4
    fstrRESULT%elem_label(i) = 'ASPECT_RATIO'
    do ie = 1, hecMESH%n_elem
      fstrRESULT%elem_val_item(n_eitem*(ie-1) + i) = elem_asp(ie)
    enddo

  end subroutine inputcheck_make_result

  !> Write check results to result file via hecmw_result_add/write interface
  subroutine inputcheck_write_result(hecMESH, elem_vol, elem_asp)
    use hecmw
    use m_fstr
    use hecmw_result
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in)         :: elem_vol(:)
    real(kind=kreal), intent(in)         :: elem_asp(:)

    integer(kind=kint) :: ig, j, jS, jE, idx, ie, isect, cid
    integer(kind=kint) :: n_ngrp, n_egrp, n_sgrp
    character(len=HECMW_NAME_LEN) :: header, comment, label, nameID
    real(kind=kreal), allocatable :: work_n(:), work_e(:)

    allocate(work_n(hecMESH%n_node))
    allocate(work_e(hecMESH%n_elem))

    header  = '*fstrresult'
    comment = 'elemcheck_result'
    call hecmw_result_init(hecMESH, 1, header, comment)

    n_ngrp = hecMESH%node_group%n_grp
    n_egrp = hecMESH%elem_group%n_grp
    n_sgrp = hecMESH%surf_group%n_grp

    ! --- Node groups ---
    do ig = 1, n_ngrp
      work_n = 0.0d0
      jS = hecMESH%node_group%grp_index(ig-1)
      jE = hecMESH%node_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%node_group%grp_item(j)
        if(idx >= 1 .and. idx <= hecMESH%n_node) work_n(idx) = 1.0d0
      enddo
      label = 'NGRP_' // trim(hecMESH%node_group%grp_name(ig))
      call hecmw_result_add(HECMW_RESULT_DTYPE_NODE, 1, label, work_n)
    enddo

    ! --- Element groups ---
    do ig = 1, n_egrp
      work_e = 0.0d0
      jS = hecMESH%elem_group%grp_index(ig-1)
      jE = hecMESH%elem_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%elem_group%grp_item(j)
        if(idx >= 1 .and. idx <= hecMESH%n_elem) work_e(idx) = 1.0d0
      enddo
      label = 'EGRP_' // trim(hecMESH%elem_group%grp_name(ig))
      call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, label, work_e)
    enddo

    ! --- Surface groups ---
    do ig = 1, n_sgrp
      work_e = 0.0d0
      jS = hecMESH%surf_group%grp_index(ig-1)
      jE = hecMESH%surf_group%grp_index(ig)
      do j = jS+1, jE
        idx = hecMESH%surf_group%grp_item(2*j-1)
        if(idx >= 1 .and. idx <= hecMESH%n_elem) then
          work_e(idx) = dble(hecMESH%surf_group%grp_item(2*j))
        endif
      enddo
      label = 'SGRP_' // trim(hecMESH%surf_group%grp_name(ig))
      call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, label, work_e)
    enddo

    ! --- MATERIAL_ID ---
    do ie = 1, hecMESH%n_elem
      isect = hecMESH%section_ID(ie)
      if(isect >= 1 .and. isect <= hecMESH%section%n_sect) then
        cid = hecMESH%section%sect_mat_ID_item(isect)
        work_e(ie) = dble(cid)
      else
        work_e(ie) = 0.0d0
      endif
    enddo
    call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, 'MATERIAL_ID', work_e)

    ! --- SECTION_ID ---
    do ie = 1, hecMESH%n_elem
      work_e(ie) = dble(hecMESH%section_ID(ie))
    enddo
    call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, 'SECTION_ID', work_e)

    ! --- VOLUME ---
    call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, 'VOLUME', elem_vol)

    ! --- ASPECT_RATIO ---
    call hecmw_result_add(HECMW_RESULT_DTYPE_ELEM, 1, 'ASPECT_RATIO', elem_asp)

    ! --- WRITE ---
    nameID = 'fstrRES'
    call hecmw_result_write_by_name(nameID)
    call hecmw_result_finalize

    deallocate(work_n)
    deallocate(work_e)

  end subroutine inputcheck_write_result

end module m_inputcheck_make_result
