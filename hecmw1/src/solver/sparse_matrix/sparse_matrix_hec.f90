!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides conversion routines between HEC data structure
!! and DOF based sparse matrix structure (CSR/COO)
module m_sparse_matrix_hec
  use hecmw_util
  use m_hecmw_comm_f
  use m_sparse_matrix
  implicit none

  private

  public :: sparse_matrix_hec_init_prof
  public :: sparse_matrix_hec_set_prof
  public :: sparse_matrix_hec_set_conv_ext
  public :: sparse_matrix_hec_set_vals
  public :: sparse_matrix_hec_set_rhs
  public :: sparse_matrix_hec_get_rhs

contains

  !!! public subroutines

  subroutine sparse_matrix_hec_init_prof(spMAT, hecMAT, hecMESH)
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint) :: ndof, ndof2, N_loc, NL, NU, NZ
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    N_loc=hecMAT%N*ndof
    if (sparse_matrix_is_sym(spMAT)) then
      NU=hecMAT%indexU(hecMAT%N)
      NZ=hecMAT%N*(ndof2+ndof)/2+NU*ndof2
    else
      NL=hecMAT%indexL(hecMAT%N)
      NU=hecMAT%indexU(hecMAT%N)
      NZ=(hecMAT%N+NU+NL)*ndof2
    endif
    call sparse_matrix_init(spMAT, N_loc, NZ)
    call sparse_matrix_hec_set_conv_ext(spMAT, hecMESH, hecMAT%NDOF)
    spMAT%iterlog = hecMAT%Iarray(21)
    spMAT%timelog = hecMAT%Iarray(22)
    call sparse_matrix_hec_set_prof(spMAT, hecMAT)
  end subroutine sparse_matrix_hec_init_prof

  subroutine sparse_matrix_hec_set_conv_ext(spMAT, hecMESH, ndof)
    type(sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    integer(kind=kint) :: i,j,nn_external,id,ierr
    integer(kint)   :: pid,lid,j0
    if (hecMESH%n_neighbor_pe==0) return
    ! create conversion list
    nn_external = hecMESH%n_node - hecMESH%nn_internal
    allocate(spMAT%conv_ext(nn_external*ndof), stat=ierr)
    if (ierr /= 0) then
      write(*,*) " Allocation error, spMAT%conv_ext"
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    spMAT%conv_ext(:) = -1
    do i=1,nn_external
      id = i + hecMESH%nn_internal
      pid = hecMESH%node_ID(id*2)
      lid = hecMESH%node_ID(id*2-1)
      j0 = spMAT%DISPLS(pid+1) + (lid-1)*ndof
      do j=1,ndof
        spMAT%conv_ext((i-1)*ndof+j) = j0+j
      enddo
    enddo
  end subroutine sparse_matrix_hec_set_conv_ext

  subroutine sparse_matrix_hec_set_prof(spMAT, hecMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: ndof, ndof2
    integer(kind=kint) :: m, i, idof, i0, ii, ls, le, l, j, j0, jdof, jdofs
    !integer(kind=kint) :: offset_l, offset_d, offset_u
    ! CONVERT TO CSR or COO STYLE
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    m=1
    ii = 0
    do i=1,hecMAT%N
      do idof=1,ndof
        i0=spMAT%offset+ndof*(i-1)
        ii=i0+idof
        if (spMAT%type==SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii-spMAT%offset)=m
        ! Lower
        if (.not. sparse_matrix_is_sym(spMAT)) then
          ls=hecMAT%indexL(i-1)+1
          le=hecMAT%indexL(i)
          do l=ls,le
            j=hecMAT%itemL(l)
            !if (j <= hecMAT%N) then
            j0=spMAT%offset+ndof*(j-1)
            !else
            !   j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
            !endif
            !offset_l=ndof2*(l-1)+ndof*(idof-1)
            do jdof=1,ndof
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
              spMAT%JCN(m)=j0+jdof
              !spMAT%A(m)=hecMAT%AL(offset_l+jdof)
              m=m+1
            enddo
          enddo
        endif
        ! Diag
        !offset_d=ndof2*(i-1)+ndof*(idof-1)
        if (sparse_matrix_is_sym(spMAT)) then; jdofs=idof; else; jdofs=1; endif
        do jdof=jdofs,ndof
          if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
          spMAT%JCN(m)=i0+jdof
          !spMAT%A(m)=hecMAT%D(offset_d+jdof)
          m=m+1
        enddo
        ! Upper
        ls=hecMAT%indexU(i-1)+1
        le=hecMAT%indexU(i)
        do l=ls,le
          j=hecMAT%itemU(l)
          if (j <= hecMAT%N) then
            j0=spMAT%offset+ndof*(j-1)
          else
            j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
            if (sparse_matrix_is_sym(spMAT) .and. j0 < i0) cycle
          endif
          !offset_u=ndof2*(l-1)+ndof*(idof-1)
          do jdof=1,ndof
            if (spMAT%type==SPARSE_MATRIX_TYPE_COO) spMAT%IRN(m)=ii
            spMAT%JCN(m)=j0+jdof
            !spMAT%A(m)=hecMAT%AU(offset_u+jdof)
            m=m+1
          enddo
        enddo
      enddo
    enddo
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) spMAT%IRN(ii+1-spMAT%offset)=m
    if (sparse_matrix_is_sym(spMAT) .and. m-1 < spMAT%NZ) spMAT%NZ=m-1
    if (m-1 /= spMAT%NZ) then
      write(*,*) 'ERROR: sparse_matrix_set_ij on rank ',hecmw_comm_get_rank()
      write(*,*) 'm-1 = ',m-1,', NZ=',spMAT%NZ
      stop
    endif
  end subroutine sparse_matrix_hec_set_prof

  subroutine sparse_matrix_hec_set_vals(spMAT, hecMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: ndof, ndof2
    integer(kind=kint) :: m, i, idof, i0, ii, ls, le, l, j, j0, jdof, jdofs
    integer(kind=kint) :: offset_l, offset_d, offset_u
    ndof=hecMAT%NDOF; ndof2=ndof*ndof
    ii = 0
    m=1
    do i=1,hecMAT%N
      do idof=1,ndof
        i0=spMAT%offset+ndof*(i-1)
        ii=i0+idof
        if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) then
          if (spMAT%IRN(ii-spMAT%offset)/=m) stop "ERROR: sparse_matrix_set_a1"
        endif
        ! Lower
        if (.not. sparse_matrix_is_sym(spMAT)) then
          ls=hecMAT%indexL(i-1)+1
          le=hecMAT%indexL(i)
          do l=ls,le
            j=hecMAT%itemL(l)
            !if (j <= hecMAT%N) then
            j0=spMAT%offset+ndof*(j-1)
            !else
            !   j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
            !endif
            offset_l=ndof2*(l-1)+ndof*(idof-1)
            do jdof=1,ndof
              if (spMAT%type==SPARSE_MATRIX_TYPE_COO) then
                if (spMAT%IRN(m)/=ii) stop "ERROR: sparse_matrix_set_a2"
              endif
              if (spMAT%JCN(m)/=j0+jdof) stop "ERROR: sparse_matrix_set_a3"
              spMAT%A(m)=hecMAT%AL(offset_l+jdof)
              m=m+1
            enddo
          enddo
        endif
        ! Diag
        offset_d=ndof2*(i-1)+ndof*(idof-1)
        if (sparse_matrix_is_sym(spMAT)) then; jdofs=idof; else; jdofs=1; endif
        do jdof=jdofs,ndof
          if (spMAT%type==SPARSE_MATRIX_TYPE_COO) then
            if (spMAT%IRN(m)/=ii) stop "ERROR: sparse_matrix_set_a4"
          endif
          if (spMAT%JCN(m)/=i0+jdof) stop "ERROR: sparse_matrix_set_a5"
          spMAT%A(m)=hecMAT%D(offset_d+jdof)
          m=m+1
        enddo
        ! Upper
        ls=hecMAT%indexU(i-1)+1
        le=hecMAT%indexU(i)
        do l=ls,le
          j=hecMAT%itemU(l)
          if (j <= hecMAT%N) then
            j0=spMAT%offset+ndof*(j-1)
          else
            j0=spMAT%conv_ext(ndof*(j-hecMAT%N))-ndof
            if (sparse_matrix_is_sym(spMAT) .and. j0 < i0) cycle
          endif
          offset_u=ndof2*(l-1)+ndof*(idof-1)
          do jdof=1,ndof
            if (spMAT%type==SPARSE_MATRIX_TYPE_COO) then
              if (spMAT%IRN(m)/=ii) stop "ERROR: sparse_matrix_set_a6"
            endif
            if (spMAT%JCN(m)/=j0+jdof) stop "ERROR: sparse_matrix_set_a7"
            spMAT%A(m)=hecMAT%AU(offset_u+jdof)
            m=m+1
          enddo
        enddo
      enddo
    enddo
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) then
      if (spMAT%IRN(ii+1-spMAT%offset)/=m) stop "ERROR: sparse_matrix_set_a8"
    endif
    if (m-1 /= spMAT%NZ) stop "ERROR: sparse_matrix_set_a9"
  end subroutine sparse_matrix_hec_set_vals

  subroutine sparse_matrix_hec_set_rhs(spMAT, hecMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: ierr,i
    allocate(spMAT%rhs(spMAT%N_loc), stat=ierr)
    if (ierr /= 0) then
      write(*,*) " Allocation error, spMAT%rhs"
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    do i=1,spMAT%N_loc
      spMAT%rhs(i)=hecMAT%b(i)
    enddo
  end subroutine sparse_matrix_hec_set_rhs

  subroutine sparse_matrix_hec_get_rhs(spMAT, hecMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint) :: i
    do i=1,spMAT%N_loc
      hecMAT%x(i)=spMAT%rhs(i)
    enddo
    deallocate(spMAT%rhs)
  end subroutine sparse_matrix_hec_get_rhs

end module m_sparse_matrix_hec
