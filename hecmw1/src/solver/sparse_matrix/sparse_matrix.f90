!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides DOF based sparse matrix data structure (CSR and COO)
module m_sparse_matrix
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private

  public :: SPARSE_MATRIX_TYPE_CSR
  public :: SPARSE_MATRIX_TYPE_COO

  public :: SPARSE_MATRIX_SYMTYPE_ASYM
  public :: SPARSE_MATRIX_SYMTYPE_SPD
  public :: SPARSE_MATRIX_SYMTYPE_SYM

  public :: sparse_matrix
  public :: sparse_matrix_set_type
  public :: sparse_matrix_init
  public :: sparse_matrix_finalize
  public :: sparse_matrix_dump
  public :: sparse_matrix_gather_rhs
  public :: sparse_matrix_scatter_rhs
  public :: sparse_matrix_is_sym

  integer(kind=kint), parameter :: SPARSE_MATRIX_TYPE_CSR=1
  integer(kind=kint), parameter :: SPARSE_MATRIX_TYPE_COO=2

  integer(kind=kint), parameter :: SPARSE_MATRIX_SYMTYPE_ASYM=0
  integer(kind=kint), parameter :: SPARSE_MATRIX_SYMTYPE_SPD=1
  integer(kind=kint), parameter :: SPARSE_MATRIX_SYMTYPE_SYM=2

  type sparse_matrix
    integer(kind=kint) :: type
    integer(kind=kint) :: symtype
    integer(kind=kint) :: N, N_loc
    integer(kind=kint) :: NZ
    integer(kind=kint), pointer :: IRN(:) => null()
    integer(kind=kint), pointer :: JCN(:) => null()
    real(kind=kreal), pointer :: A(:) => null()
    real(kind=kreal), pointer :: rhs(:)
    integer(kind=kint) :: offset
    integer(kind=kint), pointer :: N_COUNTS(:) => null()
    integer(kind=kint), pointer :: DISPLS(:) => null()
    integer(kind=kint), pointer :: conv_ext(:) => null()
    integer(kind=kint) :: iterlog
    integer(kind=kint) :: timelog
    logical :: is_initialized = .false.
  end type sparse_matrix

contains

  !!! public subroutines

  subroutine sparse_matrix_set_type(spMAT, type, symtype)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: type
    integer(kind=kint), intent(in) :: symtype
    if (type == SPARSE_MATRIX_TYPE_CSR .or. &
        type == SPARSE_MATRIX_TYPE_COO) then
      spMAT%type = type
    else
      write(*,*) 'ERROR: unknown sparse matrix type'
      stop
    endif
    if (symtype == SPARSE_MATRIX_SYMTYPE_ASYM .or. &
        symtype == SPARSE_MATRIX_SYMTYPE_SPD .or. &
        symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
      spMAT%symtype = symtype
    else
      write(*,*) 'ERROR: unknown symmetry type for sparse matrix'
      stop
    endif
  end subroutine sparse_matrix_set_type

  subroutine sparse_matrix_init(spMAT, N_loc, NZ)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: N_loc
    integer(kind=kint), intent(in) :: NZ
    if (spMAT%is_initialized) then
      !write(*,*) "WARNING: sparse_matrix_init_prof: already initialized; freeing"
      call sparse_matrix_finalize(spMAT)
    endif
    call sparse_matrix_set_dimension(spMAT, N_loc)
    call sparse_matrix_set_offsets(spMAT)
    call sparse_matrix_allocate_arrays(spMAT, NZ)
    spMAT%is_initialized = .true.
  end subroutine sparse_matrix_init

  subroutine sparse_matrix_finalize(spMAT)
    type (sparse_matrix), intent(inout) :: spMAT
    call sparse_matrix_free_arrays(spMAT)
    call sparse_matrix_clear(spMAT)
    spMAT%is_initialized = .false.
  end subroutine sparse_matrix_finalize

  subroutine sparse_matrix_dump(spMAT)
    type(sparse_matrix), intent(in) :: spMAT
    integer(kind=kint),save :: num_call
    character(len=128) :: fname
    integer(kind=kint) :: i,myrank
    myrank=hecmw_comm_get_rank()
    write(fname,'(A,I0,A,I0)') 'sparse-matrix-',num_call,'.coo.',myrank
    num_call = num_call + 1
    write(*,*) 'Dumping sparse matrix to ', fname
    open(91,file=fname)
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) then
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_ASYM) &
        write(91,*) '%SPARSE MATRIX CSR ASYM'
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) &
        write(91,*) '%SPARSE MATRIX CSR SPD'
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) &
        write(91,*) '%SPARSE MATRIX CSR SYM'
      write(91,*) spMAT%N_loc, spMAT%N_loc, spMAT%NZ
      do i=1,spMAT%N_loc+1
        write(91,*) spMAT%IRN(i)
      enddo
      do i=1,spMAT%NZ
        write(91,*) spMAT%JCN(i), spMAT%A(i)
      enddo
    elseif (spMAT%type == SPARSE_MATRIX_TYPE_COO) then
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_ASYM) &
        write(91,*) '%SPARSE MATRIX COO ASYM'
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) &
        write(91,*) '%SPARSE MATRIX COO SPD'
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) &
        write(91,*) '%SPARSE MATRIX COO SYM'
      write(91,*) spMAT%N_loc, spMAT%N_loc, spMAT%NZ
      do i=1,spMAT%NZ
        write(91,*) spMAT%IRN(i), spMAT%JCN(i), spMAT%A(i)
      enddo
    endif
    close(91)
  end subroutine sparse_matrix_dump

  subroutine sparse_matrix_gather_rhs(spMAT, rhs_all)
    type (sparse_matrix), intent(in) :: spMAT
    real(kind=kreal), intent(out) :: rhs_all(:)
    integer(kind=kint) :: ierr,i
    if (hecmw_comm_get_size() == 1) then
      do i=1,spMAT%N_loc
        rhs_all(i) = spMAT%rhs(i)
      enddo
    else
      call HECMW_GATHERV_REAL(spMAT%rhs, spMAT%N_loc, &
        rhs_all, spMAT%N_COUNTS, spMAT%DISPLS, &
        0, hecmw_comm_get_comm())
    endif
  end subroutine sparse_matrix_gather_rhs

  subroutine sparse_matrix_scatter_rhs(spMAT, rhs_all)
    type (sparse_matrix), intent(inout) :: spMAT
    real(kind=kreal), intent(in) :: rhs_all(:)
    integer(kind=kint) :: ierr,i
    if (hecmw_comm_get_size() == 1) then
      do i=1,spMAT%N_loc
        spMAT%rhs(i) = rhs_all(i)
      enddo
    else
      call HECMW_SCATTERV_REAL( &
        rhs_all, spMAT%N_COUNTS, spMAT%DISPLS, &
        spMAT%rhs, spMAT%N_loc, &
        0, hecmw_comm_get_comm())
    endif
  end subroutine sparse_matrix_scatter_rhs

  function sparse_matrix_is_sym(spMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    logical :: sparse_matrix_is_sym
    sparse_matrix_is_sym = (spMAT%symtype > 0)
  end function sparse_matrix_is_sym

  !!! private subroutines

  subroutine sparse_matrix_set_dimension(spMAT, N_loc)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: N_loc
    integer(kind=kint) :: ierr
    spMAT%N_loc = N_loc
    if (hecmw_comm_get_size() == 1) then
      spMAT%N = spMAT%N_loc
    else
      call HECMW_ALLREDUCE_INT_1(spMAT%N_loc, spMAT%N, HECMW_SUM, &
        hecmw_comm_get_comm())
    endif
  end subroutine sparse_matrix_set_dimension

  subroutine sparse_matrix_set_offsets(spMAT)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint) :: i,ierr,nprocs,myrank
    nprocs=hecmw_comm_get_size()
    myrank=hecmw_comm_get_rank()
    ! OFFSET
    !if (myrank == 0) then
    if (associated(spMAT%N_COUNTS)) deallocate(spMAT%N_COUNTS)
    if (associated(spMAT%DISPLS)) deallocate(spMAT%DISPLS)
    allocate(spMAT%N_COUNTS(nprocs), spMAT%DISPLS(nprocs), stat=ierr)
    if (ierr /= 0) then
      write(*,*) " Allocation error, spMAT%N_COUNTS, spMAT%DISPLS"
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    !endif
    if (nprocs > 1) then
      call HECMW_ALLGATHER_INT_1(spMAT%N_loc, &
        spMAT%N_COUNTS, hecmw_comm_get_comm())
    endif
    spMAT%DISPLS(1)=0
    do i=1,nprocs-1
      spMAT%DISPLS(i+1)=spMAT%DISPLS(i)+spMAT%N_COUNTS(i)
    enddo
    spMAT%offset = spMAT%DISPLS(myrank+1)
  end subroutine sparse_matrix_set_offsets

  subroutine sparse_matrix_allocate_arrays(spMAT, NZ)
    type(sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: NZ
    integer(kind=kint) :: N_loc
    integer(kind=kint) :: ierr
    if (associated(spMAT%IRN)) deallocate(spMAT%IRN)
    if (associated(spMAT%JCN)) deallocate(spMAT%JCN)
    if (associated(spMAT%A)) deallocate(spMAT%A)
    ierr = -1
    N_loc=spMAT%N_loc
    if (spMAT%type == SPARSE_MATRIX_TYPE_CSR) then
      allocate(spMAT%IRN(N_loc+1), spMAT%JCN(NZ), spMAT%A(NZ), stat=ierr)
    elseif (spMAT%type == SPARSE_MATRIX_TYPE_COO) then
      allocate(spMAT%IRN(NZ), spMAT%JCN(NZ), spMAT%A(NZ), stat=ierr)
    endif
    if (ierr /= 0) then
      write(*,*) " Allocation error, spMAT%IRN, spMAT%JCN, spMAT%A"
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    spMAT%NZ = NZ
  end subroutine sparse_matrix_allocate_arrays

  subroutine sparse_matrix_free_arrays(spMAT)
    type (sparse_matrix), intent(inout) :: spMAT
    if (associated(spMAT%IRN)) deallocate(spMAT%IRN)
    if (associated(spMAT%JCN)) deallocate(spMAT%JCN)
    if (associated(spMAT%A)) deallocate(spMAT%A)
    if (associated(spMAT%N_COUNTS)) deallocate(spMAT%N_COUNTS)
    if (associated(spMAT%DISPLS)) deallocate(spMAT%DISPLS)
    if (associated(spMAT%conv_ext)) deallocate(spMAT%conv_ext)
  end subroutine sparse_matrix_free_arrays

  subroutine sparse_matrix_clear(spMAT)
    type(sparse_matrix), intent(inout) :: spMAT
    spMAT%type = 0
    spMAT%symtype = 0
    spMAT%N = 0
    spMAT%N_loc = 0
    spMAT%NZ = 0
    spMAT%IRN => null()
    spMAT%JCN => null()
    spMAT%A => null()
    spMAT%offset = 0
    spMAT%N_COUNTS => null()
    spMAT%DISPLS => null()
    spMAT%conv_ext => null()
    spMAT%iterlog = 0
    spMAT%timelog = 0
  end subroutine sparse_matrix_clear

end module m_sparse_matrix
