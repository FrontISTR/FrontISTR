!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides wrapper for parallel sparse direct solver HeteroSolver
module m_hecmw_HeteroSolver_wrapper
#ifdef HECMW_WITH_HETEROSOLVER
  use hecmw_util
  use m_sparse_matrix
  use m_hecmw_comm_f
  use heterosolver

  private
  public :: hecmw_HeteroSolver_wrapper_init
  public :: hecmw_HeteroSolver_wrapper_factorize
  public :: hecmw_HeteroSolver_wrapper_preprocess
  public :: hecmw_HeteroSolver_wrapper_solve
  public :: hecmw_HeteroSolver_wrapper_finalize

  integer(kind=8), save :: hnd
  integer(kind=kint), pointer, save :: iaptr(:) => null()
  integer(kind=kint), pointer, save :: iaind(:) => null()
  real(kind=kreal), pointer, save :: aval(:) => null()

contains
  subroutine hecmw_HeteroSolver_wrapper_init(spMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=8) :: i,j
    integer(kind=kint) :: istat
    integer(kind=kint)  :: isymm, iformat, idist, comm, myrank, nprocs, mingind, maxgind

    iformat=HS_DCSR
    idist=HS_DIST_ABX
    comm = hecmw_comm_get_comm()
    myrank=hecmw_comm_get_rank()
    nprocs=hecmw_comm_get_size()

    if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) then
      isymm=HS_SYMMETRIC_PDS
    else if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
      isymm=HS_SYMMETRIC
    else
      isymm=HS_UNSYMMETRIC
    endif

    iaptr=>spMAT%IRN
    iaind=>spMAT%JCN
    aval=>spMAT%A

    mingind=spMAT%DISPLS(myrank+1)+1
    if(myrank .eq. nprocs-1 ) then
      maxgind=spMAT%N
    else
      maxgind=spMAT%DISPLS(myrank+2)
    endif

    call phs_init_handle(hnd, spMAT%N, spMAT%N, isymm, iformat, idist, mingind, maxgind, comm, istat);
    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver init failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    ! option settings
    ! - refinement option is not set here but it will be performed max 3 times by default
    ! - accuracy of residual norm will be specified by solve routine's argument
    call phs_set_option(hnd, HS_OUTPUT, HS_OUTPUT_B, istat) !x overwrites b after solution phase
    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver set option failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif


  end subroutine hecmw_HeteroSolver_wrapper_init
  subroutine hecmw_HeteroSolver_wrapper_preprocess(spMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint)  :: istat

    call phs_preprocess_rd(hnd, iaptr, iaind, aval, istat)

    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver preprocess failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif
  end subroutine hecmw_HeteroSolver_wrapper_preprocess
  subroutine hecmw_HeteroSolver_wrapper_factorize(spMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint)  :: istat

    call phs_factorize_rd(hnd, iaptr, iaind, aval, istat)

    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver factorize failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif
  end subroutine hecmw_HeteroSolver_wrapper_factorize
  subroutine hecmw_HeteroSolver_wrapper_solve(spMAT)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint)  :: istat
    real(kind=kreal) :: res=1.0e-8
    real(kind=kreal), pointer :: b(:),x(:)
    b=>spMAT%rhs
    x=>spMAT%rhs ! never used

    call phs_solve_rd(hnd, iaptr, iaind, aval, 1, b, x, res, istat)
    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver solve failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif

  end subroutine hecmw_HeteroSolver_wrapper_solve
  subroutine hecmw_HeteroSolver_wrapper_finalize()
    implicit none
    integer(kind=kint)  :: istat

    call phs_finalize_handle(hnd, istat)
    if (istat /= HS_RESULT_OK) then
      write(*,*) 'HeteroSolver finalize failed: ', istat
      call hecmw_abort(hecmw_comm_get_comm())
    endif

  end subroutine hecmw_HeteroSolver_wrapper_finalize
#endif
end module m_hecmw_HeteroSolver_wrapper
