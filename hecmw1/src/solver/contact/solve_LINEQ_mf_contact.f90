!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface of DIRECTmf for
!! contact problems using Lagrange multiplier without elimination.
module m_solve_LINEQ_mf_contact
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_matrix_ass
  use hecmw_matrix_dump
  use hecmw_solver_direct_mf
  use hecmw_ebc_defer

  private
  public :: solve_LINEQ_mf_contact_init
  public :: solve_LINEQ_mf_contact

  logical, save :: INITIALIZED = .false.
  logical, save :: NEED_ANALYSIS = .true.
  logical, save :: IS_SYMMETRIC = .true.
  !> the saddle-point system as a 1-DOF-block matrix: the numeric stage of DIRECTmf loads
  !> values through uniform NDOF blocks only, and the scalar Lagrange rows cannot join the
  !> node blocks, so the whole system is expanded to DOF granularity instead
  type (hecmwST_matrix), save :: mfMAT

  !> refinement defaults of this path; the refinement runs by default because the first
  !> factorization after a contact activation can leave the residual near 1e-8
  integer(kind=kint), parameter :: MFC_IR_STEPS = 5
  real(kind=kreal), parameter :: MFC_IR_TOL = 1.0d-8

contains

  subroutine solve_LINEQ_mf_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    logical, intent(in) :: is_sym

    if (INITIALIZED) then
      call hecmw_mat_finalize(mfMAT)
      INITIALIZED = .false.
    endif
    call hecmw_mat_init(mfMAT)
    IS_SYMMETRIC = is_sym
    NEED_ANALYSIS = .true.
    INITIALIZED = .true.
  end subroutine solve_LINEQ_mf_contact_init

  subroutine solve_LINEQ_mf_contact(hecMESH,hecMAT,hecLagMAT,hecEBC,istat)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT
    type (hecmwST_ebc), intent(inout) :: hecEBC
    integer(kind=kint), intent(out) :: istat
    integer(kind=kint) :: ntdf, mpc_method, loglevel, irmax, it
    real(kind=kreal) :: bnrm, rnrm, irtol
    real(kind=kreal), allocatable :: bb(:), xx(:), rr(:)
    ! the unit hecmw_solve passes to the direct solvers for error messages
    integer(kind=kint), parameter :: imsg = 51

    mpc_method = hecmw_mat_get_mpc_method(hecMAT)
    if (mpc_method < 1 .or. 3 < mpc_method) then
      mpc_method = 1
      call hecmw_mat_set_mpc_method(hecMAT,mpc_method)
    endif
    if (mpc_method /= 1) then
      write(*,*) 'ERROR: MPCMETHOD other than penalty is not available for DIRECTmf solver', &
          ' in contact analysis without elimination'
      stop
    endif
    call hecmw_mat_ass_equation(hecMESH, hecMAT)
    call hecmw_mat_ass_equation_rhs(hecMESH, hecMAT)
    call hecmw_ebc_apply(hecMESH, hecMAT, hecEBC)

    call hecmw_mat_dump(hecMAT, hecMESH)

    if (NEED_ANALYSIS) then
      call mf_contact_set_profile(hecMAT, hecLagMAT)
      NEED_ANALYSIS = .false.
    endif

    call mf_contact_set_values(hecMAT, hecLagMAT)
    ntdf = hecMAT%NP*hecMAT%NDOF + hecLagMAT%num_lagrange
    mfMAT%B(1:ntdf) = hecMAT%B(1:ntdf)
    mfMAT%X(1:ntdf) = 0.0d0
    mfMAT%Iarray(97) = 1
    mfMAT%symmetric = IS_SYMMETRIC

    call hecmw_solve_direct_mf(hecMESH, mfMAT, imsg)
    istat = 0

    ! refine on the saddle-point system itself: the refinement of the driver measures its
    ! norms through the mesh, which does not cover the Lagrange rows
    irmax = hecMAT%Iarray(44)
    if (irmax == 0) irmax = MFC_IR_STEPS
    if (irmax > 0) then
      loglevel = hecmw_mat_get_loglevel(hecMAT)
      if (loglevel < 0) loglevel = max(hecmw_mat_get_timelog(hecMAT), hecmw_mat_get_iterlog(hecMAT))
      irtol = hecMAT%Rarray(42)
      if (.not. (irtol > 0.0d0)) irtol = MFC_IR_TOL
      allocate(bb(ntdf), xx(ntdf), rr(ntdf))
      bb(1:ntdf) = mfMAT%B(1:ntdf)
      xx(1:ntdf) = mfMAT%X(1:ntdf)
      bnrm = sqrt(dot_product(bb(1:ntdf), bb(1:ntdf)))
      do it = 1, irmax
        call mf_contact_resid(ntdf, xx, bb, rr)
        rnrm = sqrt(dot_product(rr(1:ntdf), rr(1:ntdf))) / max(bnrm, tiny(bnrm))
        if (loglevel > 0) write(*,'(a,i0,a,1pe11.4)') '[DIRECTmf]: refinement ', it - 1, ': residual = ', rnrm
        if (rnrm <= irtol) exit
        mfMAT%B(1:ntdf) = rr(1:ntdf)
        call hecmw_solve_direct_mf(hecMESH, mfMAT, imsg)
        xx(1:ntdf) = xx(1:ntdf) + mfMAT%X(1:ntdf)
      enddo
      mfMAT%X(1:ntdf) = xx(1:ntdf)
      deallocate(bb, xx, rr)
    endif

    hecMAT%X(1:ntdf) = mfMAT%X(1:ntdf)

    call hecmw_mat_dump_solution(hecMAT)
  end subroutine solve_LINEQ_mf_contact

  !> Build the DOF-granularity profile of the saddle-point system: the scalar expansion of
  !> hecMAT followed by the Lagrange rows. Both halves are always built (the factorization
  !> reads whichever half an entry lands in after permutation, and the structure must stay
  !> mirrored); the Lagrange couplings of hecLagMAT are mirrored by construction.
  subroutine mf_contact_set_profile(hecMAT, hecLagMAT)
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    integer(kind=kint) :: np, ndof, nlag, ntdf, npl, npu
    integer(kind=kint) :: i, j, k, l, cl, cu

    np = hecMAT%NP
    ndof = hecMAT%NDOF
    nlag = hecLagMAT%num_lagrange
    ntdf = np*ndof + nlag
    npl = hecMAT%NPL*ndof*ndof + np*ndof*(ndof-1)/2
    npu = hecMAT%NPU*ndof*ndof + np*ndof*(ndof-1)/2
    if (nlag > 0) then
      npl = npl + hecLagMAT%numL_lagrange*ndof
      npu = npu + hecLagMAT%numU_lagrange*ndof
    endif

    call hecmw_mat_finalize(mfMAT)
    mfMAT%N = ntdf
    mfMAT%NP = ntdf
    mfMAT%NDOF = 1
    mfMAT%NPL = npl
    mfMAT%NPU = npu
    allocate(mfMAT%indexL(0:ntdf), mfMAT%indexU(0:ntdf))
    allocate(mfMAT%itemL(npl), mfMAT%itemU(npu))
    allocate(mfMAT%D(ntdf), mfMAT%AL(npl), mfMAT%AU(npu))
    allocate(mfMAT%B(ntdf), mfMAT%X(ntdf))

    cl = 0
    cu = 0
    mfMAT%indexL(0) = 0
    mfMAT%indexU(0) = 0
    do i = 1, np
      do j = 1, ndof
        do l = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
          do k = 1, ndof
            cl = cl + 1
            mfMAT%itemL(cl) = (hecMAT%itemL(l)-1)*ndof + k
          enddo
        enddo
        do k = 1, j-1
          cl = cl + 1
          mfMAT%itemL(cl) = (i-1)*ndof + k
        enddo
        do k = j+1, ndof
          cu = cu + 1
          mfMAT%itemU(cu) = (i-1)*ndof + k
        enddo
        do l = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
          do k = 1, ndof
            cu = cu + 1
            mfMAT%itemU(cu) = (hecMAT%itemU(l)-1)*ndof + k
          enddo
        enddo
        if (nlag > 0) then
          do l = hecLagMAT%indexU_lagrange(i-1)+1, hecLagMAT%indexU_lagrange(i)
            cu = cu + 1
            mfMAT%itemU(cu) = np*ndof + hecLagMAT%itemU_lagrange(l)
          enddo
        endif
        mfMAT%indexL((i-1)*ndof+j) = cl
        mfMAT%indexU((i-1)*ndof+j) = cu
      enddo
    enddo
    do i = 1, nlag
      do l = hecLagMAT%indexL_lagrange(i-1)+1, hecLagMAT%indexL_lagrange(i)
        do k = 1, ndof
          cl = cl + 1
          mfMAT%itemL(cl) = (hecLagMAT%itemL_lagrange(l)-1)*ndof + k
        enddo
      enddo
      mfMAT%indexL(np*ndof+i) = cl
      mfMAT%indexU(np*ndof+i) = cu
    enddo

    mfMAT%Iarray = hecMAT%Iarray
    mfMAT%Rarray = hecMAT%Rarray
    if (mfMAT%Iarray(43) /= 0) then
      if (hecmw_comm_get_rank() == 0) write(*,*) &
        '[DIRECTmf]: BLR is not available in contact analysis without elimination; disabled'
    endif
    ! the refinement stays with the caller (mf_contact_resid); the driver would measure it
    ! through the mesh
    mfMAT%Iarray(43) = 0
    mfMAT%Iarray(44) = -1
    mfMAT%Iarray(48) = 0
    mfMAT%Iarray(98) = 1
  end subroutine mf_contact_set_profile

  !> Load the values in the order the profile was built. The Lagrange diagonal stays zero;
  !> the delayed pivoting of the factorization handles it.
  subroutine mf_contact_set_values(hecMAT, hecLagMAT)
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    integer(kind=kint) :: np, ndof, nlag
    integer(kind=kint) :: i, j, k, l, cl, cu

    np = hecMAT%NP
    ndof = hecMAT%NDOF
    nlag = hecLagMAT%num_lagrange

    cl = 0
    cu = 0
    do i = 1, np
      do j = 1, ndof
        do l = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
          do k = 1, ndof
            cl = cl + 1
            mfMAT%AL(cl) = hecMAT%AL(((l-1)*ndof+j-1)*ndof+k)
          enddo
        enddo
        do k = 1, j-1
          cl = cl + 1
          mfMAT%AL(cl) = hecMAT%D(((i-1)*ndof+j-1)*ndof+k)
        enddo
        mfMAT%D((i-1)*ndof+j) = hecMAT%D(((i-1)*ndof+j-1)*ndof+j)
        do k = j+1, ndof
          cu = cu + 1
          mfMAT%AU(cu) = hecMAT%D(((i-1)*ndof+j-1)*ndof+k)
        enddo
        do l = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
          do k = 1, ndof
            cu = cu + 1
            mfMAT%AU(cu) = hecMAT%AU(((l-1)*ndof+j-1)*ndof+k)
          enddo
        enddo
        if (nlag > 0) then
          do l = hecLagMAT%indexU_lagrange(i-1)+1, hecLagMAT%indexU_lagrange(i)
            cu = cu + 1
            mfMAT%AU(cu) = hecLagMAT%AU_lagrange((l-1)*ndof+j)
          enddo
        endif
      enddo
    enddo
    do i = 1, nlag
      mfMAT%D(np*ndof+i) = 0.0d0
      do l = hecLagMAT%indexL_lagrange(i-1)+1, hecLagMAT%indexL_lagrange(i)
        do k = 1, ndof
          cl = cl + 1
          mfMAT%AL(cl) = hecLagMAT%AL_lagrange((l-1)*ndof+k)
        enddo
      enddo
    enddo
  end subroutine mf_contact_set_values

  !> r = b - A x on the scalar system (both halves stored)
  subroutine mf_contact_resid(ntdf, x, b, r)
    implicit none
    integer(kind=kint), intent(in) :: ntdf
    real(kind=kreal), intent(in) :: x(:), b(:)
    real(kind=kreal), intent(out) :: r(:)
    integer(kind=kint) :: i, l
    real(kind=kreal) :: s

    do i = 1, ntdf
      s = mfMAT%D(i)*x(i)
      do l = mfMAT%indexL(i-1)+1, mfMAT%indexL(i)
        s = s + mfMAT%AL(l)*x(mfMAT%itemL(l))
      enddo
      do l = mfMAT%indexU(i-1)+1, mfMAT%indexU(i)
        s = s + mfMAT%AU(l)*x(mfMAT%itemU(l))
      enddo
      r(i) = b(i) - s
    enddo
  end subroutine mf_contact_resid

end module m_solve_LINEQ_mf_contact
