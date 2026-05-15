!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief API module providing accessor functions for type(hecmwST_matrix).
!>
!> All direct member access to type(hecmwST_matrix) from fistr1 should go
!> through these accessor procedures so that the underlying data layout
!> can be changed without touching client code.
!>
!> Naming convention:
!>   hecmw_mat_get_<FIELD>   : read access (returns value or pointer)
!>   hecmw_mat_set_<FIELD>   : write access (assigns scalar or whole array)
!>   hecmw_mat_alloc_<FIELD> : allocate a pointer-array field
!>   hecmw_mat_dealloc_<FIELD>: deallocate a pointer-array field
!>   hecmw_mat_is_assoc_<FIELD>: check pointer association
module hecmw_matrix_api
  use hecmw_util
  implicit none
  private

  ! --- scalar getters / setters
  public :: hecmw_mat_get_N,    hecmw_mat_set_N
  public :: hecmw_mat_get_NP,   hecmw_mat_set_NP
  public :: hecmw_mat_get_NPL,  hecmw_mat_set_NPL
  public :: hecmw_mat_get_NPU,  hecmw_mat_set_NPU
  public :: hecmw_mat_get_NPA,  hecmw_mat_set_NPA
  public :: hecmw_mat_get_NDOF, hecmw_mat_set_NDOF
  public :: hecmw_mat_get_symmetric, hecmw_mat_set_symmetric

  ! --- pointer-array getters (return associated pointer)
  public :: hecmw_mat_get_D,      hecmw_mat_get_B,      hecmw_mat_get_X
  public :: hecmw_mat_get_AL,     hecmw_mat_get_AU,     hecmw_mat_get_A
  public :: hecmw_mat_get_indexL, hecmw_mat_get_indexU, hecmw_mat_get_indexA
  public :: hecmw_mat_get_itemL,  hecmw_mat_get_itemU,  hecmw_mat_get_itemA

  ! --- pointer-array whole-array setters / fillers
  public :: hecmw_mat_fill_D, hecmw_mat_fill_B, hecmw_mat_fill_X
  public :: hecmw_mat_fill_AL, hecmw_mat_fill_AU, hecmw_mat_fill_A
  public :: hecmw_mat_copy_D, hecmw_mat_copy_B, hecmw_mat_copy_X
  public :: hecmw_mat_copy_AL, hecmw_mat_copy_AU, hecmw_mat_copy_A

  ! --- pointer-array lifecycle
  public :: hecmw_mat_is_assoc_AL, hecmw_mat_is_assoc_AU
  public :: hecmw_mat_is_assoc_D,  hecmw_mat_is_assoc_B,  hecmw_mat_is_assoc_X
  public :: hecmw_mat_is_assoc_A
  public :: hecmw_mat_is_assoc_indexL, hecmw_mat_is_assoc_indexU
  public :: hecmw_mat_is_assoc_itemL,  hecmw_mat_is_assoc_itemU
  public :: hecmw_mat_alloc_AL, hecmw_mat_alloc_AU
  public :: hecmw_mat_alloc_D,  hecmw_mat_alloc_B,  hecmw_mat_alloc_X
  public :: hecmw_mat_alloc_A
  public :: hecmw_mat_dealloc_AL, hecmw_mat_dealloc_AU
  public :: hecmw_mat_dealloc_D,  hecmw_mat_dealloc_B,  hecmw_mat_dealloc_X
  public :: hecmw_mat_dealloc_A

  ! --- per-element get/set for pointer arrays (avoid f()(idx) syntax)
  public :: hecmw_mat_get_D_i,  hecmw_mat_set_D_i
  public :: hecmw_mat_get_B_i,  hecmw_mat_set_B_i
  public :: hecmw_mat_get_X_i,  hecmw_mat_set_X_i
  public :: hecmw_mat_get_AL_i, hecmw_mat_set_AL_i
  public :: hecmw_mat_get_AU_i, hecmw_mat_set_AU_i
  public :: hecmw_mat_get_A_i,  hecmw_mat_set_A_i
  public :: hecmw_mat_get_indexL_i, hecmw_mat_set_indexL_i
  public :: hecmw_mat_get_indexU_i, hecmw_mat_set_indexU_i
  public :: hecmw_mat_get_indexA_i, hecmw_mat_set_indexA_i
  public :: hecmw_mat_get_itemL_i,  hecmw_mat_set_itemL_i
  public :: hecmw_mat_get_itemU_i,  hecmw_mat_set_itemU_i
  public :: hecmw_mat_get_itemA_i,  hecmw_mat_set_itemA_i

  ! --- Iarray / Rarray (size 100, fixed)
  public :: hecmw_mat_get_Iarray,     hecmw_mat_set_Iarray
  public :: hecmw_mat_get_Iarray_all, hecmw_mat_set_Iarray_all
  public :: hecmw_mat_get_Rarray,     hecmw_mat_set_Rarray
  public :: hecmw_mat_get_Rarray_all, hecmw_mat_set_Rarray_all

contains

  ! ===== scalars =====================================================
  pure function hecmw_mat_get_N(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%N
  end function

  subroutine hecmw_mat_set_N(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%N = v
  end subroutine

  pure function hecmw_mat_get_NP(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%NP
  end function

  subroutine hecmw_mat_set_NP(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%NP = v
  end subroutine

  pure function hecmw_mat_get_NPL(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%NPL
  end function

  subroutine hecmw_mat_set_NPL(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%NPL = v
  end subroutine

  pure function hecmw_mat_get_NPU(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%NPU
  end function

  subroutine hecmw_mat_set_NPU(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%NPU = v
  end subroutine

  pure function hecmw_mat_get_NPA(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%NPA
  end function

  subroutine hecmw_mat_set_NPA(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%NPA = v
  end subroutine

  pure function hecmw_mat_get_NDOF(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v
    v = hecMAT%NDOF
  end function

  subroutine hecmw_mat_set_NDOF(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v
    hecMAT%NDOF = v
  end subroutine

  pure function hecmw_mat_get_symmetric(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = hecMAT%symmetric
  end function

  subroutine hecmw_mat_set_symmetric(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    logical, intent(in) :: v
    hecMAT%symmetric = v
  end subroutine

  ! ===== pointer-array getters =======================================
  function hecmw_mat_get_D(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%D
  end function

  function hecmw_mat_get_B(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%B
  end function

  function hecmw_mat_get_X(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%X
  end function

  function hecmw_mat_get_AL(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%AL
  end function

  function hecmw_mat_get_AU(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%AU
  end function

  function hecmw_mat_get_A(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), pointer :: p(:)
    p => hecMAT%A
  end function

  function hecmw_mat_get_indexL(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%indexL
  end function

  function hecmw_mat_get_indexU(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%indexU
  end function

  function hecmw_mat_get_indexA(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%indexA
  end function

  function hecmw_mat_get_itemL(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%itemL
  end function

  function hecmw_mat_get_itemU(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%itemU
  end function

  function hecmw_mat_get_itemA(hecMAT) result(p)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), pointer :: p(:)
    p => hecMAT%itemA
  end function

  ! ===== whole-array fill (scalar broadcast) =========================
  subroutine hecmw_mat_fill_D(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%D = v
  end subroutine

  subroutine hecmw_mat_fill_B(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%B = v
  end subroutine

  subroutine hecmw_mat_fill_X(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%X = v
  end subroutine

  subroutine hecmw_mat_fill_AL(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%AL = v
  end subroutine

  subroutine hecmw_mat_fill_AU(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%AU = v
  end subroutine

  subroutine hecmw_mat_fill_A(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v
    hecMAT%A = v
  end subroutine

  ! ===== whole-array copy (array RHS) ================================
  subroutine hecmw_mat_copy_D(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%D = v
  end subroutine
  subroutine hecmw_mat_copy_B(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%B = v
  end subroutine
  subroutine hecmw_mat_copy_X(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%X = v
  end subroutine
  subroutine hecmw_mat_copy_AL(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%AL = v
  end subroutine
  subroutine hecmw_mat_copy_AU(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%AU = v
  end subroutine
  subroutine hecmw_mat_copy_A(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(:)
    hecMAT%A = v
  end subroutine

  ! ===== associated() wrappers =======================================
  function hecmw_mat_is_assoc_AL(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%AL)
  end function

  function hecmw_mat_is_assoc_AU(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%AU)
  end function

  function hecmw_mat_is_assoc_D(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%D)
  end function

  function hecmw_mat_is_assoc_B(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%B)
  end function

  function hecmw_mat_is_assoc_X(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%X)
  end function

  function hecmw_mat_is_assoc_A(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%A)
  end function

  function hecmw_mat_is_assoc_indexL(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%indexL)
  end function

  function hecmw_mat_is_assoc_indexU(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%indexU)
  end function

  function hecmw_mat_is_assoc_itemL(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%itemL)
  end function

  function hecmw_mat_is_assoc_itemU(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    logical :: v
    v = associated(hecMAT%itemU)
  end function

  ! ===== allocate / deallocate =======================================
  subroutine hecmw_mat_alloc_AL(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%AL(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_alloc_AU(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%AU(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_alloc_D(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%D(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_alloc_B(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%B(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_alloc_X(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%X(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_alloc_A(hecMAT, n, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(out) :: ierr
    allocate(hecMAT%A(n), stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_AL(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%AL, stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_AU(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%AU, stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_D(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%D, stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_B(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%B, stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_X(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%X, stat=ierr)
  end subroutine

  subroutine hecmw_mat_dealloc_A(hecMAT, ierr)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(out) :: ierr
    deallocate(hecMAT%A, stat=ierr)
  end subroutine

  ! ===== Iarray / Rarray =============================================
  pure function hecmw_mat_get_Iarray(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%Iarray(i)
  end function

  subroutine hecmw_mat_set_Iarray(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i, v
    hecMAT%Iarray(i) = v
  end subroutine

  pure function hecmw_mat_get_Iarray_all(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint) :: v(100)
    v = hecMAT%Iarray
  end function

  subroutine hecmw_mat_set_Iarray_all(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: v(100)
    hecMAT%Iarray = v
  end subroutine

  pure function hecmw_mat_get_Rarray(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%Rarray(i)
  end function

  subroutine hecmw_mat_set_Rarray(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%Rarray(i) = v
  end subroutine

  pure function hecmw_mat_get_Rarray_all(hecMAT) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal) :: v(100)
    v = hecMAT%Rarray
  end function

  subroutine hecmw_mat_set_Rarray_all(hecMAT, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    real(kind=kreal), intent(in) :: v(100)
    hecMAT%Rarray = v
  end subroutine

  ! ===== element-wise accessors for pointer arrays ====================
  ! Provided to avoid the non-portable `func()(idx)` syntax on the LHS
  ! / RHS in client code; gfortran does not accept chained subscripting
  ! of a function reference even when the function returns a pointer.

  function hecmw_mat_get_D_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%D(i)
  end function
  subroutine hecmw_mat_set_D_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%D(i) = v
  end subroutine

  function hecmw_mat_get_B_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%B(i)
  end function
  subroutine hecmw_mat_set_B_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%B(i) = v
  end subroutine

  function hecmw_mat_get_X_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%X(i)
  end function
  subroutine hecmw_mat_set_X_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%X(i) = v
  end subroutine

  function hecmw_mat_get_AL_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%AL(i)
  end function
  subroutine hecmw_mat_set_AL_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%AL(i) = v
  end subroutine

  function hecmw_mat_get_AU_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%AU(i)
  end function
  subroutine hecmw_mat_set_AU_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%AU(i) = v
  end subroutine

  function hecmw_mat_get_A_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal) :: v
    v = hecMAT%A(i)
  end function
  subroutine hecmw_mat_set_A_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    real(kind=kreal), intent(in) :: v
    hecMAT%A(i) = v
  end subroutine

  function hecmw_mat_get_indexL_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%indexL(i)
  end function
  subroutine hecmw_mat_set_indexL_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%indexL(i) = v
  end subroutine

  function hecmw_mat_get_indexU_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%indexU(i)
  end function
  subroutine hecmw_mat_set_indexU_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%indexU(i) = v
  end subroutine

  function hecmw_mat_get_indexA_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%indexA(i)
  end function
  subroutine hecmw_mat_set_indexA_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%indexA(i) = v
  end subroutine

  function hecmw_mat_get_itemL_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%itemL(i)
  end function
  subroutine hecmw_mat_set_itemL_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%itemL(i) = v
  end subroutine

  function hecmw_mat_get_itemU_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%itemU(i)
  end function
  subroutine hecmw_mat_set_itemU_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%itemU(i) = v
  end subroutine

  function hecmw_mat_get_itemA_i(hecMAT, i) result(v)
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint) :: v
    v = hecMAT%itemA(i)
  end function
  subroutine hecmw_mat_set_itemA_i(hecMAT, i, v)
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: i
    integer(kind=kint), intent(in) :: v
    hecMAT%itemA(i) = v
  end subroutine

end module hecmw_matrix_api
