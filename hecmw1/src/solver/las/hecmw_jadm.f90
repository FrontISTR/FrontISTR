!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Jagged Diagonal Matrix storage for vector processors.
!> Original code was provided by JAMSTEC.

module hecmw_JAD_TYPE
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_JAD_TYPE_33
  use hecmw_JAD_TYPE_44
  use hecmw_JAD_TYPE_nn
  implicit none

  private

  public :: hecmw_JAD_INIT
  public :: hecmw_JAD_FINALIZE
  public :: hecmw_JAD_IS_INITIALIZED
  public :: hecmw_JAD_MATVEC

  !C---------------------- AU&AL
  real(kind=kreal), allocatable      :: AJAD(:)
  integer(kind=kint), allocatable    :: JAJAD(:)
  integer(kind=kint), allocatable    :: JADORD(:)
  integer(kind=kint), allocatable    :: IAJAD(:)
  integer(kind=kint) :: MJAD
  real(kind=kreal), allocatable  :: WP(:,:)
  integer(kind=kint) :: INITIALIZED = 0

contains

  subroutine hecmw_JAD_INIT(hecMAT)
    type(hecmwST_matrix) :: hecMAT
    select case(hecMAT%NDOF)
      case(3)
        call hecmw_JAD_INIT_33(hecMAT)
      case(4)
        call hecmw_JAD_INIT_44(hecMAT)
      case default
        call hecmw_JAD_INIT_nn(hecMAT)
    end select
    INITIALIZED = 1
  end subroutine hecmw_JAD_INIT

  subroutine hecmw_JAD_FINALIZE(hecMAT)
    type(hecmwST_matrix) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_JAD_FINALIZE_33()
      case(4)
        call hecmw_JAD_FINALIZE_44()
      case default
        call hecmw_JAD_FINALIZE_nn()
    end select
    INITIALIZED = 0
  end subroutine hecmw_JAD_FINALIZE

  function hecmw_JAD_IS_INITIALIZED()
    integer(kind=kint) :: hecmw_JAD_IS_INITIALIZED
    hecmw_JAD_IS_INITIALIZED = INITIALIZED
  end function hecmw_JAD_IS_INITIALIZED

  subroutine hecmw_JAD_MATVEC(hecMESH, hecMAT, X, Y, COMMtime)

    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime
    real(kind=kreal) :: START_TIME, END_TIME
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: i,idof,jdof,NDOF,NDOF2
    select case(hecMAT%NDOF)
      case(3)
        call hecmw_JAD_MATVEC_33(hecMESH, hecMAT, X, Y, COMMtime)
      case(4)
        call hecmw_JAD_MATVEC_44(hecMESH, hecMAT, X, Y, COMMtime)
      case default
        call hecmw_JAD_MATVEC_nn(hecMESH, hecMAT, X, Y, COMMtime)
    end select
  end subroutine hecmw_JAD_MATVEC

end module hecmw_JAD_TYPE
