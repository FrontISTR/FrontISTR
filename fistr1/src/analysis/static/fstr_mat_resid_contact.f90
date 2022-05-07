!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_fstr_mat_resid_contact
  use m_fstr
  use fstr_matrix_con_contact
  implicit none

  private
  public :: fstr_get_residual_contact
  public :: fstr_get_rel_resid_L2_contact
  public :: fstr_get_resid_max_contact

contains

  subroutine fstr_get_residual_contact(hecMESH, hecMAT, hecLagMAT, r, rlag)
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type fstrS  end subroutine fstr_get_residual_contact
    real(kind=kreal), intent(out) :: r(:)
    real(kind=kreal), intent(out) :: rlag(:)
    real(kind=kreal)   :: Tcomm, sum
    integer(kind=kint) :: i, idof, ii, i0, ls, le, l, loc0, ll, j, loc, l0, k
    integer(kind=kint) :: ndof, npndof
    !! r    = b     - K    x  - Ulag lag
    !! rlag = blag  - Llag x
    ndof=hecMAT%NDOF
    npndof=hecMAT%NP*ndof
    do i=1,npndof
      r(i)=0.d0
    enddo
    ! r = b - K x
    Tcomm = 0.d0
    call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, Tcomm)
    if (hecLagMAT%num_lagrange > 0) then
      ! r = r - Ulag lag
      do i=1,hecMAT%N
        i0=(i-1)*ndof
        do idof=1,ndof
          ii=i0+idof
          ls=hecLagMAT%indexU_lagrange(i-1)+1
          le=hecLagMAT%indexU_lagrange(i)
          sum=0.d0
          do l=ls,le
            loc=(l-1)*ndof+idof
            ll=hecLagMAT%itemU_lagrange(l)
            !sum=sum+hecLagMAT%AU_lagrange(loc)*hecLagMAT%Lagrange(ll)
            sum=sum+hecLagMAT%AU_lagrange(loc)*hecMAT%X(npndof+ll)
          enddo
          r(ii)=r(ii)-sum
        enddo
      enddo
      ! rlag = blag - Llag x
      do i=1,hecLagMAT%num_lagrange
        ls=hecLagMAT%indexL_lagrange(i-1)+1
        le=hecLagMAT%indexL_lagrange(i)
        sum=0.d0
        do l=ls,le
          loc0=(l-1)*ndof
          j=hecLagMAT%itemL_lagrange(l)
          l0=(j-1)*ndof
          do k=1,ndof
            loc=loc0+k
            ll=l0+k
            sum=sum+hecLagMAT%AL_lagrange(loc)*hecMAT%X(ll)
          enddo
        enddo
        rlag(i)=hecMAT%B(npndof+i)-sum
      enddo
    end if
  end subroutine fstr_get_residual_contact

  !> \brief This function calculates relative L2 residual
  function fstr_get_rel_resid_L2_contact(hecMESH, hecMAT, hecLagMAT)
    implicit none
    real(kind=kreal)         :: fstr_get_rel_resid_L2_contact
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type hecmwST_matrix_lagrange
    real(kind=kreal), allocatable :: r(:)
    real(kind=kreal), allocatable :: rlag(:)
    real(kind=kreal)   :: bnorm2, rnorm2
    real(kind=kreal)   :: rlagnorm2
    real(kind=kreal)   :: Tcomm
    integer(kind=kint) :: i
    allocate(r(hecMAT%NDOF*hecMAT%NP),rlag(hecLagMAT%num_lagrange))
    ! |b|^2
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, Tcomm)
    call fstr_get_residual_contact(hecMESH, hecMAT, hecLagMAT, r, rlag)
    ! |r|^2
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    ! |rlag|^2
    rlagnorm2=0.d0
    do i=1,hecLagMAT%num_lagrange
      rlagnorm2=rlagnorm2+rlag(i)*rlag(i)
    enddo
    deallocate(r,rlag)
    call hecmw_allreduce_R1(hecMESH, rlagnorm2, HECMW_SUM)
    ! |r_total|^2 = |r|^2 + |rlag|^2
    if (hecLagMAT%num_lagrange > 0) rnorm2=rnorm2+rlagnorm2
    ! |r_total| / |b|
    fstr_get_rel_resid_L2_contact = sqrt(rnorm2 / bnorm2)
  end function fstr_get_rel_resid_L2_contact

  !> \brief This function calculates maximum residual
  function fstr_get_resid_max_contact(hecMESH, hecMAT, hecLagMAT)
    implicit none
    real(kind=kreal)         :: fstr_get_resid_max_contact
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type hecmwST_matrix_lagrange
    real(kind=kreal), allocatable :: r(:)
    real(kind=kreal), allocatable :: rlag(:)
    real(kind=kreal) :: rmax, rlagmax
    real(kind=kreal) :: Tcomm
    allocate(r(hecMAT%NDOF*hecMAT%NP),rlag(hecLagMAT%num_lagrange))
    call fstr_get_residual_contact(hecMESH, hecMAT, hecLagMAT, r, rlag)
    rmax = maxval(dabs(r))
    if (hecLagMAT%num_lagrange > 0) then
      rlagmax = maxval(dabs(rlag))
      if (rlagmax > rmax) rmax = rlagmax
    endif
    deallocate(r,rlag)
    call hecmw_allreduce_R1(hecMESH, rmax, HECMW_MAX)
    fstr_get_resid_max_contact = rmax
  end function fstr_get_resid_max_contact

end module m_fstr_mat_resid_contact
