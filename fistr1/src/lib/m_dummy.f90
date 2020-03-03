!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module defined dummy data and function
module m_dummy
  use hecmw
  use mMechGauss

  implicit none

  type tDummy
    !> DUMMY
    integer(kind=kint) :: DUMMY_egrp_tot
    integer(kind=kint), pointer :: DUMMY_egrp_GRPID  (:)  =>null()
    integer(kind=kint), pointer :: DUMMY_egrp_ID     (:)  =>null()
    integer(kind=kint), pointer :: DUMMY_egrp_amp    (:)  =>null()
    real(kind=kreal), pointer   :: DUMMY_egrp_eps    (:)  =>null()
  end type

  integer, parameter :: kDUM_UNDEFINED = -1
  integer, parameter :: kDUM_INACTIVE  = 0
  integer, parameter :: kDUM_ACTIVE    = 1

contains

  !< print dummy info
  subroutine print_dummy_info( dum )
    type(tDummy), intent(in) :: dum  !< dummy info

<<<<<<< HEAD
    integer(kind=kint) :: i
=======
    integer(kind=kint) :: i, j
>>>>>>> 6b55241c ([dummy element] Added dummy data and reading !DUMMY keyword)

    write(*,'(A,I0)') 'DUMMY_egrp_tot: ', dum%DUMMY_egrp_tot
    do i=1,dum%DUMMY_egrp_tot
      write(*,*) 'DUMMY : ',i
      write(*,*) 'GRPID    : ',dum%DUMMY_egrp_GRPID
      write(*,*) 'EGRPID   : ',dum%DUMMY_egrp_ID
      write(*,*) 'AMPLITUDE: ',dum%DUMMY_egrp_amp
      write(*,*) 'EPSILON  : ',dum%DUMMY_egrp_eps
    end do
  end subroutine

<<<<<<< HEAD
  !< stiff matrix of a dummy element
  subroutine STF_DUMMY( ndof, nn, ecoord, u, stiff, element )
    use mMechGauss
    integer(kind=kint), intent(in)  :: ndof                   !< degree of freedum
    integer(kind=kint), intent(in)  :: nn                     !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(in)    :: u(3,nn)                !< nodal displacemwent
    type(tElement), intent(inout)   :: element                !< status of element
    real(kind=kreal),   intent(out) :: stiff(:,:)             !< stiff matrix

    integer(kind=kint) :: i,j,k,m,n
    real(kind=kreal)   :: dnn, coeff, xmax(3), xmin(3), dL

    !get average element size
    xmax(1:3) = ecoord(1:3,1)
    xmin(1:3) = ecoord(1:3,1)
    do i=2,nn
      do k=1,ndof
        if( ecoord(k,i) > xmax(k) ) xmax(k) = ecoord(k,i)
        if( ecoord(k,i) < xmin(k) ) xmin(k) = ecoord(k,i)
      end do
    end do
    dL = 0.33333333333d0*((xmax(1)-xmin(1))+(xmax(2)-xmin(2))+(xmax(3)-xmin(3)))

    coeff = element%dummy_coeff/dL
    dnn = 1.d0/dble(nn)
    stiff(:,:) = 0.d0
    do i=1,nn
      m = ndof*(i-1)

      do j=1,nn
        n = ndof*(j-1)
        do k=1,ndof
          stiff(m+k,n+k) = stiff(m+k,n+k)-coeff*dnn
        end do
      end do

      do k=1,ndof
        stiff(m+k,m+k) = stiff(m+k,m+k)+coeff
      end do
    end do

  end subroutine

  !< nodal force vector of a dummy element
  subroutine UPDATE_DUMMY( ndof, nn, ecoord, u, du, qf, element )
    use mMechGauss
    integer(kind=kint), intent(in)  :: ndof                   !< degree of freedum
    integer(kind=kint), intent(in)  :: nn                     !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    real(kind=kreal), intent(in)    :: u(3,nn)                !< nodal displacemwent
    real(kind=kreal), intent(in)    :: du(3,nn)               !< nodal displacemwent increment
    type(tElement), intent(inout)   :: element                !< status of element
    real(kind=kreal),   intent(out) :: qf(:)                  !< Internal Force

    integer(kind=kint) :: i,k,m
    real(kind=kreal)   :: dnn, coeff, aveu(3), xmax(3), xmin(3), dL

    !get average element size
    xmax(1:3) = ecoord(1:3,1)
    xmin(1:3) = ecoord(1:3,1)
    do i=2,nn
      do k=1,ndof
        if( ecoord(k,i) > xmax(k) ) xmax(k) = ecoord(k,i)
        if( ecoord(k,i) < xmin(k) ) xmin(k) = ecoord(k,i)
      end do
    end do
    dL = 0.33333333333d0*((xmax(1)-xmin(1))+(xmax(2)-xmin(2))+(xmax(3)-xmin(3)))

    coeff = element%dummy_coeff/dL
    dnn = 1.d0/dble(nn)
    aveu(:) = 0.d0
    do i=1,nn
      do k=1,ndof
        aveu(k) = aveu(k) + u(k,i) + du(k,i)
      end do
    end do
    aveu(:) = aveu(:)*dnn

    qf(:) = 0.d0
    do i=1,nn
      m = ndof*(i-1)

      do k=1,ndof
        qf(m+k) = coeff*(u(k,i)+du(k,i)-aveu(k))
      end do
    end do

    do i=1,size(element%gausses)
      element%gausses(i)%strain = 0.d0
      element%gausses(i)%stress = 0.d0
      element%gausses(i)%strain_out = 0.d0
      element%gausses(i)%stress_out = 0.d0
      if( associated(element%gausses(i)%istatus) ) element%gausses(i)%istatus = 0
      if( associated(element%gausses(i)%fstatus) ) element%gausses(i)%fstatus = 0.d0
    end do

  end subroutine


=======
>>>>>>> 6b55241c ([dummy element] Added dummy data and reading !DUMMY keyword)
end module m_dummy
