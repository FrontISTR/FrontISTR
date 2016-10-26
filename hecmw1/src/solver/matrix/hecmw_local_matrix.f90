!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2014/07/09                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  FrontISTR Forum,The University of Tokyo       !
!                                                                      !
!======================================================================!

module hecmw_local_matrix
  use hecmw_util

  private
  public :: hecmwST_local_matrix
  public :: hecmw_localmat_write
  public :: hecmw_localmat_blocking
  public :: hecmw_localmat_free
  public :: hecmw_localmat_mulvec
  public :: hecmw_trimatmul_TtKT
  public :: hecmw_trimatmul_TtKT_mpc
  public :: hecmw_localmat_transpose

  type hecmwST_local_matrix
    integer :: nr, nc, nnz, ndof
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: item(:)
    real(kind=kreal), pointer :: A(:)
  end type hecmwST_local_matrix

contains

  subroutine hecmw_localmat_write(Tmat,iunit)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer(kind=kint), intent(in) :: iunit
    integer(kind=kint) :: nr, ndof, ndof2, i, js, je, j, jj

    nr=Tmat%nr
    ndof=Tmat%ndof
    ndof2=ndof*ndof

    write(iunit,*) 'i, j, A'
    do i=1,nr
      js=Tmat%index(i-1)+1
      je=Tmat%index(i)
      do j=js,je
        jj=Tmat%item(j)
        if (ndof==1) then
          write(iunit,*) i, jj, Tmat%A(j)
        else
          write(iunit,*) i, jj
          write(iunit,*) Tmat%A((j-1)*ndof2+1:j*ndof2)
        endif
      enddo
    enddo
  end subroutine hecmw_localmat_write

  subroutine hecmw_localmat_blocking(Tmat, ndof, BTmat)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer, intent(in) :: ndof
    type (hecmwST_local_matrix), intent(out) :: BTmat
    integer, allocatable :: iw(:)
    integer :: ndof2, i, icnt, idof, idx, ls, le, l, j, jb, k, lb0, jdof, ks, ke
    ndof2=ndof*ndof

    if (mod(Tmat%nr, ndof) /= 0 .or. mod(Tmat%nc, ndof) /= 0) then
      write(0,*) Tmat%nr, Tmat%nc, ndof
      stop 'ERROR: blocking_Tmat failed'
    endif
    BTmat%nr=Tmat%nr/ndof
    BTmat%nc=Tmat%nc/ndof
    BTmat%ndof=ndof

    allocate(iw(BTmat%nc))
    allocate(BTmat%index(0:BTmat%nr))

    BTmat%index(0)=0
    do i=1,BTmat%nr
      icnt=0
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        ls=Tmat%index(idx-1)+1
        le=Tmat%index(idx)
        lcol: do l=ls,le
          j=Tmat%item(l)
          jb=(j-1)/ndof+1
          do k=1,icnt
            if (iw(k)==jb) cycle lcol
          enddo
          icnt=icnt+1
          iw(icnt)=jb
        enddo lcol
      enddo
      BTmat%index(i)=BTmat%index(i-1)+icnt
    enddo

    BTmat%nnz=BTmat%index(BTmat%nr)
    allocate(BTmat%item(BTmat%nnz))
    allocate(BTmat%A(BTmat%nnz*ndof2))
    BTmat%A=0.d0

    do i=1,BTmat%nr
      icnt=0
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        ls=Tmat%index(idx-1)+1
        le=Tmat%index(idx)
        lcol2: do l=ls,le
          j=Tmat%item(l)
          jb=(j-1)/ndof+1
          do k=1,icnt
            if (iw(k)==jb) cycle lcol2
          enddo
          icnt=icnt+1
          iw(icnt)=jb
        enddo lcol2
      enddo
      ! if (icnt /= BTmat%index(i)-BTmat%index(i-1)) stop 'ERROR: blocking Tmat'
      ! ! call qsort(iw, 1, icnt)
      lb0=BTmat%index(i-1)
      do k=1,icnt
        BTmat%item(lb0+k)=iw(k)
      enddo
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        ls=Tmat%index(idx-1)+1
        le=Tmat%index(idx)
        lcol3: do l=ls,le
          j=Tmat%item(l)
          jb=(j-1)/ndof+1
          jdof=mod((j-1), ndof)+1
          ks=BTmat%index(i-1)+1
          ke=BTmat%index(i)
          do k=ks,ke
            if (BTmat%item(k)==jb) then
              BTmat%A((k-1)*ndof2+(idof-1)*ndof+jdof)=Tmat%A(l)
              cycle lcol3
            endif
          enddo
          stop 'ERROR: something wrong in blocking Tmat'
        enddo lcol3
      enddo
    enddo
  end subroutine hecmw_localmat_blocking

  subroutine hecmw_localmat_free(Tmat)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: Tmat
    deallocate(Tmat%index)
    deallocate(Tmat%item)
    deallocate(Tmat%A)
    Tmat%nr=0
    Tmat%nc=0
    Tmat%nnz=0
    Tmat%ndof=0
  end subroutine hecmw_localmat_free

  subroutine hecmw_trimatmul_TtKT(BTtmat, hecMAT, BTmat, &
       iwS, num_lagrange)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_local_matrix) :: BTtKT

    ! perform three matrices multiplication for elimination
    call trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: BTtKT(MPC)'
    !call hecmw_localmat_write(BTtKT, 700+hecmw_comm_get_rank())

    ! place 1s where the DOF is eliminated
    call place_one_on_diag(BTtKT, iwS, num_lagrange)
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: BTtKT(MPC)'
    !call hecmw_localmat_write(BTtKT, 700+hecmw_comm_get_rank())

    ! set new values to HECMW matrix
    call set_to_hecmat(hecMAT, BTtKT)
    call hecmw_localmat_free(BTtKT)
  end subroutine hecmw_trimatmul_TtKT

  subroutine trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_matrix), intent(out) :: BTtKT
    integer :: nr, nc, ndof, ndof2, i, icnt, js, je, j, jj, ks, ke, k, kk
    integer :: ls, le, l, ll, m, ms, me, mm
    integer, allocatable :: iw(:)
    real(kind=kreal), pointer :: Ttp(:), Kp(:), Tp(:), TtKTp(:)
    ! real(kind=kreal) :: tsym_s, tsym_e, tnum_s, tnum_e

    nr=BTtmat%nr
    nc=BTmat%nc
    ndof=BTtmat%ndof
    ndof2=ndof*ndof

    BTtKT%nr=nr
    BTtKT%nc=nc
    BTtKT%ndof=ndof
    allocate(BTtKT%index(0:nr))

    ! tsym_s = hecmw_wtime()

!$omp parallel default(none), &
!$omp&         private(iw,i,icnt,js,je,j,jj,ks,ke,k,kk,ls,le,l,ll,m), &
!$omp&         shared(nr,nc,BTtmat,hecMAT,BTmat,BTtKT)
    allocate(iw(nc))
!$omp do
    do i=1,nr
      icnt=0
      js=BTtmat%index(i-1)+1
      je=BTtmat%index(i)
      do j=js,je
        jj=BTtmat%item(j)
        ! lower
        ks=hecMAT%indexL(jj-1)+1
        ke=hecMAT%indexL(jj)
        do k=ks,ke
          kk=hecMAT%itemL(k)
          ls=BTmat%index(kk-1)+1
          le=BTmat%index(kk)
          ll1: do l=ls,le
            ll=BTmat%item(l)
            do m=1,icnt
              if (iw(m)==ll) cycle ll1
            enddo
            icnt=icnt+1
            iw(icnt)=ll
            !if (i==1) write(0,*) 'l', icnt, jj, kk, ll
          enddo ll1
        enddo
        ! diagonal
        ls=BTmat%index(jj-1)+1
        le=BTmat%index(jj)
        ll2: do l=ls,le
          ll=BTmat%item(l)
          do m=1,icnt
            if (iw(m)==ll) cycle ll2
          enddo
          icnt=icnt+1
          iw(icnt)=ll
          !if (i==1) write(0,*) 'd', icnt, jj, kk, ll
        enddo ll2
        ! upper
        ks=hecMAT%indexU(jj-1)+1
        ke=hecMAT%indexU(jj)
        do k=ks,ke
          kk=hecMAT%itemU(k)
          ls=BTmat%index(kk-1)+1
          le=BTmat%index(kk)
          ll3: do l=ls,le
            ll=BTmat%item(l)
            do m=1,icnt
              if (iw(m)==ll) cycle ll3
            enddo
            icnt=icnt+1
            iw(icnt)=ll
            !if (i==1) write(0,*) 'u', icnt, jj, kk, ll
          enddo ll3
        enddo
      enddo
      if (icnt == 0) icnt=1
      !if (i==1) write(0,*) iw(1:icnt)
      BTtKT%index(i)=icnt
    enddo
!$omp end do
    deallocate(iw)
!$omp end parallel

    ! tsym_e = hecmw_wtime()
    ! write(0,*) 'tsym:',tsym_e-tsym_s

    BTtKT%index(0)=0
    do i=1,nr
      BTtKT%index(i)=BTtKT%index(i-1)+BTtKT%index(i)
    enddo
    !write(0,*) BTtKT%index(1:n)-BTtKT%index(0:n-1)

    BTtKT%nnz=BTtKT%index(nr)
    allocate(BTtKT%item(BTtKT%nnz))
    allocate(BTtKT%A(BTtKT%nnz*ndof2))
    BTtKT%item=0
    BTtKT%A=0.d0

    ! tnum_s = hecmw_wtime()

!$omp parallel default(none), &
!$omp&         private(i,icnt,js,je,j,jj,ks,ke,k,kk,ls,le,l,ll,m, &
!$omp&                 ms,me,mm,Ttp,Kp,Tp,TtKTp), &
!$omp&         shared(nr,nc,BTtmat,hecMAT,BTmat,BTtKT,ndof,ndof2)
!$omp do
    do i=1,nr
      icnt=0
      ms=BTtKT%index(i-1)+1
      !me=BTtKT%index(i)
      js=BTtmat%index(i-1)+1
      je=BTtmat%index(i)
      do j=js,je
        jj=BTtmat%item(j)
        Ttp=>BTtmat%A((j-1)*ndof2+1:j*ndof2)
        ! lower
        ks=hecMAT%indexL(jj-1)+1
        ke=hecMAT%indexL(jj)
        do k=ks,ke
          kk=hecMAT%itemL(k)
          Kp=>hecMAT%AL((k-1)*ndof2+1:k*ndof2)
          ls=BTmat%index(kk-1)+1
          le=BTmat%index(kk)
          do l=ls,le
            ll=BTmat%item(l)
            Tp=>BTmat%A((l-1)*ndof2+1:l*ndof2)
            me=ms-1+icnt
            mm=-1
            do m=ms,me
              if (BTtKT%item(m)==ll) mm=m
            enddo
            if (mm<0) then
              icnt=icnt+1
              mm=me+1
              BTtKT%item(mm)=ll
              !if (i==1) write(0,*) 'l', mm, jj, kk, ll
            endif
            TtKTp=>BTtKT%A((mm-1)*ndof2+1:mm*ndof2)
            call blk_trimatmul_add(ndof, Ttp, Kp, Tp, TtKTp)
          enddo
        enddo
        ! diagonal
        Kp=>hecMAT%D((jj-1)*ndof2+1:jj*ndof2)
        ls=BTmat%index(jj-1)+1
        le=BTmat%index(jj)
        do l=ls,le
          ll=BTmat%item(l)
          Tp=>BTmat%A((l-1)*ndof2+1:l*ndof2)
          me=ms-1+icnt
          mm=-1
          do m=ms,me
            if (BTtKT%item(m)==ll) mm=m
          enddo
          if (mm<0) then
            icnt=icnt+1
            mm=me+1
            BTtKT%item(mm)=ll
            !if (i==1) write(0,*) 'd', mm, jj, kk, ll
          endif
          TtKTp=>BTtKT%A((mm-1)*ndof2+1:mm*ndof2)
          call blk_trimatmul_add(ndof, Ttp, Kp, Tp, TtKTp)
        enddo
        ! upper
        ks=hecMAT%indexU(jj-1)+1
        ke=hecMAT%indexU(jj)
        do k=ks,ke
          kk=hecMAT%itemU(k)
          Kp=>hecMAT%AU((k-1)*ndof2+1:k*ndof2)
          ls=BTmat%index(kk-1)+1
          le=BTmat%index(kk)
          do l=ls,le
            ll=BTmat%item(l)
            Tp=>BTmat%A((l-1)*ndof2+1:l*ndof2)
            me=ms-1+icnt
            mm=-1
            do m=ms,me
              if (BTtKT%item(m)==ll) mm=m
            enddo
            if (mm<0) then
              icnt=icnt+1
              mm=me+1
              BTtKT%item(mm)=ll
              !if (i==1) write(0,*) 'u', mm, jj, kk, ll
            endif
            TtKTp=>BTtKT%A((mm-1)*ndof2+1:mm*ndof2)
            call blk_trimatmul_add(ndof, Ttp, Kp, Tp, TtKTp)
          enddo
        enddo
      enddo
      if (icnt == 0) then
        icnt=1
        BTtKT%item(ms)=i
      endif
      ! error check!
      !write(0,*) BTtKT%item(ms:ms-1+icnt)
      !write(0,*) BTtKT%index(i)-BTtKT%index(i-1), icnt
      if (ms-1+icnt /= BTtKT%index(i)) stop 'ERROR: trimatmul'
    enddo
!$omp end do
!$omp end parallel

    ! tnum_e = hecmw_wtime()
    ! write(0,*) 'tnum:',tnum_e-tnum_s
  end subroutine trimatmul_TtKT

  subroutine blk_trimatmul_add(ndof, A, B, C, ABC)
    implicit none
    integer, intent(in) :: ndof
    real(kind=kreal), intent(in) :: A(:), B(:), C(:)
    real(kind=kreal), intent(inout) :: ABC(:)
    real(kind=kreal), allocatable :: AB(:)
    integer :: ndof2, i, j, k, i0, j0, ij, ik, jk

    ndof2=ndof*ndof
    allocate(AB(ndof2))
    AB=0.d0

    do i=1,ndof
      i0=(i-1)*ndof
      do j=1,ndof
        ij=i0+j
        j0=(j-1)*ndof
        do k=1,ndof
          ik=i0+k
          jk=j0+k
          AB(ik)=AB(ik)+A(ij)*B(jk)
        enddo
      enddo
    enddo

    do i=1,ndof
      i0=(i-1)*ndof
      do j=1,ndof
        ij=i0+j
        j0=(j-1)*ndof
        do k=1,ndof
          ik=i0+k
          jk=j0+k
          ABC(ik)=ABC(ik)+AB(ij)*C(jk)
        enddo
      enddo
    enddo

    deallocate(AB)
  end subroutine blk_trimatmul_add

  subroutine place_one_on_diag(BTtKT, iwS, num_lagrange)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTtKT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint) :: ndof, ndof2, ilag, i, idof, js, je, j, jj

    ndof=BTtKT%ndof
    ndof2=ndof*ndof

    outer: do ilag=1,num_lagrange
      i=(iwS(ilag)-1)/ndof+1
      idof=mod(iwS(ilag)-1, ndof)+1
      js=BTtKT%index(i-1)+1
      je=BTtKT%index(i)
      do j=js,je
        jj=BTtKT%item(j)
        if (jj==i) then
          !write(0,*) ilag, i, idof
          BTtKT%A((j-1)*ndof2+(idof-1)*ndof+idof)=1.d0
          cycle outer
        endif
      enddo
    enddo outer
  end subroutine place_one_on_diag

  subroutine replace_hecmat(hecMAT, BTtKT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_matrix), intent(in) :: BTtKT
    integer :: nr, nc, ndof, ndof2, i, nl, nu, js, je, j, jj
    integer :: ksl, ksu, k, idof, idx

    nr=BTtKT%nr
    nc=BTtKT%nc
    ndof=hecMAT%NDOF
    ndof2=ndof*ndof

    ! free old hecMAT
    if (associated(hecMAT%AL)) deallocate(hecMAT%AL)
    if (associated(hecMAT%AU)) deallocate(hecMAT%AU)
    if (associated(hecMAT%itemL)) deallocate(hecMAT%itemL)
    if (associated(hecMAT%itemU)) deallocate(hecMAT%itemU)
    hecMAT%indexL=0
    hecMAT%indexU=0

    ! count NPL, NPU
!$omp parallel default(none),private(i,nl,nu,js,je,j,jj), &
!$omp&         shared(nr,BTtKT,hecMAT)
!$omp do
    do i=1,nr
      nl=0
      nu=0
      js=BTtKT%index(i-1)+1
      je=BTtKT%index(i)
      do j=js,je
        jj=BTtKT%item(j)
        if (jj < i) then
          nl=nl+1
        elseif (i < jj) then
          nu=nu+1
        else
          ! diagonal
        endif
      enddo
      hecMAT%indexL(i)=nl
      hecMAT%indexU(i)=nu
    enddo
!$omp end do
!$omp end parallel

    hecMAT%indexL(0)=0
    hecMAT%indexU(0)=0
    do i=1,nc
      hecMAT%indexL(i)=hecMAT%indexL(i-1)+hecMAT%indexL(i)
      hecMAT%indexU(i)=hecMAT%indexU(i-1)+hecMAT%indexU(i)
    enddo
    hecMAT%NPL=hecMAT%indexL(nc)
    hecMAT%NPU=hecMAT%indexU(nc)

    ! allocate new hecMAT
    allocate(hecMAT%itemL(hecMAT%NPL), hecMAT%itemU(hecMAT%NPU))
    allocate(hecMAT%AL(hecMAT%NPL*ndof2), hecMAT%AU(hecMAT%NPU*ndof2))
    hecMAT%itemL=0
    hecMAT%itemU=0
    hecMAT%D=0.d0
    hecMAT%AL=0.d0
    hecMAT%AU=0.d0

    ! copy from BTtKT to hecMAT
!$omp parallel default(none),private(i,nl,nu,js,je,ksl,ksu,j,jj,k), &
!$omp&  shared(nr,BTtKT,hecMAT,ndof2)
!$omp do
    do i=1,nr
      nl=0
      nu=0
      js=BTtKT%index(i-1)+1
      je=BTtKT%index(i)
      ksl=hecMAT%indexL(i-1)+1
      ksu=hecMAT%indexU(i-1)+1
      do j=js,je
        jj=BTtKT%item(j)
        if (jj < i) then
          k=ksl+nl
          hecMAT%itemL(k)=jj
          hecMAT%AL((k-1)*ndof2+1:k*ndof2)=BTtKT%A((j-1)*ndof2+1:j*ndof2)
          nl=nl+1
        elseif (i < jj) then
          k=ksu+nu
          hecMAT%itemU(k)=jj
          hecMAT%AU((k-1)*ndof2+1:k*ndof2)=BTtKT%A((j-1)*ndof2+1:j*ndof2)
          nu=nu+1
        else
          hecMAT%D((i-1)*ndof2+1:i*ndof2)=BTtKT%A((j-1)*ndof2+1:j*ndof2)
        endif
      enddo
      ! if (ksl+nl /= hecMAT%indexL(i)+1) stop 'ERROR: indexL'
      ! if (ksu+nu /= hecMAT%indexU(i)+1) stop 'ERROR: indexU'
    enddo
!$omp end do
!$omp end parallel

    ! do i=1,hecMAT%NPL
    !   if (hecMAT%itemL(i) <= 0) stop 'ERROR: negative itemL'
    !   if (hecMAT%itemL(i) > nc) stop 'ERROR: too big itemL'
    ! enddo
    ! do i=1,hecMAT%NPU
    !   if (hecMAT%itemU(i) <= 0) stop 'ERROR: negative itemU'
    !   if (hecMAT%itemU(i) > nc) stop 'ERROR: too big itemU'
    ! enddo
  end subroutine replace_hecmat

  subroutine make_new_hecmat(hecMAT, BTtKT, hecTKT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_local_matrix), intent(in) :: BTtKT
    type(hecmwST_matrix), intent(inout) :: hecTKT
    integer(kind=kint) :: nr, nc, ndof, ndof2, i

    nr=BTtKT%nr
    nc=BTtKT%nc
    ndof=BTtKT%ndof
    ndof2=ndof*ndof

    !write(0,*) 'DEBUG: nr, nc =',nr,nc

    ! if (nr /= nc) then
    !   stop 'ERROR: nr /= nc'
    ! endif
    hecTKT%N =hecMAT%N
    hecTKT%NP=nc
    hecTKT%NDOF=ndof

    allocate(hecTKT%D(nc*ndof2))
    allocate(hecTKT%B(nc*ndof))
    allocate(hecTKT%X(nc*ndof))

    allocate(hecTKT%indexL(0:nc))
    allocate(hecTKT%indexU(0:nc))

    hecTKT%Iarray=hecMAT%Iarray
    hecTKT%Rarray=hecMAT%Rarray

    call replace_hecmat(hecTKT, BTtKT)

    ! copy internal part only
    do i=1,hecMAT%N*ndof
      hecTKT%B(i)=hecMAT%B(i)
      hecTKT%X(i)=hecMAT%X(i)
    enddo
  end subroutine make_new_hecmat

  subroutine set_to_hecmat(hecMAT, BTtKT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    type(hecmwST_local_matrix), intent(in) :: BTtKT
    integer(kind=kint) :: nr, nc, ndof, ndof2, i

    if (associated(hecMAT%indexL_org2)) then
      stop 'ERROR: too many matrix backups'
    endif
    if (associated(hecMAT%indexL_org1)) then
      hecMAT%indexL_org2 => hecMAT%indexL_org1
      hecMAT%indexU_org2 => hecMAT%indexU_org1
      hecMAT%itemL_org2 => hecMAT%itemL_org1
      hecMAT%itemU_org2 => hecMAT%itemU_org1
    endif
    hecMAT%indexL_org1 => hecMAT%indexL
    hecMAT%indexU_org1 => hecMAT%indexU
    hecMAT%itemL_org1 => hecMAT%itemL
    hecMAT%itemU_org1 => hecMAT%itemU

    allocate(hecMAT%indexL(0:hecMAT%NP))
    allocate(hecMAT%indexU(0:hecMAT%NP))
    hecMAT%itemL => null()
    hecMAT%itemU => null()

    call replace_hecmat(hecMAT, BTtKT)
  end subroutine set_to_hecmat

  subroutine hecmw_localmat_mulvec(BTmat, V, TV)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    real(kind=kreal), intent(in), target :: V(:)
    real(kind=kreal), intent(out), target :: TV(:)
    real(kind=kreal), pointer :: TVp(:), Tp(:), Vp(:)
    integer :: nr, ndof, ndof2, i, js, je, j, jj, k, kl0, l
!!$    real(kind=kreal) :: vnorm

    nr=BTmat%nr
    ndof=BTmat%ndof
    ndof2=ndof*ndof

    TV=0.d0

!!$    vnorm=0.d0
!!$    do i=1,nr*ndof
!!$      vnorm=vnorm+V(i)**2
!!$    enddo
!!$    write(0,*) 'vnorm:', sqrt(vnorm)

!$omp parallel default(none),private(i,TVp,js,je,j,jj,Tp,Vp,k,kl0,l), &
!$omp&  shared(nr,TV,ndof,BTmat,ndof2,V)
!$omp do
    do i=1,nr
      TVp=>TV((i-1)*ndof+1:i*ndof)
      js=BTmat%index(i-1)+1
      je=BTmat%index(i)
      do j=js,je
        jj=BTmat%item(j)
        Tp=>BTmat%A((j-1)*ndof2+1:j*ndof2)
        Vp=>V((jj-1)*ndof+1:jj*ndof)
        do k=1,ndof
          kl0=(k-1)*ndof
          do l=1,ndof
            TVp(k)=TVp(k)+Tp(kl0+l)*Vp(l)
          enddo
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel
  end subroutine hecmw_localmat_mulvec

  subroutine hecmw_trimatmul_TtKT_mpc(hecMESH, hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_matrix) :: BTmat, BTtmat
    integer(kind=kint), allocatable :: iwS(:)
    integer(kind=kint) :: n_mpc, ndof
    integer(kind=kint) :: i, k, kk
    n_mpc=hecMESH%mpc%n_mpc
    ndof=hecMAT%NDOF
    allocate(iwS(n_mpc))
    do i=1,n_mpc
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      iwS(i)=kk
    enddo
    call make_BTmat_mpc(hecMESH, BTmat)
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: BTmat(MPC)'
    !call hecmw_localmat_write(BTmat,700+hecmw_comm_get_rank())
    ! call make_BTtmat_mpc(hecMESH, BTtmat)
    call hecmw_localmat_transpose(BTmat, BTtmat)
    ! if (hecmw_localmat_equal(BTtmat, BTtmat2) == 0) then
    !   write(0,*) 'ERROR: BTtmat2 is incorrect!!!'
    ! else
    !   write(0,*) 'DEBUG: BTtmat2 is correct'
    ! endif
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: BTtmat(MPC)'
    !call hecmw_localmat_write(BTtmat,700+hecmw_comm_get_rank())
    call hecmw_trimatmul_TtKT(BTtmat, hecMAT, BTmat, iwS, n_mpc)
    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(BTtmat)
    ! call hecmw_localmat_free(BTtmat2)
    deallocate(iwS)
  end subroutine hecmw_trimatmul_TtKT_mpc

  subroutine make_BTmat_mpc(hecMESH, BTmat)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(out) :: BTmat
    type (hecmwST_local_matrix) :: Tmat
    integer(kind=kint) :: n_mpc, ndof
    integer(kind=kint) :: i,j,k,js,jj,kk
    n_mpc=hecMESH%mpc%n_mpc
    ndof=hecMESH%n_dof
    Tmat%nr=hecMESH%n_node*ndof
    Tmat%nc=Tmat%nr
    if (n_mpc > 0) then
      Tmat%nnz=Tmat%nr+hecMESH%mpc%mpc_index(n_mpc)-2*n_mpc
    else
      Tmat%nnz=Tmat%nr
    endif
    Tmat%ndof=1
    allocate(Tmat%index(0:Tmat%nr))
    allocate(Tmat%item(Tmat%nnz), Tmat%A(Tmat%nnz))
    ! count nonzero in each row
    Tmat%index(1:Tmat%nr)=1
    do i=1,n_mpc
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      Tmat%index(kk)=hecMESH%mpc%mpc_index(i)-hecMESH%mpc%mpc_index(i-1)-1
    enddo
    ! index
    Tmat%index(0)=0
    do i=1,Tmat%nr
      Tmat%index(i)=Tmat%index(i-1)+Tmat%index(i)
    enddo
    ! diag
    do i=1,Tmat%nr
      js=Tmat%index(i-1)+1
      Tmat%item(js)=i
      Tmat%A(js)=1.d0
    enddo
    ! others
    do i=1,n_mpc
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      js=Tmat%index(kk-1)+1
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Tmat%item(js)=jj
        Tmat%A(js)=-hecMESH%mpc%mpc_val(j)
        js=js+1
      enddo
    enddo
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: Tmat(MPC)'
    !call hecmw_localmat_write(Tmat,700+hecmw_comm_get_rank())
    call hecmw_localmat_blocking(Tmat, ndof, BTmat)
    call hecmw_localmat_free(Tmat)
  end subroutine make_BTmat_mpc

  subroutine make_BTtmat_mpc(hecMESH, BTtmat)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    type (hecmwST_local_matrix) :: Ttmat
    integer(kind=kint) :: n_mpc, ndof
    integer(kind=kint) :: i,j,k,js,je,jj,kk
    integer(kind=kint), allocatable :: iw(:)
    n_mpc=hecMESH%mpc%n_mpc
    ndof=hecMESH%n_dof
    Ttmat%nr=hecMESH%n_node*ndof
    Ttmat%nc=Ttmat%nr
    if (n_mpc > 0) then
      Ttmat%nnz=Ttmat%nr+hecMESH%mpc%mpc_index(n_mpc)-2*n_mpc
    else
      Ttmat%nnz=Ttmat%nr
    endif
    Ttmat%ndof=1
    allocate(Ttmat%index(0:Ttmat%nr))
    allocate(Ttmat%item(Ttmat%nnz), Ttmat%A(Ttmat%nnz))
    ! count nonzero in each row
    Ttmat%index(1:Ttmat%nr)=1
    do i=1,n_mpc
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      Ttmat%index(kk)=0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Ttmat%index(jj)=Ttmat%index(jj)+1
      enddo
    enddo
    ! index
    Ttmat%index(0)=0
    do i=1,Ttmat%nr
      Ttmat%index(i)=Ttmat%index(i-1)+Ttmat%index(i)
    enddo
    ! diag
    do i=1,Ttmat%nr
      js=Ttmat%index(i-1)+1
      je=Ttmat%index(i)
      if (js <= je) then
        Ttmat%item(js)=i
        Ttmat%A(js)=1.d0
      endif
    enddo
    ! others
    allocate(iw(Ttmat%nr))
    iw(:)=1
    do i=1,n_mpc
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        js=Ttmat%index(jj-1)+1+iw(jj)
        Ttmat%item(js)=kk
        Ttmat%A(js)=-hecMESH%mpc%mpc_val(j)
        iw(jj)=iw(jj)+1
      enddo
    enddo
    deallocate(iw)
    !write(700+hecmw_comm_get_rank(),*) 'DEBUG: Ttmat(MPC)'
    !call hecmw_localmat_write(Ttmat,700+hecmw_comm_get_rank())
    call hecmw_localmat_blocking(Ttmat, ndof, BTtmat)
    call hecmw_localmat_free(Ttmat)
  end subroutine make_BTtmat_mpc

  subroutine hecmw_localmat_transpose(Tmat, Ttmat)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    type (hecmwST_local_matrix), intent(out) :: Ttmat
    integer(kind=kint), allocatable :: iw(:)
    integer(kind=kint) :: i, j, jj, ndof, ndof2, k, idof, jdof
    allocate(iw(Tmat%nc))
    iw = 0
    do i = 1, Tmat%nr
      do j = Tmat%index(i-1)+1, Tmat%index(i)
        jj = Tmat%item(j)
        iw(jj) = iw(jj) + 1
      enddo
    enddo
    Ttmat%nr = Tmat%nc
    Ttmat%nc = Tmat%nr
    Ttmat%nnz = Tmat%nnz
    Ttmat%ndof = Tmat%ndof
    ndof = Tmat%ndof
    ndof2 = ndof * ndof
    allocate(Ttmat%index(0:Ttmat%nr))
    allocate(Ttmat%item(Ttmat%nnz))
    allocate(Ttmat%A(Ttmat%nnz*ndof2))
    Ttmat%index(0) = 0
    do i = 1, Ttmat%nr
      Ttmat%index(i) = Ttmat%index(i-1) + iw(i)
      iw(i) = Ttmat%index(i-1) + 1
    enddo
    do i = 1, Tmat%nr
      do j = Tmat%index(i-1)+1, Tmat%index(i)
        jj = Tmat%item(j)
        k = iw(jj)
        Ttmat%item( k ) = i
        do idof = 1, ndof
          do jdof = 1, ndof
            Ttmat%A((k-1)*ndof2+(idof-1)*ndof+jdof) = &
                 Tmat%A((j-1)*ndof2+(jdof-1)*ndof+idof)
          enddo
        enddo
        iw(jj) = k + 1
      enddo
    enddo
  end subroutine hecmw_localmat_transpose

  function hecmw_localmat_equal(Tmat1, Tmat2)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat1, Tmat2
    integer(kind=kint) :: hecmw_localmat_equal
    integer(kind=kint) :: i, j, k0, k, ndof, ndof2
    hecmw_localmat_equal = 0
    if (Tmat1%nr /= Tmat2%nr) return
    if (Tmat1%nc /= Tmat2%nc) return
    if (Tmat1%nnz /= Tmat2%nnz) return
    if (Tmat1%ndof /= Tmat2%ndof) return
    ndof = Tmat1%ndof
    ndof2 = ndof * ndof
    do i = 1, Tmat1%nr
      if (Tmat1%index(i) /= Tmat2%index(i)) return
      do j = Tmat1%index(i-1)+1, Tmat1%index(i)
        if (Tmat1%item(j) /= Tmat2%item(j)) return
        k0 = (j-1)*ndof2
        do k = 1, ndof2
          if (Tmat1%A(k0+k) /= Tmat2%A(k0+k)) return
        enddo
      enddo
    enddo
    hecmw_localmat_equal = 1
  end function hecmw_localmat_equal

end module hecmw_local_matrix
