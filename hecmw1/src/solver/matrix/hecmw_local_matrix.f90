module hecmw_local_matrix
  use hecmw_util

  private
  public :: hecmwST_local_matrix
  public :: hecmw_localmat_write
  public :: hecmw_localmat_blocking
  public :: hecmw_localmat_free
  public :: hecmw_localmat_mulvec
  public :: hecmw_trimatmul_TtKT

  type hecmwST_local_matrix
    integer :: n, nnz, ndof
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: item(:)
    real(kind=kreal), pointer :: A(:)
  end type hecmwST_local_matrix

contains

  subroutine hecmw_localmat_write(Tmat,iunit)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer(kind=kint), intent(in) :: iunit
    integer(kind=kint) :: n, ndof, ndof2, i, js, je, j, jj

    n=Tmat%n
    ndof=Tmat%ndof
    ndof2=ndof*ndof

    write(iunit,*) 'i, j, A'
    do i=1,n
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

    if (mod(Tmat%n, ndof) /= 0) then
      write(0,*) Tmat%n, ndof
      stop 'ERROR: blocking_Tmat failed'
    endif
    BTmat%n=Tmat%n/ndof
    BTmat%ndof=ndof

    allocate(iw(BTmat%n))
    allocate(BTmat%index(0:BTmat%n))

    BTmat%index(0)=0
    do i=1,BTmat%n
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

    BTmat%nnz=BTmat%index(BTmat%n)
    allocate(BTmat%item(BTmat%nnz))
    allocate(BTmat%A(BTmat%nnz*ndof2))
    BTmat%A=0.d0

    do i=1,BTmat%n
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
      if (icnt /= BTmat%index(i)-BTmat%index(i-1)) stop 'ERROR: blocking Tmat'
      ! call qsort(iw, 1, icnt)
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
    Tmat%n=0
    Tmat%nnz=0
    Tmat%ndof=0
  end subroutine hecmw_localmat_free

  subroutine hecmw_trimatmul_TtKT(BTtmat, hecMAT, BTmat, iwS, num_lagrange, hecTKT)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_matrix), intent(out) :: hecTKT
    type (hecmwST_local_matrix) :: BTtKT

    ! perform three matrices multiplication for elimination
    call trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    !call hecmw_localmat_write(BTtKT, 0)

    ! place 1s where the DOF is eliminated
    call place_one_on_diag(BTtKT, iwS, num_lagrange)

    ! make new HECMW matrix
    call make_new_hecmat(hecMAT, BTtKT, hecTKT)
    call hecmw_localmat_free(BTtKT)
  end subroutine hecmw_trimatmul_TtKT

  subroutine trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_matrix), intent(out) :: BTtKT
    integer :: n, ndof, ndof2, i, icnt, js, je, j, jj, ks, ke, k, kk
    integer :: ls, le, l, ll, m, ms, me, mm
    integer, allocatable :: iw(:)
    real(kind=kreal), pointer :: Ttp(:), Kp(:), Tp(:), TtKTp(:)

    n=hecMAT%N
    ndof=hecMAT%NDOF
    ndof2=ndof*ndof

    BTtKT%n=n
    BTtKT%ndof=ndof
    allocate(BTtKT%index(0:BTtKT%n))

    allocate(iw(n))

    BTtKT%index(0)=0
    do i=1,n
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
      BTtKT%index(i)=BTtKT%index(i-1)+icnt
    enddo
    !write(0,*) BTtKT%index(1:n)-BTtKT%index(0:n-1)

    BTtKT%nnz=BTtKT%index(n)
    allocate(BTtKT%item(BTtKT%nnz))
    allocate(BTtKT%A(BTtKT%nnz*ndof2))
    BTtKT%item=0
    BTtKT%A=0.d0

    do i=1,n
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
      if (icnt == 0) icnt=1
      ! error check!
      !write(0,*) BTtKT%item(ms:ms-1+icnt)
      !write(0,*) BTtKT%index(i)-BTtKT%index(i-1), icnt
      if (ms-1+icnt /= BTtKT%index(i)) stop 'ERROR: trimatmul'
    enddo

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
    integer(kind=kint) :: n, ndof, ndof2, ilag, i, idof, js, je, j, jj

    n=BTtKT%n
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
    integer :: n, ndof, ndof2, i, nl, nu, js, je, j, jj
    integer :: ksl, ksu, k, idof, idx

    n=hecMAT%N
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
    hecMAT%indexL(0)=0
    hecMAT%indexU(0)=0
    do i=1,n
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
      hecMAT%indexL(i)=hecMAT%indexL(i-1)+nl
      hecMAT%indexU(i)=hecMAT%indexU(i-1)+nu
    enddo
    hecMAT%NPL=hecMAT%indexL(n)
    hecMAT%NPU=hecMAT%indexU(n)

    ! allocate new hecMAT
    allocate(hecMAT%itemL(hecMAT%NPL), hecMAT%itemU(hecMAT%NPU))
    allocate(hecMAT%AL(hecMAT%NPL*ndof2), hecMAT%AU(hecMAT%NPU*ndof2))
    hecMAT%itemL=0
    hecMAT%itemU=0
    hecMAT%D=0.d0
    hecMAT%AL=0.d0
    hecMAT%AU=0.d0

    ! copy from BTtKT to hecMAT
    do i=1,n
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
      if (ksl+nl /= hecMAT%indexL(i)+1) stop 'ERROR: indexL'
      if (ksu+nu /= hecMAT%indexU(i)+1) stop 'ERROR: indexU'
    enddo

    do i=1,hecMAT%NPL
      if (hecMAT%itemL(i) <= 0) stop 'ERROR: negative itemL'
      if (hecMAT%itemL(i) > n) stop 'ERROR: too big itemL'
    enddo
    do i=1,hecMAT%NPU
      if (hecMAT%itemU(i) <= 0) stop 'ERROR: negative itemU'
      if (hecMAT%itemU(i) > n) stop 'ERROR: too big itemU'
    enddo
  end subroutine replace_hecmat

  subroutine make_new_hecmat(hecMAT, BTtKT, hecTKT)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_local_matrix), intent(in) :: BTtKT
    type(hecmwST_matrix), intent(out) :: hecTKT
    integer(kind=kint) :: n, ndof, ndof2

    n=hecMAT%N
    ndof=hecMAT%NDOF
    ndof2=ndof*ndof

    hecTKT%N   =hecMAT%N
    hecTKT%NP  =hecMAT%NP
    hecTKT%NDOF=hecMAT%NDOF

    allocate(hecTKT%D(n*ndof2))
    allocate(hecTKT%B(n*ndof))
    allocate(hecTKT%X(n*ndof))

    allocate(hecTKT%indexL(0:n))
    allocate(hecTKT%indexU(0:n))

    hecTKT%Iarray=hecMAT%Iarray
    hecTKT%Rarray=hecMAT%Rarray

    call replace_hecmat(hecTKT, BTtKT)
  end subroutine make_new_hecmat

  subroutine hecmw_localmat_mulvec(BTmat, V, TV)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    real(kind=kreal), intent(in), target :: V(:)
    real(kind=kreal), intent(out), target :: TV(:)
    real(kind=kreal), pointer :: TVp(:), Tp(:), Vp(:)
    integer :: n, ndof, ndof2, i, js, je, j, jj, k, kl0, l
!!$    real(kind=kreal) :: vnorm

    n=BTmat%n
    ndof=BTmat%ndof
    ndof2=ndof*ndof

    TV=0.d0

!!$    vnorm=0.d0
!!$    do i=1,n*ndof
!!$      vnorm=vnorm+V(i)**2
!!$    enddo
!!$    write(0,*) 'vnorm:', sqrt(vnorm)

    do i=1,n
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
  end subroutine hecmw_localmat_mulvec

end module hecmw_local_matrix
