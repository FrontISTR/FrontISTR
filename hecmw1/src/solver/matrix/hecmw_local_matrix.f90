!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_local_matrix
  use hecmw_util
  use hecmw_pair_array

  private
  public :: hecmwST_local_matrix
  public :: hecmw_localmat_write
  public :: hecmw_localmat_blocking
  public :: hecmw_localmat_free
  public :: hecmw_localmat_mulvec
  public :: hecmw_trimatmul_TtKT
  public :: hecmw_trimatmul_TtKT_serial
  public :: hecmw_trimatmul_TtKT_mpc
  public :: hecmw_localmat_transpose
  public :: hecmw_localmat_assemble
  public :: hecmw_localmat_add
  public :: hecmw_localmat_init_with_hecmat
  public :: hecmw_localmat_add_hecmat
  public :: hecmw_localmat_multmat
  public :: hecmw_localmat_make_hecmat
  public :: hecmw_localmat_shrink_comm_table

  type hecmwST_local_matrix
    integer :: nr, nc, nnz, ndof
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: item(:)
    real(kind=kreal), pointer :: A(:)
  end type hecmwST_local_matrix

  integer(kind=kint), parameter :: cNCOL_ITEM = 3 !< num of column items to be migrated (2 or 3)
  integer(kind=kint), parameter :: cLID = 1       !< index for local ID in belonging rank (node_ID(2*i-1))
  integer(kind=kint), parameter :: cRANK = 2      !< index for belonging rank (node_ID(2*i))
  integer(kind=kint), parameter :: cGID = 3       !< index for global ID (used only when cNCOL_ITEM==3)

  integer(kind=kint), parameter :: DEBUG = 0
  integer(kind=kint), parameter :: DEBUG_MATRIX = 0
  integer(kind=kint), parameter :: TIMER = 0

contains

  subroutine hecmw_localmat_write(Tmat,iunit)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer(kind=kint), intent(in) :: iunit
    integer(kind=kint) :: nr, nc, nnz, ndof, ndof2, i, js, je, j, jj
    character(len=64) :: fmt
    nr=Tmat%nr
    nc=Tmat%nc
    nnz=Tmat%nnz
    ndof=Tmat%ndof
    ndof2=ndof*ndof
    write(iunit,'(a,4i10)') 'nr, nc, nnz, ndof', nr, nc, nnz, ndof
    write(iunit,'(a)') 'i, j, A'
    write(fmt,'(a,i0,a)') '(',ndof,'f12.3)'
    do i=1,nr
      js=Tmat%index(i-1)+1
      je=Tmat%index(i)
      do j=js,je
        jj=Tmat%item(j)
        if (ndof==1) then
          write(iunit,'(2i10,f12.3)') i, jj, Tmat%A(j)
        else
          write(iunit,'(2i10)') i, jj
          write(iunit,fmt) Tmat%A((j-1)*ndof2+1:j*ndof2)
        endif
      enddo
    enddo
  end subroutine hecmw_localmat_write

  subroutine hecmw_localmat_write_size(Tmat,iunit)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer(kind=kint), intent(in) :: iunit
    integer(kind=kint) :: nr, nc, nnz, ndof
    nr=Tmat%nr
    nc=Tmat%nc
    nnz=Tmat%nnz
    ndof=Tmat%ndof
    write(iunit,'(a,4i10)') 'nr, nc, nnz, ndof', nr, nc, nnz, ndof
  end subroutine hecmw_localmat_write_size

  subroutine hecmw_localmat_write_ij(Tmat,iunit)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Tmat
    integer(kind=kint), intent(in) :: iunit
    integer(kind=kint) :: nr, nc, nnz, ndof, i, js, je, j, jj
    nr=Tmat%nr
    nc=Tmat%nc
    nnz=Tmat%nnz
    ndof=Tmat%ndof
    write(iunit,'(a,4i10)') 'nr, nc, nnz, ndof', nr, nc, nnz, ndof
    write(iunit,'(a)') 'i, j'
    do i=1,nr
      js=Tmat%index(i-1)+1
      je=Tmat%index(i)
      do j=js,je
        jj=Tmat%item(j)
        write(iunit,'(2i10)') i, jj
      enddo
    enddo
  end subroutine hecmw_localmat_write_ij

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
    if (associated(Tmat%item)) deallocate(Tmat%item)
    if (associated(Tmat%A)) deallocate(Tmat%A)
    Tmat%nr=0
    Tmat%nc=0
    Tmat%nnz=0
    Tmat%ndof=0
  end subroutine hecmw_localmat_free

  subroutine hecmw_trimatmul_TtKT(hecMESH, BTtmat, hecMAT, BTmat, &
      iwS, num_lagrange, hecTKT)
    use hecmw_matrix_misc
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_local_matrix), intent(inout) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_matrix), intent(inout) :: hecTKT
    if (hecMESH%n_neighbor_pe == 0) then
      call hecmw_trimatmul_TtKT_serial(hecMESH, BTtmat, hecMAT, BTmat, &
           iwS, num_lagrange, hecTKT)
    else
      call hecmw_trimatmul_TtKT_parallel(hecMESH, BTtmat, hecMAT, BTmat, &
           iwS, num_lagrange, hecTKT)
    endif
  end subroutine hecmw_trimatmul_TtKT

  subroutine hecmw_trimatmul_TtKT_serial(hecMESH, BTtmat, hecMAT, BTmat, &
      iwS, num_lagrange, hecTKT)
    use hecmw_matrix_misc
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_matrix), intent(inout) :: hecTKT
    type (hecmwST_local_matrix) :: BTtKT
    real(kind=kreal) :: num

    ! perform three matrices multiplication for elimination
    call trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    call debug_write_matrix(BTtKT, 'BTtKT(MPC)', DEBUG_MATRIX)

    ! place small numbers where the DOF is eliminated
    !num = hecmw_mat_diag_max(hecMAT, hecMESH) * 1.0d-10
    num = 1.d0
    call place_num_on_diag(BTtKT, iwS, num_lagrange, num)
    call debug_write_matrix(BTtKT, 'BTtKT(MPC) (place 1.0 on slave diag)', DEBUG_MATRIX)

    ! make_new HECMW matrix
    call make_new_hecmat(hecMAT, BTtKT, hecTKT)
    call hecmw_localmat_free(BTtKT)
  end subroutine hecmw_trimatmul_TtKT_serial

  subroutine hecmw_trimatmul_TtKT_parallel(hecMESH, BTtmat, hecMAT, BTmat, &
      iwS, num_lagrange, hecTKT)
    use hecmw_matrix_misc
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_local_matrix), intent(inout) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (hecmwST_matrix), intent(inout) :: hecTKT
    type (hecmwST_local_matrix) :: BKmat, BTtKmat, BTtKTmat
    real(kind=kreal) :: num
    real(kind=kreal) :: t0, t1

    ! perform three matrices multiplication for elimination
    t0 = hecmw_wtime()
    call hecmw_localmat_init_with_hecmat(BKmat, hecMAT)
    call debug_write_matrix(BKmat, 'BKmat (hecMAT)', DEBUG_MATRIX)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "#### hecmw_trimatmul_TtKT_parallel (1) : ",t1-t0

    t0 = hecmw_wtime()
    call hecmw_localmat_multmat(BTtmat, BKmat, hecMESH, BTtKmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply Tt and K done'
    call debug_write_matrix(BTtKmat, 'BTtKmat', DEBUG_MATRIX)
    call hecmw_localmat_free(BKmat)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "#### hecmw_trimatmul_TtKT_parallel (2) : ",t1-t0

    t0 = hecmw_wtime()
    call hecmw_localmat_multmat(BTtKmat, BTmat, hecMESH, BTtKTmat)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply TtK and T done'
    call debug_write_matrix(BTtKTmat, 'BTtKTmat', DEBUG_MATRIX)
    call hecmw_localmat_free(BTtKmat)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "#### hecmw_trimatmul_TtKT_parallel (3) : ",t1-t0

    t0 = hecmw_wtime()
    ! place small numbers where the DOF is eliminated
    !num = hecmw_mat_diag_max(hecMAT, hecMESH) * 1.0d-10
    num = 1.d0
    call place_num_on_diag(BTtKTmat, iwS, num_lagrange, num)
    if (DEBUG >= 2) then
      write(700+hecmw_comm_get_rank(),*) 'num_lagrange =', num_lagrange
      if (DEBUG >= 3) then
        write(700+hecmw_comm_get_rank(),*) 'iwS(1:num_lagrange)'
        write(700+hecmw_comm_get_rank(),*) iwS(1:num_lagrange)
      endif
    endif
    call debug_write_matrix(BTtKTmat, 'BTtKTmat (place 1.0 on slave diag)', DEBUG_MATRIX)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "#### hecmw_trimatmul_TtKT_parallel (4) : ",t1-t0

    t0 = hecmw_wtime()
    ! make_new HECMW matrix
    call make_new_hecmat(hecMAT, BTtKTmat, hecTKT)
    call hecmw_localmat_free(BTtKTmat)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "#### hecmw_trimatmul_TtKT_parallel (5) : ",t1-t0
  end subroutine hecmw_trimatmul_TtKT_parallel

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

  subroutine place_num_on_diag(BTtKT, iwS, num_lagrange, num)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTtKT
    integer(kind=kint), intent(in) :: iwS(:)
    integer(kind=kint), intent(in) :: num_lagrange
    real(kind=kreal), intent(in) :: num
    integer(kind=kint) :: ndof, ndof2, ilag, i, idof, js, je, j, jj
    integer(kind=kint) :: nmissing, k, ks, ke
    integer(kind=kint), allocatable :: missing(:), cnt(:)
    integer(kind=kint), pointer :: index(:), item(:)
    real(kind=kreal), pointer :: A(:)

    ndof=BTtKT%ndof
    ndof2=ndof*ndof

    ! check if there are places
    allocate(missing(num_lagrange))
    nmissing = 0
    outer1: do ilag=1,num_lagrange
      i=(iwS(ilag)-1)/ndof+1
      idof=mod(iwS(ilag)-1, ndof)+1
      js=BTtKT%index(i-1)+1
      je=BTtKT%index(i)
      do j=js,je
        jj=BTtKT%item(j)
        if (jj==i) cycle outer1  ! found place
      enddo
      ! not found
      do k=1,nmissing
        if (missing(k) == i) cycle outer1 ! already marked as missing
      enddo
      nmissing = nmissing + 1
      missing(nmissing) = i
    enddo outer1

    ! if not, reallocate
    if (nmissing > 0) then
      allocate(cnt(BTtKT%nr))
      allocate(index(0:BTtKT%nr))
      do i=1,BTtKT%nr
        cnt(i) = BTtKT%index(i) - BTtKT%index(i-1)
      enddo
      do i=1,nmissing
        cnt(missing(i)) = cnt(missing(i)) + 1
      enddo
      call make_index(BTtKT%nr, cnt, index)
      allocate(item(BTtKT%nnz + nmissing))
      allocate(A(ndof2 * (BTtKT%nnz + nmissing)))
      do i=1,BTtKT%nr
        ks=index(i-1)+1
        js=BTtKT%index(i-1)+1
        je=BTtKT%index(i)
        item(ks:ks+(je-js))=BTtKT%item(js:je)
        A(ndof2*(ks-1)+1:ndof2*(ks+(je-js)))=BTtKT%A(ndof2*(js-1)+1:ndof2*je)
      enddo
      do i=1,nmissing
        ke=index(missing(i))
        item(ke)=missing(i)
        A(ndof2*(ke-1)+1:ndof2*ke)=0.d0
      enddo
      deallocate(BTtKT%index)
      deallocate(BTtKT%item)
      deallocate(BTtKT%A)
      BTtKT%index => index
      BTtKT%item => item
      BTtKT%A => A
      BTtKT%nnz = index(BTtKT%nr)
      deallocate(cnt)
    endif
    deallocate(missing)

    ! place num
    outer: do ilag=1,num_lagrange
      i=(iwS(ilag)-1)/ndof+1
      idof=mod(iwS(ilag)-1, ndof)+1
      js=BTtKT%index(i-1)+1
      je=BTtKT%index(i)
      do j=js,je
        jj=BTtKT%item(j)
        if (jj==i) then
          !write(0,*) ilag, i, idof
          BTtKT%A((j-1)*ndof2+(idof-1)*ndof+idof)=num
          cycle outer
        endif
      enddo
    enddo outer
  end subroutine place_num_on_diag

  subroutine replace_hecmat(hecMAT, BTtKT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_matrix), intent(in) :: BTtKT
    integer :: nr, nc, ndof, ndof2, i, nl, nu, js, je, j, jj
    integer :: ksl, ksu, k

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
    integer(kind=kint) :: nr, nc, ndof, ndof2

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

    if (associated(hecTKT%D)) deallocate(hecTKT%D)
    allocate(hecTKT%D(nc*ndof2))

    if (associated(hecTKT%indexL)) deallocate(hecTKT%indexL)
    if (associated(hecTKT%indexU)) deallocate(hecTKT%indexU)
    allocate(hecTKT%indexL(0:nc))
    allocate(hecTKT%indexU(0:nc))

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

  subroutine hecmw_trimatmul_TtKT_mpc(hecMESH, hecMAT, hecTKT)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_matrix), intent(inout) :: hecTKT
    type (hecmwST_local_matrix) :: BTmat, BTtmat
    integer(kind=kint), allocatable :: iwS(:)
    integer(kind=kint) :: ndof, n_mpc, i_mpc
    integer(kind=kint) :: i, j, k, kk, ilag
    integer(kind=kint) :: num_lagrange
    real(kind=kreal) :: t0, t1
    t0 = hecmw_wtime()
    ndof=hecMAT%NDOF
    n_mpc=0
    OUTER: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      n_mpc=n_mpc+1
    enddo OUTER
    allocate(iwS(n_mpc))
    i_mpc=0
    OUTER2: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER2
      enddo
      i_mpc=i_mpc+1
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      iwS(i_mpc)=kk
    enddo OUTER2
    if (DEBUG >= 2) then
      write(700+hecmw_comm_get_rank(),*) 'DEBUG: n_mpc, slaves',n_mpc,iwS(1:n_mpc)
    endif
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "### hecmw_trimatmul_TtKT_mpc (1) : ",t1-t0
    t0 = hecmw_wtime()
    call make_BTmat_mpc(hecMESH, ndof, BTmat)
    call debug_write_matrix(BTmat, 'BTmat(MPC)', DEBUG_MATRIX)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "### hecmw_trimatmul_TtKT_mpc (2) : ",t1-t0
    t0 = hecmw_wtime()
    ! call make_BTtmat_mpc(hecMESH, ndof, BTtmat)
    call hecmw_localmat_transpose(BTmat, BTtmat)
    ! if (hecmw_localmat_equal(BTtmat, BTtmat2) == 0) then
    !   write(0,*) 'ERROR: BTtmat2 is incorrect!!!'
    ! else
    !   write(0,*) 'DEBUG: BTtmat2 is correct'
    ! endif
    call debug_write_matrix(BTtmat, 'BTtmat(MPC)', DEBUG_MATRIX)

    if (DEBUG >= 3) then
      write(700+hecmw_comm_get_rank(),*) 'hecMESH%node_ID before trimatmul_TtKT'
      do i=hecMESH%nn_internal+1, hecMESH%n_node
        write(700+hecmw_comm_get_rank(),*) i,hecMESH%node_ID(2*i-1),hecMESH%node_ID(2*i),hecMESH%global_node_ID(i)
      enddo
    endif
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "### hecmw_trimatmul_TtKT_mpc (3) : ",t1-t0
    t0 = hecmw_wtime()
    call hecmw_trimatmul_TtKT(hecMESH, BTtmat, hecMAT, BTmat, iwS, n_mpc, hecTKT)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "### hecmw_trimatmul_TtKT_mpc (4) : ",t1-t0
    t0 = hecmw_wtime()
    if (DEBUG >= 3) then
      write(700+hecmw_comm_get_rank(),*) 'hecMESH%node_ID after trimatmul_TtKT'
      do i=hecMESH%nn_internal+1, hecMESH%n_node
        write(700+hecmw_comm_get_rank(),*) i,hecMESH%node_ID(2*i-1),hecMESH%node_ID(2*i),hecMESH%global_node_ID(i)
      enddo
    endif

    if (associated(hecTKT%B)) deallocate(hecTKT%B)
    if (associated(hecTKT%X)) deallocate(hecTKT%X)
    num_lagrange = size(hecMAT%B) - hecMAT%NP*ndof
    allocate(hecTKT%B(ndof*hecTKT%NP + num_lagrange))
    allocate(hecTKT%X(ndof*hecTKT%NP + num_lagrange))
    hecTKT%B(:) = 0.d0
    hecTKT%X(:) = 0.d0
    do i=1, ndof*hecMAT%NP
      hecTKT%B(i) = hecMAT%B(i)
      hecTKT%X(i) = hecMAT%X(i)
    enddo
    do i=1, num_lagrange
      hecTKT%B(ndof*hecTKT%NP+i) = hecMAT%B(ndof*hecMAT%NP+i)
      hecTKT%X(ndof*hecTKT%NP+i) = hecMAT%X(ndof*hecMAT%NP+i)
    enddo
    do ilag=1,n_mpc
      hecTKT%X(iwS(ilag)) = 0.d0
    enddo

    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(BTtmat)
    ! call hecmw_localmat_free(BTtmat2)
    deallocate(iwS)
    t1 = hecmw_wtime()
    if (TIMER >= 1) write(0, '(A,f10.4)') "### hecmw_trimatmul_TtKT_mpc (5) : ",t1-t0
  end subroutine hecmw_trimatmul_TtKT_mpc

  subroutine make_BTmat_mpc(hecMESH, ndof, BTmat)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    type (hecmwST_local_matrix), intent(out) :: BTmat
    type (hecmwST_local_matrix) :: Tmat
    integer(kind=kint) :: n_mpc
    integer(kind=kint) :: i,j,k,js,jj,kk
    n_mpc=0
    OUTER: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      n_mpc=n_mpc+1
    enddo OUTER
    Tmat%nr=hecMESH%n_node*ndof
    Tmat%nc=Tmat%nr
    Tmat%ndof=1
    allocate(Tmat%index(0:Tmat%nr))
    ! count nonzero in each row
    Tmat%index(1:Tmat%nr)=1
    OUTER2: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER2
      enddo
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      Tmat%index(kk)=hecMESH%mpc%mpc_index(i)-hecMESH%mpc%mpc_index(i-1)-1
    enddo OUTER2
    ! index
    Tmat%index(0)=0
    do i=1,Tmat%nr
      Tmat%index(i)=Tmat%index(i-1)+Tmat%index(i)
    enddo
    Tmat%nnz=Tmat%index(Tmat%nr)
    allocate(Tmat%item(Tmat%nnz), Tmat%A(Tmat%nnz))
    ! diag
    do i=1,Tmat%nr
      js=Tmat%index(i-1)+1
      Tmat%item(js)=i
      Tmat%A(js)=1.d0
    enddo
    ! others
    OUTER3: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER3
      enddo
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      js=Tmat%index(kk-1)+1
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Tmat%item(js)=jj
        Tmat%A(js)=-hecMESH%mpc%mpc_val(j)
        js=js+1
      enddo
    enddo OUTER3
    !call debug_write_matrix(Tmat, 'Tmat(MPC)', DEBUG_MATRIX)
    call hecmw_localmat_blocking(Tmat, ndof, BTmat)
    call hecmw_localmat_free(Tmat)
  end subroutine make_BTmat_mpc

  subroutine make_BTtmat_mpc(hecMESH, ndof, BTtmat)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    type (hecmwST_local_matrix) :: Ttmat
    integer(kind=kint) :: n_mpc
    integer(kind=kint) :: i,j,k,js,je,jj,kk
    integer(kind=kint), allocatable :: iw(:)
    n_mpc=0
    OUTER: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER
      enddo
      n_mpc=n_mpc+1
    enddo OUTER
    Ttmat%nr=hecMESH%n_node*ndof
    Ttmat%nc=Ttmat%nr
    Ttmat%ndof=1
    allocate(Ttmat%index(0:Ttmat%nr))
    ! count nonzero in each row
    Ttmat%index(1:Ttmat%nr)=1
    OUTER2: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER2
      enddo
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      Ttmat%index(kk)=0
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        Ttmat%index(jj)=Ttmat%index(jj)+1
      enddo
    enddo OUTER2
    ! index
    Ttmat%index(0)=0
    do i=1,Ttmat%nr
      Ttmat%index(i)=Ttmat%index(i-1)+Ttmat%index(i)
    enddo
    Ttmat%nnz=Ttmat%index(Ttmat%nr)
    allocate(Ttmat%item(Ttmat%nnz), Ttmat%A(Ttmat%nnz))
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
    OUTER3: do i=1,hecMESH%mpc%n_mpc
      do j= hecMESH%mpc%mpc_index(i-1) + 1, hecMESH%mpc%mpc_index(i)
        if (hecMESH%mpc%mpc_dof(j) > ndof) cycle OUTER3
      enddo
      k=hecMESH%mpc%mpc_index(i-1)+1
      kk=ndof*(hecMESH%mpc%mpc_item(k)-1)+hecMESH%mpc%mpc_dof(k)
      do j= hecMESH%mpc%mpc_index(i-1) + 2, hecMESH%mpc%mpc_index(i)
        jj = ndof * (hecMESH%mpc%mpc_item(j) - 1) + hecMESH%mpc%mpc_dof(j)
        js=Ttmat%index(jj-1)+1+iw(jj)
        Ttmat%item(js)=kk
        Ttmat%A(js)=-hecMESH%mpc%mpc_val(j)
        iw(jj)=iw(jj)+1
      enddo
    enddo OUTER3
    deallocate(iw)
    !call debug_write_matrix(Ttmat, 'Ttmat(MPC)', DEBUG_MATRIX)
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

!!!
!!! Subroutines for parallel contact analysis with iterative linear solver
!!!

  subroutine hecmw_localmat_assemble(BTmat, hecMESH, hecMESHnew)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTmat
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint) :: nn_int, np, ndof, ndof2, nr_ext, nnz_ext
    integer(kind=kint), allocatable :: exp_rows_index(:), exp_cols_index(:)
    integer(kind=kint), allocatable :: exp_rows_item(:,:), exp_cols_item(:,:)
    type (hecmwST_local_matrix), allocatable :: BT_ext(:)
    type (hecmwST_local_matrix) :: BT_int
    type (hecmwST_local_matrix) :: BTnew
    ! some checks
    if (DEBUG >= 1) write(0,*) 'DEBUG: nr,nc,nnz,ndof',BTmat%nr,BTmat%nc,BTmat%nnz,BTmat%ndof
    if (BTmat%nr /= hecMESH%n_node) stop 'ERROR: invalid size in hecmw_localmat_assemble'
    !
    nn_int = hecMESH%nn_internal
    np = hecMESH%n_node
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    !
    nr_ext = np - nn_int
    nnz_ext = BTmat%index(np) - BTmat%index(nn_int)
    !
    call prepare_BT_ext(BTmat, hecMESH, exp_rows_index, exp_rows_item, BT_ext)
    if (DEBUG >= 1) write(0,*) 'DEBUG: prepare_BT_ext done'
    !
    call prepare_column_info(hecMESH, BT_ext, exp_cols_index, exp_cols_item)
    if (DEBUG >= 1) write(0,*) 'DEBUG: prepare_column info done'
    !
    call send_BT_ext_and_recv_BT_int(hecMESH, exp_rows_index, exp_rows_item, BT_ext, &
         exp_cols_index, exp_cols_item, BT_int, hecMESHnew)
    if (DEBUG >= 1) write(0,*) 'DEBUG: send BT_ext and recv BT_int done'
    !
    !write(0,*) 'BTmat%ndof,BT_int%ndof',BTmat%ndof,BT_int%ndof
    call hecmw_localmat_add(BTmat, BT_int, BTnew)
    if (DEBUG >= 1) write(0,*) 'DEBUG: localmat_add done'
    !
    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(BT_int)
    !
    BTmat%nr = BTnew%nr
    BTmat%nc = BTnew%nc
    BTmat%nnz = BTnew%nnz
    BTmat%ndof = BTnew%ndof
    BTmat%index => BTnew%index
    BTmat%item => BTnew%item
    BTmat%A => BTnew%A
    !
    ! hecMESH%n_node = hecMESHnew%n_node
    ! hecMESH%n_neighbor_pe = hecMESHnew%n_neighbor_pe
    ! deallocate(hecMESH%neighbor_pe)
    ! deallocate(hecMESH%import_index)
    ! deallocate(hecMESH%export_index)
    ! deallocate(hecMESH%import_item)
    ! deallocate(hecMESH%export_item)
    ! deallocate(hecMESH%node_ID)
    ! deallocate(hecMESH%global_node_ID)
    ! hecMESH%neighbor_pe => hecMESHnew%neighbor_pe
    ! hecMESH%import_index => hecMESHnew%import_index
    ! hecMESH%export_index => hecMESHnew%export_index
    ! hecMESH%import_item => hecMESHnew%import_item
    ! hecMESH%export_item => hecMESHnew%export_item
    ! hecMESH%node_ID => hecMESHnew%node_ID
    ! hecMESH%global_node_ID => hecMESHnew%global_node_ID
    !
    if (DEBUG >= 1) write(0,*) 'DEBUG: update BTmat and hecMESH done'
  end subroutine hecmw_localmat_assemble

  subroutine prepare_BT_ext(BTmat, hecMESH, exp_rows_index, exp_rows_item, BT_ext)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), allocatable, intent(out) :: exp_rows_index(:)
    integer(kind=kint), allocatable, intent(out) :: exp_rows_item(:,:)
    type (hecmwST_local_matrix), allocatable, intent(out) :: BT_ext(:)
    integer(kind=kint), allocatable :: incl_nz(:), exp_cols_per_row(:), exp_rows_per_rank(:)
    integer(kind=kint) :: nn_int
    logical, parameter :: FLG_CHECK_NONZERO_NUMERICALLY = .true.
    nn_int = hecMESH%nn_internal
    !
    if (FLG_CHECK_NONZERO_NUMERICALLY) then
      ! efficient for assembling conMAT which has same non-zero profile as hecMAT
      ! but only dofs related to cntact are actually non-zero
      call check_external_nz_blocks(BTmat, nn_int, incl_nz)
    else
      ! probably good enough to assemble T or T^t
      call incl_all_external_nz_blocks(BTmat, nn_int, incl_nz)
    endif
    !
    call count_ext_rows_with_nz(BTmat, nn_int, incl_nz, exp_cols_per_row)
    !
    call count_exp_rows_per_rank(hecMESH, exp_cols_per_row, exp_rows_per_rank)
    !
    allocate(exp_rows_index(0:hecMESH%n_neighbor_pe))
    call make_index(hecMESH%n_neighbor_pe, exp_rows_per_rank, exp_rows_index)
    !write(0,*) 'exp_rows_index',exp_rows_index(:)
    !
    deallocate(exp_rows_per_rank)
    !
    call make_exp_rows_item(hecMESH, exp_cols_per_row, exp_rows_index, exp_rows_item)
    !
    deallocate(exp_cols_per_row)
    !
    allocate(BT_ext(hecMESH%n_neighbor_pe))
    call extract_BT_ext(hecMESH, BTmat, incl_nz, exp_rows_index, exp_rows_item, BT_ext)
    !
    deallocate(incl_nz)
  end subroutine prepare_BT_ext

  subroutine check_external_nz_blocks(BTmat, nn_internal, incl_nz)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: nn_internal
    integer(kind=kint), allocatable, intent(out) :: incl_nz(:)
    integer(kind=kint) :: ndof2, i0, nnz_ext, i, k, nnz_blk
    if (nn_internal > BTmat%nr) stop 'ERROR: invalid nn_internal'
    ndof2 = BTmat%ndof ** 2
    i0 = BTmat%index(nn_internal)
    nnz_ext = BTmat%index(BTmat%nr) - i0
    allocate(incl_nz(nnz_ext))
    nnz_blk = 0
    do i = 1, nnz_ext
      incl_nz(i) = 0
      do k = 1, ndof2
        if (BTmat%A(ndof2*(i0+i-1)+k) /= 0.0d0) then
          incl_nz(i) = 1
          nnz_blk = nnz_blk + 1
          exit
        endif
      enddo
    enddo
    if (DEBUG >= 1) write(0,*) 'DEBUG: nnz_blk',nnz_blk
  end subroutine check_external_nz_blocks

  subroutine incl_all_external_nz_blocks(BTmat, nn_internal, incl_nz)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: nn_internal
    integer(kind=kint), allocatable, intent(out) :: incl_nz(:)
    integer(kind=kint) :: i0, nnz_ext
    if (nn_internal > BTmat%nr) stop 'ERROR: invalid nn_internal'
    i0 = BTmat%index(nn_internal)
    nnz_ext = BTmat%index(BTmat%nr) - i0
    allocate(incl_nz(nnz_ext))
    incl_nz(1:nnz_ext) = 1
  end subroutine incl_all_external_nz_blocks

  subroutine count_ext_rows_with_nz(BTmat, nn_internal, incl_nz, exp_cols_per_row)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: nn_internal
    integer(kind=kint), intent(in) :: incl_nz(:)
    integer(kind=kint), allocatable, intent(out) :: exp_cols_per_row(:)
    integer(kind=kint) :: nr_ext, nnz_int, i, irow, js, je, j, jcol
    nr_ext = BTmat%nr - nn_internal
    nnz_int = BTmat%index(nn_internal)
    allocate(exp_cols_per_row(nr_ext))
    exp_cols_per_row(:) = 0
    do i = 1, nr_ext
      irow = nn_internal+i
      js = BTmat%index(irow-1)+1
      je = BTmat%index(irow)
      do j = js, je
        jcol = BTmat%item(j)
        if (incl_nz(j-nnz_int) == 1) exp_cols_per_row(i) = exp_cols_per_row(i) + 1
      enddo
    enddo
    !write(0,*) 'exp_cols_per_row',exp_cols_per_row(:)
  end subroutine count_ext_rows_with_nz

  subroutine count_exp_rows_per_rank(hecMESH, exp_cols_per_row, exp_rows_per_rank)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: exp_cols_per_row(:)
    integer(kind=kint), allocatable, intent(out) :: exp_rows_per_rank(:)
    integer(kind=kint) :: nn_int, np, nr_ext, i, irow, exp_rank, idom
    allocate(exp_rows_per_rank(hecMESH%n_neighbor_pe))
    exp_rows_per_rank(1:hecMESH%n_neighbor_pe) = 0
    nn_int = hecMESH%nn_internal
    np = hecMESH%n_node
    nr_ext = np - nn_int
    do i = 1, nr_ext
      if (exp_cols_per_row(i) > 0) then
        irow = nn_int + i
        exp_rank = hecMESH%node_ID(2*irow)
        call rank_to_idom(hecMESH, exp_rank, idom)
        exp_rows_per_rank(idom) = exp_rows_per_rank(idom) + 1
      endif
    enddo
    !write(0,*) 'exp_rows_per_rank',exp_rows_per_rank(:)
  end subroutine count_exp_rows_per_rank

  subroutine rank_to_idom(hecMESH, rank, idom)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: rank
    integer(kind=kint), intent(out) :: idom
    integer(kind=kint) :: i
    do i = 1, hecMESH%n_neighbor_pe
      if (hecMESH%neighbor_pe(i) == rank) then
        idom = i
        return
      endif
    enddo
    stop 'ERROR: exp_rank not found in neighbor_pe'
  end subroutine rank_to_idom

  subroutine make_index(len, cnt, index)
    implicit none
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(in) :: cnt(len)
    integer(kind=kint), intent(out) :: index(0:)
    integer(kind=kint) :: i
    ! write(0,*) 'make_index: len',len
    index(0) = 0
    do i = 1,len
      index(i) = index(i-1) + cnt(i)
    enddo
  end subroutine make_index

  subroutine make_exp_rows_item(hecMESH, exp_cols_per_row, exp_rows_index, exp_rows_item)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: exp_cols_per_row(:)
    integer(kind=kint), allocatable, intent(in) :: exp_rows_index(:)
    integer(kind=kint), allocatable, intent(out) :: exp_rows_item(:,:)
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: nn_int, np, nr_ext, i, irow, exp_rank, idom, idx
    allocate(exp_rows_item(2,exp_rows_index(hecMESH%n_neighbor_pe)))
    allocate(cnt(hecMESH%n_neighbor_pe))
    cnt(:) = 0
    nn_int = hecMESH%nn_internal
    np = hecMESH%n_node
    nr_ext = np - nn_int
    do i = 1, nr_ext
      if (exp_cols_per_row(i) > 0) then
        irow = nn_int + i
        exp_rank = hecMESH%node_ID(2*irow)
        call rank_to_idom(hecMESH, exp_rank, idom)
        cnt(idom) = cnt(idom) + 1
        idx = exp_rows_index(idom-1) + cnt(idom)
        exp_rows_item(1,idx) = irow
        exp_rows_item(2,idx) = exp_cols_per_row(i)
      endif
    enddo
    !write(0,*) 'cnt',cnt(:)
    do idom = 1, hecMESH%n_neighbor_pe
      if (cnt(idom) /= exp_rows_index(idom)-exp_rows_index(idom-1)) stop 'ERROR: make exp_rows_item'
    enddo
    !write(0,*) 'exp_rows_item(1,:)',exp_rows_item(1,:)
    !write(0,*) 'exp_rows_item(2,:)',exp_rows_item(2,:)
  end subroutine make_exp_rows_item

  subroutine extract_BT_ext(hecMESH, BTmat, incl_nz, exp_rows_index, exp_rows_item, BT_ext)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: incl_nz(:)
    integer(kind=kint), allocatable, intent(in) :: exp_rows_index(:)
    integer(kind=kint), intent(in) :: exp_rows_item(:,:)
    type (hecmwST_local_matrix), allocatable, intent(out) :: BT_ext(:)
    integer(kind=kint) :: ndof, ndof2, nn_int, nnz_int, idom, j, idx, ncol, cnt, jrow, ks, ke, k, kcol
    allocate(BT_ext(hecMESH%n_neighbor_pe))
    ndof = BTmat%ndof
    ndof2 = ndof * ndof
    nn_int = hecMESH%nn_internal
    nnz_int = BTmat%index(nn_int)
    do idom = 1, hecMESH%n_neighbor_pe
      BT_ext(idom)%nr = exp_rows_index(idom) - exp_rows_index(idom-1)
      BT_ext(idom)%nc = BTmat%nc
      BT_ext(idom)%nnz = 0
      BT_ext(idom)%ndof = ndof
      allocate(BT_ext(idom)%index(0:BT_ext(idom)%nr))
      BT_ext(idom)%index(0) = 0
      do j = 1, BT_ext(idom)%nr
        idx = exp_rows_index(idom-1) + j
        ncol = exp_rows_item(2,idx)
        BT_ext(idom)%index(j) = BT_ext(idom)%index(j-1) + ncol
      enddo
      BT_ext(idom)%nnz = BT_ext(idom)%index(BT_ext(idom)%nr)
      if (DEBUG >= 1) write(0,*) 'DEBUG: idom,nr,nc,nnz,ndof', &
           idom,BT_ext(idom)%nr,BT_ext(idom)%nc,BT_ext(idom)%nnz,BT_ext(idom)%ndof
      allocate(BT_ext(idom)%item(BT_ext(idom)%nnz))
      allocate(BT_ext(idom)%A(BT_ext(idom)%nnz * ndof2))
      cnt = 0
      do j = 1, BT_ext(idom)%nr
        idx = exp_rows_index(idom-1) + j
        jrow = exp_rows_item(1,idx)
        if (jrow < 1 .or. BTmat%nr < jrow) stop 'ERROR: extract BT_ext: jrow'
        ks = BTmat%index(jrow-1)+1
        ke = BTmat%index(jrow)
        do k = ks, ke
          kcol = BTmat%item(k)
          if (incl_nz(k-nnz_int) == 0) cycle
          cnt = cnt + 1
          BT_ext(idom)%item(cnt) = kcol
          BT_ext(idom)%A(ndof2*(cnt-1)+1:ndof2*cnt) = BTmat%A(ndof2*(k-1)+1:ndof2*k)
        enddo
        if (cnt /= BT_ext(idom)%index(j)) stop 'ERROR: extract BT_ext'
      enddo
      ! write(label,'(a,i0,a)') 'BT_ext(',idom,')'
      ! call debug_write_matrix(Bt_ext(idom), label, DEBUG_MATRIX)
    enddo
  end subroutine extract_BT_ext

  subroutine prepare_column_info(hecMESH, BT_ext, exp_cols_index, exp_cols_item)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BT_ext(:)
    integer(kind=kint), allocatable, intent(out) :: exp_cols_index(:)
    integer(kind=kint), allocatable, intent(out) :: exp_cols_item(:,:)
    !
    call make_exp_cols_index(hecMESH%n_neighbor_pe, BT_ext, exp_cols_index)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make exp_cols_index done'
    if (DEBUG >= 3) write(0,*) '    DEBUG3: exp_cols_index', exp_cols_index(0:hecMESH%n_neighbor_pe)
    !
    ! (col ID, rank, global ID)
    !
    call make_exp_cols_item(hecMESH, BT_ext, exp_cols_index, exp_cols_item)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make exp_cols_item done'
    ! if (DEBUG >= 3) write(0,*) '    DEBUG3: exp_cols_item', exp_cols_item(1:cNCOL_ITEM,1:exp_cols_index(hecMESH%n_neighbor_pe))
  end subroutine prepare_column_info

  subroutine make_exp_cols_index(nnb, BT_ext, exp_cols_index)
    implicit none
    integer(kind=kint), intent(in) :: nnb
    type (hecmwST_local_matrix), intent(in) :: BT_ext(:)
    integer(kind=kint), allocatable, intent(out) :: exp_cols_index(:)
    integer(kind=kint) :: idom
    allocate(exp_cols_index(0:nnb))
    exp_cols_index(0) = 0
    do idom = 1, nnb
      exp_cols_index(idom) = exp_cols_index(idom-1) + BT_ext(idom)%nnz
    enddo
  end subroutine make_exp_cols_index

  subroutine make_exp_cols_item(hecMESH, BT_ext, exp_cols_index, exp_cols_item)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BT_ext(:)
    integer(kind=kint), allocatable, intent(in) :: exp_cols_index(:)
    integer(kind=kint), allocatable, intent(out) :: exp_cols_item(:,:)
    integer(kind=kint) :: cnt, idom, j, jcol
    allocate(exp_cols_item(cNCOL_ITEM,exp_cols_index(hecMESH%n_neighbor_pe)))
    cnt = 0
    do idom = 1, hecMESH%n_neighbor_pe
      do j = 1, BT_ext(idom)%nnz
        cnt = cnt + 1
        jcol = BT_ext(idom)%item(j)
        ! if (DEBUG >= 3) write(0,*) '    DEBUG3: idom,j,cnt,jcol,nn_internal,n_node',&
        !      idom,j,cnt,jcol,hecMESH%nn_internal,hecMESH%n_node
        ! if (DEBUG >= 3) write(0,*) '    DEBUG3: size of exp_cols_item',size(exp_cols_item)
        ! if (DEBUG >= 3) write(0,*) '    DEBUG3: size of node_ID',size(hecMESH%node_ID)
        ! if (DEBUG >= 3) write(0,*) '    DEBUG3: size of global_node_ID',size(hecMESH%global_node_ID)
        exp_cols_item(cLID,cnt) = hecMESH%node_ID(2*jcol-1)
        exp_cols_item(cRANK,cnt) = hecMESH%node_ID(2*jcol)
        if (cNCOL_ITEM >= 3) exp_cols_item(cGID,cnt) = hecMESH%global_node_ID(jcol)
        ! if (DEBUG >= 3) write(0,*) '    DEBUG3: lid,rank(,gid)',exp_cols_item(1:cNCOL_ITEM,cnt)
      enddo
      if (cnt /= exp_cols_index(idom)) stop 'ERROR: make exp_cols_item'
    enddo
  end subroutine make_exp_cols_item

  subroutine send_BT_ext_and_recv_BT_int(hecMESH, exp_rows_index, exp_rows_item, BT_ext, &
       exp_cols_index, exp_cols_item, BT_int, hecMESHnew)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), allocatable, intent(inout) :: exp_rows_index(:), exp_cols_index(:)
    integer(kind=kint), allocatable, intent(inout) :: exp_rows_item(:,:), exp_cols_item(:,:)
    type (hecmwST_local_matrix), allocatable, intent(inout) :: BT_ext(:)
    type (hecmwST_local_matrix), intent(out) :: BT_int
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint), allocatable :: imp_rows_index(:), imp_cols_index(:)
    integer(kind=kint), allocatable :: imp_rows_item(:,:), imp_cols_item(:,:)
    real(kind=kreal), allocatable :: imp_vals_item(:)
    integer(kind=kint), allocatable :: map(:), add_nodes(:,:)
    integer(kind=kint) :: ndof, ndof2, idom, n_add_node, i0
    if (hecMESH%n_neighbor_pe == 0) return
    ndof = BT_ext(1)%ndof
    ndof2 = ndof*ndof
    !
    call convert_rowID_to_remote_localID(hecMESH, exp_rows_index(hecMESH%n_neighbor_pe), exp_rows_item)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: convert rowID to remote localID done'
    !
    call send_recv_BT_ext_nr_nnz(hecMESH, BT_ext, imp_rows_index, imp_cols_index)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: send recv BT_ext nr and nnz done'
    !
    call send_recv_BT_ext_contents(hecMESH, BT_ext, &
         exp_rows_index, exp_cols_index, exp_rows_item, exp_cols_item, &
         imp_rows_index, imp_cols_index, &
         imp_rows_item, imp_cols_item, imp_vals_item)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: send recv BT_ext contents done'
    !
    do idom = 1, hecMESH%n_neighbor_pe
      call hecmw_localmat_free(BT_ext(idom))
    enddo
    deallocate(BT_ext)
    !
    call allocate_BT_int(hecMESH, ndof, imp_rows_index, imp_rows_item, BT_int)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: allocate BT_int done'
    !
    ! call copy_mesh(hecMESH, hecMESHnew)
    ! if (DEBUG >= 2) write(0,*) '  DEBUG2: copy mesh done'
    !
    call map_imported_cols(hecMESHnew, imp_cols_index(hecMESH%n_neighbor_pe), &
         imp_cols_item, n_add_node, add_nodes, map, i0)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: map imported cols done'
    !
    call update_comm_table(hecMESHnew, n_add_node, add_nodes, i0)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: update comm_table done'
    !
    BT_int%nc = hecMESHnew%n_node
    !
    call copy_vals_to_BT_int(hecMESH%n_neighbor_pe, imp_rows_index, imp_cols_index, &
         imp_rows_item, map, ndof2, imp_vals_item, BT_int)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: copy vals to BT_int done'
    !
    deallocate(imp_rows_index)
    deallocate(imp_cols_index)
    deallocate(imp_rows_item)
    deallocate(imp_cols_item)
    deallocate(imp_vals_item)
    deallocate(map)
    !
    call sort_and_uniq_rows(BT_int)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: sort and uniq rows of BT_int done'
  end subroutine send_BT_ext_and_recv_BT_int

  subroutine convert_rowID_to_remote_localID(hecMESH, len, exp_rows_item)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(out) :: exp_rows_item(:,:)
    integer(kind=kint) :: i
    do i = 1, len
      exp_rows_item(1,i) = hecMESH%node_ID(2 * exp_rows_item(1,i) - 1)
    enddo
  end subroutine convert_rowID_to_remote_localID

  subroutine send_recv_BT_ext_nr_nnz(hecMESH, BT_ext, imp_rows_index, imp_cols_index)
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BT_ext(:)
    integer(kind=kint), allocatable, intent(out) :: imp_rows_index(:), imp_cols_index(:)
    integer(kind=kint) :: nnb, idom, irank, tag, recvbuf(2)
    integer(kind=kint), allocatable :: sendbuf(:,:)
    integer(kind=kint), allocatable :: requests(:)
    integer(kind=kint), allocatable :: statuses(:,:)
    nnb = hecMESH%n_neighbor_pe
    allocate(imp_rows_index(0:nnb))
    allocate(imp_cols_index(0:nnb))
    allocate(requests(nnb))
    allocate(statuses(HECMW_STATUS_SIZE, nnb))
    allocate(sendbuf(2,nnb))
    do idom = 1, nnb
      irank = hecMESH%neighbor_pe(idom)
      ! nr = exp_rows_per_rank(idom)
      sendbuf(1,idom) = BT_ext(idom)%nr
      sendbuf(2,idom) = BT_ext(idom)%nnz
      tag=2001
      call HECMW_ISEND_INT(sendbuf(1,idom), 2, irank, tag, hecMESH%MPI_COMM, &
           requests(idom))
    enddo
    imp_rows_index(0) = 0
    imp_cols_index(0) = 0
    do idom = 1, nnb
      irank = hecMESH%neighbor_pe(idom)
      tag = 2001
      call HECMW_RECV_INT(recvbuf, 2, irank, tag, &
           hecMESH%MPI_COMM, statuses(:,1))
      imp_rows_index(idom) = imp_rows_index(idom-1) + recvbuf(1)
      imp_cols_index(idom) = imp_cols_index(idom-1) + recvbuf(2)
    enddo
    call HECMW_Waitall(nnb, requests, statuses)
    deallocate(requests)
    deallocate(statuses)
    deallocate(sendbuf)
  end subroutine send_recv_BT_ext_nr_nnz

  subroutine send_recv_BT_ext_contents(hecMESH, BT_ext, &
       exp_rows_index, exp_cols_index, exp_rows_item, exp_cols_item, &
       imp_rows_index, imp_cols_index, &
       imp_rows_item, imp_cols_item, imp_vals_item)
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_local_matrix), intent(in) :: BT_ext(:)
    integer(kind=kint), allocatable, intent(inout) :: exp_rows_index(:), exp_cols_index(:)
    integer(kind=kint), allocatable, intent(inout) :: exp_rows_item(:,:), exp_cols_item(:,:)
    integer(kind=kint), allocatable, intent(in) :: imp_rows_index(:), imp_cols_index(:)
    integer(kind=kint), allocatable, intent(out) :: imp_rows_item(:,:), imp_cols_item(:,:)
    real(kind=kreal), allocatable, intent(out) :: imp_vals_item(:)
    integer(kind=kint) :: nnb, ndof2, n_send, idom, irank, tag, nr, nnz
    integer(kind=kint), allocatable :: requests(:)
    integer(kind=kint), allocatable :: statuses(:,:)
    nnb = hecMESH%n_neighbor_pe
    if (nnb < 1) return
    ndof2 = BT_ext(1)%ndof ** 2
    allocate(imp_rows_item(2,imp_rows_index(nnb)))
    allocate(imp_cols_item(cNCOL_ITEM,imp_cols_index(nnb)))
    allocate(imp_vals_item(ndof2*imp_cols_index(nnb)))
    allocate(requests(3*nnb))
    allocate(statuses(HECMW_STATUS_SIZE, 3*nnb))
    n_send = 0
    do idom = 1, nnb
      irank = hecMESH%neighbor_pe(idom)
      if (BT_ext(idom)%nr > 0) then
        n_send = n_send + 1
        tag = 2002
        call HECMW_ISEND_INT(exp_rows_item(1,exp_rows_index(idom-1)+1), &
             2*BT_ext(idom)%nr, irank, tag, hecMESH%MPI_COMM, &
             requests(n_send))
        n_send = n_send + 1
        tag = 2003
        call HECMW_ISEND_INT(exp_cols_item(1,exp_cols_index(idom-1)+1), &
             cNCOL_ITEM*BT_ext(idom)%nnz, irank, tag, hecMESH%MPI_COMM, &
             requests(n_send))
        n_send = n_send + 1
        tag = 2004
        call HECMW_ISEND_R(BT_ext(idom)%A, ndof2*BT_ext(idom)%nnz, irank, &
             tag, hecMESH%MPI_COMM, requests(n_send))
      endif
    enddo
    do idom = 1, nnb
      irank = hecMESH%neighbor_pe(idom)
      nr = imp_rows_index(idom) - imp_rows_index(idom-1)
      nnz = imp_cols_index(idom) - imp_cols_index(idom-1)
      if (nr > 0) then
        tag = 2002
        call HECMW_RECV_INT(imp_rows_item(1,imp_rows_index(idom-1)+1), &
             2*nr, irank, tag, hecMESH%MPI_COMM, statuses(:,1))
        tag = 2003
        call HECMW_RECV_INT(imp_cols_item(1,imp_cols_index(idom-1)+1), &
             cNCOL_ITEM*nnz, irank, tag, hecMESH%MPI_COMM, statuses(:,1))
        tag = 2004
        call HECMW_RECV_R(imp_vals_item(ndof2*imp_cols_index(idom-1)+1), &
             ndof2*nnz, irank, tag, hecMESH%MPI_COMM, statuses(:,1))
      endif
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
    deallocate(exp_rows_index)
    deallocate(exp_rows_item)
    deallocate(exp_cols_index)
    deallocate(exp_cols_item)
  end subroutine send_recv_BT_ext_contents

  subroutine allocate_BT_int(hecMESH, ndof, imp_rows_index, imp_rows_item, BT_int)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    integer(kind=kint), allocatable, intent(in) :: imp_rows_index(:), imp_rows_item(:,:)
    type (hecmwST_local_matrix), intent(out) :: BT_int
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: idom, is, ie, i, irow, ncol, ndof2
    ndof2 = ndof*ndof
    BT_int%nr = hecMESH%nn_internal
    BT_int%nc = hecMESH%n_node
    BT_int%nnz = 0
    BT_int%ndof = ndof
    allocate(cnt(BT_int%nr))
    cnt(:) = 0
    do idom = 1, hecMESH%n_neighbor_pe
      is = imp_rows_index(idom-1)+1
      ie = imp_rows_index(idom)
      do i = is, ie
        irow = imp_rows_item(1,i)
        ncol = imp_rows_item(2,i)
        if (irow < 1 .or. BT_int%nr < irow) stop 'ERROR: allocate BT_int'
        cnt(irow) = cnt(irow) + ncol   !!! might include duplicate cols
      enddo
    enddo
    !
    allocate(BT_int%index(0:BT_int%nr))
    call make_index(BT_int%nr, cnt, BT_int%index)
    !
    BT_int%nnz = BT_int%index(BT_int%nr)
    allocate(BT_int%item(BT_int%nnz))
    allocate(BT_int%A(BT_int%nnz * ndof2))
    BT_int%A(:) = 0.d0
  end subroutine allocate_BT_int

  subroutine copy_mesh(src, dst)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: src
    type (hecmwST_local_mesh), intent(out) :: dst
    dst%zero = src%zero
    dst%MPI_COMM = src%MPI_COMM
    dst%PETOT = src%PETOT
    dst%PEsmpTOT = src%PEsmpTOT
    dst%my_rank = src%my_rank
    dst%n_subdomain = src%n_subdomain
    dst%n_node = src%n_node
    dst%nn_internal = src%nn_internal
    dst%n_dof = src%n_dof
    dst%n_neighbor_pe = src%n_neighbor_pe
    allocate(dst%neighbor_pe(dst%n_neighbor_pe))
    dst%neighbor_pe(:) = src%neighbor_pe(:)
    allocate(dst%import_index(0:dst%n_neighbor_pe))
    allocate(dst%export_index(0:dst%n_neighbor_pe))
    dst%import_index(:)= src%import_index(:)
    dst%export_index(:)= src%export_index(:)
    allocate(dst%import_item(dst%import_index(dst%n_neighbor_pe)))
    dst%import_item(:) = src%import_item(:)
    allocate(dst%export_item(dst%export_index(dst%n_neighbor_pe)))
    dst%export_item(:) = src%export_item(:)
    allocate(dst%node_ID(2*dst%n_node))
    dst%node_ID(1:2*dst%n_node) = src%node_ID(1:2*src%n_node)
    allocate(dst%global_node_ID(dst%n_node))
    dst%global_node_ID(1:dst%n_node) = src%global_node_ID(1:src%n_node)
    dst%mpc%n_mpc = 0
    dst%node => src%node
  end subroutine copy_mesh

  subroutine map_imported_cols(hecMESHnew, ncols, cols, n_add_node, add_nodes, map, i0)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint), intent(in) :: ncols
    integer(kind=kint), intent(in) :: cols(cNCOL_ITEM,ncols)
    integer(kind=kint), allocatable, intent(out) :: map(:)
    integer(kind=kint), intent(out) :: n_add_node
    integer(kind=kint), allocatable, intent(out) :: add_nodes(:,:)
    integer(kind=kint), intent(out) :: i0
    allocate(map(ncols))
    !
    call map_present_nodes(hecMESHnew, ncols, cols, map, n_add_node)
    !
    ! add nodes == unmapped nodes
    !
    call extract_add_nodes(ncols, cols, map, n_add_node, add_nodes)
    !
    call append_nodes(hecMESHnew, n_add_node, add_nodes, i0)
    !
    call map_additional_nodes(ncols, cols, n_add_node, add_nodes, i0, map)
  end subroutine map_imported_cols

  subroutine map_present_nodes(hecMESH, ncols, cols, map, n_add_node)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ncols
    integer(kind=kint), intent(in) :: cols(cNCOL_ITEM,ncols)
    integer(kind=kint), intent(out) :: map(ncols)
    integer(kind=kint), intent(out) :: n_add_node
    integer(kind=kint) :: i, j, lid, rank, llid, n_ext_node, idx
    integer(kind=kint), allocatable :: ext_node(:)
    type (hecmwST_pair_array) :: parray
    !
    call hecmw_pair_array_init(parray, hecMESH%n_node - hecMESH%nn_internal)
    do i = hecMESH%nn_internal + 1, hecMESH%n_node
      call hecmw_pair_array_append(parray, i, hecMESH%node_ID(2*i-1), hecMESH%node_ID(2*i))
    enddo
    call hecmw_pair_array_sort(parray)
    !
    n_add_node = 0
    n_ext_node = 0
    allocate(ext_node(ncols))
    !$omp parallel default(none), &
      !$omp& private(i,lid,rank,llid,idx,j), &
      !$omp& shared(ncols,hecMESH,cols,map,n_ext_node,ext_node,parray), &
      !$omp& reduction(+:n_add_node)
    !$omp do
    do i = 1, ncols
      lid = cols(cLID,i)
      rank = cols(cRANK,i)
      !   check rank
      if (rank == hecMESH%my_rank) then  !   internal: set mapping
        map(i) = lid
      else                               !   external
        !$omp atomic capture
        n_ext_node = n_ext_node + 1
        idx = n_ext_node
        !$omp end atomic
        ext_node(idx) = i
      endif
    enddo
    !$omp end do
    !$omp do
    do j = 1, n_ext_node
      i = ext_node(j)
      lid = cols(cLID,i)
      rank = cols(cRANK,i)
      !     search node_ID in external nodes
      llid = hecmw_pair_array_find_id(parray, lid, rank)
      if (llid > 0) then  !     found: set mapping
        map(i) = llid
      else                !     not found
        map(i) = -1
        n_add_node = n_add_node + 1
      endif
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(ext_node)
    !
    call hecmw_pair_array_finalize(parray)
  end subroutine map_present_nodes

  subroutine extract_add_nodes(ncols, cols, map, n_add_node, add_nodes)
    implicit none
    integer(kind=kint), intent(in) :: ncols
    integer(kind=kint), intent(in) :: cols(cNCOL_ITEM,ncols), map(ncols)
    integer(kind=kint), intent(inout) :: n_add_node
    integer(kind=kint), allocatable, intent(out) :: add_nodes(:,:)
    integer(kind=kint) :: cnt, i
    allocate(add_nodes(cNCOL_ITEM,n_add_node))
    cnt = 0
    do i = 1, ncols
      if (map(i) == -1) then
        cnt = cnt + 1
        add_nodes(1:cNCOL_ITEM,cnt) = cols(1:cNCOL_ITEM,i)
      endif
    enddo
    if (cnt /= n_add_node) stop 'ERROR: extract add_nodes'
    call sort_and_uniq_add_nodes(n_add_node, add_nodes)
  end subroutine extract_add_nodes

  subroutine sort_and_uniq_add_nodes(n_add_node, add_nodes)
    implicit none
    integer(kind=kint), intent(inout) :: n_add_node
    integer(kind=kint), intent(inout) :: add_nodes(cNCOL_ITEM,n_add_node)
    integer(kind=kint) :: ndup
    call sort_add_nodes(add_nodes, 1, n_add_node)
    call uniq_add_nodes(add_nodes, n_add_node, ndup)
    n_add_node = n_add_node - ndup
  end subroutine sort_and_uniq_add_nodes

  recursive subroutine sort_add_nodes(add_nodes, id1, id2)
    implicit none
    integer(kind=kint), intent(inout) :: add_nodes(:,:)
    integer(kind=kint), intent(in) :: id1, id2
    integer(kind=kint) :: center, left, right
    integer(kind=kint) :: pivot(cNCOL_ITEM), tmp(cNCOL_ITEM)
    if (id1 >= id2) return
    center = (id1 + id2) / 2
    pivot(1:cNCOL_ITEM) = add_nodes(1:cNCOL_ITEM,center)
    left = id1
    right = id2
    do
      do while ((add_nodes(cRANK,left) < pivot(cRANK)) .or. &
           (add_nodes(cRANK,left) == pivot(cRANK) .and. add_nodes(cLID,left) < pivot(cLID)))
        left = left + 1
      enddo
      do while ((pivot(cRANK) < add_nodes(cRANK,right)) .or. &
           (pivot(cRANK) == add_nodes(cRANK,right) .and. pivot(cLID) < add_nodes(cLID,right)))
        right = right - 1
      enddo
      if (left >= right) exit
      tmp(1:cNCOL_ITEM) = add_nodes(1:cNCOL_ITEM,left)
      add_nodes(1:cNCOL_ITEM,left) = add_nodes(1:cNCOL_ITEM,right)
      add_nodes(1:cNCOL_ITEM,right) = tmp(1:cNCOL_ITEM)
      left = left + 1
      right = right - 1
    enddo
    if (id1 < left-1) call sort_add_nodes(add_nodes, id1, left-1)
    if (right+1 < id2) call sort_add_nodes(add_nodes, right+1, id2)
    return
  end subroutine sort_add_nodes

  subroutine uniq_add_nodes(add_nodes, len, ndup)
    implicit none
    integer(kind=kint), intent(inout) :: add_nodes(:,:)
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(out) :: ndup
    integer(kind=kint) :: i
    ndup = 0
    do i = 2,len
      if (add_nodes(cLID,i) == add_nodes(cLID,i-1-ndup) .and. &
           add_nodes(cRANK,i) == add_nodes(cRANK,i-1-ndup)) then
        ndup = ndup + 1
      else if (ndup > 0) then
        add_nodes(1:cNCOL_ITEM,i-ndup) = add_nodes(1:cNCOL_ITEM,i)
      endif
    enddo
  end subroutine uniq_add_nodes

  subroutine search_add_nodes(n_add_node, add_nodes, rank, lid, idx)
    implicit none
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), intent(in) :: add_nodes(cNCOL_ITEM,n_add_node)
    integer(kind=kint), intent(in) :: rank
    integer(kind=kint), intent(in) :: lid
    integer(kind=kint), intent(out) :: idx
    integer(kind=kint) :: left, right, center
    left = 1
    right = n_add_node
    do while (left <= right)
      center = (left + right) / 2
      if ((rank == add_nodes(cRANK,center)) .and. (lid == add_nodes(cLID,center))) then
        idx = center
        return
      else if ((rank < add_nodes(cRANK,center)) .or. &
           (rank == add_nodes(cRANK,center) .and. lid < add_nodes(cLID,center))) then
        right = center - 1
      else if ((add_nodes(cRANK,center) < rank) .or. &
           (add_nodes(cRANK,center) == rank .and. add_nodes(cLID,center) < lid)) then
        left = center + 1
      endif
    end do
    idx = -1
  end subroutine search_add_nodes

  subroutine append_nodes(hecMESHnew, n_add_node, add_nodes, i0)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), intent(in) :: add_nodes(:,:)
    integer(kind=kint), intent(out) :: i0
    integer(kind=kint) :: n_node, i, ii
    integer(kind=kint), pointer :: node_ID(:), global_node_ID(:)
    i0 = hecMESHnew%n_node
    n_node = hecMESHnew%n_node + n_add_node
    allocate(node_ID(2*n_node))
    allocate(global_node_ID(n_node))
    do i = 1, hecMESHnew%n_node
      node_ID(2*i-1) = hecMESHnew%node_ID(2*i-1)
      node_ID(2*i  ) = hecMESHnew%node_ID(2*i  )
      global_node_ID(i) = hecMESHnew%global_node_ID(i)
    enddo
    do i = 1, n_add_node
      ii = hecMESHnew%n_node + i
      node_ID(2*ii-1) = add_nodes(cLID,i)
      node_ID(2*ii  ) = add_nodes(cRANK,i)
      if (cNCOL_ITEM >= 3) then
        global_node_ID(ii) = add_nodes(cGID,i)
      else
        global_node_ID(ii) = -1
      endif
    enddo
    deallocate(hecMESHnew%node_ID)
    deallocate(hecMESHnew%global_node_ID)
    hecMESHnew%n_node = n_node
    hecMESHnew%node_ID => node_ID
    hecMESHnew%global_node_ID => global_node_ID
  end subroutine append_nodes

  subroutine map_additional_nodes(ncols, cols, n_add_node, add_nodes, i0, map)
    implicit none
    integer(kind=kint), intent(in) :: ncols
    integer(kind=kint), intent(in) :: cols(cNCOL_ITEM,ncols)
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), intent(in) :: add_nodes(cNCOL_ITEM,n_add_node)
    integer(kind=kint), intent(in) :: i0
    integer(kind=kint), intent(inout) :: map(ncols)
    integer(kind=kint) :: i, j
    do i = 1, ncols
      if (map(i) > 0) cycle
      call search_add_nodes(n_add_node, add_nodes, cols(cRANK,i), cols(cLID,i), j)
      if (j == -1) stop 'ERROR: map_additional_nodes'
      map(i) = i0 + j
    enddo
  end subroutine map_additional_nodes

  subroutine update_comm_table(hecMESHnew, n_add_node, add_nodes, i0)
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), allocatable, intent(inout) :: add_nodes(:,:)
    integer(kind=kint), intent(in) :: i0
    integer(kind=kint), allocatable :: n_add_imp(:), add_imp_index(:)
    integer(kind=kint), allocatable :: add_imp_item_remote(:), add_imp_item_local(:)
    integer(kind=kint), allocatable :: n_add_exp(:), add_exp_index(:), add_exp_item(:)
    integer(kind=kint), allocatable :: n_new_imp(:), n_new_exp(:)
    integer(kind=kint) :: npe, nnb, comm, new_nnb
    integer(kind=kint), pointer :: nbpe(:), new_nbpe(:)
    integer(kind=kint), pointer :: import_index(:), export_index(:), import_item(:), export_item(:)
    integer(kind=kint), pointer :: new_import_index(:), new_export_index(:)
    integer(kind=kint), pointer :: new_import_item(:), new_export_item(:)
    npe = hecMESHnew%PETOT
    nnb = hecMESHnew%n_neighbor_pe
    comm = hecMESHnew%MPI_COMM
    nbpe => hecMESHnew%neighbor_pe
    import_index => hecMESHnew%import_index
    export_index => hecMESHnew%export_index
    import_item => hecMESHnew%import_item
    export_item => hecMESHnew%export_item
    !
    call count_add_imp_per_rank(n_add_node, add_nodes, npe, n_add_imp)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: count add_imp per rank done'
    !
    allocate(add_imp_index(0:npe))
    call make_index(npe, n_add_imp, add_imp_index)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: make add_imp_index done'
    !
    call make_add_imp_item(n_add_node, add_nodes, npe, i0, add_imp_index, &
         add_imp_item_remote, add_imp_item_local)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: make add_imp_item done'
    !
    deallocate(add_nodes)
    !
    ! all_to_all n_add_imp -> n_add_exp
    !
    allocate(n_add_exp(npe))
    call HECMW_ALLTOALL_INT(n_add_imp, 1, n_add_exp, 1, comm)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: alltoall n_add_imp to n_add_exp done'
    !
    allocate(add_exp_index(0:npe))
    call make_index(npe, n_add_exp, add_exp_index)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: make add_exp_index done'
    !
    call send_recv_add_imp_exp_item(npe, add_imp_index, add_imp_item_remote, &
         add_exp_index, add_exp_item, comm)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: send recv add_imp/exp_item done'
    !
    ! count new import
    !
    call count_new_comm_nodes(npe, nnb, nbpe, import_index, n_add_imp, n_new_imp)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: count new comm_nodes (import) done'
    !
    ! count new export
    !
    call count_new_comm_nodes(npe, nnb, nbpe, export_index, n_add_exp, n_new_exp)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: count new comm_nodes (export) done'
    !
    call update_neighbor_pe(npe, n_new_imp, n_new_exp, new_nnb, new_nbpe)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: update neighbor_pe done'
    !
    ! merge import table: import
    !
    call merge_comm_table(npe, nnb, nbpe, import_index, import_item, &
         new_nnb, new_nbpe, add_imp_index, add_imp_item_local, n_add_imp, n_new_imp, &
         new_import_index, new_import_item)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: merge comm_table (import) done'
    !
    deallocate(n_add_imp)
    deallocate(add_imp_index)
    deallocate(add_imp_item_remote, add_imp_item_local)
    deallocate(n_new_imp)
    !
    ! merge export table: export
    !
    call merge_comm_table(npe, nnb, nbpe, export_index, export_item, &
         new_nnb, new_nbpe, add_exp_index, add_exp_item, n_add_exp, n_new_exp, &
         new_export_index, new_export_item)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: merge comm_table (export) done'
    !
    deallocate(n_add_exp)
    deallocate(add_exp_index)
    deallocate(add_exp_item)
    deallocate(n_new_exp)
    !
    deallocate(nbpe)
    deallocate(import_index,import_item)
    deallocate(export_index,export_item)
    hecMESHnew%n_neighbor_pe = new_nnb
    hecMESHnew%neighbor_pe => new_nbpe
    hecMESHnew%import_index => new_import_index
    hecMESHnew%export_index => new_export_index
    hecMESHnew%import_item => new_import_item
    hecMESHnew%export_item => new_export_item
  end subroutine update_comm_table

  subroutine count_add_imp_per_rank(n_add_node, add_nodes, npe, n_add_imp)
    implicit none
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), intent(in) :: add_nodes(cNCOL_ITEM,n_add_node)
    integer(kind=kint), intent(in) :: npe
    integer(kind=kint), allocatable, intent(out) :: n_add_imp(:)
    integer(kind=kint) :: i, rank
    allocate(n_add_imp(npe))
    n_add_imp(:) = 0
    do i = 1, n_add_node
      rank = add_nodes(cRANK,i)
      n_add_imp(rank+1) = n_add_imp(rank+1) + 1
    enddo
  end subroutine count_add_imp_per_rank

  subroutine make_add_imp_item(n_add_node, add_nodes, npe, i0, add_imp_index, &
       add_imp_item_remote, add_imp_item_local)
    implicit none
    integer(kind=kint), intent(in) :: n_add_node
    integer(kind=kint), intent(in) :: add_nodes(cNCOL_ITEM,n_add_node)
    integer(kind=kint), intent(in) :: npe, i0
    integer(kind=kint), allocatable, intent(in) :: add_imp_index(:)
    integer(kind=kint), allocatable, intent(out) :: add_imp_item_remote(:), add_imp_item_local(:)
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: i, lid, rank, ipe
    allocate(add_imp_item_remote(add_imp_index(npe)))
    allocate(add_imp_item_local(add_imp_index(npe)))
    allocate(cnt(npe))
    cnt(:) = 0
    do i = 1, n_add_node
      lid = add_nodes(cLID,i)
      rank = add_nodes(cRANK,i)
      ipe = rank + 1
      cnt(ipe) = cnt(ipe) + 1
      add_imp_item_remote(add_imp_index(ipe-1) + cnt(ipe)) = lid
      add_imp_item_local(add_imp_index(ipe-1) + cnt(ipe)) = i0 + i
    enddo
    deallocate(cnt)
  end subroutine make_add_imp_item

  subroutine send_recv_add_imp_exp_item(npe, add_imp_index, add_imp_item_remote, &
       add_exp_index, add_exp_item, mpi_comm)
    use m_hecmw_comm_f
    implicit none
    integer(kind=kint), intent(in) :: npe
    integer(kind=kint), allocatable, intent(in) :: add_imp_index(:), add_imp_item_remote(:)
    integer(kind=kint), allocatable, intent(in) :: add_exp_index(:)
    integer(kind=kint), allocatable, intent(out) :: add_exp_item(:)
    integer(kind=kint), intent(in) :: mpi_comm
    integer(kind=kint) :: n_send, i, irank, is, ie, len, tag
    integer(kind=kint), allocatable :: requests(:)
    integer(kind=kint), allocatable :: statuses(:,:)
    allocate(add_exp_item(add_exp_index(npe)))
    allocate(requests(npe))
    allocate(statuses(HECMW_STATUS_SIZE, npe))
    n_send = 0
    do i = 1, npe
      irank = i-1
      is = add_imp_index(i-1)+1
      ie = add_imp_index(i)
      len = ie - is + 1
      if (len == 0) cycle
      tag = 4001
      n_send = n_send + 1
      call HECMW_ISEND_INT(add_imp_item_remote(is:ie), len, irank, tag, &
           mpi_comm, requests(n_send))
    enddo
    !
    do i = 1, npe
      irank = i-1
      is = add_exp_index(i-1)+1
      ie = add_exp_index(i)
      len = ie - is + 1
      if (len == 0) cycle
      tag = 4001
      call HECMW_RECV_INT(add_exp_item(is:ie), len, irank, tag, &
           mpi_comm, statuses(:,1))
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
  end subroutine send_recv_add_imp_exp_item

  subroutine count_new_comm_nodes(npe, org_nnb, org_nbpe, org_index, n_add, n_new)
    implicit none
    integer(kind=kint), intent(in) :: npe, org_nnb
    !integer(kind=kint), intent(in) :: org_nbpe(org_nnb), org_index(0:org_nnb), n_add(npe)
    integer(kind=kint), pointer, intent(in) :: org_nbpe(:), org_index(:)
    integer(kind=kint), intent(in) ::  n_add(:)
    integer(kind=kint), allocatable, intent(out) :: n_new(:)
    integer(kind=kint) :: i, irank, n_org
    allocate(n_new(npe))
    n_new(:) = n_add(:)
    do i = 1, org_nnb
      irank = org_nbpe(i)
      n_org = org_index(i) - org_index(i-1)
      n_new(irank+1) = n_new(irank+1) + n_org
    enddo
  end subroutine count_new_comm_nodes

  subroutine update_neighbor_pe(npe, n_new_imp, n_new_exp, &
       new_nnb, new_nbpe)
    implicit none
    integer(kind=kint), intent(in) :: npe
    integer(kind=kint), intent(in) :: n_new_imp(npe), n_new_exp(npe)
    integer(kind=kint), intent(out) :: new_nnb
    integer(kind=kint), pointer, intent(out) :: new_nbpe(:)
    integer(kind=kint) :: i
    new_nnb = 0
    do i = 1, npe
      if (n_new_imp(i) > 0 .or. n_new_exp(i) > 0) new_nnb = new_nnb+1
    enddo
    allocate(new_nbpe(new_nnb))
    new_nnb = 0
    do i = 1, npe
      if (n_new_imp(i) > 0 .or. n_new_exp(i) > 0) then
        new_nnb = new_nnb+1
        new_nbpe(new_nnb) = i-1
      endif
    enddo
  end subroutine update_neighbor_pe

  subroutine merge_comm_table(npe, org_nnb, org_nbpe, org_index, org_item, &
       new_nnb, new_nbpe, add_index, add_item, n_add, n_new, new_index, new_item)
    implicit none
    integer(kind=kint), intent(in) :: npe, org_nnb
    !integer(kind=kint), intent(in) :: org_nbpe(org_nnb), org_index(0:org_nnb), org_item(:)
    integer(kind=kint), pointer, intent(in) :: org_nbpe(:), org_index(:), org_item(:)
    integer(kind=kint), intent(in) :: new_nnb
    !integer(kind=kint), intent(in) :: new_nbpe(new_nnb), add_index(0:npe), add_item(:)
    integer(kind=kint), pointer, intent(in) :: new_nbpe(:)
    integer(kind=kint), allocatable, intent(in) :: add_index(:), add_item(:)
    integer(kind=kint), intent(in) :: n_add(npe), n_new(npe)
    integer(kind=kint), pointer, intent(out) :: new_index(:), new_item(:)
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: i, irank, j, jrank, i0, j0, len
    ! if (associated(new_index)) deallocate(new_index)
    ! if (associated(new_item)) deallocate(new_item)
    allocate(new_index(0:new_nnb))
    new_index(0) = 0
    do i = 1, new_nnb
      irank = new_nbpe(i)
      new_index(i) = new_index(i-1) + n_new(irank+1)
    enddo
    allocate(new_item(new_index(new_nnb)))
    allocate(cnt(npe))
    cnt(:) = 0
    j = 1
    jrank = new_nbpe(j)
    do i = 1, org_nnb
      if (org_index(i) - org_index(i-1) == 0) cycle
      irank = org_nbpe(i)
      do while (jrank < irank)
        j = j + 1
        if (j > new_nnb) exit
        jrank = new_nbpe(j)
      enddo
      if (jrank /= irank) stop 'ERROR: merging comm table: org into new'
      i0 = org_index(i-1)
      len = org_index(i) - i0
      j0 = new_index(j-1)
      new_item(j0+1:j0+len) = org_item(i0+1:i0+len)
      cnt(jrank+1) = len
    enddo
    j = 1
    jrank = new_nbpe(j)
    do i = 1, npe
      if (n_add(i) == 0) cycle
      irank = i-1
      do while (jrank < irank)
        j = j + 1
        jrank = new_nbpe(j)
      enddo
      if (jrank /= irank) stop 'ERROR: merging comm table: add into new'
      i0 = add_index(i-1)
      len = add_index(i) - i0
      j0 = new_index(j-1) + cnt(jrank+1)
      new_item(j0+1:j0+len) = add_item(i0+1:i0+len)
      cnt(jrank+1) = cnt(jrank+1) + len
      if (cnt(jrank+1) /= new_index(j)-new_index(j-1)) stop 'ERROR: merging comm table'
    enddo
    deallocate(cnt)
  end subroutine merge_comm_table

  subroutine copy_vals_to_BT_int(nnb, imp_rows_index, imp_cols_index, &
       imp_rows_item, map, ndof2, imp_vals_item, BT_int)
    implicit none
    integer(kind=kint), intent(in) :: nnb
    integer(kind=kint), allocatable, intent(in) :: imp_rows_index(:), imp_cols_index(:)
    integer(kind=kint), intent(in) :: imp_rows_item(:,:), map(:)
    integer(kind=kint), intent(in) :: ndof2
    real(kind=kreal), intent(in) :: imp_vals_item(:)
    type (hecmwST_local_matrix), intent(inout) :: BT_int
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: idom, is, ie, ic0, i, irow, ncol, j0, j
    allocate(cnt(BT_int%nr))
    cnt(:) = 0
    do idom = 1, nnb
      is = imp_rows_index(idom-1)+1
      ie = imp_rows_index(idom)
      ic0 = imp_cols_index(idom-1)
      do i = is, ie
        irow = imp_rows_item(1,i)
        ncol = imp_rows_item(2,i)
        if (irow < 1 .or. BT_int%nr < irow) stop 'ERROR: copy vals to BT_int: irow'
        j0 = BT_int%index(irow-1) + cnt(irow)
        do j = 1, ncol
          BT_int%item(j0+j) = map(ic0+j)
          BT_int%A(ndof2*(j0+j-1)+1:ndof2*(j0+j)) = imp_vals_item(ndof2*(ic0+j-1)+1:ndof2*(ic0+j))
        enddo
        cnt(irow) = cnt(irow) + ncol
        ic0 = ic0 + ncol
      enddo
      if (ic0 /= imp_cols_index(idom)) stop 'ERROR: copy vals to BT_int: ic0'
    enddo
    deallocate(cnt)
  end subroutine copy_vals_to_BT_int

  subroutine sort_and_uniq_rows(BTmat)
    use hecmw_array_util
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTmat
    integer(kind=kint) :: nr, ndof, ndof2
    integer(kind=kint) :: irow, is, ie, is_new, ie_new, i, i_new
    integer(kind=kint) :: ndup, ndup_tot
    integer(kind=kint) :: js, je, js_new, je_new
    integer(kind=kint) :: new_nnz
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint), pointer :: sort_item(:), new_index(:), new_item(:)
    real(kind=kreal), pointer :: new_A(:)
    logical :: sorted
    real(kind=kreal) :: t0, t1
    t0 = hecmw_wtime()
    nr = BTmat%nr
    ! check if already sorted
    sorted = .true.
    OUTER: do irow = 1, nr
      is = BTmat%index(irow-1)+1
      ie = BTmat%index(irow)
      do i = is, ie-1
        if (BTmat%item(i) >= BTmat%item(i+1)) then
          sorted = .false.
          exit OUTER
        endif
      enddo
    end do OUTER
    t1 = hecmw_wtime()
    if (TIMER >= 4) write(0, '(A,f10.4,L2)') "####### sort_and_uniq_rows (1) : ",t1-t0,sorted
    t0 = hecmw_wtime()
    if (sorted) return
    ! perform sort
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    ! duplicate item array (sort_item)
    allocate(sort_item(BTmat%nnz))
    do i = 1, BTmat%nnz
      sort_item(i) = BTmat%item(i)
    enddo
    ! sort and uniq item for each row
    allocate(cnt(nr))
    ndup_tot = 0
    !$omp parallel do default(none), &
      !$omp& schedule(dynamic,1), &
      !$omp& private(irow,is,ie,ndup), &
      !$omp& shared(nr,BTmat,sort_item,cnt), &
      !$omp& reduction(+:ndup_tot)
    do irow = 1, nr
      is = BTmat%index(irow-1)+1
      ie = BTmat%index(irow)
      call hecmw_qsort_int_array(sort_item, is, ie)
      call hecmw_uniq_int_array(sort_item, is, ie, ndup)
      cnt(irow) = (ie-is+1) - ndup
      ndup_tot = ndup_tot + ndup
    enddo
    !$omp end parallel do
    t1 = hecmw_wtime()
    if (TIMER >= 4) write(0, '(A,f10.4,I5)') "####### sort_and_uniq_rows (2) : ",t1-t0,ndup_tot
    t0 = hecmw_wtime()
    ! make new index and item array (new_index, new_item)
    if (ndup_tot == 0) then
      new_index => BTmat%index
      new_nnz = BTmat%nnz
      new_item => sort_item
    else
      allocate(new_index(0:nr))
      call make_index(nr, cnt, new_index)
      new_nnz = new_index(nr)
      allocate(new_item(new_nnz))
      do irow = 1, nr
        is = BTmat%index(irow-1)+1
        ie = is+cnt(irow)-1
        is_new = new_index(irow-1)+1
        ie_new = is_new+cnt(irow)-1
        new_item(is_new:ie_new) = sort_item(is:ie)
      enddo
      deallocate(sort_item)
    endif
    deallocate(cnt)
    t1 = hecmw_wtime()
    if (TIMER >= 4) write(0, '(A,f10.4)') "####### sort_and_uniq_rows (3) : ",t1-t0
    t0 = hecmw_wtime()
    ! allocate and clear value array (new_A)
    allocate(new_A(ndof2*new_nnz))
    new_A(:) = 0.d0
    ! copy/add value from old A to new A
    !$omp parallel do default(none), &
      !$omp& schedule(dynamic,1), &
      !$omp& private(irow,is,ie,is_new,ie_new,i,i_new,js,je,js_new,je_new), &
      !$omp& shared(nr,BTmat,new_index,new_item,ndof2,new_A)
    do irow = 1, nr
      is = BTmat%index(irow-1)+1
      ie = BTmat%index(irow)
      is_new = new_index(irow-1)+1
      ie_new = new_index(irow)
      ! for each item in row
      do i = is, ie
        ! find place in new item
        call hecmw_bsearch_int_array(new_item, is_new, ie_new, BTmat%item(i), i_new)
        if (i_new == -1) stop 'ERROR: sort_and_uniq_rows'
        js = ndof2*(i-1)+1
        je = ndof2*i
        js_new = ndof2*(i_new-1)+1
        je_new = ndof2*i_new
        new_A(js_new:je_new) = new_A(js_new:je_new) + BTmat%A(js:je)
      enddo
    enddo
    !$omp end parallel do
    t1 = hecmw_wtime()
    if (TIMER >= 4) write(0, '(A,f10.4)') "####### sort_and_uniq_rows (4) : ",t1-t0
    t0 = hecmw_wtime()
    ! deallocate/update nnz, index, item, A
    if (ndup_tot == 0) then
      deallocate(BTmat%item)
      BTmat%item => new_item
      deallocate(BTmat%A)
      BTmat%A => new_A
    else
      BTmat%nnz = new_nnz
      deallocate(BTmat%index)
      BTmat%index => new_index
      deallocate(BTmat%item)
      BTmat%item => new_item
      deallocate(BTmat%A)
      BTmat%A => new_A
    endif
  end subroutine sort_and_uniq_rows

  subroutine hecmw_localmat_add(Amat, Bmat, Cmat)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Amat
    type (hecmwST_local_matrix), intent(in) :: Bmat
    type (hecmwST_local_matrix), intent(out) :: Cmat
    integer(kind=kint) :: ndof, ndof2, nr, nc, i, icnt, js, je, j, jcol, idx, i0, k
    integer(kind=kint), allocatable :: iw(:)
    if (Amat%ndof /= Bmat%ndof) stop 'ERROR: hecmw_localmat_add: non-matching ndof'
    ndof = Amat%ndof
    ndof2 = ndof*ndof
    nr = min(Amat%nr, Bmat%nr)
    nc = max(Amat%nc, Bmat%nc)
    Cmat%ndof = ndof
    Cmat%nr = nr
    Cmat%nc = nc
    Cmat%nnz = 0
    allocate(Cmat%index(0:nr))
    Cmat%index(0) = 0
    allocate(iw(nc))
    do i = 1, nr
      icnt = 0
      ! Amat
      js = Amat%index(i-1)+1
      je = Amat%index(i)
      do j = js, je
        jcol = Amat%item(j)
        icnt = icnt + 1
        iw(icnt) = jcol
      enddo
      ! Bmat
      js = Bmat%index(i-1)+1
      je = Bmat%index(i)
      lj1: do j = js, je
        jcol = Bmat%item(j)
        do k = 1, icnt
          if (iw(k) == jcol) cycle lj1
        enddo
        icnt = icnt + 1
        iw(icnt) = jcol
      enddo lj1
      Cmat%index(i) = Cmat%index(i-1) + icnt
    enddo
    Cmat%nnz = Cmat%index(nr)
    allocate(Cmat%item(Cmat%nnz))
    allocate(Cmat%A(ndof2*Cmat%nnz))
    do i = 1, nr
      i0 = Cmat%index(i-1)
      icnt = 0
      ! Amat
      js = Amat%index(i-1)+1
      je = Amat%index(i)
      do j = js, je
        jcol = Amat%item(j)
        icnt = icnt + 1
        idx = i0 + icnt
        Cmat%item(idx) = jcol
        Cmat%A(ndof2*(idx-1)+1:ndof2*idx) = Amat%A(ndof2*(j-1)+1:ndof2*j)
      enddo
      ! Bmat
      js = Bmat%index(i-1)+1
      je = Bmat%index(i)
      lj2: do j = js, je
        jcol = Bmat%item(j)
        do k = 1, icnt
          idx = i0 + k
          if (Cmat%item(idx) == jcol) then
            Cmat%A(ndof2*(idx-1)+1:ndof2*idx) = &
                 Cmat%A(ndof2*(idx-1)+1:ndof2*idx) + Bmat%A(ndof2*(j-1)+1:ndof2*j)
            cycle lj2
          endif
        enddo
        icnt = icnt + 1
        idx = i0 + icnt
        Cmat%item(idx) = jcol
        Cmat%A(ndof2*(idx-1)+1:ndof2*idx) = Bmat%A(ndof2*(j-1)+1:ndof2*j)
      enddo lj2
      if (i0 + icnt /= Cmat%index(i)) stop 'ERROR: merge localmat'
    enddo
    call sort_and_uniq_rows(Cmat)
  end subroutine hecmw_localmat_add

  ! subroutine hecmw_localmat_add(Amat, Bmat, Cmat)
  !   implicit none
  !   type (hecmwST_local_matrix), intent(in) :: Amat
  !   type (hecmwST_local_matrix), intent(in) :: Bmat
  !   type (hecmwST_local_matrix), intent(out) :: Cmat
  !   integer(kind=kint) :: ndof, ndof2, nr, nc, i, js, je, j, jcol, nnz_row, idx, ks, ke, k, kcol
  !   if (Amat%ndof /= Bmat%ndof) stop 'ERROR: hecmw_localmat_add: non-matching ndof'
  !   ndof = Amat%ndof
  !   ndof2 = ndof*ndof
  !   nr = min(Amat%nr, Bmat%nr)
  !   nc = max(Amat%nc, Bmat%nc)
  !   Cmat%ndof = ndof
  !   Cmat%nr = nr
  !   Cmat%nc = nc
  !   Cmat%nnz = Amat%index(nr) + Bmat%index(nr)
  !   allocate(Cmat%index(0:nr))
  !   allocate(Cmat%item(Cmat%nnz))
  !   allocate(Cmat%A(ndof2 * Cmat%nnz))
  !   Cmat%index(0) = 0
  !   idx = 0
  !   do i = 1, nr
  !     ! Amat
  !     js = Amat%index(i-1)+1
  !     je = Amat%index(i)
  !     do j = js, je
  !       idx = idx + 1
  !       Cmat%item(idx) = Amat%item(j)
  !       Cmat%A(ndof2*(idx-1)+1:ndof2*idx) = Amat%A(ndof2*(j-1)+1:ndof2*j)
  !     enddo
  !     ! Bmat
  !     js = Bmat%index(i-1)+1
  !     je = Bmat%index(i)
  !     do j = js, je
  !       idx = idx + 1
  !       Cmat%item(idx) = Bmat%item(j)
  !       Cmat%A(ndof2*(idx-1)+1:ndof2*idx) = Bmat%A(ndof2*(j-1)+1:ndof2*j)
  !     enddo
  !     Cmat%index(i) = idx
  !   enddo
  !   if (Cmat%index(nr) /= Cmat%nnz) stop 'ERROR: merge localmat'
  !   call sort_and_uniq_rows(Cmat)
  ! end subroutine hecmw_localmat_add

  subroutine hecmw_localmat_init_with_hecmat(BKmat, hecMAT, num_lagrange)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BKmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    integer(kind=kint), optional, intent(in) :: num_lagrange
    integer(kind=kint) :: ndof, ndof2, i, idx, idx2, js, je, j, k
    integer(kind=kint), allocatable :: incl_nz(:), cnt(:)
    logical :: check_nonzero
    check_nonzero = .false.
    !check_nonzero = .true.  !!! always checking nonzero seems to be faster
    !
    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    ! nr, nc, nnz
    BKmat%nr = hecMAT%NP
    BKmat%nc = hecMAT%NP
    BKmat%ndof = ndof
    !
    if (present(num_lagrange)) then       !!! TEMPORARY (DUE TO WRONG conMAT WHEN num_lagrange==0) !!!
      check_nonzero = .true.
    endif
    !
    if (check_nonzero) then
      allocate(incl_nz(hecMAT%NPL + hecMAT%NPU + hecMAT%NP))
      allocate(cnt(BKmat%nr))
      incl_nz(:) = 0
      !$omp parallel default(none), &
        !$omp& private(i,idx,js,je,j,k), &
        !$omp& shared(BKmat,hecMAT,cnt,ndof2,incl_nz)
      !$omp do
      do i = 1, BKmat%nr
        idx = hecMAT%indexL(i-1) + (i-1) + hecMAT%indexU(i-1)
        cnt(i) = 0
        ! lower
        js = hecMAT%indexL(i-1)+1
        je = hecMAT%indexL(i)
        do j = js, je
          idx = idx + 1
          do k = 1, ndof2
            if (hecMAT%AL(ndof2*(j-1)+k) /= 0.0d0) then
              incl_nz(idx) = 1
              cnt(i) = cnt(i) + 1
              exit
            endif
          enddo
        enddo
        ! diag
        idx = idx + 1
        do k = 1, ndof2
          if (hecMAT%D(ndof2*(i-1)+k) /= 0.0d0) then
            incl_nz(idx) = 1
            cnt(i) = cnt(i) + 1
            exit
          endif
        enddo
        ! upper
        js = hecMAT%indexU(i-1)+1
        je = hecMAT%indexU(i)
        do j = js, je
          idx = idx + 1
          do k = 1, ndof2
            if (hecMAT%AU(ndof2*(j-1)+k) /= 0.0d0) then
              incl_nz(idx) = 1
              cnt(i) = cnt(i) + 1
              exit
            endif
          enddo
        enddo
        if (idx /= hecMAT%indexL(i) + i + hecMAT%indexU(i)) stop 'ERROR: hecmw_localmat_init_with_hecmat: count'
      enddo
      !$omp end do
      !$omp end parallel
      ! index
      allocate(BKmat%index(0:BKmat%nr))
      call make_index(BKmat%nr, cnt, BKmat%index)
      deallocate(cnt)
      BKmat%nnz = BKmat%index(BKmat%nr)
      ! item, A
      allocate(BKmat%item(BKmat%nnz))
      allocate(BKmat%A(ndof2 * BKmat%nnz))
      !$omp parallel default(none), &
        !$omp& private(i,idx,idx2,js,je,j), &
        !$omp& shared(BKmat,hecMAT,ndof2,incl_nz)
      !$omp do
      do i = 1, BKmat%nr
        idx = hecMAT%indexL(i-1) + (i-1) + hecMAT%indexU(i-1)
        idx2 = BKmat%index(i-1)
        ! lower
        js = hecMAT%indexL(i-1)+1
        je = hecMAT%indexL(i)
        do j = js, je
          idx = idx + 1
          if (incl_nz(idx) == 1) then
            idx2 = idx2 + 1
            BKmat%item(idx2) = hecMAT%itemL(j)
            BKmat%A(ndof2*(idx2-1)+1:ndof2*idx2) = hecMAT%AL(ndof2*(j-1)+1:ndof2*j)
          endif
        enddo
        ! diag
        idx = idx + 1
        if (incl_nz(idx) == 1) then
          idx2 = idx2 + 1
          BKmat%item(idx2) = i
          BKmat%A(ndof2*(idx2-1)+1:ndof2*idx2) = hecMAT%D(ndof2*(i-1)+1:ndof2*i)
        endif
        ! upper
        js = hecMAT%indexU(i-1)+1
        je = hecMAT%indexU(i)
        do j = js, je
          idx = idx + 1
          if (incl_nz(idx) == 1) then
            idx2 = idx2 + 1
            BKmat%item(idx2) = hecMAT%itemU(j)
            BKmat%A(ndof2*(idx2-1)+1:ndof2*idx2) = hecMAT%AU(ndof2*(j-1)+1:ndof2*j)
          endif
        enddo
        if (idx /= hecMAT%indexL(i) + i + hecMAT%indexU(i)) stop 'ERROR: hecmw_localmat_init_with_hecmat: copy'
        if (idx2 /= BKmat%index(i)) stop 'ERROR: hecmw_localmat_init_with_hecmat: index'
      enddo
      !$omp end do
      !$omp end parallel
      deallocate(incl_nz)
    else
      BKmat%nnz = hecMAT%NPL + hecMAT%NP + hecMAT%NPU
      allocate(BKmat%index(0:BKmat%nr))
      allocate(BKmat%item(BKmat%nnz))
      allocate(BKmat%A(ndof2 * BKmat%nnz))
      BKmat%index(0) = 0
      !$omp parallel do default(none), &
        !$omp& private(i,idx,js,je,j), &
        !$omp& shared(BKmat,hecMAT,ndof2)
      do i = 1, BKmat%nr
        idx = hecMAT%indexL(i-1) + (i-1) + hecMAT%indexU(i-1)
        ! lower
        js = hecMAT%indexL(i-1)+1
        je = hecMAT%indexL(i)
        do j = js, je
          idx = idx + 1
          BKmat%item(idx) = hecMAT%itemL(j)
          BKmat%A(ndof2*(idx-1)+1:ndof2*idx) = hecMAT%AL(ndof2*(j-1)+1:ndof2*j)
        enddo
        ! diag
        idx = idx + 1
        BKmat%item(idx) = i
        BKmat%A(ndof2*(idx-1)+1:ndof2*idx) = hecMAT%D(ndof2*(i-1)+1:ndof2*i)
        ! upper
        js = hecMAT%indexU(i-1)+1
        je = hecMAT%indexU(i)
        do j = js, je
          idx = idx + 1
          BKmat%item(idx) = hecMAT%itemU(j)
          BKmat%A(ndof2*(idx-1)+1:ndof2*idx) = hecMAT%AU(ndof2*(j-1)+1:ndof2*j)
        enddo
        BKmat%index(i) = idx
        if (idx /= hecMAT%indexL(i) + i + hecMAT%indexU(i)) stop 'ERROR: hecmw_localmat_init_with_hecmat: copy'
      enddo
      !$omp end parallel do
    endif
  end subroutine hecmw_localmat_init_with_hecmat

  subroutine hecmw_localmat_add_hecmat(BKmat, hecMAT)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BKmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_matrix) :: W1mat, W2mat
    !! Should Be Simple If Non-Zero Profile Is Kept !!
    call hecmw_localmat_init_with_hecmat(W1mat, hecMAT)
    call debug_write_matrix(W1mat, 'BKmat (hecMAT)', DEBUG_MATRIX)
    call hecmw_localmat_add(BKmat, W1mat, W2mat)
    call hecmw_localmat_free(BKmat)
    call hecmw_localmat_free(W1mat)
    BKmat%nr = W2mat%nr
    BKmat%nc = W2mat%nc
    BKmat%nnz = W2mat%nnz
    BKmat%ndof = W2mat%ndof
    BKmat%index => W2mat%index
    BKmat%item => W2mat%item
    BKmat%A => W2mat%A
  end subroutine hecmw_localmat_add_hecmat

  subroutine hecmw_localmat_multmat(BKmat, BTmat, hecMESH, BKTmat)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BKmat
    type (hecmwST_local_matrix), intent(inout) :: BTmat
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_local_matrix), intent(out) :: BKTmat
    type (hecmwST_matrix_comm) :: hecCOMM
    type (hecmwST_local_mesh) :: hecMESHnew
    type (hecmwST_local_matrix), allocatable :: BT_exp(:)
    type (hecmwST_local_matrix) :: BT_imp, BT_all
    integer(kind=kint), allocatable :: exp_cols_index(:)
    integer(kind=kint), allocatable :: exp_cols_item(:,:)
    real(kind=kreal) :: t0, t1
    t0 = hecmw_wtime()
    !
    if (hecMESH%PETOT > 1) then
      call make_comm_table(BKmat, hecMESH, hecCOMM)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: make_comm_table done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (1) : ',t1-t0
      t0 = hecmw_wtime()
      !
      if (BTmat%nr > hecMESH%nn_internal) then
        ! consider only internal part of BTmat
        if (DEBUG >= 1) write(0,'(A)') 'DEBUG: hecmw_localmat_multmat: ignore external part of BTmat'
        BTmat%nr = hecMESH%nn_internal
        BTmat%nnz = BTmat%index(BTmat%nr)
      endif
      !
      call extract_BT_exp(BTmat, hecCOMM, BT_exp)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: extract_BT_exp done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (2) : ',t1-t0
      t0 = hecmw_wtime()
      !
      call prepare_column_info(hecMESH, BT_exp, exp_cols_index, exp_cols_item)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: prepare column info done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (3) : ',t1-t0
      t0 = hecmw_wtime()
      !
      call send_BT_exp_and_recv_BT_imp(hecMESH, hecCOMM, BT_exp, exp_cols_index, exp_cols_item, BT_imp, hecMESHnew)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: send BT_exp and recv BT_imp done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (4) : ',t1-t0
      t0 = hecmw_wtime()
      call free_comm_table(hecCOMM)
      !
      call concat_BTmat_and_BT_imp(BTmat, BT_imp, BT_all)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: concat BTmat and BT_imp into BT_all done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (5) : ',t1-t0
      t0 = hecmw_wtime()
      call hecmw_localmat_free(BT_imp)
      !
      call multiply_mat_mat(BKmat, BT_all, BKTmat)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: multiply BKmat and BT_all into BKTmat done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (6) : ',t1-t0
      t0 = hecmw_wtime()
      call hecmw_localmat_free(BT_all)
      !
      if (hecMESH%n_neighbor_pe > 0) then
        hecMESH%n_node = hecMESHnew%n_node
        hecMESH%n_neighbor_pe = hecMESHnew%n_neighbor_pe
        deallocate(hecMESH%neighbor_pe)
        deallocate(hecMESH%import_index)
        deallocate(hecMESH%export_index)
        deallocate(hecMESH%import_item)
        deallocate(hecMESH%export_item)
        deallocate(hecMESH%node_ID)
        deallocate(hecMESH%global_node_ID)
        hecMESH%neighbor_pe => hecMESHnew%neighbor_pe
        hecMESH%import_index => hecMESHnew%import_index
        hecMESH%export_index => hecMESHnew%export_index
        hecMESH%import_item => hecMESHnew%import_item
        hecMESH%export_item => hecMESHnew%export_item
        hecMESH%node_ID => hecMESHnew%node_ID
        hecMESH%global_node_ID => hecMESHnew%global_node_ID
        if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: update hecMESH done'
        t1 = hecmw_wtime()
        if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat (7) : ',t1-t0
      endif
    else
      call multiply_mat_mat(BKmat, BTmat, BKTmat)
      if (DEBUG >= 1) write(0,*) 'DEBUG: hecmw_localmat_multmat: multiply BKmat and BTmat into BKTmat done'
      t1 = hecmw_wtime()
      if (TIMER >= 2) write(0,'(A,f10.4)') '##### hecmw_localmat_multmat : ',t1-t0
    endif
  end subroutine hecmw_localmat_multmat

  subroutine make_comm_table(BKmat, hecMESH, hecCOMM)
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BKmat
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix_comm), intent(out) :: hecCOMM
    integer(kind=kint) :: nn_int, nn_ext, nnb, i, icol, irank, idom, idx, n_send, tag, js, je, len
    integer(kind=kint), allocatable :: is_nz_col(:), imp_cnt(:), exp_cnt(:), import_item_remote(:)
    integer(kind=kint), allocatable :: requests(:), statuses(:,:)
    hecCOMM%zero = hecMESH%zero
    hecCOMM%HECMW_COMM = hecMESH%MPI_COMM
    hecCOMM%PETOT = hecMESH%PETOT
    hecCOMM%PEsmpTOT = hecMESH%PEsmpTOT
    hecCOMM%my_rank = hecMESH%my_rank
    hecCOMM%errnof = hecMESH%errnof
    hecCOMM%n_subdomain = hecMESH%n_subdomain
    hecCOMM%n_neighbor_pe = hecMESH%n_neighbor_pe
    allocate(hecCOMM%neighbor_pe(hecCOMM%n_neighbor_pe))
    hecCOMM%neighbor_pe(:) = hecMESH%neighbor_pe(:)
    !
    nn_int = hecMESH%nn_internal
    nn_ext = hecMESH%n_node - hecMESH%nn_internal
    nnb = hecCOMM%n_neighbor_pe
    !
    ! check_external_nz_cols (by profile (not number))
    allocate(is_nz_col(nn_ext))
    is_nz_col(:) = 0
    do i = 1, BKmat%index(nn_int)
      icol = BKmat%item(i)
      if (icol > nn_int) is_nz_col(icol - nn_int) = 1
    enddo
    !
    ! count_nz_cols_per_rank
    allocate(imp_cnt(nnb))
    imp_cnt(:) = 0
    do i = 1, nn_ext
      if (is_nz_col(i) == 1) then
        irank = hecMESH%node_ID(2*(nn_int+i))
        call rank_to_idom(hecMESH, irank, idom)
        imp_cnt(idom) = imp_cnt(idom) + 1
      endif
    enddo
    if (DEBUG >= 3) write(0,*) '    DEBUG3: imp_cnt',imp_cnt(:)
    !
    ! make_index
    allocate(hecCOMM%import_index(0:nnb))
    call make_index(nnb, imp_cnt, hecCOMM%import_index)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: import_index',hecCOMM%import_index(:)
    !
    ! fill item
    allocate(hecCOMM%import_item(hecCOMM%import_index(nnb)))
    imp_cnt(:) = 0
    do i = 1, nn_ext
      if (is_nz_col(i) == 1) then
        irank = hecMESH%node_ID(2*(nn_int+i))
        call rank_to_idom(hecMESH, irank, idom)
        imp_cnt(idom) = imp_cnt(idom) + 1
        idx = hecCOMM%import_index(idom-1)+imp_cnt(idom)
        hecCOMM%import_item(idx) = nn_int+i
      endif
    enddo
    if (DEBUG >= 3) write(0,*) '    DEBUG3: import_item',hecCOMM%import_item(:)
    !
    allocate(import_item_remote(hecCOMM%import_index(nnb)))
    do i = 1, hecCOMM%import_index(nnb)
      import_item_remote(i) = hecMESH%node_ID(2*hecCOMM%import_item(i)-1)
    enddo
    if (DEBUG >= 3) write(0,*) '    DEBUG3: import_item_remote',import_item_remote(:)
    !
    allocate(requests(2*nnb))
    allocate(statuses(HECMW_STATUS_SIZE, 2*nnb))
    !
    ! send/recv
    n_send = 0
    do idom = 1, nnb
      irank = hecCOMM%neighbor_pe(idom)
      n_send = n_send + 1
      tag = 6001
      call HECMW_ISEND_INT(imp_cnt(idom), 1, irank, tag, hecCOMM%HECMW_COMM, requests(n_send))
      if (imp_cnt(idom) > 0) then
        js = hecCOMM%import_index(idom-1)+1
        je = hecCOMM%import_index(idom)
        len = je-js+1
        n_send = n_send + 1
        tag = 6002
        call HECMW_ISEND_INT(import_item_remote(js:je), len, irank, tag, &
             hecCOMM%HECMW_COMM, requests(n_send))
      endif
    enddo
    !
    ! index
    allocate(exp_cnt(nnb))
    do idom = 1, nnb
      irank = hecCOMM%neighbor_pe(idom)
      tag = 6001
      call HECMW_RECV_INT(exp_cnt(idom), 1, irank, tag, hecCOMM%HECMW_COMM, statuses(:,1))
    enddo
    allocate(hecCOMM%export_index(0:nnb))
    call make_index(nnb, exp_cnt, hecCOMM%export_index)
    if (DEBUG >= 3) write(0,*) '    DEBUG3: export_index',hecCOMM%export_index(:)
    !
    ! item
    allocate(hecCOMM%export_item(hecCOMM%export_index(nnb)))
    do idom = 1, nnb
      if (exp_cnt(idom) <= 0) cycle
      irank = hecCOMM%neighbor_pe(idom)
      js = hecCOMM%export_index(idom-1)+1
      je = hecCOMM%export_index(idom)
      len = je-js+1
      tag = 6002
      call HECMW_RECV_INT(hecCOMM%export_item(js:je), len, irank, tag, &
           hecCOMM%HECMW_COMM, statuses(:,1))
    enddo
    if (DEBUG >= 3) write(0,*) '    DEBUG3: export_item',hecCOMM%export_item(:)
    call HECMW_Waitall(n_send, requests, statuses)
    !
    deallocate(imp_cnt)
    deallocate(exp_cnt)
    deallocate(import_item_remote)
  end subroutine make_comm_table

  subroutine free_comm_table(hecCOMM)
    implicit none
    type (hecmwST_matrix_comm), intent(inout) :: hecCOMM
    deallocate(hecCOMM%neighbor_pe)
    deallocate(hecCOMM%import_index)
    deallocate(hecCOMM%import_item)
    deallocate(hecCOMM%export_index)
    deallocate(hecCOMM%export_item)
  end subroutine free_comm_table

  subroutine extract_BT_exp(BTmat, hecCOMM, BT_exp)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    type (hecmwST_matrix_comm), intent(in) :: hecCOMM
    type (hecmwST_local_matrix), allocatable, intent(out) :: BT_exp(:)
    integer(kind=kint) :: ndof, ndof2, idom, idx_0, idx_n, j, jrow, nnz_row, idx, ks, ke, k
    if (hecCOMM%n_neighbor_pe == 0) return
    allocate(BT_exp(hecCOMM%n_neighbor_pe))
    ndof = BTmat%ndof
    ndof2 = ndof * ndof
    do idom = 1, hecCOMM%n_neighbor_pe
      idx_0 = hecCOMM%export_index(idom-1)
      idx_n = hecCOMM%export_index(idom)
      BT_exp(idom)%nr = idx_n - idx_0
      BT_exp(idom)%nc = BTmat%nc
      BT_exp(idom)%nnz = 0
      BT_exp(idom)%ndof = ndof
      allocate(BT_exp(idom)%index(0:BT_exp(idom)%nr))
      BT_exp(idom)%index(0) = 0
      do j = 1, BT_exp(idom)%nr
        jrow = hecCOMM%export_item(idx_0 + j)
        nnz_row = BTmat%index(jrow) - BTmat%index(jrow-1)
        BT_exp(idom)%index(j) = BT_exp(idom)%index(j-1) + nnz_row
      enddo
      BT_exp(idom)%nnz = BT_exp(idom)%index(BT_exp(idom)%nr)
      allocate(BT_exp(idom)%item(BT_exp(idom)%nnz))
      allocate(BT_exp(idom)%A(ndof2 * BT_exp(idom)%nnz))
      idx = 0
      do j = 1, BT_exp(idom)%nr
        jrow = hecCOMM%export_item(idx_0 + j)
        ks = BTmat%index(jrow-1) + 1
        ke = BTmat%index(jrow)
        do k = ks, ke
          idx = idx + 1
          BT_exp(idom)%item(idx) = BTmat%item(k)
          BT_exp(idom)%A(ndof2*(idx-1)+1:ndof2*idx) = BTmat%A(ndof2*(k-1)+1:ndof2*k)
        enddo
        if (idx /= BT_exp(idom)%index(j)) stop 'ERROR: extract BT_exp'
      enddo
    enddo
  end subroutine extract_BT_exp

  subroutine send_BT_exp_and_recv_BT_imp(hecMESH, hecCOMM, BT_exp, exp_cols_index, exp_cols_item, BT_imp, hecMESHnew)
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix_comm), intent(in) :: hecCOMM
    type (hecmwST_local_matrix), allocatable, intent(inout) :: BT_exp(:)
    integer(kind=kint), allocatable, intent(inout) :: exp_cols_index(:)
    integer(kind=kint), allocatable, intent(inout) :: exp_cols_item(:,:)
    type (hecmwST_local_matrix), intent(out) :: BT_imp
    type (hecmwST_local_mesh), intent(inout) :: hecMESHnew
    integer(kind=kint), allocatable :: nnz_imp(:), cnt(:), index_imp(:)
    integer(kind=kint), allocatable :: imp_cols_index(:)
    integer(kind=kint), allocatable :: imp_cols_item(:,:)
    real(kind=kreal), allocatable :: imp_vals_item(:)
    integer(kind=kint) :: nnb, ndof, ndof2, idom, irank, nr, n_send, tag, idx_0, idx_n, j, jj, nnz
    integer(kind=kint), allocatable :: requests(:)
    integer(kind=kint), allocatable :: statuses(:,:)
    integer(kind=kint), allocatable :: map(:), add_nodes(:,:)
    integer(kind=kint) :: n_add_node, i0
    nnb = hecCOMM%n_neighbor_pe
    if (nnb == 0) then
      BT_imp%nr = 0
      BT_imp%nc = 0
      BT_imp%nnz = 0
      BT_imp%ndof = 0
      allocate(BT_imp%index(0:0))
      BT_imp%index(0) = 0
      return
    endif
    ndof = BT_exp(1)%ndof
    ndof2 = ndof*ndof
    allocate(requests(nnb*3))
    allocate(statuses(HECMW_STATUS_SIZE, nnb*3))
    n_send = 0
    do idom = 1, nnb
      irank = hecCOMM%neighbor_pe(idom)
      nr = BT_exp(idom)%nr
      if (nr == 0) cycle
      n_send = n_send + 1
      tag = 3001
      call HECMW_ISEND_INT(BT_exp(idom)%index(0:BT_exp(idom)%nr), BT_exp(idom)%nr + 1, &
           irank, tag, hecCOMM%HECMW_COMM, requests(n_send))
      if (BT_exp(idom)%nnz == 0) cycle
      n_send = n_send + 1
      tag = 3002
      call HECMW_ISEND_INT(exp_cols_item(1,exp_cols_index(idom-1)+1), &
           cNCOL_ITEM * BT_exp(idom)%nnz, irank, tag, hecCOMM%HECMW_COMM, requests(n_send))
      n_send = n_send + 1
      tag = 3003
      call HECMW_ISEND_R(BT_exp(idom)%A, ndof2 * BT_exp(idom)%nnz, &
           irank, tag, hecCOMM%HECMW_COMM, requests(n_send))
    enddo
    !
    ! BT_imp%nr = hecCOMM%import_index(nnb)
    BT_imp%nr = hecMESH%n_node - hecMESH%nn_internal
    BT_imp%nc = 0  !!! TEMPORARY
    BT_imp%nnz = 0
    BT_imp%ndof = ndof
    !
    allocate(nnz_imp(nnb))
    allocate(cnt(BT_imp%nr))
    !
    cnt(:) = 0
    do idom = 1, nnb
      irank = hecCOMM%neighbor_pe(idom)
      idx_0 = hecCOMM%import_index(idom-1)
      idx_n = hecCOMM%import_index(idom)
      nr = idx_n - idx_0
      if (nr == 0) then
        nnz_imp(idom) = 0
        cycle
      endif
      allocate(index_imp(0:nr))
      tag = 3001
      call HECMW_RECV_INT(index_imp(0:nr), nr+1, irank, tag, &
           hecCOMM%HECMW_COMM, statuses(:,1))
      nnz_imp(idom) = index_imp(nr)
      do j = 1, nr
        jj = hecCOMM%import_item(idx_0 + j) - hecMESH%nn_internal
        if (jj < 1 .or. BT_imp%nr < jj) stop 'ERROR: jj out of range'
        if (cnt(jj) /= 0) stop 'ERROR: duplicate import rows?'
        cnt(jj) = index_imp(j) - index_imp(j-1)
      enddo
      deallocate(index_imp)
    enddo
    !
    allocate(imp_cols_index(0:nnb))
    call make_index(nnb, nnz_imp, imp_cols_index)
    deallocate(nnz_imp)
    !
    allocate(BT_imp%index(0:BT_imp%nr))
    call make_index(BT_imp%nr, cnt, BT_imp%index)
    deallocate(cnt)
    !
    BT_imp%nnz = BT_imp%index(BT_imp%nr)
    if (BT_imp%nnz /= imp_cols_index(nnb)) &
         stop 'ERROR: total num of nonzero of BT_imp'
    !
    allocate(imp_cols_item(cNCOL_ITEM, BT_imp%nnz))
    allocate(imp_vals_item(ndof2 * BT_imp%nnz))
    !
    do idom = 1, nnb
      irank = hecCOMM%neighbor_pe(idom)
      idx_0 = imp_cols_index(idom-1)
      idx_n = imp_cols_index(idom)
      nnz = idx_n - idx_0
      if (nnz == 0) cycle
      tag = 3002
      call HECMW_RECV_INT(imp_cols_item(1, idx_0 + 1), cNCOL_ITEM * nnz, &
           irank, tag, hecCOMM%HECMW_COMM, statuses(:,1))
      tag = 3003
      call HECMW_RECV_R(imp_vals_item(ndof2*idx_0 + 1), ndof2 * nnz, &
           irank, tag, hecCOMM%HECMW_COMM, statuses(:,1))
    enddo
    call HECMW_Waitall(n_send, requests, statuses)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: send BT_imp and recv into temporary data done'
    !
    deallocate(requests)
    deallocate(statuses)
    !
    do idom = 1, nnb
      call hecmw_localmat_free(BT_exp(idom))
    enddo
    deallocate(BT_exp)
    deallocate(exp_cols_index)
    deallocate(exp_cols_item)
    !
    call copy_mesh(hecMESH, hecMESHnew)
    !
    call map_imported_cols(hecMESHnew, imp_cols_index(nnb), imp_cols_item, n_add_node, add_nodes, map, i0)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: map imported cols done'
    !
    call update_comm_table(hecMESHnew, n_add_node, add_nodes, i0)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: update comm_table done'
    !
    BT_imp%nc = hecMESHnew%n_node
    !
    allocate(BT_imp%item(BT_imp%nnz))
    allocate(BT_imp%A(ndof2 * BT_imp%nnz))
    call copy_vals_to_BT_imp(hecCOMM, hecMESH%nn_internal, imp_cols_index, map, imp_vals_item, BT_imp)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: copy vals to BT_imp done'
    !
    deallocate(imp_cols_index)
    deallocate(imp_cols_item)
    deallocate(imp_vals_item)
    deallocate(map)
  end subroutine send_BT_exp_and_recv_BT_imp

  subroutine copy_vals_to_BT_imp(hecCOMM, nn_internal, imp_cols_index, map, imp_vals_item, BT_imp)
    implicit none
    type (hecmwST_matrix_comm), intent(in) :: hecCOMM
    integer(kind=kint), intent(in) :: nn_internal
    integer(kind=kint), allocatable, intent(in) :: imp_cols_index(:)
    integer(kind=kint), intent(in) :: map(:)
    real(kind=kreal), intent(in) :: imp_vals_item(:)
    type (hecmwST_local_matrix), intent(inout) :: BT_imp
    integer(kind=kint) :: nnb, ndof2, idx, idom, idx_0, idx_n, nr, j, jrow, ks, ke, k
    nnb = hecCOMM%n_neighbor_pe
    ndof2 = BT_imp%ndof ** 2
    idx = 0
    do idom = 1, nnb
      idx_0 = hecCOMM%import_index(idom-1)
      idx_n = hecCOMM%import_index(idom)
      nr = idx_n - idx_0
      if (nr == 0) cycle
      do j = 1, nr
        jrow = hecCOMM%import_item(idx_0 + j) - nn_internal
        ks = BT_imp%index(jrow-1)+1
        ke = BT_imp%index(jrow)
        do k = ks, ke
          idx = idx + 1
          BT_imp%item(k) = map(idx)
          BT_imp%A(ndof2*(k-1)+1:ndof2*k) = imp_vals_item(ndof2*(idx-1)+1:ndof2*idx)
        enddo
      enddo
      if (idx /= imp_cols_index(idom)) stop 'ERROR: copy vals to BT_imp'
    enddo
  end subroutine copy_vals_to_BT_imp

  subroutine concat_BTmat_and_BT_imp(BTmat, BT_imp, BT_all)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    type (hecmwST_local_matrix), intent(in) :: BT_imp
    type (hecmwST_local_matrix), intent(out) :: BT_all
    integer(kind=kint) :: ndof, ndof2, i, ii
    ndof = BTmat%ndof
    if (BT_imp%nr > 0 .and. BT_imp%ndof /= ndof) stop 'ERROR: concat BTmat and BT_imp: ndof'
    ndof2 = ndof*ndof
    BT_all%nr = BTmat%nr + BT_imp%nr
    BT_all%nc = max(BTmat%nc, BT_imp%nc)
    BT_all%nnz = BTmat%nnz + BT_imp%nnz
    BT_all%ndof = ndof
    allocate(BT_all%index(0:BT_all%nr))
    allocate(BT_all%item(BT_all%nnz))
    allocate(BT_all%A(ndof2 * BT_all%nnz))
    BT_all%index(0) = 0
    do i = 1, BTmat%nr
      BT_all%index(i) = BTmat%index(i)
    enddo
    do i = 1, BT_imp%nr
      BT_all%index(BTmat%nr+i) = BT_all%index(BTmat%nr+i-1) + &
           BT_imp%index(i) - BT_imp%index(i-1)
    enddo
    do i = 1, BTmat%nnz
      BT_all%item(i) = BTmat%item(i)
      BT_all%A(ndof2*(i-1)+1:ndof2*i) = BTmat%A(ndof2*(i-1)+1:ndof2*i)
    enddo
    do i = 1, BT_imp%nnz
      ii = BTmat%nnz + i
      BT_all%item(ii) = BT_imp%item(i)
      BT_all%A(ndof2*(ii-1)+1:ndof2*ii) = BT_imp%A(ndof2*(i-1)+1:ndof2*i)
    enddo
  end subroutine concat_BTmat_and_BT_imp

  subroutine multiply_mat_mat(Amat, Bmat, Cmat)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: Amat
    type (hecmwST_local_matrix), intent(in) :: Bmat
    type (hecmwST_local_matrix), intent(out) :: Cmat
    integer(kind=kint) :: ndof, ndof2, nr, nc, nnz, i, icnt
    integer(kind=kint) :: js, je, j, jj, ks, ke, k, kk, l, ll, l0
    integer(kind=kint), allocatable :: iw(:)
    real(kind=kreal), pointer :: Ap(:), Bp(:), Cp(:)
    real(kind=kreal) :: t0, t1
    t0 = hecmw_wtime()
    if (Amat%ndof /= Bmat%ndof) stop 'ERROR: multiply_mat_mat: unmatching ndof'
    ndof = Amat%ndof
    ndof2 = ndof*ndof
    nr = Amat%nr
    nc = Bmat%nc
    if (Amat%nc /= Bmat%nr) then
      write(0,*) 'Amat: nr, nc = ', Amat%nr, Amat%nc
      write(0,*) 'Bmat: nr, nc = ', Bmat%nr, Bmat%nc
      stop 'ERROR: multiply_mat_mat: unmatching size'
    endif
    Cmat%ndof = ndof
    Cmat%nr = nr
    Cmat%nc = nc
    allocate(Cmat%index(0:nr))
    Cmat%index(0) = 0
    !$omp parallel default(none), &
      !$omp& private(iw,i,icnt,js,je,j,jj,ks,ke,k,kk,l), &
      !$omp& shared(nr,nc,Amat,Bmat,Cmat)
    allocate(iw(nc))
    !$omp do
    do i = 1, nr
      icnt = 0
      js = Amat%index(i-1)+1
      je = Amat%index(i)
      do j = js, je
        jj = Amat%item(j)
        ks = Bmat%index(jj-1)+1
        ke = Bmat%index(jj)
        kl1: do k = ks, ke
          kk = Bmat%item(k)
          do l = 1, icnt
            if (iw(l) == kk) cycle kl1
          enddo
          icnt = icnt + 1
          iw(icnt) = kk
        enddo kl1
      enddo
      Cmat%index(i) = icnt
    enddo
    !$omp end do
    deallocate(iw)
    !$omp end parallel
    do i = 1, nr
      Cmat%index(i) = Cmat%index(i-1) + Cmat%index(i)
    enddo
    nnz = Cmat%index(nr)
    Cmat%nnz = nnz
    !write(0,*) 'nnz',nnz
    t1 = hecmw_wtime()
    if (TIMER >= 3) write(0, '(A,f10.4)') "###### multiply_mat_mat (1) : ",t1-t0
    t0 = hecmw_wtime()
    allocate(Cmat%item(nnz))
    allocate(Cmat%A(ndof2 * nnz))
    Cmat%A(:) = 0.0d0
    !$omp parallel default(none), &
      !$omp& private(i,icnt,l0,js,je,j,jj,Ap,ks,ke,k,kk,Bp,ll,l,Cp), &
      !$omp& shared(nr,Cmat,Amat,Bmat,ndof2,ndof)
    !$omp do
    do i = 1, nr
      icnt = 0
      l0 = Cmat%index(i-1)
      ! item
      js = Amat%index(i-1)+1
      je = Amat%index(i)
      do j = js, je
        jj = Amat%item(j)
        Ap => Amat%A(ndof2*(j-1)+1:ndof2*j)
        ks = Bmat%index(jj-1)+1
        ke = Bmat%index(jj)
        do k = ks, ke
          kk = Bmat%item(k)
          Bp => Bmat%A(ndof2*(k-1)+1:ndof2*k)
          ll = -1
          do l = 1, icnt
            if (Cmat%item(l0+l) == kk) then
              ll = l0 + l
              exit
            endif
          enddo
          if (ll < 0) then
            icnt = icnt + 1
            ll = l0 + icnt
            Cmat%item(ll) = kk
          endif
          Cp => Cmat%A(ndof2*(ll-1)+1:ndof2*ll)
          call blk_matmul_add(ndof, Ap, Bp, Cp)
        enddo
      enddo
      !write(0,*) 'l0,icnt,index(i)',Cmat%index(i-1),icnt,Cmat%index(i)
      if (l0+icnt /= Cmat%index(i)) stop 'ERROR: multiply_mat_mat: unknown error'
    enddo
    !$omp end do
    !$omp end parallel
    t1 = hecmw_wtime()
    if (TIMER >= 3) write(0, '(A,f10.4)') "###### multiply_mat_mat (2) : ",t1-t0
    t0 = hecmw_wtime()
    call sort_and_uniq_rows(Cmat)
    t1 = hecmw_wtime()
    if (TIMER >= 3) write(0, '(A,f10.4)') "###### multiply_mat_mat (3) : ",t1-t0
  end subroutine multiply_mat_mat

  subroutine blk_matmul_add(ndof, A, B, AB)
    implicit none
    integer, intent(in) :: ndof
    real(kind=kreal), intent(in) :: A(:), B(:)
    real(kind=kreal), intent(inout) :: AB(:)
    integer :: ndof2, i, j, k, i0, j0, ij, ik, jk
    ndof2=ndof*ndof
    do i=1,ndof
      i0=(i-1)*ndof
      do j=1,ndof
        ij=i0+j
        j0=(j-1)*ndof
        do k=1,ndof
          ik=i0+k
          jk=j0+k
          !$omp atomic
          AB(ik)=AB(ik)+A(ij)*B(jk)
        enddo
      enddo
    enddo
  end subroutine blk_matmul_add

  subroutine hecmw_localmat_make_hecmat(hecMAT, BTtKTmat, hecTKT)
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (hecmwST_local_matrix), intent(in) :: BTtKTmat
    type (hecmwST_matrix), intent(inout) :: hecTKT
    call make_new_hecmat(hecMAT, BTtKTmat, hecTKT)
  end subroutine hecmw_localmat_make_hecmat

  subroutine hecmw_localmat_shrink_comm_table(BKmat, hecMESH)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BKmat
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix_comm) :: hecCOMM
    call make_comm_table(BKmat, hecMESH, hecCOMM)
    deallocate(hecMESH%import_index)
    deallocate(hecMESH%import_item)
    deallocate(hecMESH%export_index)
    deallocate(hecMESH%export_item)
    hecMESH%import_index => hecCOMM%import_index
    hecMESH%import_item => hecCOMM%import_item
    hecMESH%export_index => hecCOMM%export_index
    hecMESH%export_item => hecCOMM%export_item
    deallocate(hecCOMM%neighbor_pe)
  end subroutine hecmw_localmat_shrink_comm_table

  !> \brief Debug write matrix
  !>
  subroutine debug_write_matrix(Mat, label, level)
    type(hecmwST_local_matrix), intent(in) :: Mat   !< matrix
    character(len=*),           intent(in) :: label !< label for matrix
    integer(kind=kint),         intent(in) :: level !< debug level
    !
    integer(kind=kint) :: iunit

    if (level <= 0) return

    iunit = 700 + hecmw_comm_get_rank()
    write(iunit,'(a,a)') trim(label),'============================================================'
    if (level == 1) then
      call hecmw_localmat_write_size(Mat, iunit)
    else if (level == 2) then
      call hecmw_localmat_write_ij(Mat, iunit)
    else
      call hecmw_localmat_write(Mat, iunit)
    endif
  end subroutine debug_write_matrix

end module hecmw_local_matrix
