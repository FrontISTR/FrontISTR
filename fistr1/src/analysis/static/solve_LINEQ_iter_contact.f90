!> This module provides linear equation solver interface of MUMPS for
!! contact problems using Lagrange multiplier.
module m_solve_LINEQ_iter_contact
  use m_fstr
  use fstr_matrix_con_contact

  private
  public :: solve_LINEQ_iter_contact_init
  public :: solve_LINEQ_iter_contact

  logical, save :: INITIALIZED = .false.
  integer, save :: SymType = 0

  type T_matrix
    integer :: n, nnz, ndof
    integer(kind=kint), pointer :: index(:)
    integer(kind=kint), pointer :: item(:)
    real(kind=kreal), pointer :: A(:)
  end type T_matrix

contains

  subroutine solve_LINEQ_iter_contact_init(hecMESH,hecMAT,fstrMAT,is_sym)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    logical, intent(in) :: is_sym

    if (INITIALIZED) then
       INITIALIZED = .false.
    endif

    if (is_sym) then
      SymType = 1
    else
      SymType = 0
    endif

    INITIALIZED = .true.
  end subroutine solve_LINEQ_iter_contact_init

  subroutine solve_LINEQ_iter_contact(hecMESH,hecMAT,fstrMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: ndof, method_org
    integer, allocatable :: iw2(:), iwS(:)
    real(kind=kreal), allocatable :: wSL(:), wSU(:)
    type(T_matrix) :: Tmat, Ttmat
    type(T_matrix) :: BTmat, BTtmat, BTtKT
    type(hecmwST_matrix) :: hecTKT

    if (fstrMAT%num_lagrange == 0) then

      ! use CG because the matrix is symmetric
      method_org = hecmw_mat_get_method(hecMAT)
      call hecmw_mat_set_method(hecMAT, 1)

      ! solve
      call hecmw_solve_33(hecMESH, hecMAT)

      ! restore solver setting
      call hecmw_mat_set_method(hecMAT, method_org)

    else

      ndof=hecMAT%NDOF
      allocate(iw2(hecMAT%N*ndof))
      allocate(iwS(fstrMAT%num_lagrange), wSL(fstrMAT%num_lagrange), &
           wSU(fstrMAT%num_lagrange))

      call choose_slaves(hecMAT, fstrMAT, iw2, iwS, wSL)
      call make_wSU(fstrMAT, hecMAT%N, ndof, iw2, wSU)
      write(*,*) 'INFO: Slave DOFs successfully chosen'

      call make_Tmat(hecMAT, fstrMAT, iw2, wSL, Tmat)
      !call write_Tmat(Tmat)
      call make_Ttmat(hecMAT, fstrMAT, iw2, iwS, wSU, Ttmat)
      !call write_Tmat(Ttmat)
      write(*,*) 'INFO: Making Tmat and Ttmat done'

      call blocking_Tmat(Tmat, ndof, BTmat)
      !call write_Tmat(BTmat)
      call free_Tmat(Tmat)

      call blocking_Tmat(Ttmat, ndof, BTtmat)
      !call write_Tmat(BTtmat)
      call free_Tmat(Ttmat)
      write(*,*) 'INFO: Blocking Tmat and Ttmat done'

      call trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
      !call write_Tmat(BTtKT)
      write(*,*) 'INFO: trimatmul done'

      call place_one_on_diag(BTtKT, iwS, fstrMAT%num_lagrange)

      !call replace_hecmat(hecMAT, BTtKT)
      call make_new_hecmat(hecMAT, BTtKT, hecTKT)
      call free_Tmat(BTtKT)

      call make_new_b(hecMESH, hecMAT, BTtmat, iwS, wSL, &
           fstrMAT%num_lagrange, hecTKT%B)

      ! use CG when the matrix is symmetric
      if (SymType == 1) call hecmw_mat_set_method(hecTKT, 1)

      ! solve
      call hecmw_solve_33(hecMESH, hecTKT)

      ! calc u_s
      call matvec_tmat(BTmat, hecTKT%X, hecMAT%X)
      call subst_Blag(hecMAT, iwS, wSL, fstrMAT%num_lagrange)
      ! calc lambda
      call comp_lag(hecMAT, iwS, wSU, fstrMAT%num_lagrange)

      call free_Tmat(BTtmat)
      call free_Tmat(BTmat)
      call free_hecmat(hecTKT)
      deallocate(iw2, iwS)
    endif
  end subroutine solve_LINEQ_iter_contact

  subroutine choose_slaves(hecMAT, fstrMAT, iw2, iwS, wSL)
    implicit none
    type (hecmwST_matrix    ), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(out) :: iw2(:), iwS(:)
    real(kind=kreal), intent(out) :: wSL(:)
    integer :: ndof, i, j, idof, jdof, l, ls, le, idx, imax
    real(kind=kreal) :: val, vmax
    integer, allocatable :: iw1L(:), iw1U(:)

    ndof=hecMAT%NDOF
    iw2=-1
    iwS=0

    allocate(iw1L(hecMAT%N*ndof))
    allocate(iw1U(hecMAT%N*ndof))
    iw1L=0
    iw1U=0

    ! Count how many times each dof appear in Lagrange matrix
    ! lower
    do i=1,fstrMAT%num_lagrange
      ls=fstrMAT%indexL_lagrange(i-1)+1
      le=fstrMAT%indexL_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemL_lagrange(l)
        do jdof=1,ndof
          idx=(j-1)*ndof+jdof
          iw1L(idx)=iw1L(idx)+1
        enddo
      enddo
    enddo
    ! upper
    do i=1,hecMAT%N
      ls=fstrMAT%indexU_lagrange(i-1)+1
      le=fstrMAT%indexU_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemU_lagrange(l)
        do idof=1,ndof
          idx=(i-1)*ndof+idof
          iw1U(idx)=iw1U(idx)+1
        enddo
      enddo
    enddo
!!$    write(*,*) 'iw1L, iw1U:'
!!$    do i=1,hecMAT%N*ndof
!!$      if (iw1L(i) > 0 .or. iw1U(i) > 0) write(*,*) i, iw1L(i), iw1U(i)
!!$    enddo

    ! Choose dofs that
    ! - appear only onece in both lower and upper Lag. and
    ! - has greatest coefficient among them (in lower Lag.)
    do i=1,fstrMAT%num_lagrange
      ls=fstrMAT%indexL_lagrange(i-1)+1
      le=fstrMAT%indexL_lagrange(i)
      vmax = 0.d0
      imax = -1
      do l=ls,le
        j=fstrMAT%itemL_lagrange(l)
        do jdof=1,ndof
          idx=(j-1)*ndof+jdof
          val=fstrMAT%AL_lagrange((l-1)*ndof+jdof)
          if (iw1L(idx) == 1 .and. iw1U(idx) == 1 .and. abs(val) > abs(vmax)) then
            imax=idx
            vmax=val
          endif
        enddo
      enddo
      if (imax == -1) stop "ERROR: iterative solver for contact failed"
      iw2(imax)=i
      iwS(i)=imax; wSL(i)=-1.d0/vmax
    enddo
!!$    write(*,*) 'iw2:'
!!$    do i=1,hecMAT%N*ndof
!!$      if (iw2(i) > 0) write(*,*) i, iw2(i), iw1L(i), iw1U(i)
!!$    enddo
!!$    write(*,*) 'iwS:'
!!$    write(*,*) iwS(:)

    deallocate(iw1L, iw1U)
  end subroutine choose_slaves

  subroutine make_wSU(fstrMAT, n, ndof, iw2, wSU)
    implicit none
    type(fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: iw2(:)
    real(kind=kreal), intent(out) :: wSU(:)
    integer(kind=kint) :: i, idof, idx, js, je, j, k

    wSU=0.d0
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        if (iw2(idx) > 0) then
          js=fstrMAT%indexU_lagrange(i-1)+1
          je=fstrMAT%indexU_lagrange(i)
          do j=js,je
            k=fstrMAT%itemU_lagrange(j)
            if (k==iw2(idx)) then
              wSU(iw2(idx)) = -1.0/fstrMAT%AU_lagrange((j-1)*ndof+idof)
            endif
          enddo
        endif
      enddo
    enddo
    !write(*,*) wSU
  end subroutine make_wSU

  subroutine write_Tmat(Tmat)
    implicit none
    type (T_matrix), intent(in) :: Tmat
    integer(kind=kint) :: n, ndof, ndof2, i, js, je, j, jj

    n=Tmat%n
    ndof=Tmat%ndof
    ndof2=ndof*ndof

    write(*,*) 'i, j, A'
    do i=1,n
      js=Tmat%index(i-1)+1
      je=Tmat%index(i)
      do j=js,je
        jj=Tmat%item(j)
        if (ndof==1) then
          write(*,*) i, jj, Tmat%A(j)
        else
          write(*,*) i, jj
          write(*,*) Tmat%A((j-1)*ndof2+1:j*ndof2)
        endif
      enddo
    enddo
  end subroutine write_Tmat

  subroutine make_Tmat(hecMAT, fstrMAT, iw2, wSL, Tmat)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:)
    real(kind=kreal), intent(in) :: wSL(:)
    type (T_matrix), intent(out) :: Tmat
    integer :: ndof, i, nnz, l, js, je, j, k, jdof, kk, jj
    real(kind=kreal) :: factor

    ndof=hecMAT%NDOF
    Tmat%n=hecMAT%N*ndof
    Tmat%nnz=Tmat%n+fstrMAT%numL_lagrange*ndof-2*fstrMAT%num_lagrange
    Tmat%ndof=1

    allocate(Tmat%index(0:hecMAT%N*ndof))
    allocate(Tmat%item(Tmat%nnz), Tmat%A(Tmat%nnz))
    ! index
    Tmat%index(0)=0
    do i=1,Tmat%n
      if (iw2(i) > 0) then
        nnz=ndof*(fstrMAT%indexL_lagrange(iw2(i))-fstrMAT%indexL_lagrange(iw2(i)-1))-1
      else
        nnz=1
      endif
      Tmat%index(i)=Tmat%index(i-1)+nnz
    enddo
    if (Tmat%nnz /= Tmat%index(Tmat%n)) then
      write(*,*) Tmat%nnz, Tmat%index(Tmat%n)
      stop 'ERROR: Tmat%nnz wrong'
    endif
    ! item and A
    do i=1,Tmat%n
      l=Tmat%index(i-1)+1
      if (iw2(i) > 0) then
        js=fstrMAT%indexL_lagrange(iw2(i)-1)+1
        je=fstrMAT%indexL_lagrange(iw2(i))
        factor=wSL(iw2(i))
        do j=js,je
          k=fstrMAT%itemL_lagrange(j)
          do jdof=1,ndof
            kk=(k-1)*ndof+jdof
            jj=(j-1)*ndof+jdof
            if (kk==i) cycle
            Tmat%item(l)=kk
            Tmat%A(l)=fstrMAT%AL_lagrange(jj)*factor
            l=l+1
          enddo
        enddo
      else
        Tmat%item(l)=i
        Tmat%A(l)=1.0
        l=l+1
      endif
      if (l /= Tmat%index(i)+1) then
        write(*,*) l, Tmat%index(i)+1
        stop 'ERROR: Tmat%index wrong'
      endif
    enddo
  end subroutine make_Tmat

  subroutine make_Ttmat(hecMAT, fstrMAT, iw2, iwS, wSU, Ttmat)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:), iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    type (T_matrix), intent(out) :: Ttmat
    integer :: ndof, i, nnz, l, js, je, j, k, idof, jdof, idx

    ndof=hecMAT%NDOF
    Ttmat%n=hecMAT%N*ndof
    Ttmat%nnz=Ttmat%n+fstrMAT%numU_lagrange*ndof-2*fstrMAT%num_lagrange
    Ttmat%ndof=1

    allocate(Ttmat%index(0:hecMAT%N*ndof))
    allocate(Ttmat%item(Ttmat%nnz), Ttmat%A(Ttmat%nnz))
    ! index
    Ttmat%index(0)=0
    do i=1,hecMAT%N
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        if (iw2(idx) <= 0) then
          nnz=fstrMAT%indexU_lagrange(i)-fstrMAT%indexU_lagrange(i-1)+1
        else
          nnz=0
        endif
        Ttmat%index(idx)=Ttmat%index(idx-1)+nnz
      enddo
    enddo
    if (Ttmat%nnz /= Ttmat%index(Ttmat%n)) then
      write(*,*) Ttmat%nnz, Ttmat%index(Ttmat%n)
      stop 'ERROR: Ttmat%nnz wrong'
    endif
    ! item and A
    do i=1,hecMAT%N
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        l=Ttmat%index(idx-1)+1
        if (iw2(idx) <= 0) then
          ! one on diagonal
          Ttmat%item(l)=idx
          Ttmat%A(l)=1.0
          l=l+1
          ! offdiagonal
          js=fstrMAT%indexU_lagrange(i-1)+1
          je=fstrMAT%indexU_lagrange(i)
          do j=js,je
            k=fstrMAT%itemU_lagrange(j)
            Ttmat%item(l)=iwS(k)
            Ttmat%A(l)=fstrMAT%AU_lagrange((j-1)*ndof+idof)*wSU(k)
            l=l+1
          enddo
        else
          ! no element
        endif
        if (l /= Ttmat%index(idx)+1) then
          write(*,*) l, Ttmat%index(idx)+1
          stop 'ERROR: Ttmat%index wrong'
        endif
      enddo
    enddo
  end subroutine make_Ttmat

  subroutine blocking_Tmat(Tmat, ndof, BTmat)
    implicit none
    type (T_matrix), intent(in) :: Tmat
    integer, intent(in) :: ndof
    type (T_matrix), intent(out) :: BTmat
    integer, allocatable :: iw(:)
    integer :: ndof2, i, icnt, idof, idx, ls, le, l, j, jb, k, lb0, jdof, ks, ke
    ndof2=ndof*ndof

    if (mod(Tmat%n, ndof) /= 0) then
      write(*,*) Tmat%n, ndof
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
  end subroutine blocking_Tmat

  subroutine free_Tmat(Tmat)
    implicit none
    type (T_matrix), intent(inout) :: Tmat
    deallocate(Tmat%index)
    deallocate(Tmat%item)
    deallocate(Tmat%A)
    Tmat%n=0
    Tmat%nnz=0
    Tmat%ndof=0
  end subroutine free_Tmat

  subroutine trimatmul_TtKT(BTtmat, hecMAT, BTmat, BTtKT)
    implicit none
    type (T_matrix), intent(in) :: BTtmat, BTmat
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (T_matrix), intent(out) :: BTtKT
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
            !if (i==1) write(*,*) 'l', icnt, jj, kk, ll
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
          !if (i==1) write(*,*) 'd', icnt, jj, kk, ll
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
            !if (i==1) write(*,*) 'u', icnt, jj, kk, ll
          enddo ll3
        enddo
      enddo
      if (icnt == 0) icnt=1
      !if (i==1) write(*,*) iw(1:icnt)
      BTtKT%index(i)=BTtKT%index(i-1)+icnt
    enddo
    !write(*,*) BTtKT%index(1:n)-BTtKT%index(0:n-1)

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
              !if (i==1) write(*,*) 'l', mm, jj, kk, ll
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
            !if (i==1) write(*,*) 'd', mm, jj, kk, ll
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
              !if (i==1) write(*,*) 'u', mm, jj, kk, ll
            endif
            TtKTp=>BTtKT%A((mm-1)*ndof2+1:mm*ndof2)
            call blk_trimatmul_add(ndof, Ttp, Kp, Tp, TtKTp)
          enddo
        enddo
      enddo
      if (icnt == 0) icnt=1
      ! error check!
      !write(*,*) BTtKT%item(ms:ms-1+icnt)
      !write(*,*) BTtKT%index(i)-BTtKT%index(i-1), icnt
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
    type (T_matrix), intent(inout) :: BTtKT
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
          !write(*,*) ilag, i, idof
          BTtKT%A((j-1)*ndof2+(idof-1)*ndof+idof)=1.d0
          cycle outer
        endif
      enddo
    enddo outer
  end subroutine place_one_on_diag

  subroutine replace_hecmat(hecMAT, BTtKT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (T_matrix), intent(in) :: BTtKT
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
    type(T_matrix), intent(in) :: BTtKT
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

  subroutine free_hecmat(hecMAT)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    if (associated(hecMAT%D)) deallocate(hecMAT%D)
    if (associated(hecMAT%B)) deallocate(hecMAT%B)
    if (associated(hecMAT%X)) deallocate(hecMAT%X)
    if (associated(hecMAT%AL)) deallocate(hecMAT%AL)
    if (associated(hecMAT%AU)) deallocate(hecMAT%AU)
    if (associated(hecMAT%indexL)) deallocate(hecMAT%indexL)
    if (associated(hecMAT%indexU)) deallocate(hecMAT%indexU)
    if (associated(hecMAT%itemL)) deallocate(hecMAT%itemL)
    if (associated(hecMAT%itemU)) deallocate(hecMAT%itemU)
    if (associated(hecMAT%ALU)) deallocate(hecMAT%ALU)
  end subroutine free_hecmat

  subroutine matvec_tmat(BTmat, V, TV)
    implicit none
    type (T_matrix), intent(in) :: BTmat
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
!!$    write(*,*) 'vnorm:', sqrt(vnorm)

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
  end subroutine matvec_tmat

  subroutine make_new_b(hecMESH, hecMAT, BTtmat, iwS, wSL, num_lagrange, Bnew)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(T_matrix), intent(in) :: BTtmat
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSL(:)
    integer(kind=kint), intent(in) :: num_lagrange
    real(kind=kreal), intent(out) :: Bnew(:)
    real(kind=kreal), allocatable :: Btmp(:)
    integer(kind=kint) :: n, ndof, i

    n=hecMAT%N
    ndof=hecMAT%NDOF

    allocate(Btmp(n*ndof))

    !B2=-Bs^-1*Blag
    Bnew=0.d0
    !write(*,*) hecMAT%B(hecMAT%NP*ndof+1:hecMAT%NP*ndof+num_lagrange)
    do i=1,num_lagrange
      Bnew(iwS(i))=wSL(i)*hecMAT%B(hecMAT%NP*ndof+i)
    enddo
    !Btmp=B+K*B2
    call hecmw_matvec_33(hecMESH, hecMAT, Bnew, Btmp)
    Btmp=hecMAT%B+Btmp
    !B2=BTtmat*Btmp
    call matvec_tmat(BTtmat, Btmp, Bnew)

    deallocate(Btmp)
  end subroutine make_new_b

  subroutine subst_Blag(hecMAT, iwS, wSL, num_lagrange)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSL(:)
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint) :: npndof, i, ilag

    npndof=hecMAT%NP*hecMAT%NDOF
    do i=1,num_lagrange
      ilag=iwS(i)
      hecMAT%X(ilag)=hecMAT%X(ilag)-wSL(i)*hecMAT%B(npndof+i)
    enddo
  end subroutine subst_Blag

  subroutine comp_lag(hecMAT, iwS, wSU, num_lagrange)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint) :: n, ndof, ndof2, ilag, iS, i, idof
    integer(kind=kint) :: js, je, j, jj, ij0, j0, jdof
    real(kind=kreal), pointer :: xlag(:)

    n=hecMAT%N
    ndof=hecMAT%ndof
    ndof2=ndof*ndof

    xlag=>hecMAT%X(hecMAT%NP*ndof+1:hecMAT%NP*ndof+num_lagrange)

    do ilag=1,num_lagrange
      iS=iwS(ilag)
      i=(iS-1)/ndof+1
      idof=mod(iS-1, ndof)+1
      xlag(ilag)=hecMAT%B(iS)
      !lower
      js=hecMAT%indexL(i-1)+1
      je=hecMAT%indexL(i)
      do j=js,je
        jj=hecMAT%itemL(j)
        ij0=(j-1)*ndof2+(idof-1)*ndof
        j0=(jj-1)*ndof
        do jdof=1,ndof
          xlag(ilag)=xlag(ilag)-hecMAT%AL(ij0+jdof)*hecMAT%X(j0+jdof)
        enddo
      enddo
      !diag
      ij0=(i-1)*ndof2+(idof-1)*ndof
      j0=(i-1)*ndof
      do jdof=1,ndof
        xlag(ilag)=xlag(ilag)-hecMAT%D(ij0+jdof)*hecMAT%X(j0+jdof)
      enddo
      !upper
      js=hecMAT%indexU(i-1)+1
      je=hecMAT%indexU(i)
      do j=js,je
        jj=hecMAT%itemU(j)
        ij0=(j-1)*ndof2+(idof-1)*ndof
        j0=(jj-1)*ndof
        do jdof=1,ndof
          xlag(ilag)=xlag(ilag)-hecMAT%AU(ij0+jdof)*hecMAT%X(j0+jdof)
        enddo
      enddo
      xlag(ilag)=-wSU(ilag)*xlag(ilag)
    enddo
  end subroutine comp_lag

end module m_solve_LINEQ_iter_contact
