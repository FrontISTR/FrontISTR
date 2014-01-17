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
    integer :: method_org, precond_org
    logical :: fg_eliminate

    ! set if use eliminate version or not
    precond_org = hecmw_mat_get_precond(hecMAT)
    if (precond_org >= 30) then
      fg_eliminate = .false.
      call hecmw_mat_set_precond(hecMAT, precond_org - 20)
    else
      fg_eliminate = .true.
    endif

    if (fstrMAT%num_lagrange == 0) then
      write(0,*) 'INFO: no contact'
      ! use CG because the matrix is symmetric
      method_org = hecmw_mat_get_method(hecMAT)
      call hecmw_mat_set_method(hecMAT, 1)
      ! solve
      call hecmw_solve_33(hecMESH, hecMAT)
      ! restore solver setting
      call hecmw_mat_set_method(hecMAT, method_org)
    else
      write(0,*) 'INFO: with contact'
      if (fg_eliminate) then
        call solve_eliminate(hecMESH, hecMAT, fstrMAT)
      else
        call solve_no_eliminate(hecMESH, hecMAT, fstrMAT)
      endif
    endif

    if (precond_org >= 30) then
      call hecmw_mat_set_precond(hecMAT, precond_org)
    endif
  end subroutine solve_LINEQ_iter_contact

  subroutine solve_eliminate(hecMESH,hecMAT,fstrMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: ndof
    integer, allocatable :: iw2(:), iwS(:)
    real(kind=kreal), allocatable :: wSL(:), wSU(:)
    type(hecmwST_local_matrix) :: BTmat, BTtmat, BTtKT
    type(hecmwST_matrix) :: hecTKT
    real(kind=kreal) :: t1

    t1 = hecmw_wtime()
    write(0,*) 'INFO: solve_eliminate start', hecmw_wtime()-t1

    ndof=hecMAT%NDOF
    allocate(iw2(hecMAT%N*ndof))
    allocate(iwS(fstrMAT%num_lagrange), wSL(fstrMAT%num_lagrange), &
         wSU(fstrMAT%num_lagrange))

    ! choose slave DOFs to be eliminated with Lag. DOFs
    call choose_slaves(hecMAT, fstrMAT, iw2, iwS, wSL)
    call make_wSU(fstrMAT, hecMAT%N, ndof, iw2, wSU)
    write(0,*) 'INFO: Slave DOFs successfully chosen', hecmw_wtime()-t1

    ! make transformation matrix and its transpose
    call make_BTmat(hecMAT, fstrMAT, iw2, wSL, BTmat)
    !call hecmw_localmat_write(BTmat, 0)
    call make_BTtmat(hecMAT, fstrMAT, iw2, iwS, wSU, BTtmat)
    !call hecmw_localmat_write(BTtmat, 0)
    write(0,*) 'INFO: Making BTmat and BTtmat done', hecmw_wtime()-t1

    ! calc trimatmul in hecmwST_matrix data structure
    call hecmw_trimatmul_TtKT(BTtmat, hecMAT, BTmat, iwS, fstrMAT%num_lagrange, hecTKT)
    write(0,*) 'INFO: calculated hecTKT', hecmw_wtime()-t1

    ! make new RHS
    call make_new_b(hecMESH, hecMAT, BTtmat, iwS, wSL, &
         fstrMAT%num_lagrange, hecTKT%B)
    write(0,*) 'INFO: calculated RHS', hecmw_wtime()-t1

    ! use CG when the matrix is symmetric
    if (SymType == 1) call hecmw_mat_set_method(hecTKT, 1)

    ! solve
    call hecmw_solve_33(hecMESH, hecTKT)
    write(0,*) 'INFO: solver finished', hecmw_wtime()-t1

    ! calc u_s
    call hecmw_localmat_mulvec(BTmat, hecTKT%X, hecMAT%X)
    call subst_Blag(hecMAT, iwS, wSL, fstrMAT%num_lagrange)
    write(0,*) 'INFO: calculated disp', hecmw_wtime()-t1

    ! calc lambda
    call comp_lag(hecMAT, iwS, wSU, fstrMAT%num_lagrange)
    write(0,*) 'INFO: calculated lag', hecmw_wtime()-t1

    ! free matrices
    call hecmw_localmat_free(BTtmat)
    call hecmw_localmat_free(BTmat)
    call hecmw_mat_finalize(hecTKT)
    deallocate(iw2, iwS)
    write(0,*) 'INFO: solve_eliminate end', hecmw_wtime()-t1
  end subroutine solve_eliminate

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
!!$    write(0,*) 'iw1L, iw1U:'
!!$    do i=1,hecMAT%N*ndof
!!$      if (iw1L(i) > 0 .or. iw1U(i) > 0) write(0,*) i, iw1L(i), iw1U(i)
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
!!$    write(0,*) 'iw2:'
!!$    do i=1,hecMAT%N*ndof
!!$      if (iw2(i) > 0) write(0,*) i, iw2(i), iw1L(i), iw1U(i)
!!$    enddo
!!$    write(0,*) 'iwS:'
!!$    write(0,*) iwS(:)

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
    !write(0,*) wSU
  end subroutine make_wSU

  subroutine make_BTmat(hecMAT, fstrMAT, iw2, wSL, BTmat)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:)
    real(kind=kreal), intent(in) :: wSL(:)
    type (hecmwST_local_matrix), intent(out) :: BTmat
    type (hecmwST_local_matrix) :: Tmat
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
      write(0,*) Tmat%nnz, Tmat%index(Tmat%n)
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
        write(0,*) l, Tmat%index(i)+1
        stop 'ERROR: Tmat%index wrong'
      endif
    enddo
    !call hecmw_localmat_write(Tmat, 0)
    ! make 3x3-block version of Tmat
    call hecmw_localmat_blocking(Tmat, ndof, BTmat)
    call hecmw_localmat_free(Tmat)
  end subroutine make_BTmat

  subroutine make_BTtmat(hecMAT, fstrMAT, iw2, iwS, wSU, BTtmat)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:), iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    type (hecmwST_local_matrix) :: Ttmat
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
      write(0,*) Ttmat%nnz, Ttmat%index(Ttmat%n)
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
          write(0,*) l, Ttmat%index(idx)+1
          stop 'ERROR: Ttmat%index wrong'
        endif
      enddo
    enddo
    !call hecmw_localmat_write(Ttmat, 0)
    ! make 3x3-block version of Tmat
    call hecmw_localmat_blocking(Ttmat, ndof, BTtmat)
    call hecmw_localmat_free(Ttmat)
  end subroutine make_BTtmat

  subroutine make_new_b(hecMESH, hecMAT, BTtmat, iwS, wSL, num_lagrange, Bnew)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_local_matrix), intent(in) :: BTtmat
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
    !write(0,*) hecMAT%B(hecMAT%NP*ndof+1:hecMAT%NP*ndof+num_lagrange)
    do i=1,num_lagrange
      Bnew(iwS(i))=wSL(i)*hecMAT%B(hecMAT%NP*ndof+i)
    enddo
    !Btmp=B+K*B2
    call hecmw_matvec_33(hecMESH, hecMAT, Bnew, Btmp)
    Btmp=hecMAT%B+Btmp
    !B2=BTtmat*Btmp
    call hecmw_localmat_mulvec(BTtmat, Btmp, Bnew)

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

  subroutine solve_no_eliminate(hecMESH,hecMAT,fstrMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: ndof, ndof2, nb_lag, ndofextra
    integer :: i, ls, le, l, j, jb_lag, ib_lag, idof, jdof, ilag, k, ks, ke
    integer, allocatable :: iwUr(:), iwUc(:), iwLr(:), iwLc(:)
    type(hecmwST_matrix) :: hecMATLag
    real(kind=kreal) :: t1

    t1 = hecmw_wtime()
    write(0,*) 'INFO: solve_no_eliminate, start', hecmw_wtime()-t1

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    nb_lag = (fstrMAT%num_lagrange + 2)/3
    hecMATLag%NDOF = ndof
    hecMATLag%N = hecMAT%N + nb_lag
    hecMATLag%NP = hecMATLag%N
    !write(0,*) 'DEBUG: hecMAT: NDOF,N,NP=',hecMAT%NDOF,hecMAT%N,hecMAT%NP
    !write(0,*) 'DEBUG: hecMATLag: NDOF,N,NP=',hecMATLag%NDOF,hecMATLag%N,hecMATLag%NP

    ndofextra = hecMATLag%N*ndof - hecMAT%N*ndof - fstrMAT%num_lagrange
    write(0,*) 'DEBUG: num_lagrange,nb_lag,ndofextra=',fstrMAT%num_lagrange,nb_lag,ndofextra

    ! Upper: count num of blocks
    allocate(iwUr(hecMAT%N))
    allocate(iwUc(nb_lag))
    iwUr = 0
    do i = 1, hecMAT%N
      iwUc = 0
      ls=fstrMAT%indexU_lagrange(i-1)+1
      le=fstrMAT%indexU_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemU_lagrange(l)
        jb_lag = (j+2)/3
        iwUc(jb_lag) = 1
      enddo
      do j=1,nb_lag
        if (iwUc(j) > 0) iwUr(i) = iwUr(i) + 1
      enddo
      !if (iwUr(i) > 0) write(0,*) 'DEBUG: iwUr(',i,')=',iwUr(i)
    enddo

    ! Lower: count num of blocks
    allocate(iwLr(nb_lag))
    allocate(iwLc(hecMAT%N))
    iwLr = 0
    do ib_lag = 1, nb_lag
      iwLc = 0
      i = hecMAT%N + ib_lag
      do idof = 1, ndof
        ilag = (ib_lag-1)*ndof + idof
        if (ilag  > fstrMAT%num_lagrange) exit
        ls=fstrMAT%indexL_lagrange(ilag-1)+1
        le=fstrMAT%indexL_lagrange(ilag)
        do l=ls,le
          j=fstrMAT%itemL_lagrange(l)
          iwLc(j) = 1
        enddo
      enddo
      do j=1,hecMAT%N
        if (iwLc(j) > 0) iwLr(ib_lag) = iwLr(ib_lag) + 1
      enddo
      !if (iwLr(ib_lag) > 0) write(0,*) 'DEBUG: iwLr(',ib_lag,')=',iwLr(ib_lag)
    enddo

    ! Upper: indexU
    allocate(hecMATLag%indexU(0:hecMATLag%N))
    hecMATLag%indexU(0) = 0
    do i = 1, hecMAT%N
      hecMATLag%indexU(i) = hecMATLag%indexU(i-1) + &
           (hecMAT%indexU(i) - hecMAT%indexU(i-1)) + iwUr(i)
    enddo
    do i = hecMAT%N+1, hecMATLag%N
      hecMATLag%indexU(i) = hecMATLag%indexU(i-1)
    enddo
    hecMATLag%NPU = hecMATLag%indexU(hecMATLag%N)
    !write(0,*) 'DEBUG: hecMATLag%NPU=',hecMATLag%NPU

    ! Lower: indexL
    allocate(hecMATLag%indexL(0:hecMATLag%N))
    do i = 0, hecMAT%N
      hecMATLag%indexL(i) = hecMAT%indexL(i)
    enddo
    do i = hecMAT%N+1, hecMATLag%N
      hecMATLag%indexL(i) = hecMATLag%indexL(i-1) + iwLr(i-hecMAT%N)
    enddo
    hecMATLag%NPL = hecMATLag%indexL(hecMATLag%N)
    !write(0,*) 'DEBUG: hecMATLag%NPL=',hecMATLag%NPL

    ! Upper: itemU and AU
    allocate(hecMATLag%itemU(hecMATLag%NPU))
    allocate(hecMATLag%AU(hecMATLag%NPU*ndof2))
    do i = 1, hecMAT%N
      k = hecMATLag%indexU(i-1)+1
      ! copy from hecMAT
      ls = hecMAT%indexU(i-1)+1
      le = hecMAT%indexU(i)
      do l=ls,le
        hecMATLag%itemU(k) = hecMAT%itemU(l)
        hecMATLag%AU((k-1)*ndof2+1:k*ndof2)=hecMAT%AU((l-1)*ndof2+1:l*ndof2)
        k = k + 1
      enddo
      ! Lag. itemU
      iwUc = 0
      ls=fstrMAT%indexU_lagrange(i-1)+1
      le=fstrMAT%indexU_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemU_lagrange(l)
        jb_lag = (j+2)/3
        iwUc(jb_lag) = 1
      enddo
      do j=1,nb_lag
        if (iwUc(j) > 0) then
          hecMATLag%itemU(k) = hecMAT%N + j
          k = k + 1
        endif
      enddo
      if (k /= hecMATLag%indexU(i)+1) stop 'ERROR k indexU'
      ! Lag. AU
      ls=fstrMAT%indexU_lagrange(i-1)+1
      le=fstrMAT%indexU_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemU_lagrange(l)
        jb_lag = (j+2)/3
        jdof = j - (jb_lag - 1)*ndof
        ks=hecMATLag%indexU(i-1)+1
        ke=hecMATLag%indexU(i)
        do k = ks,ke
          if (hecMATLag%itemU(k) < hecMAT%N + jb_lag) cycle
          if (hecMATLag%itemU(k) > hecMAT%N + jb_lag) cycle
          !if (hecMATLag%itemU(k) /= hecMAT%N + jb_lag) stop 'ERROR itemU jb_lag'
          do idof = 1, ndof
            hecMATLag%AU((k-1)*ndof2+(idof-1)*ndof+jdof) = &
                 fstrMAT%AU_lagrange((l-1)*ndof+idof)
          enddo
        enddo
      enddo
    enddo

    ! Lower: itemL and AL
    allocate(hecMATLag%itemL(hecMATLag%NPL))
    allocate(hecMATLag%AL(hecMATLag%NPL*ndof2))
    do i = 1, hecMAT%NPL
      hecMATLag%itemL(i) = hecMAT%itemL(i)
    enddo
    do i = 1, hecMAT%NPL*ndof2
      hecMATLag%AL(i) = hecMAT%AL(i)
    enddo
    ! Lag. itemL
    k = hecMAT%NPL + 1
    do ib_lag = 1, nb_lag
      iwLc = 0
      i = hecMAT%N + ib_lag
      do idof = 1, ndof
        ilag = (ib_lag-1)*ndof + idof
        if (ilag  > fstrMAT%num_lagrange) exit
        ls=fstrMAT%indexL_lagrange(ilag-1)+1
        le=fstrMAT%indexL_lagrange(ilag)
        do l=ls,le
          j=fstrMAT%itemL_lagrange(l)
          iwLc(j) = 1
        enddo
      enddo
      do j=1,hecMAT%N
        if (iwLc(j) > 0) then
          hecMATLag%itemL(k) = j
          k = k + 1
        endif
      enddo
      if (k /= hecMATLag%indexL(hecMAT%N + ib_lag)+1) then
        stop 'ERROR k indexL'
      endif
      ! Lag. AL
      do idof = 1, ndof
        ilag = (ib_lag-1)*ndof + idof
        if (ilag  > fstrMAT%num_lagrange) exit
        ls=fstrMAT%indexL_lagrange(ilag-1)+1
        le=fstrMAT%indexL_lagrange(ilag)
        do l=ls,le
          j=fstrMAT%itemL_lagrange(l)
          ks=hecMATLag%indexL(i-1)+1
          ke=hecMATLag%indexL(i)
          do k = ks, ke
            if (hecMATLag%itemL(k) < j) cycle
            if (hecMATLag%itemL(k) > j) cycle
            !if (hecMATLag%itemL(k) /= j) stop 'ERROR itemL j'
            do jdof = 1, ndof
              hecMATLag%AL((k-1)*ndof2+(idof-1)*ndof+jdof) = &
                   fstrMAT%AL_lagrange((l-1)*ndof+jdof)
            enddo
          enddo
        enddo
      enddo
    enddo

    allocate(hecMATLag%D(hecMATLag%N*ndof2))
    hecMATLag%D = 0.d0
    do i = 1, hecMAT%N*ndof2
      hecMATLag%D(i) = hecMAT%D(i)
    enddo
    do idof = ndof - ndofextra + 1, ndof
      hecMATLag%D((hecMATLag%N-1)*ndof2 + (idof-1)*ndof + idof) = 1.0
    enddo

    allocate(hecMATLag%B(hecMATLag%N*ndof))
    allocate(hecMATLag%X(hecMATLag%N*ndof))
    hecMATLag%B = 0.d0
    hecMATLag%X = 0.d0

    do i = 1, hecMAT%N*ndof+fstrMAT%num_lagrange
      hecMATLag%B(i) = hecMAT%B(i)
    enddo

    hecMATLag%Iarray=hecMAT%Iarray
    hecMATLag%Rarray=hecMAT%Rarray

    write(0,*) 'INFO: made hecMATLag', hecmw_wtime()-t1

    call hecmw_solve_33(hecMESH, hecMATLag)

    do i = 1, hecMAT%N*ndof+fstrMAT%num_lagrange
      hecMAT%X(i) = hecMATLag%X(i)
    enddo

    call hecmw_mat_finalize(hecMATLag)

    write(0,*) 'INFO: solve_no_eliminate end', hecmw_wtime()-t1
  end subroutine solve_no_eliminate

end module m_solve_LINEQ_iter_contact
