!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides interface of iteratie linear equation solver for
!! contact problems using Lagrange multiplier.
module m_solve_LINEQ_iter_contact
  use m_fstr
  use fstr_matrix_con_contact
  use hecmw_solver

  private
  public :: solve_LINEQ_iter_contact_init
  public :: solve_LINEQ_iter_contact

  logical, save :: INITIALIZED = .false.
  integer, save :: SymType = 0

  integer, parameter :: DEBUG = 0  ! 0: no message, 1: some messages, 2: more messages, 3: vector output, 4: matrix output

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

    hecMAT%Iarray(98) = 1
    hecMAT%Iarray(97) = 1

    if (is_sym) then
      SymType = 1
    else
      SymType = 0
    endif

    INITIALIZED = .true.
  end subroutine solve_LINEQ_iter_contact_init

  subroutine solve_LINEQ_iter_contact(hecMESH,hecMAT,fstrMAT,istat,conMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer(kind=kint), intent(out) :: istat
    type (hecmwST_matrix), intent(in),optional :: conMAT
    integer :: method_org, precond_org
    logical :: fg_eliminate
    logical :: fg_amg
    integer(kind=kint) :: num_lagrange_global

    hecMAT%Iarray(97) = 1

    ! set if use eliminate version or not
    fg_amg = .false.
    precond_org = hecmw_mat_get_precond(hecMAT)
    if (precond_org >= 30 .and. precond_org <= 32) then
      fg_eliminate = .false.
      call hecmw_mat_set_precond(hecMAT, precond_org - 20)
    else
      fg_eliminate = .true.
      if (precond_org == 5) fg_amg = .true.
    endif

    num_lagrange_global = fstrMAT%num_lagrange
    call hecmw_allreduce_I1(hecMESH, num_lagrange_global, hecmw_sum)

    if (num_lagrange_global == 0) then
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: no contact'
      ! use CG because the matrix is symmetric
      method_org = hecmw_mat_get_method(hecMAT)
      call hecmw_mat_set_method(hecMAT, 1)
      ! avoid ML when no contact
      !if (fg_amg) call hecmw_mat_set_precond(hecMAT, 3) ! set diag-scaling
      ! solve
      call hecmw_solve(hecMESH, hecMAT)
      ! restore solver setting
      call hecmw_mat_set_method(hecMAT, method_org)
      !if (fg_amg) call hecmw_mat_set_precond(hecMAT, 5)
    else
      if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: with contact'
      if (fg_eliminate) then
        if(paraContactFlag.and.present(conMAT)) then
          call solve_eliminate(hecMESH, hecMAT, fstrMAT, conMAT)
        else
          call solve_eliminate(hecMESH, hecMAT, fstrMAT)
        endif
      else
        if(paraContactFlag.and.present(conMAT)) then
          write(0,*) 'ERROR: solve_no_eliminate not implemented for ParaCon'
          call hecmw_abort( hecmw_comm_get_comm())
        endif
        if (fstrMAT%num_lagrange > 0) then
          call solve_no_eliminate(hecMESH, hecMAT, fstrMAT)
        else
          call hecmw_solve(hecMESH, hecMAT)
        endif
      endif
    endif

    if (precond_org >= 30) then
      call hecmw_mat_set_precond(hecMAT, precond_org)
    endif

    istat = hecmw_mat_get_flag_diverged(hecMAT)
  end subroutine solve_LINEQ_iter_contact

  !!
  !! Solve with elimination of Lagrange-multipliers
  !!

  subroutine solve_eliminate(hecMESH,hecMAT,fstrMAT,conMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in), target :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    type (hecmwST_matrix), intent(in),optional :: conMAT
    integer(kind=kint) :: ndof
    integer(kind=kint), allocatable :: iw2(:), iwS(:)
    real(kind=kreal), allocatable :: wSL(:), wSU(:)
    type(hecmwST_local_matrix), target :: BTmat
    type(hecmwST_local_matrix) :: BTtmat
    type(hecmwST_matrix) :: hecTKT
    type(hecmwST_local_mesh), pointer :: hecMESHtmp
    type (hecmwST_local_matrix), pointer :: BT_all
    real(kind=kreal) :: t1
    integer(kind=kint) :: method_org, i, ilag
    logical :: fg_paracon
    type (hecmwST_local_matrix), pointer :: BKmat, BTtKmat, BTtKTmat, BKtmp
    integer(kind=kint), allocatable :: mark(:)
    type(fstrST_contact_comm) :: conCOMM
    integer(kind=kint) :: n_contact_dof, n_slave
    integer(kind=kint), allocatable :: contact_dofs(:), slaves(:)
    real(kind=kreal), allocatable :: Btot(:)

    t1 = hecmw_wtime()
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate start', hecmw_wtime()-t1

    fg_paracon = (paraContactFlag.and.present(conMAT))

    ndof=hecMAT%NDOF
    if(fg_paracon) then
      allocate(iw2(hecMAT%NP*ndof))
    else
      allocate(iw2(hecMAT%N*ndof))
    endif
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: num_lagrange',fstrMAT%num_lagrange
    allocate(iwS(fstrMAT%num_lagrange), wSL(fstrMAT%num_lagrange), &
      wSU(fstrMAT%num_lagrange))

    ! choose slave DOFs to be eliminated with Lag. DOFs
    call choose_slaves(hecMAT, fstrMAT, iw2, iwS, wSL, wSU, fg_paracon)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: slave DOFs chosen', hecmw_wtime()-t1

    ! make transformation matrix and its transpose
    call make_BTmat(hecMAT, fstrMAT, iw2, wSL, BTmat, fg_paracon)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make T done', hecmw_wtime()-t1
    if (DEBUG >= 4) then
      write(1000+myrank,*) 'BTmat (local)'
      call hecmw_localmat_write(BTmat, 1000+myrank)
    endif

    call make_BTtmat(hecMAT, fstrMAT, iw2, iwS, wSU, BTtmat, fg_paracon)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: make Tt done', hecmw_wtime()-t1
    if (DEBUG >= 4) then
      write(1000+myrank,*) 'BTtmat (local)'
      call hecmw_localmat_write(BTtmat, 1000+myrank)
    endif

    if(fg_paracon) then
      ! make contact dof list
      call make_contact_dof_list(hecMAT, fstrMAT, n_contact_dof, contact_dofs)

      ! make comm_table for contact dof
      ! call fstr_contact_comm_init(conCOMM, hecMESH, ndof, fstrMAT%num_lagrange, iwS)
      call fstr_contact_comm_init(conCOMM, hecMESH, ndof, n_contact_dof, contact_dofs)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: make contact comm_table done', hecmw_wtime()-t1

      ! copy hecMESH to hecMESHtmp
      allocate(hecMESHtmp)
      call copy_mesh(hecMESH, hecMESHtmp, fg_paracon)

      ! communicate and complete BTmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BTmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble T done', hecmw_wtime()-t1
      if (BTmat%nc /= hecMESH%n_node) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with T',BTmat%nc-hecMESH%n_node
      endif
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTmat (assembled)'
        call hecmw_localmat_write(BTmat, 1000+myrank)
      endif

      ! communicate and complete BTtmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BTtmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble Tt done', hecmw_wtime()-t1
      if (BTtmat%nc /= BTmat%nc) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with BTtmat',BTtmat%nc-BTmat%nc
        BTmat%nc = BTtmat%nc
      endif
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTtmat (assembled)'
        call hecmw_localmat_write(BTtmat, 1000+myrank)
      endif

      ! place 1 on diag of non-slave dofs of BTmat and BTtmat
      allocate(mark(BTmat%nr * BTmat%ndof))
      call mark_slave_dof(BTmat, mark, n_slave, slaves)
      call place_one_on_diag_of_unmarked_dof(BTmat, mark)
      call place_one_on_diag_of_unmarked_dof(BTtmat, mark)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: place 1 on diag of T and Tt done', hecmw_wtime()-t1
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTmat (1s on non-slave diag)'
        call hecmw_localmat_write(BTmat, 1000+myrank)
        write(1000+myrank,*) 'BTtmat (1s on non-slave diag)'
        call hecmw_localmat_write(BTtmat, 1000+myrank)
      endif

      ! init BKmat and substitute conMAT
      allocate(BKmat)
      call hecmw_localmat_init_with_hecmat(BKmat, conMAT, fstrMAT%num_lagrange)
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BKmat (conMAT local)'
        call hecmw_localmat_write(BKmat, 1000+myrank)
      endif

      ! communicate and complete BKmat (update hecMESHtmp)
      call hecmw_localmat_assemble(BKmat, hecMESH, hecMESHtmp)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble K (conMAT) done', hecmw_wtime()-t1
      if (BKmat%nc /= BTtmat%nc) then
        if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: node migrated with BKmat',BKmat%nc-BTtmat%nc
        BTmat%nc = BKmat%nc
        BTtmat%nc = BKmat%nc
      endif
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BKmat (conMAT assembled)'
        call hecmw_localmat_write(BKmat, 1000+myrank)
      endif

      ! add hecMAT to BKmat
      call hecmw_localmat_add_hecmat(BKmat, hecMAT)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: add hecMAT to K done', hecmw_wtime()-t1
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BKmat (hecMAT added)'
        call hecmw_localmat_write(BKmat, 1000+myrank)
      endif

      ! compute BTtKmat = BTtmat * BKmat (update hecMESHtmp)
      allocate(BTtKmat)
      call hecmw_localmat_multmat(BTtmat, BKmat, hecMESHtmp, BTtKmat)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply Tt and K done', hecmw_wtime()-t1
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTtKmat'
        call hecmw_localmat_write(BTtKmat, 1000+myrank)
      endif

      ! compute BTtKTmat = BTtKmat * BTmat (update hecMESHtmp)
      allocate(BTtKTmat)
      call hecmw_localmat_multmat(BTtKmat, BTmat, hecMESHtmp, BTtKTmat)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: multiply TtK and T done', hecmw_wtime()-t1
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTtKTmat'
        call hecmw_localmat_write(BTtKTmat, 1000+myrank)
      endif
      call hecmw_localmat_free(BTtKmat)
      deallocate(BTtKmat)

      ! shrink comm_table
      ! call hecmw_localmat_shrink_comm_table(BTtKTmat, hecMESHtmp)

      call place_num_on_diag_of_marked_dof(BTtKTmat, 1.0d0, mark)
      if (DEBUG >= 4) then
        write(1000+myrank,*) 'BTtKTmat (place 1.0 on slave diag)'
        call hecmw_localmat_write(BTtKTmat, 1000+myrank)
      endif
      call hecmw_mat_init(hecTKT)
      call hecmw_localmat_make_hecmat(hecMAT, BTtKTmat, hecTKT)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: convert TtKT to hecTKT done', hecmw_wtime()-t1
      call hecmw_localmat_free(BTtKTmat)
      deallocate(BTtKTmat)
      !
      BT_all => BTmat
    else
      if (hecMESH%n_neighbor_pe > 0) then
        ! update communication table
        allocate(hecMESHtmp, BT_all)
        call update_comm_table(hecMESH, BTmat, hecMESHtmp, BT_all)
        if (DEBUG >= 2) write(0,*) '  DEBUG2: Updating communication table done', hecmw_wtime()-t1
        call hecmw_localmat_free(BTmat)
      else
        ! in serial computation
        hecMESHtmp => hecMESH
        BT_all => BTmat
      end if

      ! calc trimatmul in hecmwST_matrix data structure
      call hecmw_mat_init(hecTKT)
      call hecmw_trimatmul_TtKT(hecMESHtmp, BTtmat, hecMAT, BT_all, iwS, fstrMAT%num_lagrange, hecTKT)
    endif
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: calculated TtKT', hecmw_wtime()-t1

    ! original RHS
    if (DEBUG >= 3) then
      write(1000+myrank,*) '======================================================================='
      write(1000+myrank,*) 'RHS(original)----------------------------------------------------------'
      write(1000+myrank,*) 'size of hecMAT%B',size(hecMAT%B)
      write(1000+myrank,*) 'hecMAT%B: 1-',hecMAT%N*ndof
      write(1000+myrank,*) hecMAT%B(1:hecMAT%N*ndof)
      if (fstrMAT%num_lagrange > 0) then
        write(1000+myrank,*) 'hecMAT%B(lag):',hecMAT%NP*ndof+1,'-',hecMAT%NP*ndof+fstrMAT%num_lagrange
        write(1000+myrank,*) hecMAT%B(hecMAT%NP*ndof+1:hecMAT%NP*ndof+fstrMAT%num_lagrange)
      endif
      if (n_slave > 0) then
        write(1000+myrank,*) 'hecMAT%B(slave):',slaves(:)
        write(1000+myrank,*) hecMAT%B(slaves(:))
      endif
    endif

    if (fg_paracon) then
      ! contact RHS
      if (DEBUG >= 3) then
        write(1000+myrank,*) 'RHS(conMAT)------------------------------------------------------------'
        write(1000+myrank,*) 'size of conMAT%B',size(conMAT%B)
        write(1000+myrank,*) 'conMAT%B: 1-',conMAT%N*ndof
        write(1000+myrank,*) conMAT%B(1:conMAT%N*ndof)
        write(1000+myrank,*) 'conMAT%B(external): ',conMAT%N*ndof+1,'-',conMAT%NP*ndof
        write(1000+myrank,*) conMAT%B(conMAT%N*ndof+1:conMAT%NP*ndof)
        if (fstrMAT%num_lagrange > 0) then
          write(1000+myrank,*) 'conMAT%B(lag):',conMAT%NP*ndof+1,'-',conMAT%NP*ndof+fstrMAT%num_lagrange
          write(1000+myrank,*) conMAT%B(conMAT%NP*ndof+1:conMAT%NP*ndof+fstrMAT%num_lagrange)
        endif
        if (n_slave > 0) then
          write(1000+myrank,*) 'conMAT%B(slave):',slaves(:)
          write(1000+myrank,*) conMAT%B(slaves(:))
        endif
      endif

      allocate(Btot(hecMAT%NP*ndof+fstrMAT%num_lagrange))
      call assemble_b_paracon(hecMESH, hecMAT, conMAT, fstrMAT%num_lagrange, conCOMM, Btot)
      if (DEBUG >= 2) write(0,*) '  DEBUG2: assemble conMAT%B done', hecmw_wtime()-t1
      if (DEBUG >= 3) then
        write(1000+myrank,*) 'RHS(conMAT assembled)--------------------------------------------------'
        write(1000+myrank,*) 'size of Btot',size(Btot)
        write(1000+myrank,*) 'Btot: 1-',conMAT%N*ndof
        write(1000+myrank,*) Btot(1:conMAT%N*ndof)
        if (fstrMAT%num_lagrange > 0) then
          write(1000+myrank,*) 'Btot(lag):',conMAT%NP*ndof+1,'-',conMAT%NP*ndof+fstrMAT%num_lagrange
          write(1000+myrank,*) Btot(conMAT%NP*ndof+1:conMAT%NP*ndof+fstrMAT%num_lagrange)
        endif
        if (n_slave > 0) then
          write(1000+myrank,*) 'Btot(slave):',slaves(:)
          write(1000+myrank,*) Btot(slaves(:))
        endif
      endif

      do i=1,hecMAT%N*ndof
        Btot(i)=Btot(i)+hecMAT%B(i)
      enddo
      ! assembled RHS
      if (DEBUG >= 3) then
        write(1000+myrank,*) 'RHS(total)-------------------------------------------------------------'
        write(1000+myrank,*) 'size of Btot',size(Btot)
        write(1000+myrank,*) 'Btot: 1-',conMAT%N*ndof
        write(1000+myrank,*) Btot(1:conMAT%N*ndof)
        if (fstrMAT%num_lagrange > 0) then
          write(1000+myrank,*) 'Btot(lag):',conMAT%NP*ndof+1,'-',conMAT%NP*ndof+fstrMAT%num_lagrange
          write(1000+myrank,*) Btot(conMAT%NP*ndof+1:conMAT%NP*ndof+fstrMAT%num_lagrange)
        endif
        if (n_slave > 0) then
          write(1000+myrank,*) 'Btot(slave):',slaves(:)
          write(1000+myrank,*) Btot(slaves(:))
        endif
      endif
    endif

    allocate(hecTKT%B(hecTKT%NP*ndof + fstrMAT%num_lagrange))
    allocate(hecTKT%X(hecTKT%NP*ndof + fstrMAT%num_lagrange))
    if (fg_paracon) then
      do i=1, hecTKT%N*ndof
        hecTKT%B(i) = Btot(i)
      enddo
    else
      do i=1, hecTKT%N*ndof
        hecTKT%B(i) = hecMAT%B(i)
      enddo
    endif
    do ilag=1, fstrMAT%num_lagrange
      hecTKT%B(hecTKT%NP*ndof + ilag) = hecMAT%B(hecMAT%NP*ndof + ilag)
    enddo
    do i=1, hecTKT%N*ndof
      hecTKT%X(i) = hecMAT%X(i)
    enddo
    do ilag=1,fstrMAT%num_lagrange
      hecTKT%X(iwS(ilag)) = 0.d0
    enddo

    ! make new RHS
    if (fg_paracon) then
      call make_new_b_paracon(hecMESH, hecMAT, conMAT, Btot, hecMESHtmp, hecTKT, BTtmat, BKmat, &
           iwS, wSL, fstrMAT%num_lagrange, conCOMM, hecTKT%B)
    else
      call make_new_b(hecMESH, hecMAT, BTtmat, iwS, wSL, &
           fstrMAT%num_lagrange, hecTKT%B)
    endif
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: converted RHS', hecmw_wtime()-t1
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'RHS(converted)---------------------------------------------------------'
      write(1000+myrank,*) 'size of hecTKT%B',size(hecTKT%B)
      write(1000+myrank,*) 'hecTKT%B: 1-',hecTKT%N*ndof
      write(1000+myrank,*) hecTKT%B(1:hecTKT%N*ndof)
    endif

    ! ! use CG when the matrix is symmetric
    ! if (SymType == 1) call hecmw_mat_set_method(hecTKT, 1)

    ! solve
    call hecmw_solve(hecMESHtmp, hecTKT)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: linear solver finished', hecmw_wtime()-t1

    ! solution in converged system
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Solution(converted)----------------------------------------------------'
      write(1000+myrank,*) 'size of hecTKT%X',size(hecTKT%X)
      write(1000+myrank,*) 'hecTKT%X: 1-',hecTKT%N*ndof
      write(1000+myrank,*) hecTKT%X(1:hecTKT%N*ndof)
    endif

    hecMAT%Iarray=hecTKT%Iarray
    hecMAT%Rarray=hecTKT%Rarray

    ! calc u_s
    if (fg_paracon) then
      call comp_x_slave_paracon(hecMESH, hecMAT, Btot, hecMESHtmp, hecTKT, BTmat, &
           fstrMAT%num_lagrange, iwS, wSL, conCOMM, n_slave, slaves)
    else
      call hecmw_localmat_mulvec(BT_all, hecTKT%X, hecMAT%X) !!!<== maybe, BT_all should be BTmat ???
      call subst_Blag(hecMAT, iwS, wSL, fstrMAT%num_lagrange)
    endif
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: recovered slave disp', hecmw_wtime()-t1

    ! calc lambda
    if (fg_paracon) then
      call comp_lag_paracon(hecMESH, hecMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
           n_slave, slaves, fstrMAT%num_lagrange, iwS, wSU, conCOMM)
      ! call comp_lag_paracon2(hecMESH, hecMAT, conMAT, fstrMAT%num_lagrange, iwS, wSU, conCOMM)
    else
      call comp_lag(hecMAT, iwS, wSU, fstrMAT%num_lagrange)
    endif
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: calculated lag', hecmw_wtime()-t1

    ! write(0,*) 'size of conMAT%X',size(conMAT%X)
    ! conMAT%X(:) = hecMAT%X(:)

    ! solution in original system
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Solution(original)-----------------------------------------------------'
      write(1000+myrank,*) 'size of hecMAT%X',size(hecMAT%X)
      write(1000+myrank,*) 'hecMAT%X: 1-',hecMAT%N*ndof
      write(1000+myrank,*) hecMAT%X(1:hecMAT%N*ndof)
      if (fstrMAT%num_lagrange > 0) then
        write(1000+myrank,*) 'hecMAT%X(lag):',hecMAT%NP*ndof+1,'-',hecMAT%NP*ndof+fstrMAT%num_lagrange
        write(1000+myrank,*) hecMAT%X(hecMAT%NP*ndof+1:hecMAT%NP*ndof+fstrMAT%num_lagrange)
      endif
      if (n_slave > 0) then
        write(1000+myrank,*) 'hecMAT%X(slave):',slaves(:)
        write(1000+myrank,*) hecMAT%X(slaves(:))
      endif
    endif

    ! check solution in the original system
    if (fg_paracon .and. DEBUG >= 2)  then
      ! call check_solution2(hecMESH, hecMAT, conMAT, fstrMAT, n_contact_dof, contact_dofs, &
      !      conCOMM, n_slave, slaves)
      call check_solution(hecMESH, hecMAT, fstrMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
           conCOMM, n_slave, slaves)
    endif

    ! free matrices
    call hecmw_localmat_free(BT_all)
    call hecmw_localmat_free(BTtmat)
    call hecmw_mat_finalize(hecTKT)
    if (fg_paracon) then
      call hecmw_localmat_free(BKmat)
      deallocate(BKmat)
      call fstr_contact_comm_finalize(conCOMM)
      call free_mesh(hecMESHtmp, fg_paracon)
      deallocate(hecMESHtmp)
    else
      if (hecMESH%n_neighbor_pe > 0) then
        call free_mesh(hecMESHtmp)
        deallocate(hecMESHtmp)
        deallocate(BT_all)
      endif
    endif
    deallocate(iw2, iwS)
    if ((DEBUG >= 1 .and. myrank==0) .or. DEBUG >= 2) write(0,*) 'DEBUG: solve_eliminate end', hecmw_wtime()-t1
  end subroutine solve_eliminate

  subroutine choose_slaves(hecMAT, fstrMAT, iw2, iwS, wSL, wSU, fg_paracon)
    implicit none
    type (hecmwST_matrix    ), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(out) :: iw2(:), iwS(:)
    real(kind=kreal), intent(out) :: wSL(:)
    real(kind=kreal), intent(out) :: wSU(:)
    logical, intent(in) :: fg_paracon
    integer :: ndof, n, i, j, idof, jdof, l, ls, le, idx, imax
    real(kind=kreal) :: val, vmax
    integer, allocatable :: iw1L(:), iw1U(:)
    integer(kind=kint) :: n_slave_in,n_slave_out

    iw2=-1

    if (fstrMAT%num_lagrange == 0) return

    ndof=hecMAT%NDOF
    iwS=0

    if (fg_paracon) then
      n = hecMAT%NP
    else
      n = hecMAT%N
    endif

    allocate(iw1L(n*ndof))
    allocate(iw1U(n*ndof))
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
    do i=1,n
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
    !!$    do i=1,n*ndof
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
    !!$    do i=1,n*ndof
    !!$      if (iw2(i) > 0) write(0,*) i, iw2(i), iw1L(i), iw1U(i)
    !!$    enddo
    !!$    write(0,*) 'iwS:'
    !!$    write(0,*) iwS(:)
    n_slave_in = 0
    do i=1,hecMAT%N*ndof
      if (iw2(i) > 0) n_slave_in = n_slave_in + 1
    enddo
    n_slave_out = 0
    do i=hecMAT%N*ndof,n*ndof
      if (iw2(i) > 0) n_slave_out = n_slave_out + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2[',myrank,']: n_slave(in,out,tot)',n_slave_in,n_slave_out,n_slave_in+n_slave_out

    deallocate(iw1L, iw1U)

    call make_wSU(fstrMAT, n, ndof, iw2, wSU)
  end subroutine choose_slaves

  subroutine make_wSU(fstrMAT, n, ndof, iw2, wSU)
    implicit none
    type(fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: iw2(:)
    real(kind=kreal), intent(out) :: wSU(:)
    integer(kind=kint) :: i, idof, idx, js, je, j, k

    if (fstrMAT%num_lagrange == 0) return

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

  subroutine make_BTmat(hecMAT, fstrMAT, iw2, wSL, BTmat, fg_paracon)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:)
    real(kind=kreal), intent(in) :: wSL(:)
    type (hecmwST_local_matrix), intent(out) :: BTmat
    logical, intent(in) :: fg_paracon
    type (hecmwST_local_matrix) :: Tmat
    integer :: ndof, n, i, nnz, l, js, je, j, k, jdof, kk, jj
    real(kind=kreal) :: factor

    ndof=hecMAT%NDOF
    if (fg_paracon) then
      n=hecMAT%NP
    else
      n=hecMAT%N
    endif
    Tmat%nr=n*ndof
    Tmat%nc=Tmat%nr
    if (fg_paracon) then
      Tmat%nnz=fstrMAT%numL_lagrange*ndof-fstrMAT%num_lagrange
    else
      Tmat%nnz=Tmat%nr+fstrMAT%numL_lagrange*ndof-2*fstrMAT%num_lagrange
    endif
    Tmat%ndof=1

    allocate(Tmat%index(0:Tmat%nr))
    allocate(Tmat%item(Tmat%nnz), Tmat%A(Tmat%nnz))
    ! index
    Tmat%index(0)=0
    do i=1,Tmat%nr
      if (iw2(i) > 0) then
        nnz=ndof*(fstrMAT%indexL_lagrange(iw2(i))-fstrMAT%indexL_lagrange(iw2(i)-1))-1
      else
        if (fg_paracon) then
          nnz = 0
        else
          nnz=1
        endif
      endif
      Tmat%index(i)=Tmat%index(i-1)+nnz
    enddo
    if (Tmat%nnz /= Tmat%index(Tmat%nr)) then
      write(0,*) Tmat%nnz, Tmat%index(Tmat%nr)
      stop 'ERROR: Tmat%nnz wrong'
    endif
    ! item and A
    do i=1,Tmat%nr
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
        if (.not. fg_paracon) then
          Tmat%item(l)=i
          Tmat%A(l)=1.0d0
          l=l+1
        endif
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

  subroutine make_BTtmat(hecMAT, fstrMAT, iw2, iwS, wSU, BTtmat, fg_paracon)
    implicit none
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer, intent(in) :: iw2(:), iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    type (hecmwST_local_matrix), intent(out) :: BTtmat
    logical, intent(in) :: fg_paracon
    type (hecmwST_local_matrix) :: Ttmat
    integer :: ndof, n, i, nnz, l, js, je, j, k, idof, jdof, idx

    ndof=hecMAT%NDOF
    if (fg_paracon) then
      n=hecMAT%NP
    else
      n=hecMAT%N
    endif
    Ttmat%nr=n*ndof
    Ttmat%nc=Ttmat%nr
    if (fg_paracon) then
      Ttmat%nnz=fstrMAT%numU_lagrange*ndof-fstrMAT%num_lagrange
    else
      Ttmat%nnz=Ttmat%nr+fstrMAT%numU_lagrange*ndof-2*fstrMAT%num_lagrange
    endif
    Ttmat%ndof=1

    allocate(Ttmat%index(0:Ttmat%nr))
    allocate(Ttmat%item(Ttmat%nnz), Ttmat%A(Ttmat%nnz))
    ! index
    Ttmat%index(0)=0
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        if (iw2(idx) <= 0) then
          if (fstrMAT%num_lagrange == 0) then
            if (fg_paracon) then
              nnz=0
            else
              nnz=1
            endif
          else
            if (fg_paracon) then
              nnz=fstrMAT%indexU_lagrange(i)-fstrMAT%indexU_lagrange(i-1)
            else
              nnz=fstrMAT%indexU_lagrange(i)-fstrMAT%indexU_lagrange(i-1)+1
            endif
          endif
        else
          nnz=0
        endif
        Ttmat%index(idx)=Ttmat%index(idx-1)+nnz
      enddo
    enddo
    if (Ttmat%nnz /= Ttmat%index(Ttmat%nr)) then
      write(0,*) Ttmat%nnz, Ttmat%index(Ttmat%nr)
      stop 'ERROR: Ttmat%nnz wrong'
    endif
    ! item and A
    do i=1,n
      do idof=1,ndof
        idx=(i-1)*ndof+idof
        l=Ttmat%index(idx-1)+1
        if (iw2(idx) <= 0) then
          if (.not. fg_paracon) then
            ! one on diagonal
            Ttmat%item(l)=idx
            Ttmat%A(l)=1.0d0
            l=l+1
          endif
          if (fstrMAT%num_lagrange > 0) then
            ! offdiagonal
            js=fstrMAT%indexU_lagrange(i-1)+1
            je=fstrMAT%indexU_lagrange(i)
            do j=js,je
              k=fstrMAT%itemU_lagrange(j)
              Ttmat%item(l)=iwS(k)
              Ttmat%A(l)=fstrMAT%AU_lagrange((j-1)*ndof+idof)*wSU(k)
              l=l+1
            enddo
          endif
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

  subroutine make_contact_dof_list(hecMAT, fstrMAT, n_contact_dof, contact_dofs)
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    integer(kind=kint), intent(out) :: n_contact_dof
    integer(kind=kint), allocatable, intent(out) :: contact_dofs(:)
    integer(kind=kint) :: ndof, icnt, i, ls, le, l, j, jdof, idx, k, idof
    integer(kind=kint), allocatable :: iw(:)
    real(kind=kreal) :: val
    if (fstrMAT%num_lagrange == 0) then
      n_contact_dof = 0
      return
    endif
    ! lower
    ndof = hecMAT%NDOF
    allocate(iw(hecMAT%NP * ndof))
    icnt = 0
    do i = 1, fstrMAT%num_lagrange
      ls = fstrMAT%indexL_lagrange(i-1)+1
      le = fstrMAT%indexL_lagrange(i)
      do l = ls, le
        j = fstrMAT%itemL_lagrange(l)
        ljdof1: do jdof = 1, ndof
          val = fstrMAT%AL_lagrange((l-1)*ndof+jdof)
          !write(0,*) 'j,jdof,val',j,jdof,val
          if (val == 0.0d0) cycle
          idx = (j-1)*ndof+jdof
          do k = 1, icnt
            if (iw(k) == idx) cycle ljdof1
          enddo
          icnt = icnt + 1
          iw(icnt) = idx
          !write(0,*) 'icnt,idx',icnt,idx
        enddo ljdof1
      enddo
    enddo
    ! upper
    do i = 1, hecMAT%NP
      ls = fstrMAT%indexU_lagrange(i-1)+1
      le = fstrMAT%indexU_lagrange(i)
      do l = ls, le
        j = fstrMAT%itemU_lagrange(l)
        lidof1: do idof = 1, ndof
          val = fstrMAT%AU_lagrange((l-1)*ndof+idof)
          !write(0,*) 'i,idof,val',i,idof,val
          if (val == 0.0d0) cycle
          idx = (i-1)*ndof+idof
          do k = 1, icnt
            if (iw(k) == idx) cycle lidof1
          enddo
          icnt = icnt + 1
          iw(icnt) = idx
          !write(0,*) 'icnt,idx',icnt,idx
        enddo lidof1
      enddo
    enddo
    call quick_sort(iw, 1, icnt)
    allocate(contact_dofs(icnt))
    contact_dofs(1:icnt) = iw(1:icnt)
    n_contact_dof = icnt
    deallocate(iw)
    if (DEBUG >= 2) write(0,*) '  DEBUG2: contact_dofs',contact_dofs(:)
  end subroutine make_contact_dof_list

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
    integer(kind=kint) :: npndof, nndof, i

    npndof=hecMAT%NP*hecMAT%NDOF
    nndof=hecMAT%N*hecMAT%NDOF

    allocate(Btmp(npndof))

    !!! BTtmat*(B+K*(-Bs^-1)*Blag)

    !B2=-Bs^-1*Blag
    Bnew=0.d0
    do i=1,num_lagrange
      Bnew(iwS(i))=wSL(i)*hecMAT%B(npndof+i)
    enddo
    !Btmp=B+K*B2
    call hecmw_matvec(hecMESH, hecMAT, Bnew, Btmp)
    do i=1,nndof
      Btmp(i)=hecMAT%B(i)+Btmp(i)
    enddo
    !B2=BTtmat*Btmp
    call hecmw_localmat_mulvec(BTtmat, Btmp, Bnew)

    deallocate(Btmp)
  end subroutine make_new_b

  subroutine assemble_b_paracon(hecMESH, hecMAT, conMAT, num_lagrange, conCOMM, Btot)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in) :: hecMAT, conMAT
    integer(kind=kint), intent(in) :: num_lagrange
    type (fstrST_contact_comm), intent(in) :: conCOMM
    real(kind=kreal), intent(out) :: Btot(:)
    integer(kind=kint) :: ndof, nndof, npndof, i
    ndof = hecMAT%NDOF
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    do i=1,npndof+num_lagrange
      Btot(i) = conMAT%B(i)
    enddo
    call fstr_contact_comm_reduce_r(conCOMM, Btot, HECMW_SUM)
  end subroutine assemble_b_paracon

  subroutine make_new_b_paracon(hecMESH, hecMAT, conMAT, Btot, hecMESHtmp, hecTKT, BTtmat, BKmat, &
       iwS, wSL, num_lagrange, conCOMM, Bnew)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH, hecMESHtmp
    type(hecmwST_matrix), intent(in) :: hecMAT, conMAT, hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type(hecmwST_local_matrix), intent(in) :: BTtmat, BKmat
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSL(:)
    integer(kind=kint), intent(in) :: num_lagrange
    type (fstrST_contact_comm), intent(in) :: conCOMM
    real(kind=kreal), intent(out) :: Bnew(:)
    real(kind=kreal), allocatable :: Btmp(:)
    integer(kind=kint) :: npndof, nndof, npndof_new, i
    ! SIZE:
    ! Btot <=> hecMAT, hecMESH
    ! Btmp <=> BKmat, hecMESHtmp
    ! Bnew <=> hecTKT, hecMESHtmp
    npndof     = hecMAT%NP*hecMAT%NDOF
    nndof      = hecMAT%N *hecMAT%NDOF
    npndof_new = hecTKT%NP*hecTKT%NDOF
    allocate(Btmp(npndof_new))
    !
    !! BTtmat*(B+K*(-Bs^-1)*Blag)
    !
    ! B2=-Bs^-1*Blag
    Bnew=0.d0
    do i=1,num_lagrange
      !Bnew(iwS(i))=wSL(i)*conMAT%B(npndof+i)
      Bnew(iwS(i))=wSL(i)*Btot(npndof+i)
    enddo
    ! send external contact dof => recv internal contact dof
    call fstr_contact_comm_reduce_r(conCOMM, Bnew, HECMW_SUM)
    ! Btmp=B+K*B2 (including update of Bnew)
    call hecmw_update_3_R(hecMESHtmp, Bnew, hecTKT%NP)
    call hecmw_localmat_mulvec(BKmat, Bnew, Btmp)
    do i=1,nndof
      Btmp(i)=Btot(i)+Btmp(i)
    enddo
    ! B2=BTtmat*Btmp
    call hecmw_update_3_R(hecMESHtmp, Btmp, hecTKT%NP)
    call hecmw_localmat_mulvec(BTtmat, Btmp, Bnew)
    deallocate(Btmp)
  end subroutine make_new_b_paracon

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

  subroutine comp_x_slave_paracon(hecMESH, hecMAT, Btot, hecMESHtmp, hecTKT, BTmat, &
       num_lagrange, iwS, wSL, conCOMM, n_slave, slaves)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH, hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), intent(in) :: hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSL(:)
    type (fstrST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint), intent(in) :: n_slave
    integer(kind=kint), intent(in) :: slaves(:)
    integer(kind=kint) :: ndof, ndof2, npndof, nndof, ilag, islave, i
    real(kind=kreal), allocatable :: Xtmp(:)
    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N  * ndof
    !!
    !! {X} = [BT] {Xp} - [-Bs^-1] {c}
    !!
    ! compute {X} = [BT] {Xp}
    call hecmw_update_3_R(hecMESHtmp, hecTKT%X, hecTKT%NP)
    call hecmw_localmat_mulvec(BTmat, hecTKT%X, hecMAT%X)
    !
    ! compute {Xtmp} = [-Bs^-1] {c}
    allocate(Xtmp(npndof))
    Xtmp(:) = 0.0d0
    do ilag = 1, num_lagrange
      islave = iwS(ilag)
      !Xtmp(islave) = wSL(ilag) * conMAT%B(npndof + ilag)
      Xtmp(islave) = wSL(ilag) * Btot(npndof + ilag)
    enddo
    !
    ! send external contact dof => recv internal contact dof
    call fstr_contact_comm_reduce_r(conCOMM, Xtmp, HECMW_SUM)
    !
    ! {X} = {X} - {Xtmp}
    do i = 1, n_slave
      islave = slaves(i)
      hecMAT%X(islave) = hecMAT%X(islave) - Xtmp(islave)
    enddo
    deallocate(Xtmp)
  end subroutine comp_x_slave_paracon

  subroutine comp_lag(hecMAT, iwS, wSU, num_lagrange)
    implicit none
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint) :: ndof, ndof2, ilag, is, i, idof
    integer(kind=kint) :: js, je, j, jj, ij0, j0, jdof
    real(kind=kreal), pointer :: xlag(:)

    ndof=hecMAT%ndof
    ndof2=ndof*ndof

    xlag=>hecMAT%X(hecMAT%NP*ndof+1:hecMAT%NP*ndof+num_lagrange)

    do ilag=1,num_lagrange
      is=iwS(ilag)
      i=(is-1)/ndof+1
      idof=mod(is-1, ndof)+1
      xlag(ilag)=hecMAT%B(is)
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

  subroutine comp_lag_paracon(hecMESH, hecMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
       n_slave, slaves, num_lagrange, iwS, wSU, conCOMM)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH, hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT, hecTKT
    real(kind=kreal), intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BKmat
    integer(kind=kint), intent(in) :: n_slave, num_lagrange
    integer(kind=kint), intent(in) :: slaves(:), iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    type (fstrST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint) :: ndof, npndof, nndof, npndof_new, i, ilag, islave
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: xlag(:)
    integer(kind=kint), allocatable :: slaves_ndup(:)
    ndof=hecMAT%ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N * ndof
    npndof_new = hecTKT%NP * ndof
    !! <SUMMARY>
    !! {lag} = [Bs^-T] ( {fs} - [Ksp Kss] {u} )
    !!
    ! 1. {Btmp} = [BKmat] {X}
    hecTKT%X(1:nndof) = hecMAT%X(1:nndof)
    call hecmw_update_3_R(hecMESHtmp, hecTKT%X, hecTKT%NP)
    allocate(Btmp(npndof))
    call hecmw_localmat_mulvec(BKmat, hecTKT%X, Btmp)
    !
    ! 2. {Btmp_s} = {fs} - {Btmp_s}
    !do i = 1, nndof
    !  Btmp(i) = Btot(i) - Btmp(i)
    !enddo
    do i = 1, n_slave
      islave = slaves(i)
      Btmp(islave) = Btot(islave) - Btmp(islave)
    enddo
    !
    ! 3. send internal contact dof => recv external contact dof
    !call hecmw_update_3_R(hecMESH, Btmp, hecMAT%NP)
    ! Btmp(nndof+1:npndof) = 0.0d0
    call fstr_contact_comm_bcast_r(conCOMM, Btmp)
    !
    ! 4. {lag} = - [-Bs^-T] {Btmp_s}
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    do ilag = 1, num_lagrange
      islave = iwS(ilag)
      !xlag(ilag)=-wSU(ilag)*(Btot(islave) - Btmp(islave))
      xlag(ilag)=-wSU(ilag)*Btmp(islave)
    enddo
    deallocate(Btmp)
  end subroutine comp_lag_paracon

  subroutine comp_lag_paracon2(hecMESH, hecMAT, conMAT, num_lagrange, iwS, wSU, conCOMM)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), intent(in) :: conMAT
    integer(kind=kint), intent(in) :: num_lagrange
    integer(kind=kint), intent(in) :: iwS(:)
    real(kind=kreal), intent(in) :: wSU(:)
    type (fstrST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint) :: ndof, ndof2, npndof, nndof, i, ilag, irow, idof, js, je, j, jcol, jdof, islave
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: xlag(:)
    ndof=hecMAT%ndof
    ndof2=ndof*ndof
    npndof = hecMAT%NP * ndof
    nndof  = hecMAT%N * ndof
    !! <SUMMARY>
    !! {lag} = [Bs^-T] ( {fs} - [Ksp Kss] {u} )
    !!
    !
    ! {fs} = {hecMAT%B} + {conMAT%B}
    ! [K] {u} = [hecMAT] {u} + [conMAT] {u}
    !
    ! {Btmp} = {hecMAT%B} - [hecMAT] {u}
    allocate(Btmp(npndof))
    call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, Btmp)
    call fstr_contact_comm_bcast_r(conCOMM, Btmp)
    !
    ! {Btmp} = {Btmp} + {conMAT%B} - [conMAT] {u}
    do ilag = 1, num_lagrange
      i = iwS(ilag)
      Btmp(i) = Btmp(i) + conMAT%B(i)
      irow = (i + ndof - 1) / ndof
      idof = i - ndof*(irow-1)
      ! lower
      js = conMAT%indexL(irow-1)+1
      je = conMAT%indexL(irow)
      do j = js, je
        jcol = conMAT%itemL(j)
        do jdof = 1, ndof
          Btmp(i) = Btmp(i) - conMAT%AL(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
        enddo
      enddo
      ! diag
      do jdof = 1, ndof
        Btmp(i) = Btmp(i) - conMAT%D(ndof2*(irow-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(irow-1)+jdof)
      enddo
      ! upper
      js = conMAT%indexU(irow-1)+1
      je = conMAT%indexU(irow)
      do j = js, je
        jcol = conMAT%itemU(j)
        do jdof = 1, ndof
          Btmp(i) = Btmp(i) - conMAT%AU(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
        enddo
      enddo
    enddo
    !
    ! {lag} = - [-Bs^-T] {Btmp_s}
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    do ilag = 1, num_lagrange
      islave = iwS(ilag)
      !xlag(ilag)=-wSU(ilag)*(conMAT%B(islave) - Btmp(islave))
      xlag(ilag)=-wSU(ilag)*Btmp(islave)
    enddo
    deallocate(Btmp)
  end subroutine comp_lag_paracon2

  subroutine check_solution(hecMESH, hecMAT, fstrMAT, Btot, hecMESHtmp, hecTKT, BKmat, &
       conCOMM, n_slave, slaves)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH, hecMESHtmp
    type (hecmwST_matrix), intent(inout) :: hecMAT, hecTKT
    type (fstrST_matrix_contact_lagrange) , intent(in) :: fstrMAT
    real(kind=kreal), target, intent(in) :: Btot(:)
    type (hecmwST_local_matrix), intent(in) :: BKmat
    type (fstrST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint), intent(in) :: n_slave
    integer(kind=kint), intent(in) :: slaves(:)
    integer(kind=kint) :: ndof, nndof, npndof, num_lagrange, i, ls, le, l, j, idof, jdof
    real(kind=kreal), allocatable, target :: r(:)
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal), pointer :: rlag(:), blag(:), xlag(:)
    real(kind=kreal) :: rnrm2, rlagnrm2
    real(kind=kreal) :: bnrm2, blagnrm2
    ndof = hecMAT%NDOF
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    num_lagrange = fstrMAT%num_lagrange
    !
    allocate(r(npndof + num_lagrange))
    r(:) = 0.0d0
    allocate(Btmp(npndof))
    !
    rlag => r(npndof+1:npndof+num_lagrange)
    blag => Btot(npndof+1:npndof+num_lagrange)
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    !
    !! {r}    = {b} - [K] {x} - [Bt] {lag}
    !! {rlag} = {c} - [B] {x}
    !
    ! {r} = {b} - [K] {x}
    hecTKT%X(1:nndof) = hecMAT%X(1:nndof)
    call hecmw_update_3_R(hecMESHtmp, hecTKT%X, hecTKT%NP)
    call hecmw_localmat_mulvec(BKmat, hecTKT%X, Btmp)
    r(1:nndof) = Btot(1:nndof) - Btmp(1:nndof)
    !
    ! {r} = {r} - [Bt] {lag}
    Btmp(:) = 0.0d0
    if (fstrMAT%num_lagrange > 0) then
      do i = 1, hecMAT%NP
        ls = fstrMAT%indexU_lagrange(i-1)+1
        le = fstrMAT%indexU_lagrange(i)
        do l = ls, le
          j = fstrMAT%itemU_lagrange(l)
          do idof = 1, ndof
            Btmp(ndof*(i-1)+idof) = Btmp(ndof*(i-1)+idof) + fstrMAT%AU_lagrange(ndof*(l-1)+idof) * xlag(j)
          enddo
        enddo
      enddo
    endif
    call fstr_contact_comm_reduce_r(conCOMM, Btmp, HECMW_SUM)
    r(1:nndof) = r(1:nndof) - Btmp(1:nndof)
    !
    ! {rlag} = {c} - [B] {x}
    call hecmw_update_3_R(hecMESH, hecMAT%X, hecMAT%NP)
    do i = 1, num_lagrange
      rlag(i) = blag(i)
      ls = fstrMAT%indexL_lagrange(i-1)+1
      le = fstrMAT%indexL_lagrange(i)
      do l = ls, le
        j = fstrMAT%itemL_lagrange(l)
        do jdof = 1, ndof
          rlag(i) = rlag(i) - fstrMAT%AL_lagrange(ndof*(l-1)+jdof) * hecMAT%X(ndof*(j-1)+jdof)
        enddo
      enddo
    enddo
    !
    ! residual in original system
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual---------------------------------------------------------------'
      write(1000+myrank,*) 'size of R',size(R)
      write(1000+myrank,*) 'R: 1-',hecMAT%N*ndof
      write(1000+myrank,*) r(1:hecMAT%N*ndof)
      if (fstrMAT%num_lagrange > 0) then
        write(1000+myrank,*) 'R(lag):',hecMAT%NP*ndof+1,'-',hecMAT%NP*ndof+fstrMAT%num_lagrange
        write(1000+myrank,*) r(hecMAT%NP*ndof+1:hecMAT%NP*ndof+fstrMAT%num_lagrange)
      endif
      if (n_slave > 0) then
        write(1000+myrank,*) 'R(slave):',slaves(:)
        write(1000+myrank,*) r(slaves(:))
      endif
    endif
    !
    call hecmw_InnerProduct_R(hecMESH, NDOF, r, r, rnrm2)
    call hecmw_InnerProduct_R(hecMESH, NDOF, Btot, Btot, bnrm2)
    rlagnrm2 = 0.0d0
    do i = 1, num_lagrange
      rlagnrm2 = rlagnrm2 + rlag(i)*rlag(i)
    enddo
    call hecmw_allreduce_R1(hecMESH, rlagnrm2, HECMW_SUM)
    blagnrm2 = 0.0d0
    do i = 1, num_lagrange
      blagnrm2 = blagnrm2 + blag(i)*blag(i)
    enddo
    call hecmw_allreduce_R1(hecMESH, blagnrm2, HECMW_SUM)
    !
    if (myrank == 0) then
      write(0,*) 'INFO: resid(x,lag,tot)',sqrt(rnrm2),sqrt(rlagnrm2),sqrt(rnrm2+rlagnrm2)
      write(0,*) 'INFO: rhs  (x,lag,tot)',sqrt(bnrm2),sqrt(blagnrm2),sqrt(bnrm2+blagnrm2)
    endif
  end subroutine check_solution

  subroutine check_solution2(hecMESH, hecMAT, conMAT, fstrMAT, n_contact_dof, contact_dofs, &
       conCOMM, n_slave, slaves)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_matrix), intent(in) :: conMAT
    type (fstrST_matrix_contact_lagrange) , intent(in) :: fstrMAT
    integer(kind=kint), intent(in) :: n_contact_dof
    integer(kind=kint), intent(in) :: contact_dofs(:)
    type (fstrST_contact_comm), intent(in) :: conCOMM
    integer(kind=kint), intent(in) :: n_slave
    integer(kind=kint), intent(in) :: slaves(:)
    integer(kind=kint) :: ndof, ndof2, nndof, npndof, num_lagrange
    integer(kind=kint) :: icon, i, irow, idof, js, je, j, jcol, jdof, ls, le, l
    real(kind=kreal), allocatable, target :: r(:)
    real(kind=kreal), allocatable :: r_con(:)
    real(kind=kreal), pointer :: rlag(:), blag(:), xlag(:)
    real(kind=kreal) :: rnrm2, rlagnrm2
    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    nndof = hecMAT%N * ndof
    npndof = hecMAT%NP * ndof
    num_lagrange = fstrMAT%num_lagrange
    !
    allocate(r(npndof + num_lagrange))
    r(:) = 0.0d0
    allocate(r_con(npndof))
    r_con(:) = 0.0d0
    !
    rlag => r(npndof+1:npndof+num_lagrange)
    blag => conMAT%B(npndof+1:npndof+num_lagrange)
    xlag => hecMAT%X(npndof+1:npndof+num_lagrange)
    !
    !! <SUMMARY>
    !! {r}    = {b} - [K] {x} - [Bt] {lag}
    !! {rlag} = {c} - [B] {x}
    !
    ! 1. {r} = {r_org} + {r_con}
    ! {r_org} = {b_org} - [hecMAT] {x}
    ! {r_con} = {b_con} - [conMAT] {x} - [Bt] {lag}
    !
    ! 1.1 {r_org}
    ! {r} = {b} - [K] {x}
    call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r)
    !
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual(original)-----------------------------------------------------'
      write(1000+myrank,*) 'size of R',size(R)
      write(1000+myrank,*) 'R: 1-',hecMAT%N*ndof
      write(1000+myrank,*) r(1:hecMAT%N*ndof)
      if (n_slave > 0) then
        write(1000+myrank,*) 'R(slave):',slaves(:)
        write(1000+myrank,*) r(slaves(:))
      endif
    endif
    !
    ! 1.2 {r_con}
    ! 1.2.1 {r_con} = {b_con} - [conMAT]{x}
    !call hecmw_update_3_R(hecMESH, hecMAT%X, hecMAT%NP)  ! X is already updated
    do i = 1, npndof
      r_con(i) = conMAT%B(i)
    enddo
    do icon = 1, n_contact_dof
      i = contact_dofs(icon)
      irow = (i + ndof - 1) / ndof
      idof = i - ndof*(irow-1)
      ! lower
      js = conMAT%indexL(irow-1)+1
      je = conMAT%indexL(irow)
      do j = js, je
        jcol = conMAT%itemL(j)
        do jdof = 1, ndof
          r_con(i) = r_con(i) - conMAT%AL(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
        enddo
      enddo
      ! diag
      do jdof = 1, ndof
        r_con(i) = r_con(i) - conMAT%D(ndof2*(irow-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(irow-1)+jdof)
      enddo
      ! upper
      js = conMAT%indexU(irow-1)+1
      je = conMAT%indexU(irow)
      do j = js, je
        jcol = conMAT%itemU(j)
        do jdof = 1, ndof
          r_con(i) = r_con(i) - conMAT%AU(ndof2*(j-1)+ndof*(idof-1)+jdof) * hecMAT%X(ndof*(jcol-1)+jdof)
        enddo
      enddo
    enddo
    !
    ! 1.2.2 {r_con} = {r_con} - [Bt] {lag}
    if (fstrMAT%num_lagrange > 0) then
      do i = 1, hecMAT%NP
        ls = fstrMAT%indexU_lagrange(i-1)+1
        le = fstrMAT%indexU_lagrange(i)
        do l = ls, le
          j = fstrMAT%itemU_lagrange(l)
          do idof = 1, ndof
            r_con(ndof*(i-1)+idof) = r_con(ndof*(i-1)+idof) - fstrMAT%AU_lagrange(ndof*(l-1)+idof) * xlag(j)
          enddo
        enddo
      enddo
    endif
    !
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual(contact,local)------------------------------------------------'
      write(1000+myrank,*) 'size of R_con',size(r_con)
      write(1000+myrank,*) 'R_con: 1-',hecMAT%N*ndof
      write(1000+myrank,*) r_con(1:hecMAT%N*ndof)
      write(1000+myrank,*) 'R_con(external): ',hecMAT%N*ndof+1,'-',hecMAT%NP*ndof
      write(1000+myrank,*) r_con(hecMAT%N*ndof+1:hecMAT%NP*ndof)
      if (n_slave > 0) then
        write(1000+myrank,*) 'R_con(slave):',slaves(:)
        write(1000+myrank,*) r_con(slaves(:))
      endif
    endif
    !
    ! 1.2.3 send external part of {r_con}
    call fstr_contact_comm_reduce_r(conCOMM, r_con, HECMW_SUM)
    !
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual(contact,assembled)--------------------------------------------'
      write(1000+myrank,*) 'size of R_con',size(r_con)
      write(1000+myrank,*) 'R_con: 1-',hecMAT%N*ndof
      write(1000+myrank,*) r_con(1:hecMAT%N*ndof)
      if (n_slave > 0) then
        write(1000+myrank,*) 'R_con(slave):',slaves(:)
        write(1000+myrank,*) r_con(slaves(:))
      endif
    endif
    !
    ! 1.3 {r} = {r_org} + {r_con}
    do i = 1,nndof
      r(i) = r(i) + r_con(i)
    enddo
    !
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual(total)--------------------------------------------------------'
      write(1000+myrank,*) 'size of R',size(r)
      write(1000+myrank,*) 'R: 1-',hecMAT%N*ndof
      write(1000+myrank,*) r(1:hecMAT%N*ndof)
      if (n_slave > 0) then
        write(1000+myrank,*) 'R(slave):',slaves(:)
        write(1000+myrank,*) r(slaves(:))
      endif
    endif
    !
    ! 2.  {rlag} = {c} - [B] {x}
    !call hecmw_update_3_R(hecMESH, hecMAT%X, hecMAT%NP)  ! X is already updated
    do i = 1, num_lagrange
      rlag(i) = blag(i)
      ls = fstrMAT%indexL_lagrange(i-1)+1
      le = fstrMAT%indexL_lagrange(i)
      do l = ls, le
        j = fstrMAT%itemL_lagrange(l)
        do jdof = 1, ndof
          rlag(i) = rlag(i) - fstrMAT%AL_lagrange(ndof*(l-1)+jdof) * hecMAT%X(ndof*(j-1)+jdof)
        enddo
      enddo
    enddo
    !
    if (DEBUG >= 3) then
      write(1000+myrank,*) 'Residual(lagrange)-----------------------------------------------------'
      if (fstrMAT%num_lagrange > 0) then
        write(1000+myrank,*) 'R(lag):',hecMAT%NP*ndof+1,'-',hecMAT%NP*ndof+fstrMAT%num_lagrange
        write(1000+myrank,*) r(hecMAT%NP*ndof+1:hecMAT%NP*ndof+fstrMAT%num_lagrange)
      endif
    endif
    !
    call hecmw_InnerProduct_R(hecMESH, NDOF, r, r, rnrm2)
    rlagnrm2 = 0.0d0
    do i = 1, num_lagrange
      rlagnrm2 = rlagnrm2 + rlag(i)*rlag(i)
    enddo
    call hecmw_allreduce_R1(hecMESH, rlagnrm2, HECMW_SUM)
    !
    if (myrank == 0) write(0,*) 'INFO: resid(x,lag,tot)',sqrt(rnrm2),sqrt(rlagnrm2),sqrt(rnrm2+rlagnrm2)
  end subroutine check_solution2

  subroutine mark_slave_dof(BTmat, mark, n_slave, slaves)
    implicit none
    type (hecmwST_local_matrix), intent(in) :: BTmat
    integer(kind=kint), intent(out) :: mark(:)
    integer(kind=kint), intent(out) :: n_slave
    integer(kind=kint), allocatable, intent(out) :: slaves(:)
    integer(kind=kint) :: ndof, ndof2, irow, js, je, j, jcol, idof, jdof, i
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    mark(:) = 0
    do irow = 1, BTmat%nr
      js = BTmat%index(irow-1)+1
      je = BTmat%index(irow)
      do j = js, je
        jcol = BTmat%item(j)
        do idof = 1, ndof
          if (mark(ndof*(irow-1)+idof) == 1) cycle
          do jdof = 1, ndof
            if (irow == jcol .and. idof == jdof) cycle
            if (BTmat%A(ndof2*(j-1)+ndof*(idof-1)+jdof) /= 0.0d0) then
              mark(ndof*(irow-1)+idof) = 1
              exit
            endif
          enddo
        enddo
      enddo
    enddo
    n_slave = 0
    do i = 1, BTmat%nr * ndof
      if (mark(i) /= 0) n_slave = n_slave + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: n_slave',n_slave
    allocate(slaves(n_slave))
    n_slave = 0
    do i = 1, BTmat%nr * ndof
      if (mark(i) /= 0) then
        n_slave = n_slave + 1
        slaves(n_slave) = i
      endif
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: slaves',slaves(:)
  end subroutine mark_slave_dof

  subroutine place_one_on_diag_of_unmarked_dof(BTmat, mark)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTmat
    integer(kind=kint), intent(in) :: mark(:)
    type (hecmwST_local_matrix) :: Imat, Wmat
    integer(kind=kint) :: ndof, ndof2, i, irow, js, je, j, jcol, idof, jdof, n_slave, n_other
    ndof = BTmat%ndof
    ndof2 = ndof*ndof
    ! Imat: unit matrix except for slave dofs
    Imat%nr = BTmat%nr
    Imat%nc = BTmat%nc
    Imat%nnz = Imat%nr
    Imat%ndof = ndof
    allocate(Imat%index(0:Imat%nr))
    allocate(Imat%item(Imat%nnz))
    Imat%index(0) = 0
    do i = 1, Imat%nr
      Imat%index(i) = i
      Imat%item(i) = i
    enddo
    allocate(Imat%A(ndof2 * Imat%nnz))
    Imat%A(:) = 0.0d0
    n_slave = 0
    n_other = 0
    do irow = 1, Imat%nr
      do idof = 1, ndof
        if (mark(ndof*(irow-1)+idof) == 0) then
          Imat%A(ndof2*(irow-1)+ndof*(idof-1)+idof) = 1.0d0
          n_other = n_other + 1
        else
          n_slave = n_slave + 1
        endif
      enddo
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG2: n_slave,n_other',n_slave,n_other
    call hecmw_localmat_add(BTmat, Imat, Wmat)
    call hecmw_localmat_free(BTmat)
    call hecmw_localmat_free(Imat)
    BTmat%nr = Wmat%nr
    BTmat%nc = Wmat%nc
    BTmat%nnz = Wmat%nnz
    BTmat%ndof = Wmat%ndof
    BTmat%index => Wmat%index
    BTmat%item => Wmat%item
    BTmat%A => Wmat%A
  end subroutine place_one_on_diag_of_unmarked_dof

  subroutine place_num_on_diag_of_marked_dof(BTtKTmat, num, mark)
    implicit none
    type (hecmwST_local_matrix), intent(inout) :: BTtKTmat
    real(kind=kreal), intent(in) :: num
    integer(kind=kint), intent(in) :: mark(:)
    integer(kind=kint) :: ndof, ndof2, irow, js, je, j, jcol, idof, jdof
    integer(kind=kint) :: n_error = 0
    ndof = BTtKTmat%ndof
    ndof2 = ndof*ndof
    do irow = 1, BTtKTmat%nr
      js = BTtKTmat%index(irow-1)+1
      je = BTtKTmat%index(irow)
      do j = js, je
        jcol = BTtKTmat%item(j)
        if (irow /= jcol) cycle
        do idof = 1, ndof
          if (mark(ndof*(irow-1)+idof) == 1) then
            if (BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof) /= 0.0d0) &
                 stop 'ERROR: nonzero diag on slave dof of BTtKTmat'
            BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof) = num
          else
            if (BTtKTmat%A(ndof2*(j-1)+ndof*(idof-1)+idof) == 0.0d0) then
              !write(0,*) 'irow,idof',irow,idof
              n_error = n_error+1
            endif
          endif
        enddo
      enddo
    enddo
    if (n_error > 0) then
      write(0,*) 'n_error',n_error
      stop 'ERROR: zero diag on non-slave dof of BTtKTmat'
    endif
  end subroutine place_num_on_diag_of_marked_dof

  subroutine update_comm_table(hecMESH, BTmat, hecMESHtmp, BT_all)
    implicit none
    type (hecmwST_local_mesh), intent(in), target :: hecMESH
    type(hecmwST_local_matrix), intent(in) :: BTmat
    type(hecmwST_local_mesh), intent(inout), target :: hecMESHtmp
    type (hecmwST_local_matrix), intent(out) :: BT_all
    type(hecmwST_local_matrix), allocatable :: BT_exp(:)
    integer(kind=kint) :: n_send, idom, irank, n_curexp, n_oldexp, n_orgexp
    integer(kind=kint) :: idx_0, idx_n, k, knod, n_newexp, j, jnod
    integer(kind=kint), pointer :: cur_export(:), org_export(:)
    integer(kind=kint), pointer :: old_export(:)
    integer(kind=kint), allocatable, target :: old_export_item(:)
    integer(kind=kint), allocatable :: new_export(:)
    integer(kind=kint) :: sendbuf(2), recvbuf(2)
    integer(kind=kint) :: n_oldimp, n_newimp, n_orgimp, i0, n_curimp
    integer(kind=kint), allocatable :: old_import(:)
    integer(kind=kint), pointer :: org_import(:), cur_import(:)
    integer(kind=kint) :: tag
    type (hecmwST_local_matrix) :: BT_imp
    integer(kind=kint) :: nnz
    integer(kind=kint),allocatable :: nnz_imp(:)
    integer(kind=kint), allocatable :: index_imp(:), item_imp(:)
    real(kind=kreal), allocatable :: val_imp(:)
    integer(kind=kint), allocatable :: requests(:)
    integer(kind=kint), allocatable :: statuses(:,:)
    integer(kind=kint) :: nr_imp, jj, ndof2, idx_0_tmp, idx_n_tmp
    integer(kind=kint) :: cnt, ks, ke, iimp, i, ii, ierror
    !!! PREPARATION FOR COMM_TABLE UPDATE
    call copy_mesh(hecMESH, hecMESHtmp)
    allocate(BT_exp(hecMESH%n_neighbor_pe))
    call extract_BT_exp(BTmat, hecMESH, BT_exp)

    !!! UPDATE COMMUNICATION TABLE for Parallel Computation
    allocate(statuses(HECMW_STATUS_SIZE,2*hecMESH%n_neighbor_pe))
    allocate(requests(2*hecMESH%n_neighbor_pe))

    allocate(old_export_item(hecMESH%export_index(hecMESH%n_neighbor_pe)))

    n_send = 0
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      allocate(cur_export(BT_exp(idom)%nnz))
      call extract_cols(BT_exp(idom), cur_export, n_curexp)
      if (DEBUG >= 1) write(0,*) 'DEBUG: extract_cols done'
      n_oldexp = 0
      idx_0 = hecMESH%export_index(idom-1)
      idx_n = hecMESH%export_index(idom)
      n_orgexp = idx_n - idx_0
      org_export => hecMESH%export_item(idx_0+1:idx_n)
      ! check location of old export nodes in original export list
      old_export => old_export_item(idx_0+1:idx_n)
      do k = 1,n_orgexp
        knod = org_export(k)
        if (.not. is_included(cur_export, n_curexp, knod)) then
          n_oldexp = n_oldexp + 1
          old_export(n_oldexp) = k
        end if
      end do
      if (DEBUG >= 1) write(0,*) 'DEBUG: making old_export done'
      ! gather new export nodes at the end of current export list
      call reorder_current_export(cur_export, n_curexp, org_export, n_orgexp, n_newexp, hecMESH%nn_internal)
      if (DEBUG >= 1) write(0,*) 'DEBUG: reorder_current_export done'
      ! check consistency
      if (n_curexp /= n_orgexp - n_oldexp + n_newexp) &
        stop 'ERROR: unknown error(num of export nodes)' !!! ASSERTION
      ! make item_exp from item of BT_exp by converting column id to place in cur_export
      call convert_BT_exp_col_id(BT_exp(idom), cur_export, n_curexp)
      if (DEBUG >= 1) write(0,*) 'DEBUG: convert_BT_expx_col_id done'
      ! add current export list to commtable
      call append_commtable(hecMESHtmp%n_neighbor_pe, hecMESHtmp%export_index, &
        hecMESHtmp%export_item, idom, cur_export, n_curexp)
      if (DEBUG >= 1) write(0,*) 'DEBUG: append_commtable (export) done'
      deallocate(cur_export)
      cur_export => hecMESHtmp%export_item(hecMESHtmp%export_index(idom-1)+1:hecMESHtmp%export_index(idom))
      ! send current export info to neighbor pe
      sendbuf(1) = n_oldexp
      sendbuf(2) = n_newexp
      tag = 1001
      call HECMW_ISEND_INT(sendbuf, 2, irank, tag, &
        hecMESH%MPI_COMM, requests(idom))
      if (n_oldexp > 0) then
        n_send = n_send + 1
        tag = 1002
        call HECMW_ISEND_INT(old_export, n_oldexp, irank, tag, &
          hecMESH%MPI_COMM, requests(hecMESH%n_neighbor_pe+n_send))
      end if
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: isend n_oldexp, n_newexp, old_export done'
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      ! receive current import info from neighbor pe
      tag = 1001
      call HECMW_RECV_INT(recvbuf, 2, irank, tag, &
        hecMESH%MPI_COMM, statuses(:,1))
      n_oldimp = recvbuf(1)
      n_newimp = recvbuf(2)
      if (n_oldimp > 0) then
        allocate(old_import(n_oldimp))
        tag = 1002
        call HECMW_RECV_INT(old_import, n_oldimp, irank, tag, &
          hecMESH%MPI_COMM, statuses(:,1))
      end if
      !
      idx_0 = hecMESH%import_index(idom-1)
      idx_n = hecMESH%import_index(idom)
      n_orgimp = idx_n - idx_0
      org_import => hecMESH%import_item(idx_0+1:idx_n)
      call append_nodes(hecMESHtmp, n_newimp, i0)
      if (DEBUG >= 1) write(0,*) 'DEBUG: append_nodes done'
      n_curimp = n_orgimp - n_oldimp + n_newimp
      allocate(cur_import(n_curimp))
      call make_cur_import(org_import, n_orgimp, old_import, n_oldimp, &
        n_newimp, i0, cur_import)
      if (n_oldimp > 0) deallocate(old_import)
      if (DEBUG >= 1) write(0,*) 'DEBUG: make_cur_import done'
      call append_commtable(hecMESHtmp%n_neighbor_pe, hecMESHtmp%import_index, &
        hecMESHtmp%import_item, idom, cur_import, n_curimp)
      if (DEBUG >= 1) write(0,*) 'DEBUG: append_commtable (import) done'
      deallocate(cur_import)
      !cur_import => hecMESHtmp%import_item(hecMESHtmp%import_index(idom-1)+1:hecMESHtmp%import_index(idom))
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: recv n_oldimp, n_newimp, old_import done'
    call HECMW_Waitall(hecMESH%n_neighbor_pe + n_send, requests, statuses)
    deallocate(old_export_item)

    !!! Send BT_exp & Recv BT_imp; nnz and index
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      sendbuf(1) = BT_exp(idom)%nr
      sendbuf(2) = BT_exp(idom)%nnz
      tag = 1003
      call HECMW_ISEND_INT(sendbuf, 2, irank, tag, &
        hecMESH%MPI_COMM, requests(2*idom-1))
      tag = 1004
      call HECMW_ISEND_INT(BT_exp(idom)%index(0:BT_exp(idom)%nr), BT_exp(idom)%nr+1, &
        irank, tag, hecMESH%MPI_COMM, requests(2*idom))
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: isend BT_exp (nnz and index) done'
    BT_imp%nr = 0
    BT_imp%nc = hecMESHtmp%n_node - hecMESHtmp%nn_internal
    BT_imp%nnz = 0
    allocate(BT_imp%index(0:hecMESH%import_index(hecMESH%n_neighbor_pe)))
    BT_imp%index(0) = 0
    allocate(nnz_imp(hecMESH%n_neighbor_pe))
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      tag = 1003
      call HECMW_RECV_INT(recvbuf, 2, irank, tag, &
        hecMESH%MPI_COMM, statuses(:,1))
      nr_imp = recvbuf(1)
      nnz_imp(idom) = recvbuf(2)
      idx_0 = hecMESH%import_index(idom-1)
      idx_n = hecMESH%import_index(idom)
      if (nr_imp /= idx_n - idx_0) &
        stop 'ERROR: num of rows of BT_imp incorrect' !!! ASSERTION
      BT_imp%nr = BT_imp%nr + nr_imp
      BT_imp%nnz = BT_imp%nnz + nnz_imp(idom)
      allocate(index_imp(0:nr_imp))
      tag = 1004
      call HECMW_RECV_INT(index_imp(0), nr_imp+1, irank, tag, &
        hecMESH%MPI_COMM, statuses(:,1))
      if (index_imp(nr_imp) /= nnz_imp(idom)) then !!! ASSERTION
        if (DEBUG >= 1) write(0,*) 'ERROR: num of nonzero of BT_imp incorrect'
        if (DEBUG >= 1) write(0,*) 'nr_imp, index_imp(nr_imp), nnz_imp', &
          nr_imp, index_imp(nr_imp), nnz_imp(idom)
        stop
      endif
      do j = 1, nr_imp
        jj = hecMESH%import_item(idx_0+j) - hecMESH%nn_internal
        BT_imp%index(jj) = index_imp(j) - index_imp(j-1)
      end do
      deallocate(index_imp)
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: recv BT_imp (nnz and index) done'
    do j = 1, hecMESH%import_index(hecMESH%n_neighbor_pe)
      BT_imp%index(j) = BT_imp%index(j-1) + BT_imp%index(j)
    end do
    if (BT_imp%index(hecMESH%import_index(hecMESH%n_neighbor_pe)) /= BT_imp%nnz) &
      stop 'ERROR: total num of nonzero of BT_imp incorrect' !!! ASSERTION
    ndof2 = BTmat%ndof ** 2
    allocate(BT_imp%item(BT_imp%nnz),BT_imp%A(BT_imp%nnz * ndof2))
    call HECMW_Waitall(hecMESH%n_neighbor_pe * 2, requests, statuses)

    !!! Send BT_exp & Recv BT_imp; item and val
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      tag = 1005
      call HECMW_Isend_INT(BT_exp(idom)%item, BT_exp(idom)%nnz, &
        irank, tag, hecMESH%MPI_COMM, requests(2*idom-1))
      tag = 1006
      call HECMW_Isend_R(BT_exp(idom)%A, BT_exp(idom)%nnz * ndof2, &
        irank, tag, hecMESH%MPI_COMM, requests(2*idom))
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: isend BT_exp (item and val) done'
    do idom = 1,hecMESH%n_neighbor_pe
      irank = hecMESH%neighbor_pe(idom)
      idx_0 = hecMESH%import_index(idom-1)
      idx_n = hecMESH%import_index(idom)
      allocate(item_imp(nnz_imp(idom)))
      tag = 1005
      call HECMW_Recv_INT(item_imp, nnz_imp(idom), &
        irank, tag, hecMESH%MPI_COMM, statuses(:,1))
      allocate(val_imp(nnz_imp(idom) * ndof2))
      tag = 1006
      call HECMW_Recv_R(val_imp, nnz_imp(idom) * ndof2, &
        irank, tag, hecMESH%MPI_COMM, statuses(:,1))

      ! convert column id of item_imp() to local id refering cur_import(:)
      idx_0_tmp = hecMESHtmp%import_index(idom-1)
      idx_n_tmp = hecMESHtmp%import_index(idom)
      cur_import => hecMESHtmp%import_item(idx_0_tmp+1:idx_n_tmp)
      n_curimp = idx_n_tmp - idx_0_tmp
      n_orgimp = idx_n - idx_0
      cnt = 0
      do j = 1, n_orgimp
        jj = hecMESH%import_item(idx_0+j) - hecMESH%nn_internal
        ks = BT_imp%index(jj-1)
        ke = BT_imp%index(jj)
        do k = ks+1, ke
          cnt = cnt + 1
          iimp = item_imp(cnt)
          if (iimp <= 0 .or. n_curimp < iimp) &
            stop 'ERROR: received column id out of range' !!! ASSERTION
          BT_imp%item(k) = cur_import(iimp)
          BT_imp%A((k-1)*ndof2+1:k*ndof2) = val_imp((cnt-1)*ndof2+1:cnt*ndof2)
        end do
      end do
      deallocate(item_imp, val_imp)
    end do
    deallocate(nnz_imp)
    if (DEBUG >= 1) write(0,*) 'DEBUG: recv BT_imp (item and val) done'
    call HECMW_Waitall(hecMESH%n_neighbor_pe * 2, requests, statuses)

    deallocate(statuses)
    deallocate(requests)

    ! make BT_all by combining BTmat and BT_exp
    BT_all%nr = BTmat%nr + BT_imp%nr
    BT_all%nc = BTmat%nc + BT_imp%nc
    BT_all%nnz = BTmat%nnz + BT_imp%nnz
    BT_all%ndof = BTmat%ndof
    allocate(BT_all%index(0:BT_all%nr))
    allocate(BT_all%item(BT_all%nnz))
    allocate(BT_all%A(BT_all%nnz * ndof2))
    BT_all%index(0) = 0
    do i = 1, BTmat%nr
      BT_all%index(i) = BTmat%index(i)
    end do
    do i = 1, BT_imp%nr
      BT_all%index(BTmat%nr+i) = BT_all%index(BTmat%nr+i-1) + &
        BT_imp%index(i) - BT_imp%index(i-1)
    end do
    do i = 1, BTmat%nnz
      BT_all%item(i) = BTmat%item(i)
      BT_all%A((i-1)*ndof2+1:i*ndof2) = BTmat%A((i-1)*ndof2+1:i*ndof2)
    end do
    do i = 1, BT_imp%nnz
      ii = BTmat%nnz + i
      BT_all%item(ii) = BT_imp%item(i)
      BT_all%A((ii-1)*ndof2+1:ii*ndof2) = BT_imp%A((i-1)*ndof2+1:i*ndof2)
    end do
    if (DEBUG >= 1) write(0,*) 'DEBUG: making BT_all done'

    ! free BT_exp(:)
    do idom=1,hecMESH%n_neighbor_pe
      call hecmw_localmat_free(BT_exp(idom))
    end do
    deallocate(BT_exp)
  end subroutine update_comm_table

  subroutine copy_mesh(src, dst, fg_paracon)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: src
    type (hecmwST_local_mesh), intent(out) :: dst
    logical, intent(in), optional :: fg_paracon
    dst%zero          = src%zero
    dst%MPI_COMM      = src%MPI_COMM
    dst%PETOT         = src%PETOT
    dst%PEsmpTOT      = src%PEsmpTOT
    dst%my_rank       = src%my_rank
    dst%n_subdomain   = src%n_subdomain
    dst%n_node        = src%n_node
    dst%nn_internal   = src%nn_internal
    dst%n_elem        = src%n_elem
    dst%ne_internal   = src%ne_internal
    dst%n_elem_type   = src%n_elem_type
    dst%n_dof         = src%n_dof
    dst%n_neighbor_pe = src%n_neighbor_pe
    allocate(dst%neighbor_pe(dst%n_neighbor_pe))
    dst%neighbor_pe(:) = src%neighbor_pe(:)
    allocate(dst%import_index(0:dst%n_neighbor_pe))
    allocate(dst%export_index(0:dst%n_neighbor_pe))
    if (present(fg_paracon) .and. fg_paracon) then
      dst%import_index(:)= src%import_index(:)
      dst%export_index(:)= src%export_index(:)
      allocate(dst%import_item(dst%import_index(dst%n_neighbor_pe)))
      dst%import_item(:) = src%import_item(:)
      allocate(dst%export_item(dst%export_index(dst%n_neighbor_pe)))
      dst%export_item(:) = src%export_item(:)
      allocate(dst%global_node_ID(dst%n_node))
      dst%global_node_ID(1:dst%n_node) = src%global_node_ID(1:dst%n_node)
    else
      dst%import_index(:)= 0
      dst%export_index(:)= 0
      dst%import_item => null()
      dst%export_item => null()
    endif
    allocate(dst%node_ID(2*dst%n_node))
    dst%node_ID(1:2*dst%n_node) = src%node_ID(1:2*dst%n_node)
    allocate(dst%elem_type_item(dst%n_elem_type))
    dst%elem_type_item(:) = src%elem_type_item(:)
    !dst%mpc            = src%mpc
    ! MPC is already set outside of here
    dst%mpc%n_mpc = 0
    dst%node => src%node
  end subroutine copy_mesh

  subroutine free_mesh(hecMESH, fg_paracon)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    logical, intent(in), optional :: fg_paracon
    deallocate(hecMESH%neighbor_pe)
    deallocate(hecMESH%import_index)
    deallocate(hecMESH%export_index)
    if (present(fg_paracon) .and. fg_paracon) then
      deallocate(hecMESH%import_item)
      deallocate(hecMESH%export_item)
      deallocate(hecMESH%global_node_ID)
    else
      if (associated(hecMESH%import_item)) deallocate(hecMESH%import_item)
      if (associated(hecMESH%export_item)) deallocate(hecMESH%export_item)
    endif
    deallocate(hecMESH%node_ID)
    deallocate(hecMESH%elem_type_item)
    !hecMESH%node => null()
  end subroutine free_mesh

  subroutine extract_BT_exp(BTmat, hecMESH, BT_exp)
    implicit none
    type(hecmwST_local_matrix), intent(in) :: BTmat
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_local_matrix), intent(out) :: BT_exp(:)
    integer(kind=kint) :: i, j, k, n, idx_0, idx_n, jrow, ndof2
    ndof2 = BTmat%ndof ** 2
    do i = 1,hecMESH%n_neighbor_pe
      idx_0 = hecMESH%export_index(i-1)
      idx_n = hecMESH%export_index(i)
      BT_exp(i)%nr = idx_n - idx_0
      BT_exp(i)%nc = BTmat%nc
      BT_exp(i)%nnz = 0
      BT_exp(i)%ndof = BTmat%ndof
      allocate(BT_exp(i)%index(0:BT_exp(i)%nr))
      BT_exp(i)%index(0) = 0
      do j = 1,BT_exp(i)%nr
        jrow = hecMESH%export_item(j + idx_0)
        n = BTmat%index(jrow) - BTmat%index(jrow-1)
        BT_exp(i)%nnz = BT_exp(i)%nnz + n
        BT_exp(i)%index(j) = BT_exp(i)%index(j-1) + n
      end do
      allocate(BT_exp(i)%item(BT_exp(i)%nnz))
      allocate(BT_exp(i)%A(BT_exp(i)%nnz * ndof2))
      n = 0
      do j = 1,BT_exp(i)%nr
        jrow = hecMESH%export_item(j + idx_0)
        do k = BTmat%index(jrow-1)+1,BTmat%index(jrow)
          n = n + 1
          !write(0,*) j, jrow, k, n
          BT_exp(i)%item(n) = BTmat%item(k)
          BT_exp(i)%A(ndof2*(n-1)+1:ndof2*n) = BTmat%A(ndof2*(k-1)+1:ndof2*k)
        end do
      end do
    end do
  end subroutine extract_BT_exp

  subroutine extract_cols(BT_exp, cur_export, n_curexp)
    implicit none
    type(hecmwST_local_matrix), intent(in) :: BT_exp
    integer(kind=kint), intent(out) :: cur_export(:)
    integer(kind=kint), intent(out) :: n_curexp
    ! write(0,*) 'BT_exp%item(1:',BT_exp%nnz,')'
    ! write(0,*) BT_exp%item(1:BT_exp%nnz)
    cur_export(1:BT_exp%nnz) = BT_exp%item(1:BT_exp%nnz)
    call quick_sort(cur_export, 1, BT_exp%nnz)
    call unique(cur_export, BT_exp%nnz, n_curexp)
    ! write(0,*) 'cur_export(1:',n_curexp,')'
    ! write(0,*) cur_export(1:n_curexp)
  end subroutine extract_cols

  subroutine reorder_current_export(cur_export, n_curexp, org_export, n_orgexp, n_newexp, nn_internal)
    implicit none
    integer(kind=kint), intent(inout) :: cur_export(:)
    integer(kind=kint), intent(in) :: n_curexp
    integer(kind=kint), intent(in) :: org_export(:)
    integer(kind=kint), intent(in) :: n_orgexp
    integer(kind=kint), intent(out) :: n_newexp
    integer(kind=kint), intent(in) :: nn_internal
    integer(kind=kint), allocatable :: new_export(:)
    integer(kind=kint) :: j, jnod
    n_newexp = 0
    allocate(new_export(n_curexp))
    do j = 1,n_curexp
      jnod = cur_export(j)
      if (jnod > nn_internal) &
        stop 'ERROR: unknown error (jnod)'  !!! ASSERTION
      if (.not. is_included(org_export, n_orgexp, jnod)) then
        n_newexp = n_newexp + 1
        new_export(n_newexp) = jnod
        !write(0,*) 'found new export', jnod
      else if (n_newexp > 0) then
        cur_export(j - n_newexp) = jnod
      end if
    end do
    do j = 1,n_newexp
      cur_export(n_curexp - n_newexp + j) = new_export(j)
    end do
    deallocate(new_export)
    ! write(0,*) 'reordered cur_export(1:',n_curexp,')'
    ! write(0,*) cur_export(1:n_curexp)
  end subroutine reorder_current_export

  subroutine convert_BT_exp_col_id(BT_exp, cur_export, n_curexp)
    implicit none
    type(hecmwST_local_matrix), intent(inout) :: BT_exp
    integer(kind=kint), intent(in) :: cur_export(:)
    integer(kind=kint), intent(in) :: n_curexp
    integer(kind=kint) :: i, icol, j
    logical :: found
    ! make item_exp from item of BT_exp by converting column id to place in cur_export
    do i = 1, BT_exp%nnz
      icol = BT_exp%item(i)
      found = .false.
      do j = 1, n_curexp
        if (icol == cur_export(j)) then
          BT_exp%item(i) = j
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        write(0,*) icol
        stop 'ERROR: unknown error (item not found in cur_export)' !!! ASSERTION
      end if
    end do
  end subroutine convert_BT_exp_col_id

  subroutine append_commtable(n, index, item, idom, cur, ncur)
    implicit none
    integer(kind=kint), intent(in) :: n, idom, ncur
    integer(kind=kint), pointer :: index(:), item(:)
    integer(kind=kint), pointer :: cur(:)
    integer(kind=kint), allocatable :: tmp_index(:), tmp_item(:)
    integer(kind=kint) :: norg, j
    allocate(tmp_index(0:n))
    tmp_index(:) = index(:)
    norg = index(n)
    allocate(tmp_item(norg))
    if (norg > 0) then
      tmp_item(:) = item(:)
      if (associated(item)) deallocate(item)
    end if
    allocate(item(norg + ncur))
    do j = idom,n
      index(j) = index(j) + ncur
    end do
    do j = 1,tmp_index(idom)
      item(j) = tmp_item(j)
    end do
    do j = 1,ncur
      item(tmp_index(idom)+j) = cur(j)
    end do
    do j = tmp_index(idom)+1,tmp_index(n)
      item(j+ncur) = tmp_item(j)
    end do
    deallocate(tmp_index, tmp_item)
  end subroutine append_commtable

  subroutine append_nodes(hecMESHtmp, n_newimp, i0)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESHtmp
    integer(kind=kint), intent(in) :: n_newimp
    integer(kind=kint), intent(out) :: i0
    i0 = hecMESHtmp%n_node
    hecMESHtmp%n_node = hecMESHtmp%n_node + n_newimp
  end subroutine append_nodes

  subroutine make_cur_import(org_import, n_orgimp, old_import, n_oldimp, &
      n_newimp, i0, cur_import)
    implicit none
    integer(kind=kint), intent(in) :: org_import(:), old_import(:)
    integer(kind=kint), intent(in) :: n_orgimp, n_oldimp, n_newimp, i0
    integer(kind=kint), intent(out) :: cur_import(:)
    ! integer(kind=kint), intent(out) :: n_curimp
    integer(kind=kint) :: ndel, i, j
    ndel = 0
    i = 1
    do while (i <= n_orgimp .and. ndel < n_oldimp)
      if (org_import(i) == old_import(ndel+1)) then
        ndel = ndel + 1
      else
        cur_import(i-ndel) = org_import(i)
      endif
      i = i + 1
    enddo
    if (ndel /= n_oldimp) stop 'ERROR: unknown error (ndel)' !!! ASSERTION
    do j = i, n_orgimp
      cur_import(j-ndel) = org_import(j)
    enddo
    i = n_orgimp - ndel
    do j = 1, n_newimp
      cur_import(i + j) = i0+j
    end do
  end subroutine make_cur_import

  recursive subroutine quick_sort(array, id1, id2)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: id1, id2
    integer(kind=kint) :: pivot, center, left, right, tmp
    if (id1 >= id2) return
    center = (id1 + id2) / 2
    pivot = array(center)
    left = id1
    right = id2
    do
      do while (array(left) < pivot)
        left = left + 1
      end do
      do while (pivot < array(right))
        right = right - 1
      end do
      if (left >= right) exit
      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp
      left = left + 1
      right = right - 1
    end do
    if (id1 < left-1) call quick_sort(array, id1, left-1)
    if (right+1 < id2) call quick_sort(array, right+1, id2)
    return
  end subroutine quick_sort

  subroutine unique(array, len, newlen)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(out) :: newlen
    integer(kind=kint) :: i, ndup
    ndup = 0
    do i=2,len
      if (array(i) == array(i - 1 - ndup)) then
        ndup = ndup + 1
      else if (ndup > 0) then
        array(i - ndup) = array(i)
      endif
    end do
    newlen = len - ndup
  end subroutine unique

  function is_included(array, len, ival)
    implicit none
    logical :: is_included
    integer(kind=kint), intent(in) :: array(:)
    integer(kind=kint), intent(in) :: len
    integer(kind=kint), intent(in) :: ival
    integer(kind=kint) :: i
    is_included = .false.
    do i=1,len
      if (array(i) == ival) then
        is_included = .true.
        exit
      end if
    end do
  end function is_included

  !!
  !! Solve without elimination of Lagrange-multipliers
  !!

  subroutine solve_no_eliminate(hecMESH,hecMAT,fstrMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer :: ndof, ndof2, nb_lag, ndofextra
    integer :: i, ls, le, l, j, jb_lag, ib_lag, idof, jdof, ilag, k
    integer :: idx, idx_lag_s, idx_lag_e, ll
    integer, allocatable :: iwUr(:), iwUc(:), iwLr(:), iwLc(:)
    type(hecmwST_matrix) :: hecMATLag
    real(kind=kreal) :: t1

    t1 = hecmw_wtime()
    if (DEBUG >= 1) write(0,*) 'DEBUG: solve_no_eliminate, start', hecmw_wtime()-t1

    call hecmw_mat_init(hecMATLag)

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    nb_lag = (fstrMAT%num_lagrange + 2)/3
    hecMATLag%NDOF = ndof
    hecMATLag%N = hecMAT%N + nb_lag
    hecMATLag%NP = hecMAT%NP + nb_lag
    !write(0,*) 'DEBUG: hecMAT: NDOF,N,NP=',hecMAT%NDOF,hecMAT%N,hecMAT%NP
    !write(0,*) 'DEBUG: hecMATLag: NDOF,N,NP=',hecMATLag%NDOF,hecMATLag%N,hecMATLag%NP

    ndofextra = hecMATLag%N*ndof - hecMAT%N*ndof - fstrMAT%num_lagrange
    if (DEBUG >= 1) write(0,*) 'DEBUG: num_lagrange,nb_lag,ndofextra=',fstrMAT%num_lagrange,nb_lag,ndofextra

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
    allocate(hecMATLag%indexU(0:hecMATLag%NP))
    hecMATLag%indexU(0) = 0
    do i = 1, hecMAT%N
      hecMATLag%indexU(i) = hecMATLag%indexU(i-1) + &
        (hecMAT%indexU(i) - hecMAT%indexU(i-1)) + iwUr(i)
    enddo
    do i = hecMAT%N+1, hecMATLag%N
      hecMATLag%indexU(i) = hecMATLag%indexU(i-1)
    enddo
    do i = hecMATLag%N+1, hecMATLag%NP
      hecMATLag%indexU(i) = hecMATLag%indexU(i-1) + &
        (hecMAT%indexU(i-nb_lag) - hecMAT%indexU(i-1-nb_lag))
    enddo
    hecMATLag%NPU = hecMATLag%indexU(hecMATLag%NP)
    !write(0,*) 'DEBUG: hecMATLag%NPU=',hecMATLag%NPU

    ! Lower: indexL
    allocate(hecMATLag%indexL(0:hecMATLag%NP))
    do i = 0, hecMAT%N
      hecMATLag%indexL(i) = hecMAT%indexL(i)
    enddo
    do i = hecMAT%N+1, hecMATLag%N
      hecMATLag%indexL(i) = hecMATLag%indexL(i-1) + iwLr(i-hecMAT%N)
    enddo
    do i = hecMATLag%N+1, hecMATLag%NP
      hecMATLag%indexL(i) = hecMATLag%indexL(i-1) + &
        (hecMAT%indexL(i-nb_lag) - hecMAT%indexL(i-1-nb_lag))
    enddo
    hecMATLag%NPL = hecMATLag%indexL(hecMATLag%NP)
    !write(0,*) 'DEBUG: hecMATLag%NPL=',hecMATLag%NPL

    ! Upper: itemU and AU
    allocate(hecMATLag%itemU(hecMATLag%NPU))
    allocate(hecMATLag%AU(hecMATLag%NPU*ndof2))
    hecMATLag%AU = 0.d0
    do i = 1, hecMAT%N
      idx = hecMATLag%indexU(i-1)+1
      ! copy from hecMAT; internal
      ls = hecMAT%indexU(i-1)+1
      le = hecMAT%indexU(i)
      do l=ls,le
        ll = hecMAT%itemU(l)
        if (ll > hecMAT%N) cycle  ! skip external
        hecMATLag%itemU(idx) = ll
        hecMATLag%AU((idx-1)*ndof2+1:idx*ndof2)=hecMAT%AU((l-1)*ndof2+1:l*ndof2)
        idx = idx + 1
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
      idx_lag_s = idx
      do j=1,nb_lag
        if (iwUc(j) > 0) then
          hecMATLag%itemU(idx) = hecMAT%N + j
          idx = idx + 1
        endif
      enddo
      idx_lag_e = idx - 1
      ! Lag. AU
      ls=fstrMAT%indexU_lagrange(i-1)+1
      le=fstrMAT%indexU_lagrange(i)
      do l=ls,le
        j=fstrMAT%itemU_lagrange(l)
        jb_lag = (j+2)/3
        jdof = j - (jb_lag - 1)*ndof
        do k = idx_lag_s, idx_lag_e
          if (hecMATLag%itemU(k) < hecMAT%N + jb_lag) cycle
          if (hecMATLag%itemU(k) > hecMAT%N + jb_lag) cycle
          !if (hecMATLag%itemU(k) /= hecMAT%N + jb_lag) stop 'ERROR itemU jb_lag'
          do idof = 1, ndof
            hecMATLag%AU((k-1)*ndof2+(idof-1)*ndof+jdof) = &
              fstrMAT%AU_lagrange((l-1)*ndof+idof)
          enddo
        enddo
      enddo
      ! copy from hecMAT; externl
      ls = hecMAT%indexU(i-1)+1
      le = hecMAT%indexU(i)
      do l=ls,le
        ll = hecMAT%itemU(l)
        if (ll <= hecMAT%N) cycle  ! skip internal
        hecMATLag%itemU(idx) = ll + nb_lag
        hecMATLag%AU((idx-1)*ndof2+1:idx*ndof2)=hecMAT%AU((l-1)*ndof2+1:l*ndof2)
        idx = idx + 1
      enddo
      if (idx /= hecMATLag%indexU(i)+1) stop 'ERROR idx indexU'
    enddo
    do i = hecMAT%N, hecMAT%NP
      idx = hecMATLag%indexU(i+nb_lag-1)+1
      ls=hecMAT%indexU(i-1)+1
      le=hecMAT%indexU(i)
      do l=ls,le
        ll = hecMAT%itemU(l)
        hecMATLag%itemU(idx) = ll + nb_lag
        hecMATLag%AU((idx-1)*ndof2+1:idx*ndof2)=hecMAT%AU((l-1)*ndof2+1:l*ndof2)
        idx = idx + 1
      enddo
    enddo

    ! Lower: itemL and AL
    allocate(hecMATLag%itemL(hecMATLag%NPL))
    allocate(hecMATLag%AL(hecMATLag%NPL*ndof2))
    hecMATLag%AL = 0.d0
    do i = 1, hecMAT%indexL(hecMAT%N)
      hecMATLag%itemL(i) = hecMAT%itemL(i)
    enddo
    do i = 1, hecMAT%indexL(hecMAT%N)*ndof2
      hecMATLag%AL(i) = hecMAT%AL(i)
    enddo
    ! Lag. itemL
    idx = hecMAT%indexL(hecMAT%N) + 1
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
      idx_lag_s = idx
      do j=1,hecMAT%N
        if (iwLc(j) > 0) then
          hecMATLag%itemL(idx) = j
          idx = idx + 1
        endif
      enddo
      idx_lag_e = idx - 1
      if (idx /= hecMATLag%indexL(hecMAT%N + ib_lag)+1) then
        stop 'ERROR idx indexL'
      endif
      ! Lag. AL
      do idof = 1, ndof
        ilag = (ib_lag-1)*ndof + idof
        if (ilag  > fstrMAT%num_lagrange) exit
        ls=fstrMAT%indexL_lagrange(ilag-1)+1
        le=fstrMAT%indexL_lagrange(ilag)
        do l=ls,le
          j=fstrMAT%itemL_lagrange(l)
          do k = idx_lag_s, idx_lag_e
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
    do i = hecMAT%N+1, hecMAT%NP
      idx = hecMATLag%indexL(i+nb_lag-1)+1
      ls=hecMAT%indexL(i-1)+1
      le=hecMAT%indexL(i)
      do l=ls,le
        ll = hecMAT%itemL(l)
        if (ll <= hecMAT%N) then
          hecMATLag%itemL(idx) = ll
        else
          hecMATLag%itemL(idx) = ll + nb_lag
        endif
        hecMATLag%AL((idx-1)*ndof2+1:idx*ndof2)=hecMAT%AL((l-1)*ndof2+1:l*ndof2)
        idx = idx + 1
      enddo
    enddo

    deallocate(iwUr, iwUc, iwLr, iwLc)

    allocate(hecMATLag%D(hecMATLag%NP*ndof2))
    hecMATLag%D = 0.d0
    ! internal
    do i = 1, hecMAT%N*ndof2
      hecMATLag%D(i) = hecMAT%D(i)
    enddo
    ! Lag.
    do idof = ndof - ndofextra + 1, ndof
      hecMATLag%D((hecMATLag%N-1)*ndof2 + (idof-1)*ndof + idof) = 1.d0
    enddo
    ! external
    do i = hecMAT%N*ndof2+1, hecMAT%NP*ndof2
      hecMATLag%D(i + nb_lag*ndof2) = hecMAT%D(i)
    enddo

    allocate(hecMATLag%B(hecMATLag%NP*ndof))
    allocate(hecMATLag%X(hecMATLag%NP*ndof))
    hecMATLag%B = 0.d0
    hecMATLag%X = 0.d0

    ! internal
    do i = 1, hecMAT%N*ndof
      hecMATLag%B(i) = hecMAT%B(i)
    enddo
    ! Lag.
    do i = 1, fstrMAT%num_lagrange
      hecMATLag%B(hecMAT%N*ndof + i) = hecMAT%B(hecMAT%NP*ndof + i)
    enddo
    ! external
    do i = hecMAT%N*ndof+1, hecMAT%NP*ndof
      hecMATlag%B(i + nb_lag*ndof) = hecMAT%B(i)
    enddo

    hecMATLag%Iarray=hecMAT%Iarray
    hecMATLag%Rarray=hecMAT%Rarray

    if (DEBUG >= 1) write(0,*) 'DEBUG: made hecMATLag', hecmw_wtime()-t1

    if (hecMESH%n_neighbor_pe > 0) then
      do i = 1, hecMESH%import_index(hecMESH%n_neighbor_pe)
        hecMESH%import_item(i) = hecMESH%import_item(i) + nb_lag
      enddo
    endif

    call hecmw_solve(hecMESH, hecMATLag)

    hecMAT%Iarray=hecMATLag%Iarray
    hecMAT%Rarray=hecMATLag%Rarray

    if (hecMESH%n_neighbor_pe > 0) then
      do i = 1, hecMESH%import_index(hecMESH%n_neighbor_pe)
        hecMESH%import_item(i) = hecMESH%import_item(i) - nb_lag
      enddo
    endif

    hecMAT%X = 0.d0
    ! internal
    do i = 1, hecMAT%N*ndof
      hecMAT%X(i) = hecMATLag%X(i)
    enddo
    ! external
    do i = hecMAT%N*ndof+1, hecMAT%NP*ndof
      hecMAT%X(i) = hecMATLag%X(i + nb_lag*ndof)
    enddo
    ! Lag.
    do i = 1, fstrMAT%num_lagrange
      hecMAT%X(hecMAT%NP*ndof + i) = hecMATLag%X(hecMAT%N*ndof + i)
    enddo

    call hecmw_mat_finalize(hecMATLag)

    if (DEBUG >= 1) write(0,*) 'DEBUG: solve_no_eliminate end', hecmw_wtime()-t1
  end subroutine solve_no_eliminate

end module m_solve_LINEQ_iter_contact
