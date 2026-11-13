!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to solve sparse system of
!> \linear equitions in the case of contact analysis using standard
!> \Lagrange multiplier algorithm

module m_solve_LINEQ_contact

  use hecmw_util
  use m_solve_LINEQ_MKL_contact
  use m_solve_LINEQ_direct_serial_lag
  use m_solve_LINEQ_MUMPS_contact
  use m_solve_LINEQ_mf_contact
  use m_solve_LINEQ_contact_elim
  use hecmw_matrix_misc
  use m_hecmw_comm_f
  use hecmw_ebc_defer
  use hecmw_solver_las

  implicit none

  private
  public :: solve_LINEQ_contact_init
  public :: solve_LINEQ_contact

contains

  !> \brief This subroutine
  subroutine solve_LINEQ_contact_init(hecMESH,hecMAT,is_sym)
    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    logical :: is_sym
    integer(kind=kint) :: contact_elim, solver_type, loglevel, myrank
    logical :: writelog

    contact_elim = hecmw_mat_get_contact_elim(hecMAT)
    solver_type = hecmw_mat_get_solver_type(hecMAT)
    loglevel = max(0, hecmw_mat_get_loglevel(hecMAT)) ! LOGLEVEL is -1 when not given
    myrank = hecmw_comm_get_rank()
    writelog = (myrank==0 .and. loglevel>0)

    if( contact_elim==0 )then ! auto
      if( solver_type==1 )then ! iterative
        contact_elim=1
      else ! direct
        contact_elim=-1
      endif
      call hecmw_mat_set_contact_elim(hecMAT,contact_elim)
    endif

    ! without a Lagrange multiplier anywhere (ALAGRANGE, or SLAGRANGE before any activation)
    ! there is no saddle point part to keep, so the single-process limitation of the
    ! no-elimination DIRECTmf path does not apply: the elimination path solves the plain
    ! distributed matrix as is. The stored CONTACT_ELIM is not changed, so the limitation
    ! still guards the solves after multipliers appear.
    if( contact_elim/=1 .and. solver_type==7 .and. hecmw_comm_get_size() > 1 )then
      if( .not. any_lagrange(hecMESH, hecMAT) ) contact_elim = 1
    endif

    if( contact_elim==1 )then
      if( writelog ) write(*,*) 'solve contact with elimination'
      call solve_LINEQ_contact_elim_init(hecMESH,hecMAT,is_sym)
    else
      if( writelog ) write(*,*) 'solve contact without elimination'
      if( solver_type==1 )then
        write(*,*) 'ERROR: iterative solver without elimination not available in contact analysis'
        call hecmw_abort(hecmw_comm_get_comm())
      elseif( solver_type==2 )then
        call solve_LINEQ_serial_lag_hecmw_init(hecMAT,is_sym)
      else if( solver_type==3 )then
        call solve_LINEQ_MKL_contact_init(hecMESH,is_sym)
      elseif( solver_type==5 ) then
        call solve_LINEQ_mumps_contact_init(hecMESH,hecMAT,is_sym)
      elseif( solver_type==7 ) then
        if( hecmw_comm_get_size() > 1) then
          write(*,*) 'ERROR: !SOLVER,METHOD=DIRECTmf not available in parallel contact analysis',&
              ' without elimination; specify CONTACT_ELIM=1 or use MUMPS or DIRECTmkl instead'
          call hecmw_abort(hecmw_comm_get_comm())
        endif
        call solve_LINEQ_mf_contact_init(hecMESH,hecMAT,is_sym)
      else
        write(*,*) 'ERROR: specified solver not available in contact analysis without elimination;',&
            ' please use MUMPS or DIRECTmkl instead'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
    endif
  end subroutine solve_LINEQ_contact_init


  !> \brief This subroutine
  subroutine solve_LINEQ_contact(hecMESH,hecMAT,hecEBC,istat,rf,is_contact_active)

    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_ebc)                       :: hecEBC         !< prescribed displacements imposed after the MPC processing
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), optional               :: rf
    logical                                  :: is_contact_active

    real(kind=kreal)                         :: factor
    real(kind=kreal) :: t1, t2, resid
    integer(kind=kint) :: ndof
    integer(kind=kint) :: contact_elim, solver_type
    logical :: conMAT_added

    contact_elim = hecmw_mat_get_contact_elim(hecMAT)
    solver_type = hecmw_mat_get_solver_type(hecMAT)

    ! the same bypass as the init: no multiplier anywhere means no saddle point part
    if( contact_elim/=1 .and. solver_type==7 .and. hecmw_comm_get_size() > 1 )then
      if( .not. any_lagrange(hecMESH, hecMAT) ) contact_elim = 1
    endif

    factor = 1.0d0
    if( present(rf) )factor = rf

    t1 = hecmw_wtime()

    istat = 0
    conMAT_added = .false.
    ! the direct-solver branches below report failure through istat only and never
    ! touch these flags; clear them so a stale verdict from an earlier solve is not
    ! mistaken for this one's. The elimination branch overwrites them with the real result.
    call hecmw_mat_set_flag_converged(hecMAT, 1)
    call hecmw_mat_set_flag_diverged(hecMAT, 0)
    if( contact_elim==1 )then
      call solve_LINEQ_contact_elim(hecMESH,hecMAT,hecEBC,istat,is_contact_active)
    else
      if( solver_type==1 )then
        write(*,*) 'ERROR: iterative solver without elimination not available in contact analysis'
        call hecmw_abort(hecmw_comm_get_comm())
      elseif( solver_type==2 )then
        if( hecmw_comm_get_size() > 1) then
          write(*,*) 'ERROR: !SOLVER,METHOD=DIRECT not available in parallel contact analysis;',&
              ' please use MUMPS or DIRECTmkl instead'
          call hecmw_abort(hecmw_comm_get_comm())
        else
          call add_conMAT_to_hecMAT(hecMAT)
          conMAT_added = .true.
          call solve_LINEQ_serial_lag_hecmw(hecMESH,hecMAT,hecEBC)
        endif
      elseif( solver_type==3 )then
        if( hecmw_comm_get_size() == 1 )then
          call add_conMAT_to_hecMAT(hecMAT)
          conMAT_added = .true.
        endif
        call solve_LINEQ_MKL_contact(hecMESH,hecMAT,hecEBC,istat)
      elseif( solver_type==5 ) then
        call solve_LINEQ_mumps_contact(hecMESH,hecMAT,hecEBC,istat)
      elseif( solver_type==7 ) then
        if( hecmw_comm_get_size() > 1) then
          write(*,*) 'ERROR: !SOLVER,METHOD=DIRECTmf not available in parallel contact analysis',&
              ' without elimination; specify CONTACT_ELIM=1 or use MUMPS or DIRECTmkl instead'
          call hecmw_abort(hecmw_comm_get_comm())
        else
          call add_conMAT_to_hecMAT(hecMAT)
          conMAT_added = .true.
          call solve_LINEQ_mf_contact(hecMESH,hecMAT,hecEBC,istat)
        endif
      else
        write(*,*) 'ERROR: specified solver not available in contact analysis without elimination;',&
            ' please use MUMPS or DIRECTmkl instead'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
    endif

    ndof = hecMAT%NDOF
    call hecmw_update_R(hecMESH,hecMAT%X,hecMESH%n_node, ndof)

    ! the final residual of every direct solver is printed here, as hecmw_solve does on the
    ! other paths; the solvers print at most the history of their own refinement
    if( contact_elim/=1 .and. hecmw_mat_get_iterlog(hecMAT) > 0 )then
      ! hecMAT%conmat keeps its values after it is added to hecMAT, so it is left out then
      resid = rel_resid_lagrange(hecMESH,hecMAT,.not. conMAT_added)
      if( hecmw_comm_get_rank()==0 ) write(*,"(a,1pe12.5)") '### Relative residual =', resid
    endif

    t2 = hecmw_wtime()
    if (hecmw_mat_get_timelog(hecMAT) .ge. 1) then
      if ( hecmw_comm_get_rank() ==0) write(*,*) ' solve time :', t2 - t1
    endif

    hecMAT%X=factor*hecMAT%X

  end subroutine solve_LINEQ_contact


  !> a Lagrange multiplier exists on some rank (the saddle point part is present)
  logical function any_lagrange(hecMESH, hecMAT)
    type (hecmwST_local_mesh)      :: hecMESH
    type (hecmwST_matrix)          :: hecMAT
    integer(kind=kint) :: n

    n = hecMAT%lagmat%num_lagrange
    call hecmw_allreduce_I1(hecMESH, n, hecmw_max)
    any_lagrange = n > 0
  end function any_lagrange


  subroutine add_conMAT_to_hecMAT(hecMAT)
    type (hecmwST_matrix),          intent(inout) :: hecMAT


    integer(kind=kint) :: ndof,ndof2,i

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof

    do i=1,hecMAT%NP*ndof + hecMAT%lagmat%num_lagrange
      hecMAT%B(i) = hecMAT%B(i) + hecMAT%conmat%B(i)
    enddo

    do i=1,hecMAT%NP*ndof2
      hecMAT%D(i) = hecMAT%D(i) + hecMAT%conmat%D(i)
    enddo

    do i=1,hecMAT%NPL*ndof2
      hecMAT%AL(i) = hecMAT%AL(i) + hecMAT%conmat%AL(i)
    enddo

    do i=1,hecMAT%NPU*ndof2
      hecMAT%AU(i) = hecMAT%AU(i) + hecMAT%conmat%AU(i)
    enddo
  end subroutine add_conMAT_to_hecMAT


  !> relative residual of the system with the Lagrange multiplier rows, measured on the matrices
  !> the solver input is built from; with_conmat when hecMAT%conmat has not been added to hecMAT
  function rel_resid_lagrange(hecMESH, hecMAT, with_conmat) result(resid)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix),     intent(in) :: hecMAT
    logical,                   intent(in) :: with_conmat
    real(kind=kreal) :: resid

    integer(kind=kint) :: ndof, ndof2, nndof, npndof, i, idof, jdof, j, k, l, irow
    real(kind=kreal), allocatable :: y(:), yc(:), bc(:)
    real(kind=kreal), pointer :: blag(:)
    real(kind=kreal) :: bi, rl, sums(2)

    ndof = hecMAT%NDOF
    ndof2 = ndof*ndof
    nndof = hecMAT%N*ndof
    npndof = hecMAT%NP*ndof

    allocate(y(npndof), yc(npndof), bc(npndof))
    call hecmw_matvec(hecMESH, hecMAT, hecMAT%X, y)

    ! the contact part and the multiplier columns also have entries on the external rows,
    ! which are this rank's share of the rows of the owner rank
    yc(:) = 0.0d0
    bc(:) = 0.0d0
    if( with_conmat )then
      do i = 1, hecMAT%NP
        do idof = 1, ndof
          irow = ndof*(i-1)+idof
          bc(irow) = hecMAT%conmat%B(irow)
          do jdof = 1, ndof
            yc(irow) = yc(irow) + hecMAT%conmat%D(ndof2*(i-1)+ndof*(idof-1)+jdof)*hecMAT%X(ndof*(i-1)+jdof)
          enddo
          do k = hecMAT%conmat%indexL(i-1)+1, hecMAT%conmat%indexL(i)
            j = hecMAT%conmat%itemL(k)
            do jdof = 1, ndof
              yc(irow) = yc(irow) + hecMAT%conmat%AL(ndof2*(k-1)+ndof*(idof-1)+jdof)*hecMAT%X(ndof*(j-1)+jdof)
            enddo
          enddo
          do k = hecMAT%conmat%indexU(i-1)+1, hecMAT%conmat%indexU(i)
            j = hecMAT%conmat%itemU(k)
            do jdof = 1, ndof
              yc(irow) = yc(irow) + hecMAT%conmat%AU(ndof2*(k-1)+ndof*(idof-1)+jdof)*hecMAT%X(ndof*(j-1)+jdof)
            enddo
          enddo
        enddo
      enddo
      call hecmw_assemble_R(hecMESH, bc, hecMAT%NP, ndof)
    endif
    if( hecMAT%lagmat%num_lagrange > 0 )then
      do i = 1, hecMAT%NP
        do l = hecMAT%lagmat%indexU_lagrange(i-1)+1, hecMAT%lagmat%indexU_lagrange(i)
          j = hecMAT%lagmat%itemU_lagrange(l)
          do idof = 1, ndof
            irow = ndof*(i-1)+idof
            yc(irow) = yc(irow) + hecMAT%lagmat%AU_lagrange(ndof*(l-1)+idof)*hecMAT%X(npndof+j)
          enddo
        enddo
      enddo
    endif
    call hecmw_assemble_R(hecMESH, yc, hecMAT%NP, ndof)

    ! near convergence the parts of the right-hand side in hecMAT and conmat nearly cancel, so
    ! they are summed before the products are subtracted, as in the solver input
    sums(:) = 0.0d0
    do i = 1, nndof
      bi = hecMAT%B(i) + bc(i)
      sums(1) = sums(1) + (bi - (y(i) + yc(i)))**2
      sums(2) = sums(2) + bi**2
    enddo

    ! the solvers take the right-hand side of the multiplier rows from conmat when it is separate
    if( with_conmat )then
      blag => hecMAT%conmat%B
    else
      blag => hecMAT%B
    endif
    do i = 1, hecMAT%lagmat%num_lagrange
      rl = blag(npndof+i)
      do l = hecMAT%lagmat%indexL_lagrange(i-1)+1, hecMAT%lagmat%indexL_lagrange(i)
        j = hecMAT%lagmat%itemL_lagrange(l)
        do jdof = 1, ndof
          rl = rl - hecMAT%lagmat%AL_lagrange(ndof*(l-1)+jdof)*hecMAT%X(ndof*(j-1)+jdof)
        enddo
      enddo
      sums(1) = sums(1) + rl**2
      sums(2) = sums(2) + blag(npndof+i)**2
    enddo

    call hecmw_allreduce_R(hecMESH, sums, 2, hecmw_sum)
    if( .not. (sums(2) > 0.0d0) ) sums(2) = 1.0d0
    resid = sqrt(sums(1)/sums(2))

    deallocate(y, yc, bc)
  end function rel_resid_lagrange


end module m_solve_LINEQ_contact
