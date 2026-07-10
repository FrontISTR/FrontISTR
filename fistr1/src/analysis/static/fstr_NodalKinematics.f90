!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Finite-rotation nodal kinematics for NLGEOM.
!>
!> This module manages the per-node geometric state (nodal frame / triad plus a
!> drilling scalar) needed when the Newton solution increment cannot be applied as
!> a simple vector addition to the nodal degrees of freedom.
!>
!> This path is used by elastic MITC4 shell (741) under Total/Updated
!> Lagrangian kinematics.  Element-level rotation algebra lives in
!> m_static_LIB_shell; this module stores and advances the nodal frame state.
module m_fstr_NodalKinematics
  use m_fstr
  use m_fstr_FiniteRotationKinematics, only: fstr_uses_finite_rotation_kinematics
  implicit none

  private

  public :: fstr_ensure_finite_rotation_state
  public :: fstr_begin_nodal_kinematics_step
  public :: fstr_apply_solution_increment
  public :: fstr_commit_solution_increment
  public :: fstr_get_shell_trial_directors
  public :: fstr_get_shell_trial_drilling
  public :: fstr_get_shell_current_directors
  public :: fstr_get_shell_reference_directors

contains

  !> Build the per-node reference frames once, by averaging element shell triads
  !> at shared nodes. Already-initialized nodes are left untouched, so repeated
  !> calls are idempotent.
  subroutine fstr_ensure_finite_rotation_state( hecMESH, fstrSOLID, ndof )
    use elementInfo, only: fe_mitc4_shell
    use m_static_LIB_shell, only: ShellOrthonormalizeTriad
    implicit none

    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    integer(kind=kint), intent(in)        :: ndof

    integer(kind=kint) :: itype, is, iE, ic_type, icel, iiS, nn, j, node_id
    integer(kind=kint), allocatable :: node_mode(:), node_count(:)
    real(kind=kreal), allocatable :: director_sum(:,:), tangent_sum(:,:)
    real(kind=kreal) :: ecoord(3, 8), triad(3, 3), trial(3, 3)
    real(kind=kreal) :: director(3), tangent(3), ref_axis(3), normv, proj

    if( .not. fstrSOLID%has_finite_rotation_kinematics ) return
    if( fstrSOLID%finite_rotation_state_ready ) return
    if( .not. associated(fstrSOLID%shell_rot_state) ) return

    allocate( director_sum(3, hecMESH%n_node) )
    allocate( tangent_sum(3, hecMESH%n_node) )
    allocate( node_mode(hecMESH%n_node) )
    allocate( node_count(hecMESH%n_node) )
    director_sum(:, :) = 0.0D0
    tangent_sum(:, :) = 0.0D0
    node_mode(:) = 0
    node_count(:) = 0

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype)
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= fe_mitc4_shell ) cycle
      do icel = is, iE
        iiS = hecMESH%elem_node_index(icel-1)
        nn = hecMESH%elem_node_index(icel) - iiS
        if( .not. associated( fstrSOLID%elements(icel)%gausses ) ) cycle
        if( .not. fstr_uses_finite_rotation_kinematics( ic_type, nn, &
            fstrSOLID%elements(icel)%gausses(1)%pMaterial ) ) cycle

        do j = 1, min(nn, 8)
          node_id = hecMESH%elem_node_item(iiS+j)
          ecoord(1:3, j) = hecMESH%node(3*node_id-2:3*node_id)
        end do

        if( ndof >= 6 ) then
          do j = 1, nn
            call fstr_reference_shell_triad( nn, ecoord(1:3, 1:nn), j, triad )
            node_id = hecMESH%elem_node_item(iiS+j)
            if( node_id <= 0 .or. node_id > hecMESH%n_node ) cycle
            if( fstrSOLID%shell_rot_state(node_id) /= 0 ) cycle
            ! average shell triads at shared nodes
            director_sum(1:3, node_id) = director_sum(1:3, node_id) + triad(1:3, 3)
            tangent_sum(1:3, node_id) = tangent_sum(1:3, node_id) + triad(1:3, 1)
            node_count(node_id) = node_count(node_id) + 1
            node_mode(node_id) = 1
          end do
        endif
      end do
    end do

    do node_id = 1, hecMESH%n_node
      if( node_count(node_id) <= 0 ) cycle
      if( fstrSOLID%shell_rot_state(node_id) /= 0 ) cycle

      director(1:3) = director_sum(1:3, node_id)
      normv = dsqrt( dot_product( director(1:3), director(1:3) ) )
      if( normv <= 1.0D-14 ) then
        call fstr_set_identity_triad( triad )
        call fstr_store_shell_triad_node( fstrSOLID, node_id, triad, node_mode(node_id) )
        cycle
      endif
      director(1:3) = director(1:3)/normv

      tangent(1:3) = tangent_sum(1:3, node_id)
      proj = dot_product( tangent(1:3), director(1:3) )
      tangent(1:3) = tangent(1:3) - proj*director(1:3)
      normv = dsqrt( dot_product( tangent(1:3), tangent(1:3) ) )
      if( normv <= 1.0D-14 ) then
        ref_axis(1:3) = (/ 1.0D0, 0.0D0, 0.0D0 /)
        if( dabs(director(1)) > 0.9D0 ) ref_axis(1:3) = (/ 0.0D0, 1.0D0, 0.0D0 /)
        proj = dot_product( ref_axis(1:3), director(1:3) )
        tangent(1:3) = ref_axis(1:3) - proj*director(1:3)
        normv = dsqrt( dot_product( tangent(1:3), tangent(1:3) ) )
      endif
      tangent(1:3) = tangent(1:3)/normv

      call fstr_set_identity_triad( trial )
      trial(1:3, 1) = tangent(1:3)
      trial(1:3, 3) = director(1:3)
      call ShellOrthonormalizeTriad( trial, triad )
      call fstr_store_shell_triad_node( fstrSOLID, node_id, triad, node_mode(node_id) )
    end do

    deallocate( director_sum )
    deallocate( tangent_sum )
    deallocate( node_mode )
    deallocate( node_count )
    fstrSOLID%finite_rotation_state_ready = .true.

  end subroutine fstr_ensure_finite_rotation_state

  !> Snapshot the converged rotation state at the start of a load step and reset
  !> the Newton trial state to it.
  subroutine fstr_begin_nodal_kinematics_step( hecMESH, fstrSOLID, ndof )
    implicit none

    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    integer(kind=kint), intent(in)        :: ndof

    call fstr_ensure_finite_rotation_state( hecMESH, fstrSOLID, ndof )
    if( .not. associated(fstrSOLID%shell_triad) ) return

    fstrSOLID%shell_triad_bak(:) = fstrSOLID%shell_triad(:)
    fstrSOLID%shell_drill_bak(:) = fstrSOLID%shell_drill(:)
    fstrSOLID%shell_dtriad(:) = fstrSOLID%shell_triad(:)
    fstrSOLID%shell_ddrill(:) = fstrSOLID%shell_drill(:)

  end subroutine fstr_begin_nodal_kinematics_step

  !> Apply the linear-solver solution increment x to the step displacement dunode.
  !>
  !> For ordinary nodes this is the usual dunode += x. For finite-rotation shell
  !> nodes the translational part is added directly while the rotational part is
  !> composed onto the trial nodal frame (dtriad) and drilling scalar.
  subroutine fstr_apply_solution_increment( hecMESH, fstrSOLID, ndof, x )
    use m_static_LIB_shell, only: ShellUpdateTriadWithIncrement, ShellComposeRotationVector
    implicit none

    type (hecmwST_local_mesh), intent(in)    :: hecMESH
    type (fstr_solid), intent(inout)         :: fstrSOLID
    integer(kind=kint), intent(in)           :: ndof
    real(kind=kreal), intent(in)             :: x(:)

    integer(kind=kint) :: node_id, idx, base
    real(kind=kreal) :: theta_inc(3), theta_compat(3)
    real(kind=kreal) :: triad_old(3, 3), triad_new(3, 3), drill_new

    if( .not. fstrSOLID%has_finite_rotation_kinematics ) then
      do node_id = 1, hecMESH%n_node
        idx = ndof*(node_id-1)
        fstrSOLID%dunode(idx+1:idx+ndof) = fstrSOLID%dunode(idx+1:idx+ndof) + x(idx+1:idx+ndof)
      end do
      return
    endif
    call fstr_ensure_finite_rotation_state( hecMESH, fstrSOLID, ndof )
    if( .not. associated(fstrSOLID%shell_node_mode) ) then
      do node_id = 1, hecMESH%n_node
        idx = ndof*(node_id-1)
        fstrSOLID%dunode(idx+1:idx+ndof) = fstrSOLID%dunode(idx+1:idx+ndof) + x(idx+1:idx+ndof)
      end do
      return
    endif

    do node_id = 1, hecMESH%n_node
      idx = ndof*(node_id-1)
      if( fstrSOLID%shell_node_mode(node_id) == 1 ) then
        fstrSOLID%dunode(idx+1:idx+3) = fstrSOLID%dunode(idx+1:idx+3) + x(idx+1:idx+3)
        theta_inc(1:3) = x(idx+4:idx+6)
        base = 9*(node_id-1)
        triad_old(1:3, 1) = fstrSOLID%shell_dtriad(base+1:base+3)
        triad_old(1:3, 2) = fstrSOLID%shell_dtriad(base+4:base+6)
        triad_old(1:3, 3) = fstrSOLID%shell_dtriad(base+7:base+9)
        call ShellUpdateTriadWithIncrement( triad_old, fstrSOLID%shell_ddrill(node_id), &
          theta_inc, triad_new, drill_new )
        fstrSOLID%shell_dtriad(base+1:base+3) = triad_new(1:3, 1)
        fstrSOLID%shell_dtriad(base+4:base+6) = triad_new(1:3, 2)
        fstrSOLID%shell_dtriad(base+7:base+9) = triad_new(1:3, 3)
        fstrSOLID%shell_ddrill(node_id) = drill_new
        ! update nodal rotation vector for output
        call ShellComposeRotationVector( fstrSOLID%dunode(idx+4:idx+6), x(idx+4:idx+6), theta_compat )
        fstrSOLID%dunode(idx+4:idx+6) = theta_compat(1:3)
        if( ndof > 6 ) then
          fstrSOLID%dunode(idx+7:idx+ndof) = fstrSOLID%dunode(idx+7:idx+ndof) + x(idx+7:idx+ndof)
        endif
      else
        fstrSOLID%dunode(idx+1:idx+ndof) = fstrSOLID%dunode(idx+1:idx+ndof) + x(idx+1:idx+ndof)
      endif
    end do

  end subroutine fstr_apply_solution_increment

  !> Commit the converged step increment dunode into the total displacement unode.
  !>
  !> For finite-rotation shell nodes the converged trial frame (dtriad) and
  !> drilling scalar become the new reference state.
  subroutine fstr_commit_solution_increment( hecMESH, fstrSOLID, ndof )
    use m_static_LIB_shell, only: ShellComposeRotationVector
    implicit none

    type (hecmwST_local_mesh), intent(in)    :: hecMESH
    type (fstr_solid), intent(inout)         :: fstrSOLID
    integer(kind=kint), intent(in)           :: ndof

    integer(kind=kint) :: node_id, idx, base
    real(kind=kreal) :: theta_compat(3)

    if( .not. fstrSOLID%has_finite_rotation_kinematics ) then
      do node_id = 1, hecMESH%n_node
        idx = ndof*(node_id-1)
        fstrSOLID%unode(idx+1:idx+ndof) = fstrSOLID%unode(idx+1:idx+ndof) + fstrSOLID%dunode(idx+1:idx+ndof)
      end do
      return
    endif
    call fstr_ensure_finite_rotation_state( hecMESH, fstrSOLID, ndof )
    if( .not. associated(fstrSOLID%shell_node_mode) ) then
      do node_id = 1, hecMESH%n_node
        idx = ndof*(node_id-1)
        fstrSOLID%unode(idx+1:idx+ndof) = fstrSOLID%unode(idx+1:idx+ndof) + fstrSOLID%dunode(idx+1:idx+ndof)
      end do
      return
    endif

    do node_id = 1, hecMESH%n_node
      idx = ndof*(node_id-1)
      if( fstrSOLID%shell_node_mode(node_id) == 1 ) then
        fstrSOLID%unode(idx+1:idx+3) = fstrSOLID%unode(idx+1:idx+3) + fstrSOLID%dunode(idx+1:idx+3)
        base = 9*(node_id-1)
        fstrSOLID%shell_triad(base+1:base+9) = fstrSOLID%shell_dtriad(base+1:base+9)
        fstrSOLID%shell_drill(node_id) = fstrSOLID%shell_ddrill(node_id)
        ! update nodal rotation vector for output
        call ShellComposeRotationVector( fstrSOLID%unode(idx+4:idx+6), fstrSOLID%dunode(idx+4:idx+6), theta_compat )
        fstrSOLID%unode(idx+4:idx+6) = theta_compat(1:3)
        if( ndof > 6 ) then
          fstrSOLID%unode(idx+7:idx+ndof) = fstrSOLID%unode(idx+7:idx+ndof) + fstrSOLID%dunode(idx+7:idx+ndof)
        endif
      else
        fstrSOLID%unode(idx+1:idx+ndof) = fstrSOLID%unode(idx+1:idx+ndof) + fstrSOLID%dunode(idx+1:idx+ndof)
      endif
    end do

  end subroutine fstr_commit_solution_increment

  !> Half-thickness director from the current Newton trial frame (dtriad).
  subroutine fstr_get_shell_trial_directors( fstrSOLID, thick, nn, nodLOCAL, directors )
    implicit none

    type (fstr_solid), intent(in)      :: fstrSOLID
    real(kind=kreal), intent(in)       :: thick
    integer(kind=kint), intent(in)     :: nn
    integer(kind=kint), intent(in)     :: nodLOCAL(:)
    real(kind=kreal), intent(out)      :: directors(3, nn)

    integer(kind=kint) :: j, node_id, base

    directors(:, :) = 0.0D0
    if( .not. associated(fstrSOLID%shell_dtriad) ) return

    do j = 1, nn
      node_id = nodLOCAL(j)
      if( node_id <= 0 .or. node_id > size(fstrSOLID%shell_rot_state) ) cycle
      base = 9*(node_id-1)
      directors(1:3, j) = 0.5D0*thick*fstrSOLID%shell_dtriad(base+7:base+9)
    end do

  end subroutine fstr_get_shell_trial_directors

  !> Drilling scalar from the current Newton trial frame.
  subroutine fstr_get_shell_trial_drilling( fstrSOLID, nn, nodLOCAL, drilling )
    implicit none

    type (fstr_solid), intent(in)      :: fstrSOLID
    integer(kind=kint), intent(in)     :: nn
    integer(kind=kint), intent(in)     :: nodLOCAL(:)
    real(kind=kreal), intent(out)      :: drilling(nn)

    integer(kind=kint) :: j, node_id

    drilling(:) = 0.0D0
    if( .not. associated(fstrSOLID%shell_ddrill) ) return

    do j = 1, nn
      node_id = nodLOCAL(j)
      if( node_id <= 0 .or. node_id > size(fstrSOLID%shell_ddrill) ) cycle
      drilling(j) = fstrSOLID%shell_ddrill(node_id)
    end do

  end subroutine fstr_get_shell_trial_drilling

  !> Half-thickness director from the converged frame (triad).
  subroutine fstr_get_shell_current_directors( fstrSOLID, thick, nn, nodLOCAL, directors )
    implicit none

    type (fstr_solid), intent(in)      :: fstrSOLID
    real(kind=kreal), intent(in)       :: thick
    integer(kind=kint), intent(in)     :: nn
    integer(kind=kint), intent(in)     :: nodLOCAL(:)
    real(kind=kreal), intent(out)      :: directors(3, nn)

    integer(kind=kint) :: j, node_id, base

    directors(:, :) = 0.0D0
    if( .not. associated(fstrSOLID%shell_triad) ) return

    do j = 1, nn
      node_id = nodLOCAL(j)
      if( node_id <= 0 .or. node_id > size(fstrSOLID%shell_rot_state) ) cycle
      base = 9*(node_id-1)
      directors(1:3, j) = 0.5D0*thick*fstrSOLID%shell_triad(base+7:base+9)
    end do

  end subroutine fstr_get_shell_current_directors

  !> Half-thickness director from the fixed reference frame (ref_triad).
  subroutine fstr_get_shell_reference_directors( fstrSOLID, thick, nn, nodLOCAL, directors )
    implicit none

    type (fstr_solid), intent(in)      :: fstrSOLID
    real(kind=kreal), intent(in)       :: thick
    integer(kind=kint), intent(in)     :: nn
    integer(kind=kint), intent(in)     :: nodLOCAL(:)
    real(kind=kreal), intent(out)      :: directors(3, nn)

    integer(kind=kint) :: j, node_id, base

    directors(:, :) = 0.0D0
    if( .not. associated(fstrSOLID%shell_ref_triad) ) return

    do j = 1, nn
      node_id = nodLOCAL(j)
      if( node_id <= 0 .or. node_id > size(fstrSOLID%shell_rot_state) ) cycle
      base = 9*(node_id-1)
      directors(1:3, j) = 0.5D0*thick*fstrSOLID%shell_ref_triad(base+7:base+9)
    end do

  end subroutine fstr_get_shell_reference_directors

  ! ---------------------------------------------------------------------------
  ! Private helpers: reference-frame construction at element nodes.
  ! ---------------------------------------------------------------------------

  subroutine fstr_set_identity_triad( triad )
    implicit none

    real(kind=kreal), intent(out) :: triad(3, 3)

    triad(:, :) = 0.0D0
    triad(1, 1) = 1.0D0
    triad(2, 2) = 1.0D0
    triad(3, 3) = 1.0D0

  end subroutine fstr_set_identity_triad

  !> Reference nodal frame at corner inode of a MITC4 element: e3 from the surface
  !> normal (g1 x g2), e2 = e3 x e0, e1 = e2 x e3, then orthonormalized.
  subroutine fstr_reference_shell_triad( nn, ecoord, inode, triad )
    use elementInfo, only: fe_mitc4_shell, getShapeDeriv
    use m_static_LIB_shell, only: ShellOrthonormalizeTriad
    implicit none

    integer(kind=kint), intent(in) :: nn
    integer(kind=kint), intent(in) :: inode
    real(kind=kreal), intent(in)  :: ecoord(3, nn)
    real(kind=kreal), intent(out) :: triad(3, 3)

    real(kind=kreal) :: xi, eta
    real(kind=kreal) :: shapederiv(4, 2)
    real(kind=kreal) :: g1(3), g2(3), e0(3), trial(3, 3), normv
    integer(kind=kint) :: i

    call fstr_set_identity_triad( trial )
    if( nn /= 4 ) then
      triad(1:3, 1:3) = trial(1:3, 1:3)
      return
    endif

    call getShapeDeriv( fe_mitc4_shell, (/ 0.0D0, 0.0D0 /), shapederiv )
    e0(1:3) = 0.0D0
    do i = 1, 4
      e0(1:3) = e0(1:3) + shapederiv(i, 1)*ecoord(1:3, i)
    end do

    select case( inode )
    case( 1 )
      xi = -1.0D0
      eta = -1.0D0
    case( 2 )
      xi =  1.0D0
      eta = -1.0D0
    case( 3 )
      xi =  1.0D0
      eta =  1.0D0
    case default
      xi = -1.0D0
      eta =  1.0D0
    end select

    call getShapeDeriv( fe_mitc4_shell, (/ xi, eta /), shapederiv )
    g1(1:3) = 0.0D0
    g2(1:3) = 0.0D0
    do i = 1, 4
      g1(1:3) = g1(1:3) + shapederiv(i, 1)*ecoord(1:3, i)
      g2(1:3) = g2(1:3) + shapederiv(i, 2)*ecoord(1:3, i)
    end do

    trial(1, 3) = g1(2)*g2(3) - g1(3)*g2(2)
    trial(2, 3) = g1(3)*g2(1) - g1(1)*g2(3)
    trial(3, 3) = g1(1)*g2(2) - g1(2)*g2(1)

    trial(1, 2) = trial(2, 3)*e0(3) - trial(3, 3)*e0(2)
    trial(2, 2) = trial(3, 3)*e0(1) - trial(1, 3)*e0(3)
    trial(3, 2) = trial(1, 3)*e0(2) - trial(2, 3)*e0(1)
    normv = dsqrt( dot_product( trial(1:3, 2), trial(1:3, 2) ) )
    if( normv > 1.0D-14 ) trial(1:3, 2) = trial(1:3, 2)/normv

    trial(1, 1) = trial(2, 2)*trial(3, 3) - trial(3, 2)*trial(2, 3)
    trial(2, 1) = trial(3, 2)*trial(1, 3) - trial(1, 2)*trial(3, 3)
    trial(3, 1) = trial(1, 2)*trial(2, 3) - trial(2, 2)*trial(1, 3)

    call ShellOrthonormalizeTriad( trial, triad )

  end subroutine fstr_reference_shell_triad

  !> Store the reference frame into all rotation-state arrays for a fresh node.
  subroutine fstr_store_shell_triad_node( fstrSOLID, node_id, triad, mode )
    implicit none

    type (fstr_solid), intent(inout)  :: fstrSOLID
    integer(kind=kint), intent(in)    :: node_id, mode
    real(kind=kreal), intent(in)      :: triad(3, 3)

    integer(kind=kint) :: base

    if( node_id <= 0 ) return
    if( node_id > size(fstrSOLID%shell_rot_state) ) return
    if( fstrSOLID%shell_rot_state(node_id) /= 0 ) return

    base = 9*(node_id-1)
    fstrSOLID%shell_rot_state(node_id) = mode
    ! initialize reference and current shell triads
    fstrSOLID%shell_ref_triad(base+1:base+3) = triad(1:3, 1)
    fstrSOLID%shell_ref_triad(base+4:base+6) = triad(1:3, 2)
    fstrSOLID%shell_ref_triad(base+7:base+9) = triad(1:3, 3)
    fstrSOLID%shell_triad(base+1:base+3) = triad(1:3, 1)
    fstrSOLID%shell_triad(base+4:base+6) = triad(1:3, 2)
    fstrSOLID%shell_triad(base+7:base+9) = triad(1:3, 3)
    fstrSOLID%shell_triad_bak(base+1:base+9) = fstrSOLID%shell_triad(base+1:base+9)
    fstrSOLID%shell_dtriad(base+1:base+9) = fstrSOLID%shell_triad(base+1:base+9)
    fstrSOLID%shell_drill(node_id) = 0.0D0
    fstrSOLID%shell_drill_bak(node_id) = 0.0D0
    fstrSOLID%shell_ddrill(node_id) = 0.0D0

  end subroutine fstr_store_shell_triad_node

end module m_fstr_NodalKinematics
