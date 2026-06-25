!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Rigid-body / near-kernel modes from mesh coordinates
!!
!! Shared by the ML and SA-AMG preconditioners (single source of truth for the
!! near-kernel used to seed smoothed-aggregation coarsening).  Builds the
!! translational + rotational rigid-body modes from the nodal coordinates of the
!! INTERNAL nodes, in the flat column-major (dof, mode) layout both callers use:
!! mode k occupies rbm((k-1)*vec_leng+1 : k*vec_leng), vec_leng = nn_internal*ndof.
!!
!!   ndof = 1 -> 1 mode  (constant; heat / scalar diffusion)
!!   ndof = 2 -> 3 modes (trans-x, trans-y, rot-z; plane elasticity)
!!   ndof = 3 -> 6 modes (trans*3, rot*3; solid, incl. 33struct shell)
!!   ndof = 6 -> 6 modes (as ndof=3 plus the rotational-DOF identity for shells)
!!
!! 33struct shell rotation nodes (the latter-half nodes of each 33struct element)
!! carry only the rotational modes.  Coordinates are shifted to the GLOBAL centroid
!! (allreduce over hecMESH) so the rotation modes stay well-conditioned far from the
!! origin; this does not change the spanned near-kernel space.
module hecmw_precond_rbm
  use hecmw_util
  use hecmw_etype
  use m_hecmw_comm_f, only: hecmw_allreduce_R
  implicit none

  private
  public :: hecmw_precond_rbm_from_mesh

contains

  subroutine hecmw_precond_rbm_from_mesh(hecMESH, ndof, rbm)
    implicit none
    type(hecmwST_local_mesh), intent(in)  :: hecMESH
    integer(kind=kint),       intent(in)  :: ndof
    real(kind=kreal),         intent(out) :: rbm(*)   !< flat (dof, mode), filled for internal nodes
    integer(kind=kint), allocatable :: mark(:)
    integer(kind=kint) :: nint, vl, nmode, node, base
    integer(kind=kint) :: itype, ic_type, nn, is, iE, icel, iiS, j, nod
    real(kind=kreal) :: cx, cy, cz, x, y, z, csum(4)

    nint = hecMESH%nn_internal
    vl   = nint * ndof
    select case (ndof)
    case (1);    nmode = 1
    case (2);    nmode = 3
    case default; nmode = 6
    end select

    ! mark 33struct shell rotation nodes (latter-half element nodes)
    allocate(mark(hecMESH%n_node)); mark = 0
    if (ndof == 3) then
      do itype = 1, hecMESH%n_elem_type
        ic_type = hecMESH%elem_type_item(itype)
        if (hecmw_is_etype_33struct(ic_type)) then
          nn = hecmw_get_max_node(ic_type)
          is = hecMESH%elem_type_index(itype-1)+1
          iE = hecMESH%elem_type_index(itype)
          do icel = is, iE
            iiS = hecMESH%elem_node_index(icel-1)
            do j = nn/2+1, nn
              nod = hecMESH%elem_node_item(iiS+j); mark(nod) = 1
            end do
          end do
        end if
      end do
    end if

    ! global centroid, for well-conditioned rotation modes far from the origin
    ! (sum of coordinates and node count reduced together in one allreduce)
    csum = 0.0d0
    do node = 1, nint
      csum(1) = csum(1) + hecMESH%node(3*node-2)
      csum(2) = csum(2) + hecMESH%node(3*node-1)
      csum(3) = csum(3) + hecMESH%node(3*node)
    end do
    csum(4) = real(nint, kreal)
    call hecmw_allreduce_R(hecMESH, csum, 4, HECMW_SUM)
    cx = 0.0d0; cy = 0.0d0; cz = 0.0d0
    if (csum(4) > 0.0d0) then
      cx = csum(1)/csum(4); cy = csum(2)/csum(4); cz = csum(3)/csum(4)
    end if

    rbm(1:vl*nmode) = 0.0d0
    if (ndof == 1) then
      do node = 1, nint
        rbm(node) = 1.0d0                                   ! constant
      end do
    else if (ndof == 2) then
      do node = 1, nint
        x = hecMESH%node(3*node-2) - cx; y = hecMESH%node(3*node-1) - cy
        base = 2*(node-1)
        rbm(        base+1) = 1.0d0                          ! trans-x
        rbm(   vl + base+2) = 1.0d0                          ! trans-y
        rbm( 2*vl + base+1) = -y                            ! rot-z
        rbm( 2*vl + base+2) =  x
      end do
    else   ! ndof = 3 or 6
      do node = 1, nint
        base = ndof*(node-1)
        if (ndof == 3 .and. mark(node) == 1) then
          ! 33struct rotation node: rotational rigid-body modes only [0 | I3]
          rbm( 3*vl + base+1) = 1.0d0
          rbm( 4*vl + base+2) = 1.0d0
          rbm( 5*vl + base+3) = 1.0d0
        else
          x = hecMESH%node(3*node-2) - cx
          y = hecMESH%node(3*node-1) - cy
          z = hecMESH%node(3*node)   - cz
          rbm(        base+1) = 1.0d0                        ! trans-x
          rbm(   vl + base+2) = 1.0d0                        ! trans-y
          rbm( 2*vl + base+3) = 1.0d0                        ! trans-z
          rbm( 3*vl + base+2) = -z                          ! rot-x
          rbm( 3*vl + base+3) =  y
          rbm( 4*vl + base+1) =  z                          ! rot-y
          rbm( 4*vl + base+3) = -x
          rbm( 5*vl + base+1) = -y                          ! rot-z
          rbm( 5*vl + base+2) =  x
          if (ndof == 6) then
            rbm( 3*vl + base+4) = 1.0d0                      ! shell rotational DOF identity
            rbm( 4*vl + base+5) = 1.0d0
            rbm( 5*vl + base+6) = 1.0d0
          end if
        end if
      end do
    end if

    deallocate(mark)
  end subroutine hecmw_precond_rbm_from_mesh

end module hecmw_precond_rbm
