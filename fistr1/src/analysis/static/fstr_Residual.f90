!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calculate residual of nodal force.

module m_fstr_Residual
  use hecmw
  implicit none

  public :: fstr_Update_NDForce
  public :: fstr_Update_NDForce_SPC
  public :: fstr_get_residual
  public :: fstr_get_norm_contact
  public :: fstr_get_norm_para_contact
  public :: fstr_get_x_norm_contact
  public :: fstr_get_potential

  private :: fstr_Update_NDForce_MPC

contains

  !C---------------------------------------------------------------------*
  subroutine fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT)
    !C---------------------------------------------------------------------*
    !> In this subroutine, nodal force arising from prescribed displacement constraints
    !> are cleared and nodal force residual is calculated.
    !> Those constraints considered here includes:
    !!-#  nodal displacement
    !!-#  equation (or mpc)
    use m_fstr
    use mULoad
    use m_fstr_spring
    integer(kind=kint), intent(in)       :: cstep !< current step
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    type(hecmwST_matrix), intent(inout)  :: hecMAT !< linear equation, its right side modified here
    type(fstr_solid), intent(inout)      :: fstrSOLID !< we need boundary conditions of curr step
    type(hecmwST_matrix), intent(inout), optional  :: conMAT
    !    Local variables
    integer(kind=kint) :: ndof, idof
    real(kind=kreal)   :: factor

    factor = fstrSOLID%factor(2)

    !    Set residual load
    do idof=1, hecMESH%n_node*  hecMESH%n_dof
      hecMAT%B(idof)=fstrSOLID%GL(idof)-fstrSOLID%QFORCE(idof)
    end do
    ndof = hecMAT%NDOF

    call fstr_Update_NDForce_spring( cstep, hecMESH, fstrSOLID, hecMAT%B )

    !    Consider Uload
    call uResidual( cstep, factor, hecMAT%B )

    !    Consider EQUATION condition
    call fstr_Update_NDForce_MPC( hecMESH, hecMAT%B )

    !    Consider SPC condition
    call fstr_Update_NDForce_SPC( cstep, hecMESH, fstrSOLID, hecMAT%B )
    if(present(conMAT)) call fstr_Update_NDForce_SPC( cstep, hecMESH, fstrSOLID, conMAT%B )

    !
    call hecmw_update_R(hecMESH,hecMAT%B,hecMESH%n_node, ndof)
  end subroutine fstr_Update_NDForce

  subroutine fstr_Update_NDForce_MPC( hecMESH, B )
    use m_fstr
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    real(kind=kreal), intent(inout)      :: B(:) !< right hand side
    !    Local variables
    integer(kind=kint) ndof, ig0, iS0, iE0, ik, in, idof
    real(kind=kreal) :: rhs, lambda

    ndof = hecMESH%n_dof
    OUTER: do ig0=1,hecMESH%mpc%n_mpc
      iS0= hecMESH%mpc%mpc_index(ig0-1)+1
      iE0= hecMESH%mpc%mpc_index(ig0)
      do ik= iS0, iE0
        if (hecMESH%mpc%mpc_dof(ik) > ndof) cycle OUTER
      enddo
      ! Suppose the lagrange multiplier= first dof of first node
      in = hecMESH%mpc%mpc_item(iS0)
      idof = hecMESH%mpc%mpc_dof(iS0)
      rhs = hecMESH%mpc%mpc_val(iS0)
      lambda = B(ndof*(in-1)+idof)/rhs
      ! update nodal residual
      do ik= iS0, iE0
        in = hecMESH%mpc%mpc_item(ik)
        idof = hecMESH%mpc%mpc_dof(ik)
        rhs = hecMESH%mpc%mpc_val(ik)
        B(ndof*(in-1)+idof) = B(ndof*(in-1)+idof) - rhs*lambda
      enddo
    enddo OUTER
  end subroutine fstr_Update_NDForce_MPC

  subroutine fstr_Update_NDForce_SPC( cstep, hecMESH, fstrSOLID, B )
    use m_fstr
    integer(kind=kint), intent(in)       :: cstep !< current step
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    type(fstr_solid), intent(in)         :: fstrSOLID !< we need boundary conditions of curr step
    real(kind=kreal), intent(inout)      :: B(:) !< right hand side
    !    Local variables
    integer(kind=kint) ndof, ig0, ig, ityp, iS0, iE0, ik, in, idof1, idof2, idof
    integer(kind=kint) :: grpid
    real(kind=kreal) :: rhs

    ndof = hecMESH%n_dof
    fstrSOLID%REACTION = 0.d0

    do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      rhs= fstrSOLID%BOUNDARY_ngrp_val(ig0)
      ityp= fstrSOLID%BOUNDARY_ngrp_type(ig0)
      iS0= hecMESH%node_group%grp_index(ig-1) + 1
      iE0= hecMESH%node_group%grp_index(ig  )
      do ik= iS0, iE0
        in   = hecMESH%node_group%grp_item(ik)
        idof1 = ityp/10
        idof2 = ityp - idof1*10
        if( fstrSOLID%BOUNDARY_ngrp_rotID(ig0) > 0 ) then
          idof1 = 1
          idof2 = ndof
        end if
        do idof=idof1,idof2
          B( ndof*(in-1) + idof ) = 0.d0
          !for output reaction force
          fstrSOLID%REACTION(ndof*(in-1)+idof) = fstrSOLID%QFORCE(ndof*(in-1)+idof)
          !count embed force as reaction force
          if( associated(fstrSOLID%EMBED_NFORCE) ) fstrSOLID%REACTION(ndof*(in-1)+idof) = &
                &  fstrSOLID%QFORCE(ndof*(in-1)+idof) - fstrSOLID%EMBED_NFORCE(ndof*(in-1)+idof)
        enddo
      enddo
    enddo
  end subroutine fstr_Update_NDForce_SPC


  !> \breaf This subroutine calculate residual vector
  subroutine fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    use m_fstr_Update
    use m_fstr_ass_load
    implicit none
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter

    ! ----- update the strain, stress, and internal force
    call fstr_UpdateNewton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter)

    ! ----- Set residual
    if( fstrSOLID%DLOAD_follow /= 0 .or. fstrSOLID%CLOAD_ngrp_rot /= 0 ) &
      & call fstr_ass_load(cstep, ctime+dtime, hecMESH, hecMAT, fstrSOLID, fstrPARAM )

    call fstr_Update_NDForce(cstep, hecMESH, hecMAT, fstrSOLID)
  end subroutine fstr_calc_residual_vector

  !> \breaf This subroutine calculate residual vector
  subroutine fstr_calc_residual_vector_with_X(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    use m_fstr_Update
    implicit none
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter

    integer :: i
    real(kind=kreal) :: dunode_bak(hecMAT%ndof*hecMESH%n_node)

    do i=1, hecMAT%ndof*hecMESH%n_node
      dunode_bak(i) = fstrSOLID%dunode(i)
      fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
    enddo
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    do i=1, hecMAT%ndof*hecMESH%n_node
      fstrSOLID%dunode(i) = dunode_bak(i)
    enddo
  end subroutine fstr_calc_residual_vector_with_X

  !> Calculate magnitude of a real vector
  real(kind=kreal) function fstr_get_residual( force, hecMESH )
    use m_fstr
    real(kind=kreal), intent(in)         :: force(:)
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    integer :: ndof
    ndof = hecMESH%n_dof
    call hecmw_innerProduct_R(hecMESH,ndof,force,force,fstr_get_residual)
  end function

  !> Calculate magnitude of a real vector
  real(kind=kreal) function fstr_get_energy( force, displacement, hecMESH )
    use m_fstr
    real(kind=kreal), intent(in)         :: force(:), displacement(:)
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    integer :: ndof
    ndof = hecMESH%n_dof
    call hecmw_innerProduct_R(hecMESH, ndof, force, displacement, fstr_get_energy)
  end function

  !> Calculate square norm
  real(kind=kreal) function fstr_get_norm_contact(flag,hecMESH,hecMAT,fstrSOLID,hecLagMAT)
    use m_fstr
    type(hecmwST_local_mesh), intent(in)             :: hecMESH !< mesh information
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(fstr_solid), intent(in)                     :: fstrSOLID
    type(hecmwST_matrix_lagrange), intent(in)        :: hecLagMAT
    character(len=13)                                :: flag
    real(kind=kreal) :: tmp1, tmp2, bi
    integer :: i, i0, ndof
    if( flag=='residualForce' )then
      ndof = hecMESH%n_dof
      call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%B,hecMAT%B,tmp1)
      tmp2 = 0.0d0
      i0 = hecMESH%n_node*ndof
      do i=1,hecLagMAT%num_lagrange
        bi = hecMAT%B(i0+i)
        tmp2 = tmp2 + bi*bi
      enddo
      call hecmw_allreduce_R1(hecMESH,tmp2,HECMW_SUM)
      fstr_get_norm_contact = tmp1 + tmp2
    elseif( flag=='        force' )then
      call hecmw_innerProduct_R(hecMESH,ndof,fstrSOLID%QFORCE,fstrSOLID%QFORCE,fstr_get_norm_contact)
    endif
  end function

  !
  function fstr_get_norm_para_contact(hecMAT,hecLagMAT,conMAT,hecMESH) result(rhsB)
    use m_fstr
    implicit none
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(hecmwST_matrix_lagrange), intent(in)        :: hecLagMAT
    type(hecmwST_matrix), intent(in)                 :: conMAT
    type(hecmwST_local_mesh), intent(in)             :: hecMESH
    !
    real(kind=kreal) ::  rhsB
    integer(kind=kint) ::  i,ndof,nndof,npndof,num_lagrange
    real(kind=kreal), allocatable   :: rhs_con(:)
    real(kind=kreal), pointer :: rhs_lag(:)

    ndof = conMAT%ndof
    nndof = conMAT%N * ndof
    npndof = conMAT%NP * ndof
    num_lagrange = hecLagMAT%num_lagrange

    allocate(rhs_con(npndof))
    do i=1,npndof
      rhs_con(i) = conMAT%B(i)
    enddo
    call hecmw_assemble_R(hecMESH, rhs_con, conMAT%NP, conMAT%NDOF)

    do i=1,nndof
      rhs_con(i) = rhs_con(i) + hecMAT%B(i)
    enddo

    rhs_lag => conMAT%B(npndof+1:npndof+num_lagrange)

    rhsB = dot_product(rhs_con(1:nndof), rhs_con(1:nndof)) + dot_product(rhs_lag(:), rhs_lag(:))
    call hecmw_allreduce_R1(hecMESH, rhsB, hecmw_sum)
    deallocate(rhs_con)

  end function fstr_get_norm_para_contact

  function fstr_get_x_norm_contact(hecMAT,hecLagMAT,hecMESH) result(rhsX)
    use m_fstr
    implicit none
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(hecmwST_matrix_lagrange), intent(in)        :: hecLagMAT
    type(hecmwST_local_mesh), intent(in)             :: hecMESH
    real(kind=kreal)   ::  rhsX
    integer(kind=kint) :: nndof, npndof, i

    nndof = hecMAT%N * hecMAT%NDOF
    npndof = hecMAT%NP * hecMAT%NDOF
    rhsX = 0.d0
    do i=1,nndof
      rhsX = rhsX + hecMAT%X(i) ** 2
    end do
    do i=1,hecLagMAT%num_lagrange
      rhsX = rhsX + hecMAT%X(npndof+i) ** 2
    end do
    call hecmw_allreduce_R1(hecMESH, rhsX, hecmw_sum)

  end function fstr_get_x_norm_contact

  !C---------------------------------------------------------------------*
  function fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,ptype) result(potential)
    use m_fstr
    use mULoad
    use m_fstr_spring
    integer(kind=kint), intent(in)       :: cstep !< current step
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    type(hecmwST_matrix), intent(inout)  :: hecMAT !< linear equation, its right side modified here
    type(fstr_solid), intent(inout)      :: fstrSOLID !< we need boundary conditions of curr step
    integer(kind=kint), intent(in)       :: ptype
    real(kind=kreal)   ::  potential
    !    Local variables
    integer(kind=kint) :: ndof, i, icel, icel0, ig
    real(kind=kreal), allocatable :: totdisp(:), totload(:)
    real(kind=kreal) :: factor
    real(kind=kreal) :: extenal_work, internal_work

    ndof = hecMAT%NDOF
    factor = fstrSOLID%factor(2)
    potential = 0.d0

    allocate(totload(ndof*hecMESH%n_node))

    ! Set external load
    do i=1,ndof*hecMESH%n_node
      totload(i) = 0.5d0*(fstrSOLID%GL(i)+fstrSOLID%GL0(i))
    end do
    !    Consiter Spring
    call fstr_Update_NDForce_spring( cstep, hecMESH, fstrSOLID, totload )
    !    Consider Uload
    call uResidual( cstep, factor, totload )
    !    Consider EQUATION condition
    call fstr_Update_NDForce_MPC( hecMESH, totload )

    if( ptype == 1 ) then !< assemble strain potential
      ! calc internal work
      internal_work = 0.d0
      do icel0=1,hecMESH%ne_internal
        icel=hecMESH%elem_internal_list(icel0)
        do ig=1,size(fstrSOLID%elements(icel)%gausses)
          internal_work = internal_work + fstrSOLID%elements(icel)%gausses(ig)%strain_energy
        enddo
      enddo
      call hecmw_allreduce_R1(hecMESH, internal_work, hecmw_sum)
      ! calc external work
      call hecmw_innerProduct_R(hecMESH,ndof,totload,fstrSOLID%dunode,extenal_work)      
      potential = internal_work - extenal_work
    else if( ptype == 2 ) then ! multiply internal force with displacement
      ! calc internal work
      call hecmw_innerProduct_R(hecMESH,ndof,0.5d0*(fstrSOLID%QFORCE+fstrSOLID%QFORCE_bak),fstrSOLID%dunode,internal_work)
      ! calc external work
      call hecmw_innerProduct_R(hecMESH,ndof,totload,fstrSOLID%dunode,extenal_work)
      potential = internal_work - extenal_work
    else if( ptype == 3 ) then ! norm of residual
      ! Consider internal force
      totload(:) = fstrSOLID%QFORCE(:) - totload(:)
      !    Consider SPC condition
      call fstr_Update_NDForce_SPC( cstep, hecMESH, fstrSOLID, totload )
      ! calc norm of residual
      call hecmw_innerProduct_R(hecMESH,ndof,totload,totload,potential)
      potential = dsqrt(potential)
    else
      stop "ptype error in fstr_get_potential"
    endif
    !write(IMSG,*) "internal_work,extenal_work",internal_work,extenal_work
  end function fstr_get_potential

  function fstr_get_potential_with_X(cstep,hecMESH,hecMAT,fstrSOLID,ptype) result(potential)
    use m_fstr
    use mULoad
    use m_fstr_spring
    implicit none
    integer(kind=kint), intent(in)       :: cstep !< current step
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< mesh information
    type(hecmwST_matrix), intent(inout)  :: hecMAT !< linear equation, its right side modified here
    type(fstr_solid), intent(inout)      :: fstrSOLID !< we need boundary conditions of curr step
    integer(kind=kint), intent(in)       :: ptype
    real(kind=kreal)   ::  potential

    integer :: i
    real(kind=kreal) :: dunode_bak(hecMAT%ndof*hecMESH%n_node)

    do i=1, hecMAT%ndof*hecMESH%n_node
      dunode_bak(i) = fstrSOLID%dunode(i)
      fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
    enddo
    potential = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,ptype)
    do i=1, hecMAT%ndof*hecMESH%n_node
      fstrSOLID%dunode(i) = dunode_bak(i)
    enddo
  end function fstr_get_potential_with_X

  subroutine plot_potential_graph( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, &
    restrt_step_num, sub_step, ctime, dtime, iter, tincr )
    use m_fstr
    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint) :: restrt_step_num

    integer(kind=kint) :: iter, i, j
    integer(kind=kint), parameter :: ntot = 30
    real(kind=kreal)   :: tincr, alpha, dulen
    real(kind=kreal) :: pot(3)
    real(kind=kreal), allocatable :: dunode_bak(:)

    dulen = dsqrt(dot_product(fstrSOLID%dunode(:),fstrSOLID%dunode(:)))
    write(IMSG,*) "dulen",dulen

    allocate(dunode_bak(size(fstrSOLID%dunode)))
    do i=1, hecMAT%ndof*hecMESH%n_node
      dunode_bak(i) = fstrSOLID%dunode(i)
    enddo

    do i=-ntot,ntot
      alpha = 5.d-3*dble(i)/dble(ntot)
      hecMAT%X(:) = alpha*fstrSOLID%dunode(:)
      do j=1, size(fstrSOLID%dunode)
        fstrSOLID%dunode(j) = dunode_bak(j)+hecMAT%X(j)
      enddo
      call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
      pot(1) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,1)
      pot(2) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,2)
      pot(3) = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,3)
      write(IMSG,'(I4,5(",",1pE16.9))') i,alpha,alpha*dulen,pot(1:3)
    enddo

    do i=1, size(fstrSOLID%dunode)
      fstrSOLID%dunode(i) = dunode_bak(i)
    enddo
    hecMAT%X(:) = 0.d0
    call fstr_calc_residual_vector_with_X(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
  end subroutine

end module m_fstr_Residual
