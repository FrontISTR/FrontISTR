!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calcualte residual of nodal force.

module m_fstr_Residual
  use hecmw
  implicit none

  public :: fstr_Update_NDForce
  public :: fstr_Update_NDForce_SPC
  public :: fstr_get_residual
  public :: fstr_get_norm_contact
  public :: fstr_get_norm_para_contact
  public :: fstr_get_x_norm_contact

  private :: fstr_Update_NDForce_MPC

contains

  !C---------------------------------------------------------------------*
  subroutine fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,conMAT)
    !C---------------------------------------------------------------------*
    !> In this subroutine, nodal force arose from prescribed displacement constarints
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
    integer(kind=kint) :: ndof, idof, num
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
    if( ndof==3 ) then
      call hecmw_update_3_R(hecMESH,hecMAT%B,hecMESH%n_node)
    else if( ndof==2 ) then
      call hecmw_update_2_R(hecMESH,hecMAT%B,hecMESH%n_node)
    else if( ndof==6 ) then
      call hecmw_update_m_R(hecMESH,hecMAT%B,hecMESH%n_node,6)
    endif

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
        enddo
      enddo
    enddo
  end subroutine fstr_Update_NDForce_SPC

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
  real(kind=kreal) function fstr_get_norm_contact(flag,hecMESH,hecMAT,fstrSOLID,fstrMAT)
    use m_fstr
    use fstr_matrix_con_contact
    type(hecmwST_local_mesh), intent(in)             :: hecMESH !< mesh information
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(fstr_solid), intent(in)                     :: fstrSOLID
    type(fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    character(len=13)                                :: flag
    real(kind=kreal) :: tmp1, tmp2, bi
    integer :: i, i0, ndof
    if( flag=='residualForce' )then
      ndof = hecMESH%n_dof
      call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%B,hecMAT%B,tmp1)
      tmp2 = 0.0d0
      i0 = hecMESH%n_node*ndof
      do i=1,fstrMAT%num_lagrange
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
  function fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH) result(rhsB)
    use m_fstr
    use fstr_matrix_con_contact
    implicit none
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    type(hecmwST_matrix), intent(in)                 :: conMAT
    type(hecmwST_local_mesh), intent(in)             :: hecMESH
    !
    real(kind=kreal) ::  rhsB
    integer(kind=kint) ::  i,j,N,i0,ndof,N_loc,nndof
    integer(kind=kint) :: offset, pid, lid
    integer(kind=kint), allocatable :: displs(:)
    real(kind=kreal), allocatable   :: rhs_con_all(:), rhs_con(:)
    !
    ndof = hecMAT%ndof
    nndof = hecMAT%N*ndof
    N_loc = nndof + fstrMAT%num_lagrange
    allocate(displs(0:nprocs))
    displs(:) = 0
    displs(myrank+1) = N_loc
    call hecmw_allreduce_I(hecMESH, displs, nprocs+1, hecmw_sum)
    do i=1,nprocs
      displs(i) = displs(i-1) + displs(i)
    end do
    offset = displs(myrank)
    N = displs(nprocs)
    allocate(rhs_con_all(N))
    do i=1,N
      rhs_con_all(i) = 0.0D0
    end do
    do i= hecMAT%N+1,hecMAT%NP
      !        i0 = getExternalGlobalIndex(i,ndof,hecMAT%N)
      pid = hecMESH%node_ID(i*2)
      lid = hecMESH%node_ID(i*2-1)
      i0 = displs(pid) + (lid-1)*ndof
      if((i0 < 0.or.i0 > N)) then
        !          print *,myrank,'ext dummy',i,i0/ndof
        do j=1,ndof
          if(conMAT%b((i-1)*ndof+j) /= 0.0D0) then
            print *,myrank,'i0',i,'conMAT%b',conMAT%b((i-1)*ndof+j)
            stop
          endif
        enddo
      else
        rhs_con_all(i0+1:i0+ndof) = conMAT%b((i-1)*ndof+1:i*ndof)
      endif
    enddo
    deallocate(displs)
    call hecmw_allreduce_R(hecMESH, rhs_con_all, N, hecmw_sum)
    allocate(rhs_con(N_loc))
    do i=1,nndof
      rhs_con(i) = rhs_con_all(offset+i) + hecMAT%B(i) + conMAT%B(i)
    end do
    deallocate(rhs_con_all)
    i0 = hecMAT%NP*ndof
    do i=1,fstrMAT%num_lagrange
      rhs_con(nndof+i) = conMAT%B(i0+i)
    enddo
    rhsB = 0.d0
    rhsB = dot_product(rhs_con, rhs_con)
    call hecmw_allreduce_R1(hecMESH, rhsB, hecmw_sum)
    deallocate(rhs_con)

  end function fstr_get_norm_para_contact

  function fstr_get_x_norm_contact(hecMAT,fstrMAT,hecMESH) result(rhsX)
    use m_fstr
    use fstr_matrix_con_contact
    implicit none
    type(hecmwST_matrix), intent(in)                 :: hecMAT
    type(fstrST_matrix_contact_lagrange), intent(in) :: fstrMAT
    type(hecmwST_local_mesh), intent(in)             :: hecMESH
    real(kind=kreal)   ::  rhsX
    integer(kind=kint) :: nndof, npndof, i

    nndof = hecMAT%N * hecMAT%NDOF
    npndof = hecMAT%NP * hecMAT%NDOF
    rhsX = 0.d0
    do i=1,nndof
      rhsX = rhsX + hecMAT%X(i) ** 2
    end do
    do i=1,fstrMAT%num_lagrange
      rhsX = rhsX + hecMAT%X(npndof+i) ** 2
    end do
    call hecmw_allreduce_R1(hecMESH, rhsX, hecmw_sum)

  end function fstr_get_x_norm_contact

end module m_fstr_Residual
