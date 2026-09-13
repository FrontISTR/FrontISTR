!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief MPC (Multi-Point Constraint) processing for contact analysis
!>
!> This module provides internal functions to convert contact constraints to MPC form.
!> All functions are private and used only within mContact module.
module m_fstr_contact_mpc
  use hecmw
  use m_fstr
  use mContactDef
  use m_fstr_contact_element
  implicit none

  private
  public :: contact2mpcval
  public :: l_contact2mpc
  public :: l_tied2mpc

contains

  subroutine contact2mpcval( cstate, etype, nnode, mpcval )
    type(tContactState), intent(in) :: cstate !< contact state
    integer, intent(in)             :: etype !< type of contacting surface
    integer, intent(in)             :: nnode !< number of elemental nodes
    real(kind=kreal), intent(out)   :: mpcval(nnode*3 + 4) !< MPC constraint

    integer          :: i,j
    real(kind=kreal) :: shapefunc(nnode)

    call getShapeFunc( etype, cstate%lpos(1:2), shapefunc )
    mpcval(1:3) = cstate%direction(1:3)
    do i=1,nnode
      do j=1,3
        mpcval( i*3+j ) = -cstate%direction(j)*shapefunc(i)
      enddo
    enddo
    mpcval( 3*nnode+4 )=cstate%distance
  end subroutine

  subroutine l_contact2mpc( contact, mpcs, nmpc )
    use fstr_ctrl_modifier
    type( tContact ), intent(in)         :: contact  !< current contact state
    type( hecmwST_mpc ), intent(inout)   :: mpcs     !< to who mpc be appended
    integer(kind=kint), intent(out)      :: nmpc     !< number of mpc conditions appended
    integer(kind=kint), parameter :: ndof = 3        ! 3D problem only, currently
    real(kind=kreal), parameter   :: tol =1.d-10
    integer(kind=kint) :: i, j, k, nn, csurf, nenode, etype, tdof
    integer(kind=kint) :: nodes(l_max_surface_node*ndof+ndof), dofs(l_max_surface_node*ndof+ndof)
    real(kind=kreal) :: values(l_max_surface_node*ndof+ndof+1),val(l_max_surface_node*ndof+ndof+1)
    nmpc=0
    do i=1,size(contact%states)
      if( contact%states(i)%state == -1 ) cycle    ! in free
      csurf = contact%states(i)%surface
      if( csurf<=0 ) stop "error in contact state"
      etype =  contact%master(csurf)%etype
      nenode = size(contact%master(csurf)%nodes)
      tdof = nenode*ndof+ndof
      call contact2mpcval( contact%states(i), etype, nenode, values(1:tdof+1) )
      tdof = 0
      do j=1,ndof
        if( dabs(values(j))<tol ) cycle
        tdof = tdof+1
        nodes(tdof) = contact%slave(i)
        dofs(tdof) = j
        val(tdof) = values(j)
      enddo
      do j=1,nenode
        nn =  contact%master(csurf)%nodes(j)
        nodes( j*ndof+1:j*ndof+ndof ) =  nn
        do k=1,ndof
          if( dabs(values(j*ndof+k)) < tol ) cycle
          tdof=tdof+1
          nodes(tdof)=nn
          dofs(tdof ) = k
          val(tdof)=values(j*ndof+k)
        enddo
      enddo
      val(tdof+1) = values(nenode*ndof+ndof+1)

      call fstr_append_mpc( tdof, nodes(1:tdof), dofs(1:tdof), val(1:tdof+1), mpcs )
      nmpc=nmpc+1
    enddo
  end subroutine l_contact2mpc

  !> Rigid connect condition to equation conditions
  subroutine l_tied2mpc( contact, mpcs, nmpc )
    use fstr_ctrl_modifier
    type( tContact ), intent(in)         :: contact  !< current contact state
    type( hecmwST_mpc ), intent(inout)   :: mpcs     !< to who mpc be appended
    integer(kind=kint), intent(out)      :: nmpc     !< number of mpc conditions appended
    integer(kind=kint) :: i, j, csurf, nenode, etype, tdof
    integer(kind=kint) :: nodes(l_max_surface_node+1), dofs(l_max_surface_node+1)
    real(kind=kreal) :: values(l_max_surface_node+2)
    nmpc=0
    do i=1,size(contact%slave)
      csurf = contact%states(i)%surface
      if( csurf<=0 ) cycle                         ! contactor not exists
      nenode = size(contact%master(csurf)%nodes)
      tdof = nenode+1
      nodes(1) = contact%slave(i)
      nodes( 2:tdof ) =  contact%master(csurf)%nodes(:)
      values(1) = -1.d0
      values(2:tdof) = 1.d0
      values(tdof+1) = 0.d0
      etype =  contact%master(csurf)%etype
      do j=1,3
        dofs(1:tdof) = j
        call fstr_append_mpc( tdof, nodes(1:tdof), dofs(1:tdof), values(1:tdof+1), mpcs )
        nmpc=nmpc+1
      enddo
    enddo
  end subroutine l_tied2mpc

end module m_fstr_contact_mpc
