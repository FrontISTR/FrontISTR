!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Contact Analysis                                  !
!                                                                      !
!            Written by X. YUAN(AdavanceSoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!>  \brief   This module provides functions to calcualte contact stiff matrix
!>

!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/10/12
!>  \version    0.00
module mContact

use mContactDef
use hecmw
use m_fstr
implicit none

  private :: l_contact2mpc, l_tied2mpc
  integer, save :: n_contact_mpc
  logical, private :: active
  
  real(kind=kreal), save :: mu=1.d10  !< penalty, default value
  real(kind=kreal), save :: mut=1.d6  !< penalty along tangent direction
  real(kind=kreal), save :: cdotp=1.d3  !< mu=cdotp*maxval
  
  real(kind=kreal), save :: cgn=1.d-5 !< convergent condition of penetration
  real(kind=kreal), save :: cgt=1.d-3 !< convergent condition of relative tangent disp
  
  real(kind=kreal), save :: gnt(2)    !< 1:current avarage penetration;
                                      !< 2:current relative tangent displacement
  real(kind=kreal), save :: bakgnt(2) !< 1:current avarage penetration;
                                      !< 2:current relative tangent displacement

contains

  !> Write out the contact definition read from mesh file
  subroutine print_contatct_pair( file, pair )
       integer, intent(in)                      :: file
       type( hecmwST_contact_pair ), intent(in) :: pair
       
       integer :: i
       write(file,*) "Number of contact pair", pair%n_pair
       do i=1,pair%n_pair
         write(file,*) trim(pair%name(i)), pair%type(i), pair%slave_grp_id(i)  &
                       ,pair%master_grp_id(i)
       enddo
  end subroutine
  
  subroutine fstr_set_contact_penalty( maxv )
       real(kind=kreal), intent(in) :: maxv
       mu = cdotp * maxv
       if( gnt(1)<1.d-3 ) mu=cdotp*10.d0* maxv
       bakgnt = 0.d0
  end subroutine
 
  logical function fstr_is_contact_active()
       fstr_is_contact_active = active
  end function
  
  subroutine fstr_set_contact_active( a )
       logical, intent(in) :: a
       active = a
  end subroutine
  
  logical function fstr_is_contact_conv()
      fstr_is_contact_conv = .false.
      if( gnt(1)<cgn .and. gnt(2)<cgt ) fstr_is_contact_conv = .true.
      write(*,'(a,2e15.7)') "--Relative displacement in contact surface",gnt
    !  if( dabs( bakgnt(1)-gnt(1) )<cgn*1.d-2 .and.    &
    !      dabs( bakgnt(2)-gnt(2) )<cgn ) fstr_is_contact_conv = .true.
      bakgnt = gnt
  end function
  
  !> Contact state to equation conditions
  subroutine l_contact2mpc( contact, mpcs, nmpc )
       use fstr_ctrl_modifier
       type( tContact ), intent(in)         :: contact  !< current contact state
       type( hecmwST_mpc ), intent(inout)   :: mpcs     !< to who mpc be appended
       integer, intent(out)                 :: nmpc     !< number of mpc conditions appended
       integer, parameter          :: ndof = 3          ! 3D problem only, currently
       real(kind=kreal), parameter :: tol =1.d-10
       integer :: i, j, k, nn, csurf, nenode, etype, tdof
       integer :: nodes(l_max_surface_node*ndof+ndof), dofs(l_max_surface_node*ndof+ndof)
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
       integer, intent(out)                 :: nmpc     !< number of mpc conditions appended
       integer :: i, j, csurf, nenode, etype, tdof
       integer :: nodes(l_max_surface_node+1), dofs(l_max_surface_node+1)
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

  !> Contact states to equation conditions
  subroutine fstr_contact2mpc( contacts, mpcs )
      type( tContact ), intent(in)       :: contacts(:)  !< current contact state
      type( hecmwST_mpc ), intent(inout) :: mpcs         !< to who mpc be appended
      integer :: i, nmpc
      n_contact_mpc = 0
      do i=1,size(contacts)
        if( contacts(i)%algtype == CONTACTUNKNOWN ) cycle     ! not initialized
        if( contacts(i)%algtype == CONTACTFSLID ) then
           print *, "Cannot deal with finite slip problems by MPC!"
           cycle
        endif
        if( contacts(i)%algtype == CONTACTSSLID ) then
          call l_contact2mpc( contacts(i), mpcs, nmpc )
          n_contact_mpc = n_contact_mpc + nmpc
        elseif( contacts(i)%algtype == CONTACTTIED ) then
          call l_tied2mpc( contacts(i), mpcs, nmpc )
          n_contact_mpc = n_contact_mpc + nmpc
        endif
      enddo
  end subroutine
  
  !> Delete mpcs derived from contact conditions
  subroutine fstr_del_contactmpc( mpcs )
      use fstr_ctrl_modifier
      type( hecmwST_mpc ), intent(inout) :: mpcs       !<  mpcs to be modified
      call fstr_delete_mpc( n_contact_mpc, mpcs )
  end subroutine
  
  !> Print out mpc conditions
  subroutine fstr_write_mpc( file, mpcs )
      integer, intent(in)  :: file                  !<  file number
      type( hecmwST_mpc ), intent(in) :: mpcs       !<  mpcs to be printed

      integer :: i,j,n0,n1
      write(file, *) "Number of equation", mpcs%n_mpc
      do i=1,mpcs%n_mpc
        write(file,*) "--Equation",i
        n0=mpcs%mpc_index(i-1)+1
        n1=mpcs%mpc_index(i)
        write(file, *) n0,n1
        write(file,'(30i5)') (mpcs%mpc_item(j),j=n0,n1)
        write(file,'(30i5)') (mpcs%mpc_dof(j),j=n0,n1)
        write(file,'(30f7.2)') (mpcs%mpc_val(j),j=n0,n1),mpcs%mpc_const(i)
      enddo
  end subroutine
  
   !> Scanning contact state
  subroutine fstr_scan_contact_state( cstep, hecMESH, fstrSOLID, changed, B )
      integer, intent(in)                    :: cstep      !< current step number
      type( hecmwST_local_mesh ), intent(in) :: hecMESH    !< type mesh
      type(fstr_solid), intent(inout)        :: fstrSOLID  !< type fstr_solid
      logical, intent(inout)                 :: changed    !< if contact state changed
      real(kind=kreal), optional             :: B(:)       !< nodal force residual
      integer :: i, free2contact, contact2free, grpid
      logical :: iactive
	 ! P.A. We redefine fstrSOLID%ddunode as current coordinate of every nodes
    !  fstrSOLID%ddunode(:) = fstrSOLID%unode(:) + fstrSOLID%dunode(:) 
      fstrSOLID%ddunode(:) = hecMESH%node(:) + fstrSOLID%unode(:) + fstrSOLID%dunode(:) 
      active = .false.
      free2contact = 0;   contact2free = 0
      do i=1,size(fstrSOLID%contacts)
        grpid = fstrSOLID%contacts(i)%group
        if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) then
          call clear_contact_state(fstrSOLID%contacts(i));  cycle
        endif
        if( present(B) ) then
          call scan_contact_state( fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%QFORCE(:),  &
            free2contact,contact2free,hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), iactive,  &
			mu, B )
        else
          call scan_contact_state( fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%QFORCE(:),  &
            free2contact,contact2free,hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), iactive, mu )
        endif
        if( .not. active ) active = iactive
      enddo
      if( free2contact+contact2free > 0 ) then
        changed = .true.
    !    write(*,*) " --Contact state changes:",free2contact, contact2free
      endif
  end subroutine
  
  !> Update lagrangian multiplier
  subroutine fstr_update_contact1( hecMESH, fstrSOLID, incdisp, B )
      type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< mesh
      type(fstr_solid), intent(inout)        :: fstrSOLID   !< fstr_solid
      real(kind=kreal), intent(in)           :: incdisp(:)  !< disp increment of curr step
      real(kind=kreal), intent(inout)        :: B(:)        !< nodal force residual
	  
      integer :: i, nc
      gnt = 0.d0
      nc = size(fstrSOLID%contacts)
      do i=1, nc
     !   if( contacts(i)%mpced ) cycle
        call calcu_contact_force1( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
            , fstrSOLID%dunode(:), incdisp, fstrSOLID%contacts(i)%fcoeff, mu, mut, B, gnt )
      enddo
      if( nc>0 ) gnt = gnt/nc
  end subroutine
  
   !> Update lagrangian multiplier
  subroutine fstr_update_contact0( hecMESH, fstrSOLID, B )
      type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
      type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
      real(kind=kreal), intent(inout)        :: B(:)        !< nodal force residual
	  
      integer :: i
      do i=1, size(fstrSOLID%contacts)
     !   if( contacts(i)%mpced ) cycle
        call calcu_contact_force0( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
            , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, B )
      enddo
  end subroutine
  
 
  !> Update lagrangian multiplier
  subroutine fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchange )
      type( hecmwST_local_mesh ), intent(in) :: hecMESH
      type(fstr_solid), intent(inout)        :: fstrSOLID
      logical, intent(out)                   :: ctchange
	  
      integer :: i, nc
      gnt = 0.d0;  ctchange = .false.
      nc = size(fstrSOLID%contacts)
      do i=1, nc
     !   if( contacts(i)%mpced ) cycle
        call update_contact_multiplier( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
            , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, gnt, ctchange )
      enddo
      if( nc>0 ) gnt = gnt/nc
  end subroutine
  
  !> Introduce contact stiff into global stiff matrix or mpc conditions into hecMESH
  subroutine fstr_contactBC( iter, hecMESH, hecMAT, fstrSOLID )
      use fstr_ctrl_modifier
      integer(kind=kint), intent(in)           :: iter      !< NR iterations
      type (hecmwST_local_mesh), intent(inout) :: hecMESH   !< type mesh
      type (hecmwST_matrix), intent(inout)     :: hecMAT    !< type matrix
      type(fstr_solid), intent(inout)          :: fstrSOLID !< type fstr_solid
	  
      integer, parameter :: NDOF=3

      integer :: i, j, k, m, nnode, nd, etype
      integer :: ctsurf, ndLocal(l_max_surface_node+1)
      real(kind=kreal) :: factor, elecoord( 3, l_max_surface_node)
      real(kind=kreal) :: stiff(l_max_surface_node*3+3, l_max_surface_node*3+3)
      real(kind=kreal) :: nrlforce, force(l_max_surface_node*3+3)
   !   if( n_contact_mpc>0 ) call fstr_delete_mpc( n_contact_mpc, hecMESH%mpc )
   !   call fstr_contact2mpc( fstrSOLID%contacts, hecMESH%mpc )
	  ! temp. need modification
   !   do i=1,size(fstrSOLID%contacts) 
!	    fstrSOLID%contacts(i)%mpced = .true.
 !     enddo	 
 !    call fstr_write_mpc( 6, hecMESH%mpc )
 !print *,"Contact to mpc ok!",n_contact_mpc,"equations generated"
      factor = fstrSOLID%FACTOR(2)
      call hecmw_cmat_clear( hecMAT%cmat )
      do i=1,size(fstrSOLID%contacts)
        do j=1, size(fstrSOLID%contacts(i)%slave)
          if( fstrSOLID%contacts(i)%states(j)%state==CONTACTFREE ) cycle   ! free
          ctsurf = fstrSOLID%contacts(i)%states(j)%surface          ! contacting surface
          etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
          nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
          ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
          do k=1,nnode
            ndLocal(k+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(k)
            elecoord(1,k)=hecMESH%node(3*ndLocal(k+1)-2)
            elecoord(2,k)=hecMESH%node(3*ndLocal(k+1)-1)
            elecoord(3,k)=hecMESH%node(3*ndLocal(k+1))
          enddo
          call contact2stiff( fstrSOLID%contacts(i)%algtype, fstrSOLID%contacts(i)%states(j),    &
                              etype, nnode, elecoord(:,:), mu, mut, fstrSOLID%contacts(i)%fcoeff,    &
                              fstrSOLID%contacts(i)%symmetric, stiff(:,:), force(:) )
          call hecmw_mat_ass_contact( hecMAT,nnode+1,ndLocal,stiff )

          if( iter>1 ) cycle
        !  if( fstrSOLID%contacts(i)%states(j)%multiplier(1)/=0.d0 ) cycle  
          ! In case of new contact nodes, add enforced disp constraint
          fstrSOLID%contacts(i)%states(j)%wkdist = fstrSOLID%contacts(i)%states(j)%distance
          nrlforce = -mu*fstrSOLID%contacts(i)%states(j)%distance
          force(1:nnode*NDOF+NDOF) = force(1:nnode*NDOF+NDOF)*nrlforce
          do m=1,nnode+1
            nd = ndLocal(m)
            do k=1,NDOF
               hecMAT%B(NDOF*(nd-1)+k)=hecMAT%B(NDOF*(nd-1)+k)-force((m-1)*NDOF+k)
            enddo
          enddo

        enddo
      enddo

  end subroutine
  
end module mContact
