!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by X. YUAN                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief  This module manage surface elements in 3D
!!
!>  It provides basic definition of surface elements (trianglar and quadrilateral) 
!!  and functions for fetch its neighborhood information

!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/07/16
!>  \version    0.00

module mSurfElement
implicit none
integer, parameter, private :: kreal = kind(0.0d0)

integer, parameter       :: l_max_surface_node =20
integer, parameter       :: l_max_elem_node = 100

    !> Structure to define surface group
    type tSurfElement
      integer              :: eid                  !< elemental index(global)
      integer              :: etype                !< type of surface element
      integer, pointer     :: nodes(:)=>null()     !< nodes index(global)
      integer              :: n_neighbor           !< size of neighber(:)
      integer, pointer     :: neighbor(:)=>null()  !< index(global) of neighborhood tSurfElement
    end type
    
contains

  !> Initializer  
  subroutine initialize_surf( eid, etype, nsurf, surf )
    use elementInfo
    integer, intent(in)               :: eid
    integer, intent(in)               :: etype
    integer, intent(in)               :: nsurf
    type(tSurfElement), intent(inout) :: surf
    integer :: n, outtype, nodes(100)
    surf%eid = eid
    call getSubFace( etype, nsurf, outtype, nodes )  
    surf%etype = outtype
    n=getNumberOfNodes( outtype )
    allocate( surf%nodes(n) )
    surf%nodes(1:n)=nodes(1:n)
    n=getNumberOfSubface( outtype )
    surf%n_neighbor = 0
  end subroutine

  !> Memeory management subroutine
  subroutine finalize_surf( surf )
    type(tSurfElement), intent(inout) :: surf
    if( associated(surf%nodes) ) deallocate( surf%nodes )
    if( associated(surf%neighbor) ) deallocate( surf%neighbor )
  end subroutine

  !> Write out elemental surface
  subroutine write_surf( file, surf )
    integer, intent(in)               :: file   !< file number
    type(tSurfElement), intent(in)    :: surf   !< elemental surface
    integer :: i
    write(file,*) "Element:",surf%eid,"Surface type=",surf%etype
    if( associated( surf%nodes ) )   &
      write(file,*) ( surf%nodes(i),i=1,size(surf%nodes) )
    if( associated( surf%neighbor ) )   &
      write(file,*) ( surf%neighbor(i),i=1,surf%n_neighbor )
  end subroutine
  
  !> Find neighboring surface elements
  subroutine find_surface_neighbor( surf )
    use m_utilities
    type(tSurfElement), intent(inout) :: surf(:)
    integer :: i, j, ii,jj, nd1, nd2, nsurf
    nsurf = size(surf)
    do i=1,nsurf
      do ii=1, size(surf(i)%nodes)
        nd1 = surf(i)%nodes(ii)
        do j=1,nsurf
          if( i==j ) cycle
          if( associated(surf(i)%neighbor) ) then
            if ( any( surf(i)%neighbor==j ) ) cycle
          endif
          do jj=1, size(surf(j)%nodes)
            nd2 = surf(i)%nodes(jj)
            if( nd1==nd2 ) then
              surf(i)%n_neighbor = surf(i)%n_neighbor+1
              call insert_int2array( j, surf(i)%neighbor )
            endif
          enddo
        enddo		
      enddo
    enddo
  end subroutine

  !> Tracing next contact position
  integer function next_position( surf, cpos )
    use elementInfo
    type(tSurfElement), intent(in) :: surf      !< current surface element
    real(kind=kreal), intent(in)   :: cpos(2)   !< current position(local coordinate)

    integer          :: i
    real(kind=kreal) :: maxv(3)
    next_position = surf%eid
    if( .not. associated(surf%neighbor) ) return   ! do nothing when not initialized
    maxv(:) = 0.d0
    select case(surf%etype)
    case( fe_tri3n, fe_tri6n )
      if( all(cpos>0.d0) .and. all(cpos<1.d0) ) return
      if( size(surf%neighbor)<3 )  return
      if( all(cpos(:)>0.d0) ) return
      do i=1,3
        if( cpos(i)< 0.d0 ) maxv(i) = -cpos(i)
      enddo
      next_position = maxloc(maxv,1)
    case( fe_quad4n, fe_quad8n )
      if( all(cpos>-1.d0) .and. all(cpos<1.d0) ) return
      if( size(surf%neighbor)<4 ) return
      if( cpos(1)<-1.d0 .and. dabs(cpos(2))<1.d0 ) then
        next_position = 1
      elseif( cpos(1)>1.d0 .and. dabs(cpos(2))<1.d0 ) then
        next_position = 3
      elseif( dabs(cpos(1))<1.d0 .and. cpos(2)<-1.d0 ) then
        next_position = 2
      elseif( dabs(cpos(1))<1.d0 .and. cpos(2)>1.d0 ) then
        next_position = 4
      elseif( cpos(1)<-1.d0 .and. cpos(2)<-1.d0 ) then
        if( cpos(1)>cpos(2) ) then
          next_position = 2
        else
          next_position = 1
        endif
      elseif( cpos(1)<-1.d0 .and. cpos(2)>1.d0 ) then
        if( dabs(cpos(1))>cpos(2) ) then
          next_position = 1
        else
          next_position = 4
        endif
      elseif( cpos(1)>1.d0 .and. cpos(2)<-1.d0 ) then
        if( cpos(1)>dabs(cpos(2)) ) then
          next_position = 3
        else
          next_position = 2
        endif
      elseif( cpos(1)>1.d0 .and. cpos(2)>1.d0 ) then
        if( cpos(1)>cpos(2) ) then
          next_position =  3
        else
          next_position = 4
        endif
      endif
    case default
      stop "type of surface element not defined"
    end select
    next_position = surf%neighbor( next_position )

  end function next_position
  
end module
