!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module manage surface elements in 3D
!>  It provides basic definition of surface elements (trianglar and quadrilateral)
!!  and functions for fetch its neighborhood information
module mSurfElement
  use hecmw, only: kint, kreal
  implicit none

  integer(kind=kint), parameter       :: l_max_surface_node =20
  integer(kind=kint), parameter       :: l_max_elem_node = 100

  integer(kind=kint), parameter       :: N_NEIGHBOR_MAX_INIT = 8

  !> Structure to define surface group
  type tSurfElement
    integer(kind=kint)              :: eid                  !< elemental index(global)
    integer(kind=kint)              :: etype                !< type of surface element
    integer(kind=kint), pointer     :: nodes(:)=>null()     !< nodes index(global)
    integer(kind=kint)              :: n_neighbor           !< number of neighborhood tSurfElement
    integer(kind=kint)              :: n_neighbor_max       !< maximum size of neighber(:)
    integer(kind=kint), pointer     :: neighbor(:)=>null()  !< index(global) of neighborhood tSurfElement
    real(kind=kreal)                :: reflen               !< reference length
    real(kind=kreal)                :: xavg(3)              !< current coordinate of element center
    real(kind=kreal)                :: dmax                 !< half length of edge of cube that include surf
  end type tSurfElement

contains

  !> Initializer
  subroutine initialize_surf( eid, etype, nsurf, surf )
    use elementInfo
    integer(kind=kint), intent(in)               :: eid
    integer(kind=kint), intent(in)               :: etype
    integer(kind=kint), intent(in)               :: nsurf
    type(tSurfElement), intent(inout) :: surf
    integer(kind=kint) :: n, outtype, nodes(100)
    surf%eid = eid
    call getSubFace( etype, nsurf, outtype, nodes )
    surf%etype = outtype
    n=getNumberOfNodes( outtype )
    allocate( surf%nodes(n) )
    surf%nodes(1:n)=nodes(1:n)
    n=getNumberOfSubface( outtype )
    surf%n_neighbor = 0
    surf%n_neighbor_max = 0
    surf%reflen  = -1.d0
    surf%xavg(:) =  0.d0
    surf%dmax    = -1.d0
  end subroutine

  subroutine initialize_surf_reflen( surf, coord )
    use elementInfo
    type(tSurfElement), intent(inout) :: surf
    real(kind=kreal), intent(in)      :: coord(:)
    real(kind=kreal) :: elem(3, l_max_surface_node), r0(2)
    integer(kind=kint) :: nn, i, iss
    nn = size(surf%nodes)
    do i=1,nn
      iss = surf%nodes(i)
      elem(1:3,i) = coord(3*iss-2:3*iss)
    enddo
    call getElementCenter( surf%etype, r0 )
    surf%reflen = getReferenceLength( surf%etype, nn, r0, elem )
  end subroutine initialize_surf_reflen

  !> Memeory management subroutine
  subroutine finalize_surf( surf )
    type(tSurfElement), intent(inout) :: surf
    if( associated(surf%nodes) ) deallocate( surf%nodes )
    if( associated(surf%neighbor) ) deallocate( surf%neighbor )
  end subroutine

  !> Write out elemental surface
  subroutine write_surf( file, surf )
    integer(kind=kint), intent(in)    :: file   !< file number
    type(tSurfElement), intent(in)    :: surf   !< elemental surface
    integer(kind=kint) :: i
    write(file,*) "Element:",surf%eid,"Surface type=",surf%etype
    if( associated( surf%nodes ) )   &
      write(file,*) ( surf%nodes(i),i=1,size(surf%nodes) )
    if( associated( surf%neighbor ) )   &
      write(file,*) ( surf%neighbor(i),i=1,surf%n_neighbor )
  end subroutine

  !> Find neighboring surface elements
  subroutine find_surface_neighbor( surf )
    use m_utilities
    use hecmw_util
    type(tSurfElement), intent(inout) :: surf(:)
    integer(kind=kint) :: i, j, ii,jj, nd1, nd2, nsurf
    integer(kind=kint) :: k, oldsize, newsize
    integer(kind=kint), pointer :: dumarray(:) => null()

    nsurf = size(surf)

    !$omp parallel do default(none),private(i,ii,nd1,j,jj,nd2,oldsize,newsize,dumarray,k), &
      !$omp&shared(nsurf,surf)
    do i=1,nsurf
      JLOOP: do j=1,nsurf
        if( i==j ) cycle
        if( associated(surf(i)%neighbor) ) then
          if ( any( surf(i)%neighbor(1:surf(i)%n_neighbor)==j ) ) cycle
        endif
        do ii=1, size(surf(i)%nodes)
          nd1 = surf(i)%nodes(ii)
          do jj=1, size(surf(j)%nodes)
            nd2 = surf(j)%nodes(jj)
            if( nd1==nd2 ) then
              surf(i)%n_neighbor = surf(i)%n_neighbor+1

              if( .not. associated(surf(i)%neighbor) ) then
                allocate( surf(i)%neighbor(N_NEIGHBOR_MAX_INIT) )
                surf(i)%n_neighbor_max = N_NEIGHBOR_MAX_INIT
                surf(i)%neighbor(1) = j
              else if( surf(i)%n_neighbor > surf(i)%n_neighbor_max ) then
                oldsize = surf(i)%n_neighbor_max
                newsize = oldsize * 2
                dumarray => surf(i)%neighbor
                allocate( surf(i)%neighbor(newsize) )
                surf(i)%n_neighbor_max = newsize
                do k=1,oldsize
                  surf(i)%neighbor(k) = dumarray(k)
                enddo
                surf(i)%neighbor(oldsize+1) = j
                deallocate( dumarray )
              else
                surf(i)%neighbor(surf(i)%n_neighbor) = j
              endif

              cycle JLOOP
            endif
          enddo
        enddo
      enddo JLOOP
    enddo
    !$omp end parallel do

  end subroutine

  !> Tracing next contact position
  integer(kind=kint) function next_position( surf, cpos )
    use elementInfo
    type(tSurfElement), intent(in) :: surf      !< current surface element
    real(kind=kreal), intent(in)   :: cpos(2)   !< current position(local coordinate)

    integer(kind=kint)          :: i
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

  subroutine update_surface_box_info( surf, currpos )
    type(tSurfElement), intent(inout) :: surf    !< current surface element
    real(kind=kreal), intent(in) :: currpos(:)   !< current coordinate of all nodes
    integer(kind=kint) :: nn, i, iss
    real(kind=kreal) :: elem(3,l_max_surface_node),xmin(3), xmax(3), xsum(3), lx(3)
    nn = size(surf%nodes)
    do i=1,nn
      iss = surf%nodes(i)
      elem(1:3,i) = currpos(3*iss-2:3*iss)
    enddo
    do i=1,3
      xmin(i) = minval(elem(i,1:nn))
      xmax(i) = maxval(elem(i,1:nn))
      xsum(i) = sum(elem(i,1:nn))
    enddo
    surf%xavg(:) = xsum(:) / nn
    lx(:) = xmax(:) - xmin(:)
    surf%dmax = maxval(lx) * 0.5d0
  end subroutine update_surface_box_info

  logical function is_in_surface_box(surf, x0, exp_rate)
    type(tSurfElement), intent(inout) :: surf  !< current surface element
    real(kind=kreal), intent(in) :: x0(3)      !< coordinate of slave node
    real(kind=kreal), intent(in) :: exp_rate   !< expansion rate (>1.0)
    real(kind=kreal) :: dist(3), er
    integer(kind=kint) :: i
    er = max(exp_rate, 1.d0)
    dist(:) = abs(x0(:) - surf%xavg(:))
    if ( maxval(dist(:)) < surf%dmax * er ) then
      is_in_surface_box = .true.
    else
      is_in_surface_box = .false.
    endif
  end function is_in_surface_box

end module mSurfElement
