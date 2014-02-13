!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Data structure                                    !
!                                                                      !
!            Written by Xi Yuan (AdvanceSoft) with reference to FLIBS  !
!                 http://fortranwiki.org/fortran/show/FLIBS            !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides data structure table which would be
!>       dictionaried afterwards
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/08/28
!>  \version    0.00
!
!======================================================================!
module m_table
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)
  
!
! The length of the keys
!
  integer, parameter :: DICT_KEY_LENGTH = 128
  integer, parameter :: MAXINDEX = 20
  integer, parameter :: MAXNVAL = 10000

!
! The data that will be stored with each key
  type tTable
    integer                     :: ndepends
    integer                     :: tbcol, tbrow
    integer                     :: tbindex(MAXINDEX)
    real(kind=kreal), pointer   :: tbval(:,:)=>null()
  end type tTable

!
! The "null" value for these data
!
  type(tTable), parameter :: DICT_NULL = tTable( 0,0,0, 0,null() )
  
! overload of =
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE TableCopy
  END INTERFACE
  
! overload of ==
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE TableCompare
  END INTERFACE
  
  contains
  
    subroutine init_table( table, ndp, col, row, tbval )
      type( tTable ), intent(inout)   :: table
      integer, intent(in)             :: ndp, col, row
      real(kind=kreal), intent(in)    :: tbval(col,row)
      integer :: i,j
      real(kind=kreal) :: tbindexval(MAXNVAL)
      if( associated( table%tbval ) ) deallocate( table%tbval )
      table%tbcol = col
      table%tbrow = row
      table%ndepends = ndp
      allocate( table%tbval( col, row ) )
      do i=1, col
      do j=1, row
        table%tbval(i,j) = tbval(i,j)
      enddo
      enddo
 	  
      do i=1,ndp
        table%tbindex(i) =1
		tbindexval(1) = table%tbval(table%tbcol-i+1, 1)
	!	print *, 0, i, table%tbindex(i), tbindexval(1)

        do j=2, table%tbrow
          if( any( tbindexval(1:table%tbindex(i))==table%tbval(table%tbcol-i+1, j) ) ) cycle
          table%tbindex(i) = table%tbindex(i)+1
          tbindexval(table%tbindex(i)) = table%tbval(table%tbcol-i+1, j)
	!	  print *, 1, i, table%tbindex(i), tbindexval(table%tbindex(i))
        enddo
      enddo
      
	  j=1
      do i=1,ndp
        j=j*table%tbindex(i)
      enddo
      if( j/= row) stop "Error in table defnition!"
    !  print *, j,row, table%tbindex(1:ndp); pause
    end subroutine
	
    subroutine finalize_table( table )
      type( tTable ), intent(inout) :: table
      if( associated( table%tbval ) ) deallocate( table%tbval )
    end subroutine
	  
    subroutine print_table( table, fname )
      type( tTable ), intent(in) :: table
      integer, intent(in)        :: fname
      integer :: i,j
      write(fname,*) table%ndepends, table%tbcol, table%tbrow
      do i=1,table%tbrow
        write(fname,*) i,(table%tbval(j,i),j=1,table%tbcol)
      enddo
    end subroutine
	
	logical function TableCompare( lhs, rhs )
      TYPE(tTable), INTENT(IN) :: lhs
      TYPE(tTable), INTENT(IN) :: rhs
      integer :: i,j
      TableCompare = .false.
      if( lhs%ndepends /= rhs%ndepends ) return
      if( lhs%tbcol /= rhs%tbcol ) return
      if( lhs%tbrow /= rhs%tbrow ) return
      do i=1, rhs%ndepends
        if( lhs%tbindex(i) /= rhs%tbindex(i) ) return
      enddo
      do i=1, rhs%tbcol
      do j=1, rhs%tbrow
        if( lhs%tbval(i,j) /= rhs%tbval(i,j) ) return
      enddo
      enddo
      TableCompare = .true.
    END function
	
	SUBROUTINE TableCopy( lhs, rhs )
      TYPE(tTable), INTENT(OUT) :: lhs
      TYPE(tTable), INTENT(IN)  :: rhs
      integer :: i,j
      lhs%ndepends = rhs%ndepends
      lhs%tbcol = rhs%tbcol
      lhs%tbrow = rhs%tbrow
      lhs%tbindex(:) = rhs%tbindex(:)
      IF( associated( lhs%tbval ) ) deallocate( lhs%tbval )
      allocate( lhs%tbval( rhs%tbcol, rhs%tbrow ) )
      do i=1, rhs%tbcol
      do j=1, rhs%tbrow
        lhs%tbval(i,j) = rhs%tbval(i,j)
      enddo
      enddo
    END SUBROUTINE

end module m_table

!======================================================================!
!
!> \brief  This module provides data structure of dictionaried table list
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/08/28
!>  \version    0.00
!
!======================================================================!

module Table_DICTS
    use m_Table
    implicit none
	integer, parameter, private :: kreal = kind(0.0d0)
	
	private :: GetTableGrad, GetTableData

    include "dictionary.f90"


	!> fetch a data table itself.
	!! P.A. it sould be deleted by users of this subroutine
    subroutine fetch_Table( key, dict, dicval, ierr )
      character(len=*), intent(in)   :: key     !< parameter key
      type(DICT_STRUCT), pointer     :: dict    !< data table
      logical, intent(out)           :: ierr
      type(tTable)                :: dicval

      dicval = dict_get_key( dict, key )
      ierr = .false.
      if( dicval==DICT_NULL ) then
        ierr=.true.;  return
      endif
      
    end subroutine
	
	
	!> fetch a data table row
    integer function fetch_TableRow( key, dict )
      character(len=*), intent(in)   :: key     !< parameter key
      type(DICT_STRUCT), pointer     :: dict    !< data table
      type(tTable)                :: dicval

      dicval = dict_get_key( dict, key )
      fetch_TableRow = -1
      if( dicval==DICT_NULL ) return
      fetch_TableRow = dicval%tbrow	
      call finalize_table( dicval )
    end function
	
	!> fetch a data by interpolation (This subroutine is specified for calculating
	!! temperature dependent hardening coefficient)	 
    subroutine fetch_TableGrad( key, a, dict, outa, ierr )
      character(len=*), intent(in)   :: key     !< parameter key
      real(kind=kreal), intent(in)   :: a(:)    !< automatic variables
      type(DICT_STRUCT), pointer     :: dict    !< data table
      real(kind=kreal), intent(out)  :: outa    !< gradient
      logical, intent(out)           :: ierr
	  
      type(tTable)                :: dicval
      integer          :: i, j, na, dd, crow, cindex
      integer          :: cindex1(MAXINDEX), cindex2(MAXINDEX)
      dicval = dict_get_key( dict, key )
      ierr = .false.
      if( dicval==DICT_NULL ) then
        ierr=.true.;  return
      endif
	  
      cindex = 1;  crow=1
      dd = dicval%tbrow	
	  na = 1
      if( size(a) > dicval%ndepends ) then   ! in case no dependent definition
        na = size(a) - dicval%ndepends+1
      endif

      call GetTableGrad( a(na:), cindex, dicval, dd, crow, outa )
      call finalize_table( dicval )
      
    end subroutine
	
	!> fetch a grad value by interpolation  
    recursive subroutine GetTableGrad( a, cindex, table, dd, crow, outa )
      real(kind=kreal), intent(in)   :: a(:)
      integer, intent(inout)         :: cindex
      type(tTable)                :: table
      integer, intent(inout)         :: dd, crow
      real(kind=kreal), intent(out)  :: outa

      integer          :: i, j, na, nn, ccol, ddd
      real(kind=kreal) :: cval, val1, val2, lambda
      logical          :: isok

      ddd = dd / table%tbindex(cindex)
      ccol = table%tbcol-cindex+1

      if( ddd==1 ) then
        if( a(cindex)<table%tbval(2, crow) ) then
           outa = 0.d0
        elseif( a(cindex)>=table%tbval(2, crow+1) ) then
           outa = 0.d0
        else
          do i=crow, crow+dd-1
            if( a(cindex)>=table%tbval(2, i) .and. a(cindex)<table%tbval(2, i+1) ) then 
              outa = (table%tbval(1, i+1)-table%tbval(1, i))/(table%tbval(2, i+1)-table%tbval(2, i))
              exit
            endif
          enddo
        endif
      else
        if( a(cindex)<=table%tbval(ccol, crow) ) then
            outa = 0.d0
        elseif( a(cindex)>=table%tbval(ccol, crow+dd-1) ) then
            outa = 0.d0
        else
          do i=crow, crow+dd-1, ddd
            if( a(cindex)==table%tbval(ccol, i) ) then
              crow = i
              cindex = cindex+1
              dd = ddd
              call GetTableGrad( a, cindex, table, dd, crow, outa )	
              exit	
            elseif( a(cindex)==table%tbval(ccol, i+ddd) ) then
              crow = i+ddd
              cindex = cindex+1
              dd = ddd
              call GetTableGrad( a, cindex, table, dd, crow, outa )	
              exit
            elseif( a(cindex)>table%tbval(ccol, i) .and. a(cindex)<table%tbval(ccol, i+ddd) ) then
		      crow=i
              cindex = cindex+1
              dd = ddd
              call GetTableGrad( a, cindex, table, dd, crow, val1 )
              crow = i+ddd
              call GetTableGrad( a, cindex, table, dd, crow, val2 )	
              lambda = (a(cindex-1)-table%tbval(ccol, i))/(table%tbval(ccol, crow)-table%tbval(ccol, i))
              outa = (1.d0-lambda)*val1+ lambda* val2
              exit
            endif
          enddo
        endif
      endif
      
    end subroutine
	
	!> fetch a data by interpolation (This subroutine is specified for calculating
	!! temperature dependent hardening coefficient)	 
    subroutine fetch_TableData( key, dict, outa, ierr, a )
      character(len=*), intent(in)   :: key     !< parameter key
      type(DICT_STRUCT), pointer     :: dict    !< data table
      real(kind=kreal), intent(out)  :: outa(:) !< output data
      logical, intent(out)           :: ierr    !< erro message
      real(kind=kreal), intent(in), optional   :: a(:)    !< automatic variables
	  
      type(tTable)                :: dicval
      integer          :: nval, na, dd, crow, cindex
	  
      dicval = dict_get_key( dict, key )
      ierr = .false.
      if( dicval==DICT_NULL ) then
        ierr=.true.;  return
      endif
	!  call print_table( dicval, 6 )
	  
      nval = dicval%tbcol-dicval%ndepends
      if( nval /= size(outa) ) then
        ierr=.true.;  return
      endif
      if( dicval%tbrow==1 ) then      
         outa(:) =dicval%tbval(1:nval, 1); return
      endif	 
      if( .not. present(a) ) then
         outa(:) =dicval%tbval(1:nval, 1); return
      endif
	  
      cindex = 1;  crow=1
      dd = dicval%tbrow	
	  na = 1
      if( size(a) > dicval%ndepends ) then   ! in case no dependent definition
        na = size(a) - dicval%ndepends+1
      endif
    !  if( size(a) < dicval%ndepends ) then   ! in case take no consider of dependent definition
    !    na = dicval%ndepends- size(a)+1
    !    cindex = na
    !  endif
	       
      call GetTableData( a(na:), cindex, dicval, dd, crow, outa )
      call finalize_table( dicval )
      
    end subroutine
	
	!> fetch a data value by interpolation  
    recursive subroutine GetTableData( a, cindex, table, dd, crow, outa )
      real(kind=kreal), intent(in)   :: a(:)
      integer, intent(inout)         :: cindex
      type(tTable)                :: table
      integer, intent(inout)         :: dd, crow
      real(kind=kreal), intent(out)  :: outa(:)
	  
      integer          :: i, j, na, nn, ccol, ddd, nval
      real(kind=kreal) :: cval, lambda, val1(MAXINDEX), val2(MAXINDEX)
      logical          :: isok

      ddd = dd / table%tbindex(cindex)
      ccol = table%tbcol-cindex+1
	  nval = size(outa)

      if( ddd==1 ) then
        if( a(cindex)<=table%tbval(ccol, crow) ) then
           outa(:) = table%tbval(1:ccol-1, crow)
        elseif( a(cindex)>=table%tbval(ccol, crow+1) ) then
           outa(:) = table%tbval(1:ccol-1, crow+dd-1)
        else
          do i=crow, crow+dd-1
            if( a(cindex)>table%tbval(ccol, i) .and. a(cindex)<table%tbval(ccol, i+1) ) then 
              lambda = (a(cindex)-table%tbval(ccol, i))/(table%tbval(ccol, i+1)-table%tbval(ccol, i))
              outa(:) = (1.d0-lambda)*table%tbval(1:ccol-1, i)+ lambda* table%tbval(1:ccol-1, i+1) 
              exit
            endif
          enddo
        endif
      else
        if( a(cindex)<=table%tbval(ccol, crow) ) then
            cindex = cindex+1
            dd = ddd
            call GetTableData( a, cindex, table, dd, crow, outa )	
        elseif( a(cindex)>=table%tbval(ccol, crow+dd-1) ) then
            crow = crow+dd-1
            cindex = cindex+1
            dd = ddd
            call GetTableData( a, cindex, table, dd, crow, outa )	
        else
          do i=crow, crow+dd-1, ddd
            if( a(cindex)==table%tbval(ccol, i) ) then
              crow = i
              cindex = cindex+1
              dd = ddd
              call GetTableData( a, cindex, table, dd, crow, outa )	
              exit	
            elseif( a(cindex)==table%tbval(ccol, i+ddd) ) then
              crow = i+ddd
              cindex = cindex+1
              dd = ddd
              call GetTableData( a, cindex, table, dd, crow, outa )	
              exit
            elseif( a(cindex)>table%tbval(ccol, i) .and. a(cindex)<table%tbval(ccol, i+ddd) ) then
		      crow=i
              cindex = cindex+1
              dd = ddd
              call GetTableData( a, cindex, table, dd, crow, val1(1:nval) )
              crow = i+ddd
              call GetTableData( a, cindex, table, dd, crow, val2(1:nval) )	
              lambda = (a(cindex-1)-table%tbval(ccol, i))/(table%tbval(ccol, crow)-table%tbval(ccol, i))
              outa(:) = (1.d0-lambda)*val1(1:nval)+ lambda* val2(1:nval)
              exit
            endif
          enddo
        endif
      endif
      
    end subroutine


	!> Print our the contents of a dictionary	
    subroutine print_TableData( dict, fname )
      type(DICT_STRUCT), pointer     :: dict
      integer, intent(in)            :: fname
	  
      type(tTable)             :: dicval
      type(LINKED_LIST), pointer  :: current
      integer :: i
      do i = 1,size(dict%table)
        if ( associated( dict%table(i)%list ) ) then
           current => dict%table(i)%list
           do while ( associated(current) )
             if( trim(current%data%key) /= 'INIT' ) then
               write( fname, * ) trim(current%data%key)
               call print_table( current%data%value, fname )  
             endif
             current => current%next
           enddo
        endif
      enddo
    end subroutine
	
end module Table_DICTS

  		
