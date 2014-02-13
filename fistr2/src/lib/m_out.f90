!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by K. SATOH (Advancesoft), X. YUAN (AdavanceSoft) !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!>   This module manages step infomation
!>
!>  \author     K. SATOH (Advancesoft), X. YUAN (AdavanceSoft)
!>  \date       2009/11/09
!>  \version    0.00
!!
!======================================================================!
module m_out
  use hecmw
  implicit none
  
  integer, parameter :: MAXOUT     = 100
  integer, parameter, private :: MAXNAMELEN = 30
  character(len=20), parameter, private :: OUTFILENAME="ifstr.out"

  include 'fstr_ctrl_util_f.inc'

  !> output information
  type output_info
    integer                   :: num_items          !< number of output items
    character(HECMW_NAME_LEN) :: grp_id_name        !< output nodal or elemental group name
    character(len=MAXNAMELEN) :: keyWord(MAXOUT)    !< name of output varibales
    integer                   :: vtype(MAXOUT)      !< 0: not defined; -1: scale -2:vector; -3:symmetric tensor -4:tensor
                                                    !< >0: user given number of viariable components
    integer                   :: location(MAXOUT)   !< -1: nodal output otherwise elemental output
                                                    !< 0: center; 1:average(default); 2:quadrature points; 3: nodes
    logical                   :: on(MAXOUT)         !< if output needed
    integer                   :: actn               !< 0: do nothing; 1: sum
  end type output_info

  !> output control such as output filename, output freqency etc  
  type t_output_ctrl
    character(len=HECMW_NAME_LEN) :: filename       !< output file name
    integer                       :: filenum        !< file number
    integer                       :: freqency       !< output frenqency
    type( output_info )           :: outinfo        !< info of output
  end type t_output_ctrl
  
contains

! ----------------------------------------------------------------------------
subroutine initOutInfo( outinfo )
  type( output_info ), intent(out) :: outinfo
  outinfo%grp_id_name      = "ALL"
  outinfo%keyWord(:)       =  ""
  outinfo%vtype(:)         =  0
  outinfo%location(:)      = -1
  outinfo%on(:)            = .false.
  
  outinfo%num_items        = 4
  ! following default value
  outinfo%keyWord(1) = "DISPLACEMENT"
  outinfo%vtype(1) = -2
  outinfo%on(1) = .true.
  
  outinfo%keyWord(2) = "REACTIONFORCE"
  outinfo%vtype(2) = -2
  outinfo%on(2) = .true.
  
  outinfo%keyWord(3) = "STRAIN"
  outinfo%vtype(3) = -3
  outinfo%location(3) = 1
  outinfo%on(3) = .false.
  
  outinfo%keyWord(4) = "STRESS"
  outinfo%vtype(4) = -3
  outinfo%location(4) = 1
  outinfo%on(4) = .false.
end subroutine initOutInfo


subroutine write_outinfo( fnum, nitem, outinfo, outdata )
  integer, intent(in)             :: fnum
  integer, intent(in)             :: nitem
  type( output_info ), intent(in) :: outinfo
  real(kind=kreal), intent(in)    :: outdata(:,:)
  
  integer  :: i,j, nsize, ncomp
  if( nitem>outinfo%num_items ) return
  if( .not. outinfo%on(nitem) ) return
  nsize = size( outdata, 1 )
  ncomp = size( outdata, 2 )
  write(fnum,'(a)') trim(outinfo%keyWord(nitem))
  do i=1,nsize
    write( fnum, * ) (outdata(i,j),j=1,ncomp)
  enddo
end subroutine 


integer function n_comp_valtype(vtype, ndim)
  integer, intent(in) :: vtype    !< variable type
  integer, intent(in) :: ndim     !< space dimension
  if( vtype>0 ) then
    n_comp_valtype = vtype
  else if( vtype==-1 ) then
    n_comp_valtype = 1
  else if( vtype==-2 ) then
    n_comp_valtype = ndim
  else if( vtype==-3 ) then
    n_comp_valtype = ndim*(ndim+1)/2
  else if( vtype==-4 ) then
    n_comp_valtype = ndim*ndim
  else 
    n_comp_valtype = 0
  endif
end function



! ----following t_output_ctrl------
subroutine fstr_init_outctrl(outctrl)
   type(t_output_ctrl), intent(out) :: outctrl
   outctrl%filename=trim(OUTFILENAME)
   outctrl%filenum = -1
   outctrl%freqency = 1
   call initOutInfo( outctrl%outinfo )
end subroutine  

subroutine fstr_copy_outctrl(outctrl1, outctrl2)
   type(t_output_ctrl), intent(out) :: outctrl1
   type(t_output_ctrl), intent(in)  :: outctrl2
   outctrl1%filename = outctrl2%filename
   outctrl1%filenum  = outctrl2%filenum
   outctrl1%freqency = outctrl2%freqency
end subroutine

logical function fstr_output_active( cstep, outctrl )
   integer, intent(in)             :: cstep  !< current step number
   type(t_output_ctrl), intent(in) :: outctrl
   fstr_output_active = .false.
   if( mod( cstep, outctrl%freqency )==0 ) fstr_output_active=.true.
end function

subroutine fstr_ctrl_get_filename( ctrl, ss )
  integer(kind=kint), intent(in)             :: ctrl
  character(len=HECMW_NAME_LEN), intent(out) :: ss

  integer(kind=kint) :: rcode

  ss = trim(OUTFILENAME)
  rcode = fstr_ctrl_get_param_ex( ctrl, 'FILENAME ',   '# ',  0, 'S', ss )
end subroutine

subroutine fstr_ctrl_get_output( ctrl, outctrl, islog, res, visual, femap )
  integer(kind=kint), intent(in)     :: ctrl
  type(t_output_ctrl), intent(inout) :: outctrl
  integer(kind=kint), intent(out)    :: islog
  integer(kind=kint), intent(out)    :: res, visual, femap

  integer(kind=kint) :: rcode, n
  character(len=HECMW_NAME_LEN) :: ss

  islog=0
  res=0;  visual=0;  femap=0
  if( fstr_ctrl_get_param_ex( ctrl, 'LOG ',   '# ',    0,   'E',   islog  )/= 0) return
  if( fstr_ctrl_get_param_ex( ctrl, 'RESULT ',  '# ',    0,   'E',   res    )/= 0) return
  if( fstr_ctrl_get_param_ex( ctrl, 'VISUAL ',  '# ',    0,   'E',   visual )/= 0) return
  if( fstr_ctrl_get_param_ex( ctrl, 'FEMAP ',   '# ',    0,   'E',   femap  )/= 0) return
  if( res==1 ) islog=1
  if( islog /= 1 ) return
  
  call fstr_init_outctrl(outctrl)
  ss =""
  rcode = fstr_ctrl_get_param_ex( ctrl, 'FILE ',   '# ',  0, 'S', ss )
  if( len(trim(ss))>0 ) then
    outctrl%filename = trim(ss)
    islog = 0
  endif
  n=-1
  rcode = fstr_ctrl_get_param_ex( ctrl, 'FREQENCY ',   '# ',  0, 'I', n )
  if( n>0 ) outctrl%freqency = n
    
end subroutine

subroutine print_output_ctrl( nfile, outctrl )
  integer, intent(in)                :: nfile
  type(t_output_ctrl), intent(inout) :: outctrl
  integer :: i
  write( nfile, *) trim(outctrl%filename),outctrl%filenum,outctrl%freqency,outctrl%outinfo%num_items
  do i=1,outctrl%outinfo%num_items
    write( nfile, *) trim(outctrl%outinfo%keyWord(i)),outctrl%outinfo%on(i),  &
                     outctrl%outinfo%vtype(i), outctrl%outinfo%location(i)
  enddo
end subroutine

end module m_out
