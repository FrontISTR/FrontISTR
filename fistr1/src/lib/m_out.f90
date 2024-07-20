!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module manages step information
module m_out
  use hecmw
  implicit none

  integer, parameter, private :: MAXOUT = 256
  integer, parameter, private :: MAXNAMELEN = 24
  character(len=20), parameter, private :: OUTFILENAME = "ifstr.out"

  include 'fstr_ctrl_util_f.inc'

  !> output information
  type output_info
    integer                   :: num_items          !< number of output items
    character(len=MAXNAMELEN) :: keyWord(MAXOUT)    !< name of output variables
    integer                   :: vtype(MAXOUT)      !< 0:not defined; -1:scaler; -2:vector; -3:symmetric tensor; -4:tensor
    !< >0: user given number of viariable components
    logical                   :: on(MAXOUT)         !< if output needed
    character(HECMW_NAME_LEN) :: grp_id_name        !< output nodal or elemental group name
    integer                   :: grp_id
    integer                   :: actn               !< 0: do nothing; 1: sum
  end type output_info

  !> output control such as output filename, output frequency etc.
  type t_output_ctrl
    character(len=HECMW_NAME_LEN) :: filename       !< output file name
    integer                       :: filenum        !< file number
    integer                       :: frequency       !< output frenqency
    type( output_info )           :: outinfo        !< info of output
  end type t_output_ctrl

contains

  ! ----------------------------------------------------------------------------
  subroutine initOutInfo( outinfo )
    type( output_info ), intent(out) :: outinfo

    outinfo%grp_id_name = "ALL"
    outinfo%grp_id      = -1
    outinfo%on(:)       = .false.
    outinfo%num_items   = 43

    outinfo%keyWord(1)  = "DISP"
    outinfo%vtype(1)    = -2
    outinfo%on(1)       = .true.

    outinfo%keyWord(2)  = "REACTION"
    outinfo%vtype(2)    = -2

    outinfo%keyWord(3)  = "NSTRAIN"
    outinfo%vtype(3)    = -3

    outinfo%keyWord(4)  = "NSTRESS"
    outinfo%vtype(4)    = -3
    outinfo%on(4)       = .true.

    outinfo%keyWord(5)  = "NMISES"
    outinfo%vtype(5)    = -1
    outinfo%on(5)       = .true.

    outinfo%keyWord(6)  = "ESTRAIN"
    outinfo%vtype(6)    = -3

    outinfo%keyWord(7)  = "ESTRESS"
    outinfo%vtype(7)    = -3

    outinfo%keyWord(8)  = "EMISES"
    outinfo%vtype(8)    = -1

    outinfo%keyWord(9)  = "ISTRAIN"
    outinfo%vtype(9)    = -3

    outinfo%keyWord(10) = "ISTRESS"
    outinfo%vtype(10)   = -3

    outinfo%keyWord(11) = "PL_ISTRAIN"
    outinfo%vtype(11)   = -1

    outinfo%keyWord(12) = "TH_NSTRAIN"
    outinfo%vtype(12)   = -3

    outinfo%keyWord(13) = "TH_ESTRAIN"
    outinfo%vtype(13)   = -3

    outinfo%keyWord(14) = "TH_ISTRAIN"
    outinfo%vtype(14)   = -3

    outinfo%keyWord(15) = "VEL"
    outinfo%vtype(15)   = -2

    outinfo%keyWord(16) = "ACC"
    outinfo%vtype(16)   = -2

    outinfo%keyWord(17) = "TEMP"
    outinfo%vtype(17)   = -1

    outinfo%keyWord(18) = "ROT"
    outinfo%vtype(18)   = -2
    outinfo%on(18)      = .true.

    outinfo%keyWord(19) = "PRINC_NSTRESS"
    outinfo%vtype(19)   = -2

    outinfo%keyWord(20) = "PRINC_ESTRESS"
    outinfo%vtype(20)   = -2

    outinfo%keyWord(21) = "PRINC_NSTRAIN"
    outinfo%vtype(21)   = -2

    outinfo%keyWord(22) = "PRINC_ESTRAIN"
    outinfo%vtype(22)   = -2

    outinfo%keyWord(23) = "PRINCV_NSTRESS"
    outinfo%vtype(23)   = -2

    outinfo%keyWord(24) = "PRINCV_ESTRESS"
    outinfo%vtype(24)   = -2

    outinfo%keyWord(25) = "PRINCV_NSTRAIN"
    outinfo%vtype(25)   = -2

    outinfo%keyWord(26) = "PRINCV_ESTRAIN"
    outinfo%vtype(26)   = -2

    outinfo%keyWord(27) = "SHELL_LAYER"
    outinfo%vtype(27)   = -1

    outinfo%keyWord(28) = "SHELL_SURFACE"
    outinfo%vtype(28)   = -1

    outinfo%keyWord(29) = "YIELD_RATIO"
    outinfo%vtype(29)   = -1

    outinfo%keyWord(30) = "CONTACT_NFORCE"
    outinfo%vtype(30)   = -2

    outinfo%keyWord(31) = "CONTACT_FRICTION"
    outinfo%vtype(31)   = -2

    outinfo%keyWord(32) = "CONTACT_RELVEL"
    outinfo%vtype(32)   = -2

    outinfo%keyWord(33) = "CONTACT_STATE"
    outinfo%vtype(33)   = -1

    outinfo%keyWord(34) = "MATERIAL_ID"
    outinfo%vtype(34)   = -1

    outinfo%keyWord(35) = "BEAM_NQM"
    outinfo%vtype(35)   = -5

    outinfo%keyWord(36) = "CONTACT_NTRACTION"
    outinfo%vtype(36)   = -2

    outinfo%keyWord(37) = "CONTACT_FTRACTION"
    outinfo%vtype(37)   = -2

    outinfo%keyWord(38) = "NODE_ID"
    outinfo%vtype(38)   = -1

    outinfo%keyWord(39) = "ELEM_ID"
    outinfo%vtype(39)   = -1

    outinfo%keyWord(40) = "SECTION_ID"
    outinfo%vtype(40)   = -1

    outinfo%keyWord(41) = "NOT ASSIGNED"
    outinfo%vtype(41)   = -1

    outinfo%keyWord(42) = "NOT ASSIGNED"
    outinfo%vtype(42)   = -1

    outinfo%keyWord(43)  = "PL_ESTRAIN"
    outinfo%vtype(43)    = -1

  end subroutine initOutInfo

  subroutine write_outinfo( fnum, nitem, outinfo, outdata )
    integer, intent(in)             :: fnum
    integer, intent(in)             :: nitem
    type( output_info ), intent(in) :: outinfo
    real(kind=kreal), intent(in)    :: outdata(:,:)

    integer  :: i, j, nsize, ncomp
    if( nitem>outinfo%num_items ) return
    if( .not. outinfo%on(nitem) ) return
    nsize = size( outdata, 1 )
    ncomp = size( outdata, 2 )
    write( fnum, '(a)' ) trim(outinfo%keyWord(nitem))
    do i=1,nsize
      write( fnum, * ) (outdata(i,j),j=1,ncomp)
    enddo
  end subroutine write_outinfo


  integer function n_comp_valtype( vtype, ndim )
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
      if(ndim == 4 .or. ndim == 6)n_comp_valtype=6
    else if( vtype==-4 ) then
      n_comp_valtype = ndim*ndim
    else if( vtype==-5 ) then
      n_comp_valtype = 12
    else
      n_comp_valtype = 0
    endif
  end function n_comp_valtype


  ! ----following t_output_ctrl------
  subroutine fstr_init_outctrl(outctrl)
    type(t_output_ctrl), intent(out) :: outctrl
    outctrl%filename= trim(OUTFILENAME)
    outctrl%filenum = -1
    outctrl%frequency= 1
    call initOutInfo( outctrl%outinfo )
  end subroutine

  subroutine fstr_copy_outctrl(outctrl1, outctrl2)
    type(t_output_ctrl), intent(out) :: outctrl1
    type(t_output_ctrl), intent(in)  :: outctrl2
    outctrl1%filename = outctrl2%filename
    outctrl1%filenum  = outctrl2%filenum
    outctrl1%frequency = outctrl2%frequency
  end subroutine

  logical function fstr_output_active( cstep, outctrl )
    integer, intent(in)             :: cstep  !< current step number
    type(t_output_ctrl), intent(in) :: outctrl
    fstr_output_active = .false.
    if( mod( cstep, outctrl%frequency )==0 ) fstr_output_active=.true.
  end function

  subroutine fstr_ctrl_get_filename( ctrl, ss )
    integer(kind=kint), intent(in)             :: ctrl
    character(len=HECMW_NAME_LEN), intent(out) :: ss
    integer(kind=kint) :: rcode
    ss = trim(OUTFILENAME)
    rcode = fstr_ctrl_get_param_ex( ctrl, 'FILENAME ', '# ', 0, 'S', ss )
  end subroutine

  subroutine fstr_ctrl_get_output( ctrl, outctrl, islog, res, visual, femap )
    integer(kind=kint), intent(in)     :: ctrl
    type(t_output_ctrl), intent(inout) :: outctrl
    integer(kind=kint), intent(out)    :: islog, res, visual, femap
    integer(kind=kint) :: rcode, n
    character(len=HECMW_NAME_LEN) :: ss
    islog=0; res=0; visual=0; femap=0
    if( fstr_ctrl_get_param_ex( ctrl, 'LOG ', '# ', 0, 'E', islog )/= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'RESULT ', '# ', 0, 'E', res )/= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'VISUAL ', '# ', 0, 'E', visual )/= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'UTABLE ', '# ', 0, 'E', femap  )/= 0 ) return

    call fstr_init_outctrl(outctrl)
    ss = ""
    rcode = fstr_ctrl_get_param_ex( ctrl, 'FILE ', '# ', 0, 'S', ss )
    if( len(trim(ss))>0 ) then
      outctrl%filename = trim(ss)
    endif
    outctrl%frequency = 1
    n = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'FREQUENCY ', '# ', 0, 'I', n )
    if( n>0 ) outctrl%frequency = n
  end subroutine

  subroutine print_output_ctrl( nfile, outctrl )
    integer, intent(in)                :: nfile
    type(t_output_ctrl), intent(inout) :: outctrl
    integer :: i
    write( nfile, *) trim(outctrl%filename),outctrl%filenum,outctrl%frequency,outctrl%outinfo%num_items
    do i=1,outctrl%outinfo%num_items
      write( nfile, *) trim(outctrl%outinfo%keyWord(i)),outctrl%outinfo%on(i),outctrl%outinfo%vtype(i)
    enddo
  end subroutine

end module m_out
