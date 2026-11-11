!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
program solver_duplicate_test
  use m_fstr, only: kint, HECMW_FILENAME_LEN, ILOG
  use fstr_setup_util, only: fstr_param_pack
  use m_fstr_param_setup, only: fstr_setup_SOLVER
  use fstr_ctrl_util_f
  implicit none
  character(len=HECMW_FILENAME_LEN) :: filename, logfile
  integer(kind=kint) :: ctrl, rcode
  type(fstr_param_pack) :: P

  call get_command_argument(1, filename)
  call get_command_argument(2, logfile)
  open(unit=ILOG, file=trim(logfile), status='replace', action='write', iostat=rcode)
  if( rcode /= 0 ) stop 3
  if( fstr_ctrl_is_inp_file(filename) ) then
    ctrl = fstr_ctrl_open_inp(filename)
  else
    ctrl = fstr_ctrl_open(filename)
  endif
  if( ctrl < 0 ) stop 3
  rcode = fstr_ctrl_seek_next_header(ctrl)
  if( rcode == 0 ) stop 3
  P%STATE%slot%first = 2
  P%STATE%slot%last = 2
  call fstr_setup_SOLVER(ctrl, P)
  stop 1
end program solver_duplicate_test
