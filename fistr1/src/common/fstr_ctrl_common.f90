!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains fstr control file data obtaining functions

module fstr_ctrl_common
  use m_fstr
  use hecmw
  use mContact
  use m_timepoint

  implicit none

  include 'fstr_ctrl_util_f.inc'

  private :: pc_strupr

contains

  subroutine pc_strupr( s )
    implicit none
    character(*) :: s
    integer :: i, n, a, da

    n = len_trim(s)
    da = iachar('a') - iachar('A')
    do i = 1, n
      a = iachar(s(i:i))
      if( a > iachar('Z')) then
        a = a - da
        s(i:i) = achar(a)
      end if
    end do
  end subroutine pc_strupr


  !> Read in !SOLUTION
  function fstr_ctrl_get_SOLUTION( ctrl, type, nlgeom )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: type
    logical            :: nlgeom
    integer(kind=kint) :: fstr_ctrl_get_SOLUTION

    integer(kind=kint) :: ipt
    character(len=80) :: s

    fstr_ctrl_get_SOLUTION = -1

    s = 'ELEMCHECK,STATIC,EIGEN,HEAT,DYNAMIC,NLSTATIC,STATICEIGEN,NZPROF '
    if( fstr_ctrl_get_param_ex( ctrl,      'TYPE ',     s,    1,   'P',  type )/= 0) return
    type = type -1

    ipt=0
    if( fstr_ctrl_get_param_ex( ctrl,    'NONLINEAR ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 .and. ( type == kstSTATIC .or. type == kstDYNAMIC )) nlgeom = .true.

    if( type == 5 ) then !if type == NLSTATIC
      type = kstSTATIC
      nlgeom = .true.
    end if
    if( type == kstSTATICEIGEN ) nlgeom = .true.

    fstr_ctrl_get_SOLUTION = 0
  end function fstr_ctrl_get_SOLUTION


  !> Read in !SOLVER
  function fstr_ctrl_get_SOLVER( ctrl, method, precond, nset, iterlog, timelog, steplog, nier, &
      iterpremax, nrest, nBFGS, scaling, &
      dumptype, dumpexit, usejad, ncolor_in, mpc_method, estcond, method2, recyclepre, &
      solver_opt, &
      resid, singma_diag, sigma, thresh, filter )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: method
    integer(kind=kint) :: precond
    integer(kind=kint) :: nset
    integer(kind=kint) :: iterlog
    integer(kind=kint) :: timelog
    integer(kind=kint) :: steplog
    integer(kind=kint) :: nier
    integer(kind=kint) :: iterpremax
    integer(kind=kint) :: nrest
    integer(kind=kint) :: nBFGS
    integer(kind=kint) :: scaling
    integer(kind=kint) :: dumptype
    integer(kind=kint) :: dumpexit
    integer(kind=kint) :: usejad
    integer(kind=kint) :: ncolor_in
    integer(kind=kint) :: mpc_method
    integer(kind=kint) :: estcond
    integer(kind=kint) :: method2
    integer(kind=kint) :: recyclepre
    integer(kind=kint) :: solver_opt(10)
    real(kind=kreal) :: resid
    real(kind=kreal) :: singma_diag
    real(kind=kreal) :: sigma
    real(kind=kreal) :: thresh
    real(kind=kreal) :: filter
    integer(kind=kint) :: fstr_ctrl_get_SOLVER

    character(92) :: mlist = '1,2,3,4,101,CG,BiCGSTAB,GMRES,GPBiCG,GMRESR,GMRESREN,DIRECT,DIRECTmkl,DIRECTlag,MUMPS,MKL ' 
    !character(92) :: mlist = '1,2,3,4,5,101,CG,BiCGSTAB,GMRES,GPBiCG,DIRECT,DIRECTmkl,DIRECTlag,MUMPS,MKL '
    character(24) :: dlist = '0,1,2,3,NONE,MM,CSR,BSR '

    integer(kind=kint) :: number_number = 5
    integer(kind=kint) :: indirect_number = 6 ! GMRESR and GMRESREN need to be added
    integer(kind=kint) :: iter, time, sclg, dmpt, dmpx, usjd, step

    fstr_ctrl_get_SOLVER = -1

    iter = iterlog+1
    time = timelog+1
    step = steplog+1
    sclg = scaling+1
    dmpt = dumptype+1
    dmpx = dumpexit+1
    usjd = usejad+1
    !* parameter in header line -----------------------------------------------------------------*!

    ! JP-0
    if( fstr_ctrl_get_param_ex( ctrl, 'METHOD ',   mlist,              1,   'P',   method  ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'PRECOND ', '1,2,3,4,5,6,7,8,9,10,11,12,20,21,30,31,32 ' ,0, 'I', precond ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'NSET ',    '0,-1,+1 ',          0,   'I',   nset    ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'ITERLOG ', 'NO,YES ',           0,   'P',   iter ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'TIMELOG ', 'NO,YES,VERBOSE ',   0,   'P',   time ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'STEPLOG ', 'NO,YES ',           0,   'P',   step ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'SCALING ', 'NO,YES ',           0,   'P',   sclg ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'DUMPTYPE ', dlist,              0,   'P',   dmpt ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'DUMPEXIT ','NO,YES ',           0,   'P',   dmpx ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'USEJAD '  ,'NO,YES ',           0,   'P',   usjd ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'MPCMETHOD ','# ',               0, 'I',mpc_method) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'ESTCOND '  ,'# ',               0,   'I',estcond ) /= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'METHOD2 ',  mlist,              0,   'P',   method2 ) /= 0) return
    ! JP-1
    if( method > number_number ) then  ! JP-2
      method = method - number_number
      if( method > indirect_number ) then
        ! JP-3
        method = method - indirect_number + 100
        if( method == 103 ) method = 101 ! DIRECTlag => DIRECT
        if( method == 105 ) method = 102 ! MKL => DIRECTmkl
      end if
    end if
    if( method2 > number_number ) then  ! JP-2
      method2 = method2 - number_number
      if( method2 > indirect_number ) then
        ! JP-3
        method2 = method2 - indirect_number + 100
      end if
    end if

    dumptype = dmpt - 1
    if( dumptype >= 4 ) then
      dumptype = dumptype - 4
    end if

    !* data --------------------------------------------------------------------------------------- *!
    ! JP-4
    if( fstr_ctrl_get_data_ex( ctrl, 1,   'iiiiii ', nier, iterpremax, nrest, ncolor_in, recyclepre, nBFGS )/= 0) return
    if( fstr_ctrl_get_data_ex( ctrl, 2,   'rrr ', resid, singma_diag, sigma )/= 0) return

    if( precond == 20 .or. precond == 21) then
      if( fstr_ctrl_get_data_ex( ctrl, 3, 'rr ', thresh, filter)/= 0) return
    else if( precond == 5 ) then
      if( fstr_ctrl_get_data_ex( ctrl, 3, 'iiiiiiiiii ', &
           solver_opt(1), solver_opt(2), solver_opt(3), solver_opt(4), solver_opt(5), &
           solver_opt(6), solver_opt(7), solver_opt(8), solver_opt(9), solver_opt(10) )/= 0) return
    else if( method == 101 ) then
      if( fstr_ctrl_get_data_ex( ctrl, 3, 'i ', solver_opt(1) )/= 0) return
    end if

    iterlog = iter -1
    timelog = time -1
    steplog = step -1
    scaling = sclg -1
    dumpexit = dmpx -1
    usejad = usjd -1

    fstr_ctrl_get_SOLVER = 0

  end function fstr_ctrl_get_SOLVER


  !> Read in !STEP
  function fstr_ctrl_get_STEP( ctrl, amp, iproc )
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    integer(kind=kint) :: iproc
    integer(kind=kint) :: fstr_ctrl_get_STEP

    integer(kind=kint) :: ipt = 0
    integer(kind=kint) :: ip = 0

    fstr_ctrl_get_STEP = -1

    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',     '# ',  0, 'S', amp )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',   'STANDARD,NLGEOM ', 0, 'P',   ipt    )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'NLGEOM ',  '# ',           0,    'E',   ip     )/= 0) return

    if( ipt == 2 .or. ip == 1 ) iproc = 1

    fstr_ctrl_get_STEP = 0

  end function fstr_ctrl_get_STEP

  !> Read in !STEP and !ISTEP
  logical function fstr_ctrl_get_ISTEP( ctrl, hecMESH, steps, tpname, apname )
    use fstr_setup_util
    use m_step
    integer(kind=kint), intent(in)        :: ctrl      !< ctrl file
    type (hecmwST_local_mesh), intent(in) :: hecMESH   !< mesh information
    type(step_info), intent(out)          :: steps     !< step control info
    character(len=*), intent(out)         :: tpname    !< name of timepoints
    character(len=*), intent(out)         :: apname    !< name of auto increment parameter

    character(len=HECMW_NAME_LEN) :: data_fmt,ss, data_fmt1
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN) :: header_name
    integer(kind=kint) :: bcid
    integer(kind=kint) :: i, n, sn, ierr
    integer(kind=kint) :: bc_n, load_n, contact_n
    real(kind=kreal) :: fn, f1, f2, f3

    fstr_ctrl_get_ISTEP = .false.

    write(ss,*)  HECMW_NAME_LEN
    write( data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'I '
    write( data_fmt1, '(a,a,a)') 'S', trim(adjustl(ss)),'rrr '

    call init_stepInfo(steps)
    steps%solution = stepStatic
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',   'STATIC,VISCO ', 0, 'P', steps%solution )/= 0) return
    steps%inc_type = stepFixedInc
    if( fstr_ctrl_get_param_ex( ctrl, 'INC_TYPE ', 'FIXED,AUTO ', 0, 'P', steps%inc_type )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'SUBSTEPS ',  '# ',  0, 'I', steps%num_substep )/= 0) return
    steps%initdt = 1.d0/steps%num_substep
    if( fstr_ctrl_get_param_ex( ctrl, 'ITMAX ',  '# ',  0, 'I', steps%max_iter )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'MAXITER ',  '# ',  0, 'I', steps%max_iter )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'MAXCONTITER ',  '# ',  0, 'I', steps%max_contiter )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'CONVERG ',  '# ',  0, 'R', steps%converg )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'CONVERG_LAG ',  '# ',  0, 'R', steps%converg_lag )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'MAXRES ',  '# ',  0, 'R', steps%maxres )/= 0) return
    amp = ""
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
    if( len( trim(amp) )>0 ) then
      call amp_name_to_id( hecMESH, '!STEP', amp, steps%amp_id )
    endif
    tpname=""
    if( fstr_ctrl_get_param_ex( ctrl, 'TIMEPOINTS ',  '# ',  0, 'S', tpname )/= 0) return
    apname=""
    if( fstr_ctrl_get_param_ex( ctrl, 'AUTOINCPARAM ',  '# ',  0, 'S', apname )/= 0) return

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) then
      fstr_ctrl_get_ISTEP = .true.;  return
    endif

    f2 = steps%mindt
    f3 = steps%maxdt
    if( fstr_ctrl_get_data_ex( ctrl, 1, data_fmt1, ss, f1, f2, f3  )/= 0) return
    read( ss, * , iostat=ierr ) fn
    sn=1
    if( ierr==0 ) then
      steps%initdt = fn
      steps%elapsetime = f1
      if( steps%inc_type == stepAutoInc ) then
        steps%mindt = min(f2,steps%initdt)
        steps%maxdt = f3
      endif
      steps%num_substep = max(int((f1+0.999999999d0*fn)/fn),steps%num_substep)
      !if( mod(f1,fn)/=0 ) steps%num_substep =steps%num_substep+1
      sn = 2
    endif

    bc_n = 0
    load_n = 0
    contact_n = 0
    do i=sn,n
      if( fstr_ctrl_get_data_ex( ctrl, i, data_fmt, header_name, bcid  )/= 0) return
      if( trim(header_name) == 'BOUNDARY' ) then
        bc_n = bc_n + 1
      else if( trim(header_name) == 'LOAD' ) then
        load_n = load_n +1
      else if( trim(header_name) == 'CONTACT' ) then
        contact_n = contact_n+1
      else if( trim(header_name) == 'TEMPERATURE' ) then
        !   steps%Temperature = .true.
      endif
    end do

    if( bc_n>0 ) allocate( steps%Boundary(bc_n) )
    if( load_n>0 ) allocate( steps%Load(load_n) )
    if( contact_n>0 ) allocate( steps%Contact(contact_n) )

    bc_n = 0
    load_n = 0
    contact_n = 0
    do i=sn,n
      if( fstr_ctrl_get_data_ex( ctrl, i, data_fmt, header_name, bcid  )/= 0) return
      if( trim(header_name) == 'BOUNDARY' ) then
        bc_n = bc_n + 1
        steps%Boundary(bc_n) = bcid
      else if( trim(header_name) == 'LOAD' ) then
        load_n = load_n +1
        steps%Load(load_n) = bcid
      else if( trim(header_name) == 'CONTACT' ) then
        contact_n = contact_n+1
        steps%Contact(contact_n) = bcid
      endif
    end do

    fstr_ctrl_get_ISTEP = .true.
  end function fstr_ctrl_get_ISTEP

  !> Read in !SECTION
  integer function fstr_ctrl_get_SECTION( ctrl, hecMESH, sections )
    use fstr_setup_util
    integer(kind=kint), intent(in)           :: ctrl
    type (hecmwST_local_mesh), intent(inout) :: hecMESH   !< mesh information
    type (tSection), pointer, intent(inout)  :: sections(:)

    integer(kind=kint)            :: j, k, sect_id, ori_id, elemopt
    integer(kind=kint),save       :: cache = 1
    character(len=HECMW_NAME_LEN) :: sect_orien
    character(19) :: form341list = 'FI,SELECTIVE_ESNS '
    character(16) :: form361list = 'FI,BBAR,IC,FBAR '

    fstr_ctrl_get_SECTION = -1

    if( fstr_ctrl_get_param_ex( ctrl, 'SECNUM ',  '# ',  1, 'I', sect_id )/= 0) return
    if( sect_id > hecMESH%section%n_sect ) return

    elemopt = 0
    if( fstr_ctrl_get_param_ex( ctrl, 'FORM341 ',   form341list, 0, 'P', elemopt )/= 0) return
    if( elemopt > 0 ) sections(sect_id)%elemopt341 = elemopt

    elemopt = 0
    if( fstr_ctrl_get_param_ex( ctrl, 'FORM361 ',   form361list, 0, 'P', elemopt )/= 0) return
    if( elemopt > 0 ) sections(sect_id)%elemopt361 = elemopt

    ! sectional orientation ID
    hecMESH%section%sect_orien_ID(sect_id) = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'ORIENTATION ',  '# ',  0, 'S', sect_orien )/= 0) return

    if( associated(g_LocalCoordSys) ) then
      call fstr_strupr(sect_orien)
      k = size(g_LocalCoordSys)

      if(cache < k)then
        if( sect_orien == g_LocalCoordSys(cache)%sys_name ) then
          hecMESH%section%sect_orien_ID(sect_id) = cache
          cache = cache + 1
          fstr_ctrl_get_SECTION = 0
          return
        endif
      endif

      do j=1, k
        if( sect_orien == g_LocalCoordSys(j)%sys_name ) then
          hecMESH%section%sect_orien_ID(sect_id) = j
          cache = j + 1
          exit
        endif
      enddo
    endif

    fstr_ctrl_get_SECTION = 0

  end function fstr_ctrl_get_SECTION


  !> Read in !WRITE
  function fstr_ctrl_get_WRITE( ctrl, res, visual, femap )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: res
    integer(kind=kint) :: visual
    integer(kind=kint) :: femap
    integer(kind=kint) :: fstr_ctrl_get_WRITE

    fstr_ctrl_get_WRITE = -1

    ! JP-6
    if( fstr_ctrl_get_param_ex( ctrl, 'RESULT ',  '# ',    0,   'E',   res    )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'VISUAL ',  '# ',    0,   'E',   visual )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'FEMAP ',   '# ',    0,   'E',   femap  )/= 0) return

    fstr_ctrl_get_WRITE = 0

  end function fstr_ctrl_get_WRITE

  !> Read in !ECHO
  function fstr_ctrl_get_ECHO( ctrl, echo )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: echo
    integer(kind=kint) :: fstr_ctrl_get_ECHO

    echo = kON;

    fstr_ctrl_get_ECHO = 0

  end function fstr_ctrl_get_ECHO

  !> Read in !COUPLE
  function fstr_ctrl_get_COUPLE( ctrl, fg_type, fg_first, fg_window, surf_id, surf_id_len )
    integer(kind=kint) :: ctrl                           !< readed data
    integer(kind=kint) :: fg_type                        !< if type
    integer(kind=kint) :: fg_first                       !< if first
    integer(kind=kint) :: fg_window                      !< if window
    character(len=HECMW_NAME_LEN),target  :: surf_id(:)  !< surface id
    character(len=HECMW_NAME_LEN),pointer :: surf_id_p   !< surface id
    integer(kind=kint) :: surf_id_len
    integer(kind=kint) :: fstr_ctrl_get_COUPLE

    character(len=HECMW_NAME_LEN) :: data_fmt,ss
    write(ss,*)  surf_id_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),' '

    fstr_ctrl_get_COUPLE = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ', '1,2,3,4,5,6 ', 0, 'I', fg_type )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'ISTEP ', '# ', 0, 'I', fg_first )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'WINDOW ', '# ', 0, 'I', fg_window )/= 0) return

    surf_id_p => surf_id(1)
    fstr_ctrl_get_COUPLE = &
      fstr_ctrl_get_data_array_ex( ctrl, data_fmt, surf_id_p )

  end function fstr_ctrl_get_COUPLE

  !> Read in !MPC
  function fstr_ctrl_get_MPC( ctrl, penalty )
    integer(kind=kint), intent(in) :: ctrl      !< readed data
    real(kind=kreal), intent(out)  :: penalty   !< penalty
    integer(kind=kint) :: fstr_ctrl_get_MPC

    fstr_ctrl_get_MPC = fstr_ctrl_get_data_ex( ctrl, 1,   'r ', penalty )
    if( penalty <= 1.0 ) then
      if (myrank == 0) then
        write(IMSG,*) "Warging : !MPC : too small penalty: ", penalty
        write(*,*) "Warging : !MPC : too small penalty: ", penalty
      endif
    endif

  end function fstr_ctrl_get_MPC

  !> Read in !OUTPUT_RES & !OUTPUT_VIS
  logical function fstr_ctrl_get_outitem( ctrl, hecMESH, outinfo )
    use fstr_setup_util
    use m_out
    integer(kind=kint), intent(in)        :: ctrl      !< readed data
    type (hecmwST_local_mesh), intent(in) :: hecMESH   !< mesh information
    type( output_info ), intent(out)      :: outinfo   !< output information

    integer(kind=kint) :: rcode, ipos
    integer(kind=kint) :: n, i, j
    character(len=HECMW_NAME_LEN) :: data_fmt, ss
    character(len=HECMW_NAME_LEN), allocatable :: header_name(:), onoff(:), vtype(:)

    write( ss, * )  HECMW_NAME_LEN
    write( data_fmt, '(a,a,a,a,a)') 'S', trim(adjustl(ss)), 'S', trim(adjustl(ss)), ' '
    !  write( data_fmt, '(a,a,a,a,a,a,a)') 'S', trim(adjustl(ss)), 'S', trim(adjustl(ss)), 'S', trim(adjustl(ss)), ' '

    fstr_ctrl_get_outitem = .false.

    outinfo%grp_id_name = "ALL"
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GROUP ', '# ', 0, 'S', outinfo%grp_id_name )
    ipos = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'ACTION ', 'SUM ', 0, 'P', ipos )
    outinfo%actn = ipos

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return
    allocate( header_name(n), onoff(n), vtype(n) )
    header_name(:) = ""; vtype(:) = "";  onoff(:) = ""
    rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, header_name, onoff )
    !  rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, header_name, onoff, vtype )

    do i = 1, n
      do j = 1, outinfo%num_items
        if( trim(header_name(i)) == outinfo%keyWord(j) ) then
          outinfo%on(j) = .true.
          if( trim(onoff(i)) == 'OFF' ) outinfo%on(j) = .false.
          if( len( trim(vtype(i)) )>0 ) then
            if( fstr_str2index( vtype(i), ipos ) ) then
              outinfo%vtype(j) = ipos
            else if( trim(vtype(i)) == "SCALER" ) then
              outinfo%vtype(j) = -1
            else if( trim(vtype(i)) == "VECTOR" ) then
              outinfo%vtype(j) = -2
            else if( trim(vtype(i)) == "SYMTENSOR" ) then
              outinfo%vtype(j) = -3
            else if( trim(vtype(i)) == "TENSOR" ) then
              outinfo%vtype(j) = -4
            endif
          endif
        endif
      enddo
    enddo

    deallocate( header_name, onoff, vtype )
    fstr_ctrl_get_outitem = .true.

  end function fstr_ctrl_get_outitem

  !> Read in !CONTACT
  function fstr_ctrl_get_CONTACTALGO( ctrl, algo )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: algo
    integer(kind=kint) :: fstr_ctrl_get_CONTACTALGO

    integer(kind=kint) :: rcode
    character(len=80) :: s
    algo = kcaSLagrange
    s = 'SLAGRANGE,ALAGRANGE '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', s, 0, 'P', algo )
    fstr_ctrl_get_CONTACTALGO = rcode
  end function fstr_ctrl_get_CONTACTALGO

  !>  Read in contact definition
  logical function fstr_ctrl_get_CONTACT( ctrl, n, contact, np, tp, ntol, ttol, ctAlgo, cpname )
    use fstr_setup_util
    integer(kind=kint), intent(in)    :: ctrl          !< ctrl file
    integer(kind=kint), intent(in)    :: n             !< number of item defined in this section
    integer(kind=kint), intent(in)    :: ctAlgo        !< contact algorithm
    type(tContact), intent(out)       :: contact(n)    !< contact definition
    real(kind=kreal), intent(out)      :: np             !< penalty along contact nomral
    real(kind=kreal), intent(out)      :: tp             !< penalty along contact tangent
    real(kind=kreal), intent(out)      :: ntol           !< tolrence along contact nomral
    real(kind=kreal), intent(out)      :: ttol           !< tolrence along contact tangent
    character(len=*), intent(out)      :: cpname         !< name of contact parameter

    integer           :: rcode, ipt
    character(len=30) :: s1 = 'TIED,GLUED,SSLID,FSLID '
    character(len=HECMW_NAME_LEN) :: data_fmt,ss
    character(len=HECMW_NAME_LEN) :: cp_name(n)
    real(kind=kreal)  :: fcoeff(n),tPenalty(n)

    tPenalty = 1.0d6

    write(ss,*)  HECMW_NAME_LEN

    fstr_ctrl_get_CONTACT = .false.
    contact(1)%ctype = 1   ! pure slave-master contact; default value
    contact(1)%algtype = CONTACTSSLID ! small sliding contact; default value
    rcode = fstr_ctrl_get_param_ex( ctrl, 'INTERACTION ', s1, 0, 'P', contact(1)%algtype )
    if( contact(1)%algtype==CONTACTGLUED ) contact(1)%algtype=CONTACTFSLID  ! not complemented yet
    if( fstr_ctrl_get_param_ex( ctrl, 'GRPID ', '# ', 1, 'I', contact(1)%group )/=0) return
    do rcode=2,n
      contact(rcode)%ctype = contact(1)%ctype
      contact(rcode)%group = contact(1)%group
      contact(rcode)%algtype = contact(1)%algtype
    end do

    if( contact(1)%algtype==CONTACTSSLID .or. contact(1)%algtype==CONTACTFSLID ) then
      write( data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)),'Rr '
      if(  fstr_ctrl_get_data_array_ex( ctrl, data_fmt, cp_name, fcoeff, tPenalty ) /= 0 ) return
      do rcode=1,n
        call fstr_strupr(cp_name(rcode))
        contact(rcode)%pair_name = cp_name(rcode)
        contact(rcode)%fcoeff = fcoeff(rcode)
        contact(rcode)%tPenalty = tPenalty(rcode)
      enddo
    else if( contact(1)%algtype==CONTACTTIED ) then
      write( data_fmt, '(a,a)') 'S', trim(adjustl(ss))
      if(  fstr_ctrl_get_data_array_ex( ctrl, data_fmt, cp_name ) /= 0 ) return
      do rcode=1,n
        call fstr_strupr(cp_name(rcode))
        contact(rcode)%pair_name = cp_name(rcode)
        contact(rcode)%fcoeff = 0.d0
        contact(rcode)%tPenalty = 1.d+4
      enddo
    endif

    np = 0.d0;  tp=0.d0
    ntol = 0.d0;  ttol=0.d0
    if( fstr_ctrl_get_param_ex( ctrl, 'NPENALTY ',  '# ',  0, 'R', np ) /= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'TPENALTY ', '# ', 0, 'R', tp ) /= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'NTOL ',  '# ',  0, 'R', ntol ) /= 0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'TTOL ', '# ', 0, 'R', ttol ) /= 0 ) return
    cpname=""
    if( fstr_ctrl_get_param_ex( ctrl, 'CONTACTPARAM ',  '# ',  0, 'S', cpname )/= 0) return
    fstr_ctrl_get_CONTACT = .true.
  end function fstr_ctrl_get_CONTACT

  !> Read in !CONTACT_PARAM                                                             !
  function fstr_ctrl_get_CONTACTPARAM( ctrl, contactparam )
    implicit none
    integer(kind=kint)    :: ctrl
    type( tContactParam ) :: contactparam !< contact parameter
    integer(kind=kint) :: fstr_ctrl_get_CONTACTPARAM

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=128) :: msg
    real(kind=kreal) :: CLEARANCE, CLR_SAME_ELEM, CLR_DIFFLPOS, CLR_CAL_NORM
    real(kind=kreal) :: DISTCLR_INIT, DISTCLR_FREE, DISTCLR_NOCHECK, TENSILE_FORCE
    real(kind=kreal) :: BOX_EXP_RATE

    fstr_ctrl_get_CONTACTPARAM = -1

    !parameters
    contactparam%name = ''
    if( fstr_ctrl_get_param_ex( ctrl, 'NAME ', '# ', 1, 'S', contactparam%name ) /=0 ) return

    !read first line
    data_fmt = 'rrrr '
    rcode = fstr_ctrl_get_data_ex( ctrl, 1, data_fmt, &
      &  CLEARANCE, CLR_SAME_ELEM, CLR_DIFFLPOS, CLR_CAL_NORM )
    if( rcode /= 0 ) return
    contactparam%CLEARANCE = CLEARANCE
    contactparam%CLR_SAME_ELEM = CLR_SAME_ELEM
    contactparam%CLR_DIFFLPOS  = CLR_DIFFLPOS
    contactparam%CLR_CAL_NORM  = CLR_CAL_NORM

    !read second line
    data_fmt = 'rrrrr '
    rcode = fstr_ctrl_get_data_ex( ctrl, 2, data_fmt, &
      &  DISTCLR_INIT, DISTCLR_FREE, DISTCLR_NOCHECK, TENSILE_FORCE, BOX_EXP_RATE )
    if( rcode /= 0 ) return
    contactparam%DISTCLR_INIT = DISTCLR_INIT
    contactparam%DISTCLR_FREE = DISTCLR_FREE
    contactparam%DISTCLR_NOCHECK = DISTCLR_NOCHECK
    contactparam%TENSILE_FORCE = TENSILE_FORCE
    contactparam%BOX_EXP_RATE = BOX_EXP_RATE

    !input check
    rcode = 1
    if( CLEARANCE<0.d0 .OR. 1.d0<CLEARANCE ) THEN
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : CLEARANCE must be 0 < CLEARANCE < 1.'
    else if( CLR_SAME_ELEM<0.d0 .or. 1.d0<CLR_SAME_ELEM ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : CLR_SAME_ELEM must be 0 < CLR_SAME_ELEM < 1.'
    else if( CLR_DIFFLPOS<0.d0 .or. 1.d0<CLR_DIFFLPOS ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : CLR_DIFFLPOS must be 0 < CLR_DIFFLPOS < 1.'
    else if( CLR_CAL_NORM<0.d0 .or. 1.d0<CLR_CAL_NORM ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : CLR_CAL_NORM must be 0 < CLR_CAL_NORM < 1.'
    else if( DISTCLR_INIT<0.d0 .or. 1.d0<DISTCLR_INIT ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : DISTCLR_INIT must be 0 < DISTCLR_INIT < 1.'
    else if( DISTCLR_FREE<-1.d0 .or. 1.d0<DISTCLR_FREE ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : DISTCLR_FREE must be -1 < DISTCLR_FREE < 1.'
    else if( DISTCLR_NOCHECK<0.5d0 ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : DISTCLR_NOCHECK must be >= 0.5.'
    else if( TENSILE_FORCE>=0.d0 ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : TENSILE_FORCE must be < 0.'
    else if( BOX_EXP_RATE<=1.d0 .or. 2.0<BOX_EXP_RATE ) then
      write(msg,*) 'fstr control file error : !CONTACT_PARAM : BOX_EXP_RATE must be 1 < BOX_EXP_RATE <= 2.'
    else
      rcode =0
    end if
    if( rcode /= 0 ) then
      write(*,*) trim(msg)
      write(ILOG,*) trim(msg)
      return
    endif

    fstr_ctrl_get_CONTACTPARAM = 0
  end function fstr_ctrl_get_CONTACTPARAM

  !> Read in !ELEMOPT
  function fstr_ctrl_get_ELEMOPT( ctrl, elemopt361 )
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: elemopt361
    integer(kind=kint) :: fstr_ctrl_get_ELEMOPT

    character(72) :: o361list = 'IC,Bbar '

    integer(kind=kint) :: o361

    fstr_ctrl_get_ELEMOPT = -1

    o361 = elemopt361 + 1

    !* parameter in header line -----------------------------------------------------------------*!
    if( fstr_ctrl_get_param_ex( ctrl, '361 ', o361list, 0, 'P', o361 ) /= 0) return

    elemopt361 = o361 - 1

    fstr_ctrl_get_ELEMOPT = 0

  end function fstr_ctrl_get_ELEMOPT


  !> Read in !AUTOINC_PARAM                                                             !
  function fstr_get_AUTOINC( ctrl, aincparam )
    implicit none
    integer(kind=kint)    :: ctrl
    type( tParamAutoInc ) :: aincparam !< auto increment parameter
    integer(kind=kint) :: fstr_get_AUTOINC

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=128) :: msg
    integer(kind=kint) :: bound_s(10), bound_l(10)
    real(kind=kreal) :: Rs, Rl

    fstr_get_AUTOINC = -1

    bound_s(:) = 0
    bound_l(:) = 0

    !parameters
    aincparam%name = ''
    if( fstr_ctrl_get_param_ex( ctrl, 'NAME ', '# ', 1, 'S', aincparam%name ) /=0 ) return

    !read first line ( decrease criteria )
    data_fmt = 'riiii '
    rcode = fstr_ctrl_get_data_ex( ctrl, 1, data_fmt, Rs, &
      &  bound_s(1), bound_s(2), bound_s(3), aincparam%NRtimes_s )
    if( rcode /= 0 ) return
    aincparam%ainc_Rs = Rs
    aincparam%NRbound_s(knstMAXIT) = bound_s(1)
    aincparam%NRbound_s(knstSUMIT) = bound_s(2)
    aincparam%NRbound_s(knstCITER) = bound_s(3)

    !read second line ( increase criteria )
    data_fmt = 'riiii '
    rcode = fstr_ctrl_get_data_ex( ctrl, 2, data_fmt, Rl, &
      &  bound_l(1), bound_l(2), bound_l(3), aincparam%NRtimes_l )
    if( rcode /= 0 ) return
    aincparam%ainc_Rl = Rl
    aincparam%NRbound_l(knstMAXIT) = bound_l(1)
    aincparam%NRbound_l(knstSUMIT) = bound_l(2)
    aincparam%NRbound_l(knstCITER) = bound_l(3)

    !read third line ( cutback criteria )
    data_fmt = 'ri '
    rcode = fstr_ctrl_get_data_ex( ctrl, 3, data_fmt, &
      &  aincparam%ainc_Rc, aincparam%CBbound )
    if( rcode /= 0 ) return

    !input check
    rcode = 1
    if( Rs<0.d0 .or. Rs>1.d0 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : decrease ratio Rs must 0 < Rs < 1.'
    else if( any(bound_s<0) ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : decrease NR bound must >= 0.'
    else if( aincparam%NRtimes_s < 1 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : # of times to decrease must > 0.'
    else if( Rl<1.d0 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : increase ratio Rl must > 1.'
    else if( any(bound_l<0) ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : increase NR bound must >= 0.'
    else if( aincparam%NRtimes_l < 1 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : # of times to increase must > 0.'
    elseif( aincparam%ainc_Rc<0.d0 .or. aincparam%ainc_Rc>1.d0 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : cutback decrease ratio Rc must 0 < Rc < 1.'
    else if( aincparam%CBbound < 1 ) then
      write(msg,*) 'fstr control file error : !AUTOINC_PARAM : maximum # of cutback times must > 0.'
    else
      rcode =0
    end if
    if( rcode /= 0 ) then
      write(*,*) trim(msg)
      write(ILOG,*) trim(msg)
      return
    endif

    fstr_get_AUTOINC = 0
  end function fstr_get_AUTOINC

  !> Read in !TIME_POINTS
  function fstr_ctrl_get_TIMEPOINTS( ctrl, tp )
    integer(kind=kint) :: ctrl
    type(time_points)  :: tp
    integer(kind=kint) :: fstr_ctrl_get_TIMEPOINTS

    integer(kind=kint) :: i, n, rcode
    logical            :: generate
    real(kind=kreal)   :: stime, etime, interval

    fstr_ctrl_get_TIMEPOINTS = -1

    tp%name = ''
    if( fstr_ctrl_get_param_ex( ctrl, 'NAME ', '# ', 1, 'S', tp%name ) /=0 ) return
    tp%range_type = 1
    if( fstr_ctrl_get_param_ex( ctrl, 'TIME ', 'STEP,TOTAL ', 0, 'P', tp%range_type ) /= 0 ) return
    generate = .false.
    if( fstr_ctrl_get_param_ex( ctrl, 'GENERATE ',  '# ', 0, 'E', generate ) /= 0) return

    if( generate ) then
      stime = 0.d0;  etime = 0.d0;  interval = 1.d0
      if( fstr_ctrl_get_data_ex( ctrl, 1, 'rrr ', stime, etime, interval ) /= 0) return
      tp%n_points = int((etime-stime)/interval)+1
      allocate(tp%points(tp%n_points))
      do i=1,tp%n_points
        tp%points(i) = stime + dble(i-1)*interval
      end do
    else
      n = fstr_ctrl_get_data_line_n( ctrl )
      if( n == 0 ) return
      tp%n_points = n
      allocate(tp%points(tp%n_points))
      if( fstr_ctrl_get_data_array_ex( ctrl, 'r ', tp%points ) /= 0 ) return
      do i=1,tp%n_points-1
        if( tp%points(i) < tp%points(i+1) ) cycle
        write(*,*) 'Error in reading !TIME_POINT: time points must be given in ascending order.'
        return
      end do
    end if

    fstr_ctrl_get_TIMEPOINTS = 0
  end function fstr_ctrl_get_TIMEPOINTS

  !> Read in !AMPLITUDE
  function fstr_ctrl_get_AMPLITUDE( ctrl, nline, name, type_def, type_time, type_val, n, val, table )
    implicit none
    integer(kind=kint),            intent(in)  :: ctrl
    integer(kind=kint),            intent(in)  :: nline
    character(len=HECMW_NAME_LEN), intent(out) :: name
    integer(kind=kint),            intent(out) :: type_def
    integer(kind=kint),            intent(out) :: type_time
    integer(kind=kint),            intent(out) :: type_val
    integer(kind=kint),            intent(out) :: n
    real(kind=kreal),              pointer     :: val(:)
    real(kind=kreal),              pointer     :: table(:)
    integer(kind=kint) :: fstr_ctrl_get_AMPLITUDE

    integer(kind=kint) :: t_def, t_time, t_val
    integer(kind=kint) :: i, j
    real(kind=kreal) :: r(4), t(4)

    fstr_ctrl_get_AMPLITUDE = -1

    name = ''
    t_def = 1
    t_time = 1
    t_val = 1

    if( fstr_ctrl_get_param_ex( ctrl, 'NAME ', '# ', 1, 'S', name )/=0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'DEFINITION ', 'TABULAR ', 0, 'P', t_def )/=0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'TIME ', 'STEP ', 0, 'P', t_time )/=0 ) return
    if( fstr_ctrl_get_param_ex( ctrl, 'VALUE ', 'RELATIVE,ABSOLUTE ', 0, 'P', t_val )/=0 ) return

    if( t_def==1 ) then
      type_def = HECMW_AMP_TYPEDEF_TABULAR
    else
      write(*,*) 'Error in reading !AMPLITUDE: invalid value for parameter DEFINITION.'
    endif
    if( t_time==1 ) then
      type_time = HECMW_AMP_TYPETIME_STEP
    else
      write(*,*) 'Error in reading !AMPLITUDE: invalid value for parameter TIME.'
    endif
    if( t_val==1 ) then
      type_val = HECMW_AMP_TYPEVAL_RELATIVE
    elseif( t_val==2 ) then
      type_val = HECMW_AMP_TYPEVAL_ABSOLUTE
    else
      write(*,*) 'Error in reading !AMPLITUDE: invalid value for parameter VALUE.'
    endif

    n = 0
    do i = 1, nline
      r(:)=huge(0.0d0); t(:)=huge(0.0d0)
      if( fstr_ctrl_get_data_ex( ctrl, 1, 'RRrrrrrr ', r(1), t(1), r(2), t(2), r(3), t(3), r(4), t(4) ) /= 0) return
      n = n+1
      val(n) = r(1)
      table(n) = t(1)
      do j = 2, 4
        if (r(j) < huge(0.0d0) .and. t(j) < huge(0.0d0)) then
          n = n+1
          val(n) = r(j)
          table(n) = t(j)
        else
          exit
        endif
      enddo
    enddo
    fstr_ctrl_get_AMPLITUDE = 0

  end function fstr_ctrl_get_AMPLITUDE

end module fstr_ctrl_common
