!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Steady-state harmonic response analysis           !
!                                                                      !
!            Written by Kuniaki Koike (ASTOM)                          !
!                                                                      !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This source file contains subroutine for reading control data for harmonic response analysis (this implementation intend that this file is included by fstr_setup.f90)


!-----------------------------------------------------------------------------!
!> Read in !FLOAD                                                             !
!-----------------------------------------------------------------------------!
  subroutine fstr_setup_FLOAD( ctrl, counter, P )
  !---- args
    integer(kind=kint)   :: ctrl
    integer(kind=kint)   :: counter
    type(fstr_param_pack) :: P  
  !---- vals
    integer(kind=kint)                  :: rcode
    character(HECMW_NAME_LEN)           :: amp
    integer(kind=kint)                  :: amp_id
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal), pointer           :: val_ptr(:)
    integer(kind=kint), pointer        :: id_ptr(:), type_ptr(:)
    integer(kind=kint)                  :: i, n, old_size, new_size
    integer(kind=kint)                  :: gid, loadcase    
  !---- body
    
    if( P%SOLID%file_type /= kbcfFSTR) return
    
    !read grpid
    gid = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',  0, 'I', gid )
    !read loadcase (real=1:default, img=2)
    loadcase = kFLOADCASE_RE
    rcode = fstr_ctrl_get_param_ex( ctrl, 'LOAD CASE ', '# ', 0, 'I', loadcase)
    !write(*,*) "loadcase=", loadcase
    !pause
    
    !read the num of dataline
    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return
    old_size = P%FREQ%FLOAD_ngrp_tot
    new_size = old_size + n
    
    !expand data array
    P%FREQ%FLOAD_ngrp_tot = new_size
    call fstr_expand_integer_array( P%FREQ%FLOAD_ngrp_GRPID, old_size, new_size )
    call fstr_expand_integer_array( P%FREQ%FLOAD_ngrp_ID,    old_size, new_size )
    call fstr_expand_integer_array( P%FREQ%FLOAD_ngrp_TYPE,  old_size, new_size )
    call fstr_expand_integer_array( P%FREQ%FLOAD_ngrp_DOF,   old_size, new_size )
    call fstr_expand_real_array   ( P%FREQ%FLOAD_ngrp_valre, old_size, new_size )
    call fstr_expand_real_array   ( P%FREQ%FLOAD_ngrp_valim, old_size, new_size )
   
    !fill bc data
    allocate( grp_id_name(n) )
    if(loadcase == kFLOADCASE_RE) then    
      val_ptr  => P%FREQ%FLOAD_ngrp_valre(old_size+1:)
    else if(loadcase == kFLOADCASE_IM) then
      val_ptr  => P%FREQ%FLOAD_ngrp_valim(old_size+1:)
    else
      !error
      write(*,*)    "Error this load set is not defined!"
      write(ilog,*) "Error this load set is not defined!"
      stop
    end if
    id_ptr   => P%FREQ%FLOAD_ngrp_DOF(old_size+1:)
    type_ptr => P%FREQ%FLOAD_ngrp_TYPE(old_size+1:)
    val_ptr = 0.0D0
    rcode = fstr_ctrl_get_FLOAD( ctrl, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr)
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    P%FREQ%FLOAD_ngrp_GRPID(old_size+1:new_size) = gid
    call nodesurf_grp_name_to_id_ex( P%MESH, '!FLOAD', n, grp_id_name, &
         P%FREQ%FLOAD_ngrp_ID(old_size+1:), P%FREQ%FLOAD_ngrp_TYPE(old_size+1:))
    
    deallocate( grp_id_name )
    return
    
    contains
    
    function fstr_ctrl_get_FLOAD(ctrl, node_id, node_id_len, dof_id, value)
!      include 'fstr_ctrl_util_f.inc'                         !Fortran->C interface for fstr ctrl API    
      integer(kind=kint)                    :: ctrl
      character(len=HECMW_NAME_LEN),target :: node_id(:)  !Node group name
      integer(kind=kint), pointer          :: dof_id(:)
      integer(kind=kint)                    :: node_id_len    
      real(kind=kreal), pointer             :: value(:)
      integer(kind=kint)                    :: fstr_ctrl_get_FLOAD !return value
      character(len=HECMW_NAME_LEN)        :: data_fmt, ss
      character(len=HECMW_NAME_LEN),pointer :: node_id_p
      
      write(ss,*) node_id_len
      write(data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IR '
      
      node_id_p => node_id(1)
      fstr_ctrl_get_FLOAD = fstr_ctrl_get_data_array_ex(ctrl, data_fmt, node_id_p, dof_id, value)
    end function
      
  end subroutine
  
!-----------------------------------------------------------------------------!
!> Read in !EIGENREAD                                                         !
!-----------------------------------------------------------------------------!
  subroutine fstr_setup_eigenread( ctrl, counter, P )
  !---- args
    integer(kind=kint)    :: ctrl
    integer(kind=kint)    :: counter
    type(fstr_param_pack) :: P
  !---- vals    
    integer(kind=kint)                :: filename_len
    character(len=HECMW_NAME_LEN) :: datafmt, ss
  !---- body
  
    filename_len = HECMW_FILENAME_LEN
    write(ss,*) filename_len
    write(datafmt, '(a,a,a)') 'S', trim(adjustl(ss)), ' '
    
    if( fstr_ctrl_get_data_ex( ctrl, 1, datafmt, P%FREQ%eigenlog_filename ) /= 0) return
    if( fstr_ctrl_get_data_ex( ctrl, 2, 'ii ', P%FREQ%start_mode, P%FREQ%end_mode ) /= 0) return
    
    return

  end subroutine
 
