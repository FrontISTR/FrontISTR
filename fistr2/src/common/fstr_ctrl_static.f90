!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : I/O and Utility                                   !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains control file data obtaining functions for static analysis
module fstr_ctrl_static
use m_fstr
use hecmw
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

!* ----------------------------------------------------------------------------------------------- *!
!> Read in !STATIC      
!* ----------------------------------------------------------------------------------------------- *!

function fstr_ctrl_get_STATIC( ctrl, &
        & dtime, etime, itime, eps, restart_nout, &
        & idx_elpl, &
        & iout_list, &
        & sig_y0, h_dash, &
        & nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1 )
        implicit none
        integer(kind=kint) :: ctrl
        real(kind=kreal)   :: dtime
        real(kind=kreal)   :: etime
        integer(kind=kint) :: itime
        real(kind=kreal)   :: eps
        integer(kind=kint) :: restart_nout
        integer(kind=kint) :: idx_elpl
        real(kind=kreal)   :: sig_y0, h_dash
        integer(kind=kint) :: nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1
        integer(kind=kint) :: iout_list(6)
        integer(kind=kint) :: fstr_ctrl_get_STATIC

        fstr_ctrl_get_STATIC = -1

        if( fstr_ctrl_get_data_ex( ctrl, 1, 'rriri ', dtime, etime, itime, eps, restart_nout ) /= 0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 2, 'i ', idx_elpl ) /= 0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 3, 'iiiiii ', &
        & iout_list(1), iout_list(2), iout_list(3), iout_list(4), iout_list(5), iout_list(6)) /= 0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 4, 'rr ', sig_y0, h_dash ) /= 0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 5, 'iiiii ', &
                & nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1 ) /= 0 ) return

        fstr_ctrl_get_STATIC = 0
end function fstr_ctrl_get_STATIC

!* ----------------------------------------------------------------------------------------------- *!
!> Read in !BOUNDARY            
!* ----------------------------------------------------------------------------------------------- *!

integer function fstr_ctrl_get_BOUNDARY ( ctrl, npos, boundary_grp )
        implicit none
        integer(kind=kint)        :: ctrl
        integer, intent(inout)    :: npos
        type( fstr_boundary_grp ) :: boundary_grp(:)

        character(len=HECMW_NAME_LEN) :: data_fmt,ss, amp
        integer :: i, j, n, gid, amp_id, pid
        character(len=HECMW_NAME_LEN),pointer :: node_id_name(:)
        integer(kind=kint),pointer :: dof_ids(:), dof_ide(:)
        real(kind=kreal),pointer   :: fvalue(:)
	
        write(ss,*)  HECMW_NAME_LEN
        write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'IIr '

        fstr_ctrl_get_BOUNDARY = -1

        amp = ' '
        if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
      !  call amp_name_to_id( P%MESH, '!CLOAD', amp, amp_id ) 
        gid = 1
        if( fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ', 0, 'I', gid  )/=0) return
        pid = 0
        if( fstr_ctrl_get_param_ex( ctrl, 'PARTID ',  '# ', 0, 'I', pid  )/=0) return 

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n<=0 ) return
        allocate( node_id_name(n))
        allocate( dof_ids(n), dof_ide(n) )
        allocate( fvalue(n))

        if( fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_name, dof_ids, dof_ide, fvalue )/= 0) return
        do i=1,n
          npos = npos+1
          boundary_grp(npos)%gid = gid
          boundary_grp(npos)%part_id = pid
          boundary_grp(npos)%grp_name = node_id_name(i)
          if( dof_ids(i)<1 .or. dof_ids(i)>6 ) stop "Error in BOUNDARY definition!"
          if( dof_ide(i)<1 .or. dof_ide(i)>6 ) stop "Error in BOUNDARY definition!"
          do j=dof_ids(i), dof_ide(i)
            boundary_grp(npos)%dof(j) = 1
          enddo
          boundary_grp(npos)%fval = fvalue(i)
        enddo
		
        deallocate( node_id_name )
        deallocate( dof_ids, dof_ide )
        deallocate( fvalue )
		
		fstr_ctrl_get_BOUNDARY = 0

end function fstr_ctrl_get_BOUNDARY


!* ----------------------------------------------------------------------------------------------- *!
!> Read in !CLOAD       
!* ----------------------------------------------------------------------------------------------- *!

integer function fstr_ctrl_get_CLOAD( ctrl, npos, cload_grp )
        implicit none
        integer(kind=kint)        :: ctrl
        integer, intent(inout)    :: npos
        type( fstr_boundary_grp ) :: cload_grp(:)

        character(len=HECMW_NAME_LEN) :: data_fmt,ss, amp
        integer :: i, n, gid, amp_id, pid
        character(len=HECMW_NAME_LEN),pointer :: node_id_name(:)
        integer(kind=kint),pointer :: dof_id(:)
        real(kind=kreal),pointer   :: fvalue(:)
	
        write(ss,*)  HECMW_NAME_LEN
        write( data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IR '

        fstr_ctrl_get_CLOAD = -1

        amp = ' '
        if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
      !  call amp_name_to_id( P%MESH, '!CLOAD', amp, amp_id ) 
        gid = 1
        if( fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ', 0, 'I', gid  )/=0) return
        pid = 0
        if( fstr_ctrl_get_param_ex( ctrl, 'PARTID ',  '# ', 0, 'I', pid  )/=0) return 
		
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n<=0 ) return
        allocate( node_id_name(n))
        allocate( dof_id(n))
        allocate( fvalue(n))
        
        if( fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_name, dof_id, fvalue )/= 0) return
        do i=1,n
          npos = npos+1
          cload_grp(npos)%gid = gid
          cload_grp(npos)%part_id = pid
          cload_grp(npos)%grp_name = node_id_name(i)
          if( dof_id(i)<1 .or. dof_id(i)>6 ) stop "Error in CLOAD definition!"
          cload_grp(npos)%dof(dof_id(i)) = 1
          cload_grp(npos)%fval = fvalue(i)
        enddo
		
        deallocate( node_id_name )
        deallocate( dof_id )
        deallocate( fvalue )
		
		fstr_ctrl_get_CLOAD = 0

end function fstr_ctrl_get_CLOAD

!* ----------------------------------------------------------------------------------------------- *!
!> Read in !DLOAD       
!* ----------------------------------------------------------------------------------------------- *!

integer function fstr_ctrl_get_DLOAD( ctrl, npos, dload_grp )
        implicit none
        integer(kind=kint)        :: ctrl
        integer, intent(inout)    :: npos
        type( fstr_dload_grp )    :: dload_grp(:)

        character(len=HECMW_NAME_LEN) :: data_fmt,s1,s2, amp
        integer :: i, n, gid, amp_id, pid, lid, rcode

        real(kind=kreal),pointer   :: params(:,:)
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        character(len=HECMW_NAME_LEN),pointer :: type_name_list(:)

	
        write(s1,*)  HECMW_NAME_LEN
        write(s2,*)  HECMW_NAME_LEN
        write( data_fmt, '(a,a,a,a,a)') 'S', trim(adjustl(s1)), 'S', trim(adjustl(s2)),'Rrrrrrr '

        fstr_ctrl_get_DLOAD = -1

        amp = ' '
        if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
      !  call amp_name_to_id( P%MESH, '!CLOAD', amp, amp_id ) 
        gid = 1
        if( fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ', 0, 'I', gid  )/=0) return
        pid = 0
        if( fstr_ctrl_get_param_ex( ctrl, 'PARTID ',  '# ', 0, 'I', pid  )/=0) return 
		
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n<=0 ) return
        allocate( grp_id_name(n))
        allocate( params(7,n))
        allocate( type_name_list(n) )
       
        rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_id_name, type_name_list, &
                params(1,:), params(2,:), params(3,:), params(4,:), params(5,:),params(6,:), &
                params(7,:) )
				
        if( rcode /= 0 ) then
            deallocate( grp_id_name )
            deallocate( params )
            deallocate( type_name_list )
            return
        end if
				
        do i=1, n
                npos = npos+1
                call pc_strupr( type_name_list(i) )
                lid = -1;
                if(      type_name_list(i)(1:2) == 'BX'  ) then; lid = 1
                else if( type_name_list(i)(1:2) == 'BY'  ) then; lid = 2
                else if( type_name_list(i)(1:2) == 'BZ'  ) then; lid = 3
                else if( type_name_list(i)(1:4) == 'GRAV') then; lid = 4
                else if( type_name_list(i)(1:4) == 'CENT') then; lid = 5
                else if( type_name_list(i)(1:2) == 'PP'  ) then; lid = 10
                else if( type_name_list(i)(1:2) == 'P0'  ) then; lid = 10
                else if( type_name_list(i)(1:2) == 'P1'  ) then; lid = 10
                else if( type_name_list(i)(1:2) == 'P2'  ) then; lid = 20
                else if( type_name_list(i)(1:2) == 'P3'  ) then; lid = 30
                else if( type_name_list(i)(1:2) == 'P4'  ) then; lid = 40
                else if( type_name_list(i)(1:2) == 'P5'  ) then; lid = 50
                else if( type_name_list(i)(1:2) == 'P6'  ) then; lid = 60
                else if( type_name_list(i)(1:1) == 'S'   ) then; lid = 100
                else
                    write(ILOG, *) 'Error : !DLOAD : Load  type ',type_name_list(i), ' is unknown'
                    deallocate( grp_id_name )
                    deallocate( params )
                    deallocate( type_name_list )
                    return
                end if
                dload_grp(npos)%gid = gid
                dload_grp(npos)%part_id = pid
                dload_grp(npos)%grp_name = grp_id_name(i)
                dload_grp(npos)%itype = lid
                dload_grp(npos)%fval = params(:,i)
        end do
		
        deallocate( grp_id_name )
        deallocate( params )
        deallocate( type_name_list )
		
		fstr_ctrl_get_DLOAD = 0

end function fstr_ctrl_get_DLOAD


!* ----------------------------------------------------------------------------------------------- *!
!> Read in !REFTEMP        
!* ----------------------------------------------------------------------------------------------- *!

function fstr_ctrl_get_REFTEMP( ctrl, value )
        implicit none
        integer(kind=kint) :: ctrl
        real(kind=kreal)   :: value
        integer(kind=kint) :: fstr_ctrl_get_REFTEMP,rcode

        rcode = fstr_ctrl_get_data_array_ex( ctrl, 'r ', value )
        fstr_ctrl_get_REFTEMP = rcode

end function fstr_ctrl_get_REFTEMP

!* ----------------------------------------------------------------------------------------------- *!
!> Read in !TEMPERATURE        
!* ----------------------------------------------------------------------------------------------- *!

integer function fstr_ctrl_get_TEMPERATURE( ctrl, npos, temp_grp, irres, tstep )
        implicit none
        integer(kind=kint)        :: ctrl
        integer, intent(inout)    :: npos
        type( fstr_ndscalar_grp ) :: temp_grp(:)
        integer                   :: irres, tstep

        character(len=HECMW_NAME_LEN) :: data_fmt,ss, amp
        integer :: i, n, gid, amp_id, pid
        character(len=HECMW_NAME_LEN),pointer :: node_id_name(:)
        real(kind=kreal),pointer   :: fvalue(:)
		
        if( fstr_ctrl_get_param_ex( ctrl, 'READRESULT ', '# ', 0, 'E', irres )/= 0) return
        if( fstr_ctrl_get_param_ex( ctrl, 'TSTEP ',      '# ', 0, 'I', tstep )/= 0) return
        if( irres == 1 ) then
          fstr_ctrl_get_TEMPERATURE = 0
          return
        endif
	
        write(ss,*)  HECMW_NAME_LEN
        write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'R '

        fstr_ctrl_get_TEMPERATURE = -1

        amp = ' '
        if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
      !  call amp_name_to_id( P%MESH, '!CLOAD', amp, amp_id ) 
        gid = 1
        if( fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ', 0, 'I', gid  )/=0) return
        pid = 0
        if( fstr_ctrl_get_param_ex( ctrl, 'PARTID ',  '# ', 0, 'I', pid  )/=0) return 
		
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n<=0 ) return
        allocate( node_id_name(n))
        allocate( fvalue(n))
        
        if( fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_name, fvalue )/=0 ) return
		
        do i=1,n
          npos = npos+1
          temp_grp(npos)%gid = gid
          temp_grp(npos)%part_id = pid
          temp_grp(npos)%grp_name = node_id_name(i)
          temp_grp(npos)%fval = fvalue(i)
        enddo
		
        deallocate( node_id_name )
        deallocate( fvalue )
		
		fstr_ctrl_get_TEMPERATURE = 0

end function fstr_ctrl_get_TEMPERATURE


!----------------------------------------------------------------------
!> Read in !ULOAD
integer function fstr_ctrl_get_USERLOAD( ctrl )
        use mULoad
        integer(kind=kint), intent(in)    :: ctrl
		
        character(len=256) :: fname

        fstr_ctrl_get_USERLOAD = -1
        if( fstr_ctrl_get_param_ex( ctrl, 'FILE  ', '# ',           0,   'S',   fname )/=0 ) return
        if( fname=="" ) STOP "You must define a file name before read in user-defined material"
        if( ureadload(fname)/=0 ) return
		
        fstr_ctrl_get_USERMATERIAL = 0
end function fstr_ctrl_get_USERLOAD

end module fstr_ctrl_static




