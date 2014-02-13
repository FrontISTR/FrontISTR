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
!> \brief This module contains control file data obtaining functions for dynamic analysis

module fstr_ctrl_dynamic
use m_fstr
use hecmw
include 'fstr_ctrl_util_f.inc'
        private :: fstr_ctrl_get_nval
contains

!* ----------------------------------------------------------------------------------------------- *!
!* fstr_ctrl_get_nval                                                                              *!
!* ----------------------------------------------------------------------------------------------- *!

function fstr_ctrl_get_nval( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )
        implicit none
        integer(kind=kint) :: ctrl
        character(len=HECMW_NAME_LEN) :: amp
        character(len=HECMW_NAME_LEN),target :: node_id(:)
        character(len=HECMW_NAME_LEN),pointer :: node_id_p
        integer(kind=kint) :: node_id_len
        integer(kind=kint),pointer :: dof_ids (:)
        integer(kind=kint),pointer :: dof_ide (:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: fstr_ctrl_get_nval

        character(len=HECMW_NAME_LEN) :: data_fmt,ss
        write(ss,*)  node_id_len
        write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'IIr '

        fstr_ctrl_get_nval = -1
        if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
        node_id_p => node_id(1)
        fstr_ctrl_get_nval = &
                fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_p, dof_ids, dof_ide, value )

end function fstr_ctrl_get_nval

!> Read in !VELOCITY                                                                                      
function fstr_ctrl_get_VELOCITY( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )
        implicit none
        integer(kind=kint) :: ctrl
        character(len=HECMW_NAME_LEN) :: amp
        character(len=HECMW_NAME_LEN) :: node_id(:)
        integer(kind=kint) :: node_id_len
        integer(kind=kint),pointer :: dof_ids (:)
        integer(kind=kint),pointer :: dof_ide (:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: fstr_ctrl_get_VELOCITY

        fstr_ctrl_get_VELOCITY = &
          fstr_ctrl_get_nval( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )

end function fstr_ctrl_get_VELOCITY

!> Read in !ACCELERATION                                                               
function fstr_ctrl_get_ACCELERATION( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )
        implicit none
        integer(kind=kint) :: ctrl
        character(len=HECMW_NAME_LEN) :: amp
        character(len=HECMW_NAME_LEN) :: node_id(:)
        integer(kind=kint) :: node_id_len
        integer(kind=kint),pointer :: dof_ids (:)
        integer(kind=kint),pointer :: dof_ide (:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: fstr_ctrl_get_ACCELERATION

        fstr_ctrl_get_ACCELERATION = &
          fstr_ctrl_get_nval( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )

end function fstr_ctrl_get_ACCELERATION

!> Read in !DYNAMIC                                            
function fstr_ctrl_get_DYNAMIC( ctrl, &
                idx_eqa, idx_resp, n_step, t_start, t_end, t_delta, restart_nout, &
                ganma, beta, idx_mas, idx_dmp, ray_m, ray_k, &
                nout, node_monit_1, nout_monit, iout_list )
        implicit none
        integer(kind=kint) :: ctrl

        ! ANALYSIS TYPE CONTROL
        integer(kind=kint) :: idx_eqa
        integer(kind=kint) :: idx_resp

        ! TIME CONTROL
        integer(kind=kint) :: n_step
        real(kind=kreal)   :: t_start
        real(kind=kreal)   :: t_end
        real(kind=kreal)   :: t_delta

        ! Newmark-beta parameter
        real(kind=kreal)   :: ganma
        real(kind=kreal)   :: beta

        ! mass matrix control
        integer(kind=kint) :: idx_mas

        ! damping control
        integer(kind=kint) :: idx_dmp
        real(kind=kreal)   :: ray_m
        real(kind=kreal)   :: ray_k

        ! OUTPUT CONTROL
        integer(kind=kint) :: restart_nout
        integer(kind=kint) :: nout
        integer(kind=kint) :: node_monit_1
        integer(kind=kint) :: nout_monit
        integer(kind=kint) :: iout_list(6)

        integer(kind=kint) :: fstr_ctrl_get_DYNAMIC

        fstr_ctrl_get_DYNAMIC = -1
        if( fstr_ctrl_get_data_ex( ctrl, 1, 'ii ',   idx_eqa, idx_resp )/=0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 2, 'rriri ', t_start, t_end, n_step, t_delta, restart_nout )/=0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 3, 'rr ',   ganma, beta )/=0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 4, 'iirr ', idx_mas, idx_dmp, ray_m, ray_k )/=0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 5, 'iii ',   nout, node_monit_1, nout_monit )/=0 ) return
        if( fstr_ctrl_get_data_ex( ctrl, 6, 'iiiiii ', &
                iout_list(1), iout_list(2), iout_list(3), iout_list(4), iout_list(5), iout_list(6) )/=0 ) return
        fstr_ctrl_get_DYNAMIC = 0

end function fstr_ctrl_get_DYNAMIC

end module fstr_ctrl_dynamic
