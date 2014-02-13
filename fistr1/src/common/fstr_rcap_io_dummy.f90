!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : common                                            !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> DUMMY ROUTINES FOR FSTR_RCAP_IO

module m_fstr_rcap_io
use m_fstr

    public :: fstr_rcap_initialize ! call after fstr_setup
    public :: fstr_rcap_finalize   ! call before hecmw_finalize
    public :: fstr_rcap_send
    public :: fstr_rcap_get

contains

!------------------------------------------------------------------------------
subroutine fstr_rcap_initialize( hecMESH, fstrPARAM, fstrCPL )
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type( fstr_param  ) :: fstrPARAM
    type( fstr_couple ) :: fstrCPL
    character( len=16)  :: portfile
    integer(kind=kint)  :: myrank
    integer(kind=kint)  :: err

    if( fstrPARAM%fg_couple == 0 ) return

    if( hecmw_comm_get_rank() == 0 ) then
        write(*,*) "##Error : REVOCAP functions are not supported"
    end if
    call hecmw_abort( hecmw_comm_get_comm() )

end subroutine fstr_rcap_initialize

!------------------------------------------------------------------------------
subroutine fstr_rcap_finalize( fstrPARAM, fstrCPL )
    implicit none
    type( fstr_param  ) :: fstrPARAM
    type( fstr_couple ) :: fstrCPL

    if( fstrPARAM%fg_couple == 0 ) return

    if( hecmw_comm_get_rank() == 0 ) then
        write(*,*) "##Error : REVOCAP functions are not supported"
    end if
    call hecmw_abort( hecmw_comm_get_comm() )

end subroutine fstr_rcap_finalize
!------------------------------------------------------------------------------
subroutine fstr_rcap_send( fstrCPL )
    implicit none
    type( fstr_couple ) :: fstrCPL

    if( hecmw_comm_get_rank() == 0 ) then
        write(*,*) "##Error : REVOCAP functions are not supported"
    end if
    call hecmw_abort( hecmw_comm_get_comm() )

end subroutine fstr_rcap_send
!------------------------------------------------------------------------------
subroutine fstr_rcap_get( fstrCPL )
    implicit none
    type( fstr_couple ) :: fstrCPL

    if( hecmw_comm_get_rank() == 0 ) then
        write(*,*) "##Error : REVOCAP functions are not supported"
    end if
    call hecmw_abort( hecmw_comm_get_comm() )

end subroutine fstr_rcap_get
!------------------------------------------------------------------------------
subroutine fstr_get_convergence( revocap_flag )
    implicit none
    integer(kind=kint)  :: revocap_flag

    if( hecmw_comm_get_rank() == 0 ) then
        write(*,*) "##Error : REVOCAP functions are not supported"
    end if
    call hecmw_abort( hecmw_comm_get_comm() )

end subroutine fstr_get_convergence
!------------------------------------------------------------------------------

end module m_fstr_rcap_io

