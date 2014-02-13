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

module fstr_ctrl_eigen
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


!> Read in !EIGEN (struct)                                                                          
function fstr_ctrl_get_EIGEN( ctrl, nget, lcztol, lczmax)
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: nget
        real(kind=kreal) :: lcztol
        integer(kind=kint) :: lczmax
        integer(kind=kint) :: fstr_ctrl_get_EIGEN

        ! JP-16
        fstr_ctrl_get_EIGEN = fstr_ctrl_get_data_ex( ctrl, 1,  'Iri ',  nget, lcztol, lczmax )

end function fstr_ctrl_get_EIGEN


!* ----------------------------------------------------------------------------------------------- *!
end module fstr_ctrl_eigen

