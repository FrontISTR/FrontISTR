!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module defines parameters for Lanczos eigenvalue solver
module lczparm
use hecmw
implicit none
public
     
      integer(kind=kint),parameter :: LENG = 256
      integer(kind=kint),parameter :: lvecq_size = 1000

	    !> Allocatable array, used or Lanczos eigenvalue analysis
        type lczvec
                 real(kind=kreal), pointer, dimension(:) :: q
        end type lczvec


		!> Package of data used by Lanczos eigenvalue solver
        type lczparam
                integer   (kind=kint)  :: eqset         ! Flag (1:eigen analysis,  0:not eigen ana.)
                integer   (kind=kint)  :: nget          ! Solved eigen value number (default:5)
                real      (kind=kreal) :: lczsgm        ! 0.0
                integer   (kind=kint)  :: lczmax        ! Max. Lcz iterations (default:60)
                real      (kind=kreal) :: lcztol        ! Lcz tolerance (default:1.0e-8)
                real      (kind=kreal) :: lczrod,lczrot ! lczrod = 1.0, lczrot = 0.0
                real      (kind=kreal) :: iluetol
                real      (kind=kreal), pointer :: mass(:)
        end type lczparam

contains

        subroutine fstr_nullify_lczparam( E )
        implicit none
        type( lczparam ) :: E
        nullify( E%mass )
        end subroutine fstr_nullify_lczparam

end module lczparm

