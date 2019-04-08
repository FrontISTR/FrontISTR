!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!
!C***
!C*** hecmw_solve_init
!C***
!

module m_hecmw_solve_init

contains

  subroutine hecmw_solve_init (hecMAT)
    use hecmw_util
    implicit real*8 (A-H,O-Z)
    type (hecmwST_matrix)     :: hecMAT

    hecMAT%Iarray= 0
    hecMAT%Rarray= 0.d0

    hecMAT%Iarray(1)= 100
    hecMAT%Iarray(6)=  10

    hecMAT%Rarray(1)= 1.d-8
    hecMAT%Rarray(2)= 1.d0
    hecMAT%Rarray(3)= 0.d0
    hecMAT%Rarray(4)= 0.10d0
    hecMAT%Rarray(5)= 0.10d0

  end subroutine hecmw_solve_init

end module m_hecmw_solve_init


