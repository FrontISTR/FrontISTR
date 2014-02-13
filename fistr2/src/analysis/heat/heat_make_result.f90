!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Heat Analysis                                     !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provide a function to prepare output data of heat analysis
module m_heat_make_RESULT
   contains
!C***
!C***  MAKE RESULT
!C***
   subroutine heat_make_RESULT( hecMESH, fstrHEAT, fstrRESULT )

      use m_fstr
      use hecmw_result
      implicit none
      integer(kind=kint) i
      type (fstr_heat          ) :: fstrHEAT
      type (hecmwST_local_mesh ) :: hecMESH
      type (hecmwST_result_data) :: fstrRESULT
!C*** 

      do i= 1, hecMESH%n_node
        fstrRESULT%node_val_item(i)=fstrHEAT%TEMP(i)
      enddo

   end subroutine heat_make_RESULT
end module m_heat_make_RESULT
