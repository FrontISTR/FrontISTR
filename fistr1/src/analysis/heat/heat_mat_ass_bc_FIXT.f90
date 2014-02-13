!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> \brief This module provides a subroutine for setting fixed temperature
!! boundary conditions
module m_heat_mat_ass_bc_FIXT
   contains
!C
!C***
!C*** MAT_ASS_BC_FIXT
!C***
!C
   subroutine heat_mat_ass_bc_FIXT( hecMAT, fstrHEAT, CTIME )

      use m_fstr
      use m_heat_get_amplitude

      implicit none
      integer(kind=kint) ib,ii,id
      real(kind=kreal)   CTIME,QQ
      type (fstr_heat         ) :: fstrHEAT
      type (hecmwST_matrix    ) :: hecMAT
      logical :: OutOfRange

      do ib= 1, fstrHEAT%T_FIX_tot
        ii= fstrHEAT%T_FIX_node(ib)

        ID = fstrHEAT%T_FIX_ampl(ib)
        call heat_get_amplitude( fstrHEAT, ID, CTIME, QQ, OutOfRange )

        if (OutOfRange) cycle

        call hecmw_mat_ass_bc(hecMAT, ii, 1, fstrHEAT%T_FIX_VAL(ib) * QQ)

      enddo

      return

   end subroutine heat_mat_ass_bc_FIXT
end module m_heat_mat_ass_bc_FIXT
