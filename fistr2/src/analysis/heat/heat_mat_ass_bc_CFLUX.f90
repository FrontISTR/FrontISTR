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
!> \brief This module provides a subroutine for setting concerntrated heat 
!! flux boundary conditions
module m_heat_mat_ass_bc_CFLUX
   contains
!C
!C***
!C*** MAT_ASS_BC_CFLUX
!C***
!C
   subroutine heat_mat_ass_bc_CFLUX( hecMAT, fstrHEAT, CTIME )
                                  
      use m_fstr
      use m_heat_get_amplitude

      implicit none
      integer(kind=kint) ib,in,ia
      real(kind=kreal)    CTIME,QQ
      type (fstr_heat         ) :: fstrHEAT
      type (hecmwST_matrix    ) :: hecMAT

      do ib = 1, fstrHEAT%Q_NOD_tot
        in = fstrHEAT%Q_NOD_node(ib)
        ia = fstrHEAT%Q_NOD_ampl(ib)
        call heat_get_amplitude ( fstrHEAT,ia,CTIME,QQ )
        hecMAT%B(in) = hecMAT%B(in) + QQ*fstrHEAT%Q_NOD_val(ib)
      enddo

      return
   end subroutine heat_mat_ass_bc_CFLUX
end module m_heat_mat_ass_bc_CFLUX
