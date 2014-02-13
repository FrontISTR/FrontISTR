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
!> \brief This module provides a subroutine for all boundary conditions
!! needed in heat anaylsis
module m_heat_mat_ass_boundary
   contains
!C***
!C*** MAT_ASS_BOUNDARY
!C***
   subroutine heat_mat_ass_boundary ( hecMESH,hecMAT,fstrHEAT,ATIME, BTIME, DTIME )

      use m_fstr
      use m_heat_mat_ass_bc_CFLUX
      use m_heat_mat_ass_bc_DFLUX
      use m_heat_mat_ass_bc_FIXT
      use m_heat_mat_ass_bc_FILM
      use m_heat_mat_ass_bc_RADIATE

      implicit none
      real(kind=kreal)   ATIME,BTIME,CTIME, DTIME
      type (fstr_heat         ) :: fstrHEAT
      type (hecmwST_matrix    ) :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH

      CTIME = ATIME + BTIME

!C
!C +---------+
!C | !CFLUX  |
!C +---------+
!C===
      call heat_mat_ass_bc_CFLUX ( hecMAT, fstrHEAT, CTIME )
!C
!C +---------+
!C | !DFLUX  |
!C +---------+
!C===
      call heat_mat_ass_bc_DFLUX ( hecMESH, hecMAT, fstrHEAT, CTIME, DTIME )
!C
!C +--------+
!C | !FILM  |
!C +--------+
!C===
      call heat_mat_ass_bc_FILM ( hecMESH, hecMAT, fstrHEAT, CTIME )
!C
!C +-----------+
!C | !RADIATE  |
!C +-----------+
!C===
      call heat_mat_ass_bc_RADIATE ( hecMESH, hecMAT, fstrHEAT, CTIME )
!C
!C +------------+
!C | !BOUNDARY  |
!C +------------+
!C===
      call heat_mat_ass_bc_FIXT ( hecMAT, fstrHEAT, CTIME )
      return

   end subroutine heat_mat_ass_boundary
end module m_heat_mat_ass_boundary
