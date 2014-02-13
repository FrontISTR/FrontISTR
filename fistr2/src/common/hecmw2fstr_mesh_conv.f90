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



!======================================================================!
!                                                                      !
!> \brief HECMW to FSTR Mesh Data Converter.                              
!! Convering Conectivity of Element Type 232, 342 and 352             
!                                                                      !
!======================================================================!


module m_hecmw2fstr_mesh_conv
use hecmw
        external hecmw2fstr_connect_conv

contains

subroutine hecmw2fstr_mesh_conv( hecMESH )
        implicit none
        type (hecmwST_local_mesh) :: hecMESH

        call hecmw2fstr_connect_conv( hecMESH%n_elem,    &
                                hecMESH%elem_type,       &
                                hecMESH%elem_node_index, &
                                hecMESH%elem_node_item )

end subroutine hecmw2fstr_mesh_conv


subroutine fstr2hecmw_mesh_conv( hecMESH )
        implicit none
        type (hecmwST_local_mesh) :: hecMESH

        call fstr2hecmw_connect_conv( hecMESH%n_elem,    &
                                hecMESH%elem_type,       &
                                hecMESH%elem_node_index, &
                                hecMESH%elem_node_item )

end subroutine fstr2hecmw_mesh_conv

end module m_hecmw2fstr_mesh_conv




