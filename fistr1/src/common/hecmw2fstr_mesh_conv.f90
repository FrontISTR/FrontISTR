!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief HECMW to FSTR Mesh Data Converter.
!! Converting Connectivity of Element Type 232, 342 and 352

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




