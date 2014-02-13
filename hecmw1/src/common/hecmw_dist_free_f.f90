!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuaki Sakane (RIST)                          !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_dist_free_f
    use hecmw_util
    implicit none

    public :: hecmw_dist_free
    
    contains

    subroutine hecmw_dist_free(mesh)
        type(hecmwST_local_mesh) :: mesh

        call free_etc(mesh)
        call free_node(mesh)
        call free_elem(mesh)
        call free_comm(mesh)
        call free_adapt(mesh)
        call free_sect(mesh%section)
        call free_mat(mesh%material)
        call free_mpc(mesh%mpc)
        call free_amp(mesh%amp)
        call free_ngrp(mesh%node_group)
        call free_egrp(mesh%elem_group)
        call free_sgrp(mesh%surf_group)
        call free_cpair(mesh%contact_pair)
        call free_reforg(mesh%refine_origin)
    end subroutine hecmw_dist_free


    subroutine free_etc(mesh)
        type(hecmwST_local_mesh) :: mesh

        if(associated(mesh%files)) deallocate(mesh%files)
    end subroutine free_etc


    subroutine free_node(mesh)
        type(hecmwST_local_mesh) :: mesh

        if(associated(mesh%node)) deallocate(mesh%node)
        if(associated(mesh%node_ID)) deallocate(mesh%node_ID)
        if(associated(mesh%global_node_ID)) deallocate(mesh%global_node_ID)
        if(associated(mesh%node_val_index)) deallocate(mesh%node_val_index)
        if(associated(mesh%node_val_item)) deallocate(mesh%node_val_item)
        if(associated(mesh%node_dof_index)) deallocate(mesh%node_dof_index)
        if(associated(mesh%node_dof_item)) deallocate(mesh%node_dof_item)
        if(associated(mesh%node_init_val_index)) deallocate(mesh%node_init_val_index)
        if(associated(mesh%node_init_val_item)) deallocate(mesh%node_init_val_item)
        if(associated(mesh%node_internal_list)) deallocate(mesh%node_internal_list)
    end subroutine free_node


    subroutine free_elem(mesh)
        type(hecmwST_local_mesh) :: mesh

        if(associated(mesh%elem_type_index)) deallocate(mesh%elem_type_index)
        if(associated(mesh%elem_type_item)) deallocate(mesh%elem_type_item)
        if(associated(mesh%elem_type)) deallocate(mesh%elem_type)
        if(associated(mesh%section_ID)) deallocate(mesh%section_ID)
        if(associated(mesh%elem_mat_ID_index)) deallocate(mesh%elem_mat_ID_index)
        if(associated(mesh%elem_mat_ID_item)) deallocate(mesh%elem_mat_ID_item)
        if(associated(mesh%elem_node_index)) deallocate(mesh%elem_node_index)
        if(associated(mesh%elem_node_item)) deallocate(mesh%elem_node_item)
        if(associated(mesh%elem_ID)) deallocate(mesh%elem_ID)
        if(associated(mesh%global_elem_ID)) deallocate(mesh%global_elem_ID)
        if(associated(mesh%elem_internal_list)) deallocate(mesh%elem_internal_list)
        if(associated(mesh%elem_mat_int_index)) deallocate(mesh%elem_mat_int_index)
        if(associated(mesh%elem_mat_int_val)) deallocate(mesh%elem_mat_int_val)
        if(associated(mesh%elem_val_index)) deallocate(mesh%elem_val_index)
        if(associated(mesh%elem_val_item)) deallocate(mesh%elem_val_item)
    end subroutine free_elem


    subroutine free_comm(mesh)
        type(hecmwST_local_mesh) :: mesh

        if(associated(mesh%neighbor_pe)) deallocate(mesh%neighbor_pe)
        if(associated(mesh%import_index)) deallocate(mesh%import_index)
        if(associated(mesh%import_item)) deallocate(mesh%import_item)
        if(associated(mesh%export_index)) deallocate(mesh%export_index)
        if(associated(mesh%export_item)) deallocate(mesh%export_item)
        if(associated(mesh%shared_index)) deallocate(mesh%shared_index)
        if(associated(mesh%shared_item)) deallocate(mesh%shared_item)
    end subroutine free_comm


    subroutine free_adapt(mesh)
        type(hecmwST_local_mesh) :: mesh

        if(associated(mesh%when_i_was_refined_node)) deallocate(mesh%when_i_was_refined_node)
        if(associated(mesh%when_i_was_refined_elem)) deallocate(mesh%when_i_was_refined_elem)
        if(associated(mesh%adapt_parent_type)) deallocate(mesh%adapt_parent_type)
        if(associated(mesh%adapt_type)) deallocate(mesh%adapt_type)
        if(associated(mesh%adapt_level)) deallocate(mesh%adapt_level)
        if(associated(mesh%adapt_parent)) deallocate(mesh%adapt_parent)
        if(associated(mesh%adapt_children_index)) deallocate(mesh%adapt_children_index)
        if(associated(mesh%adapt_children_item)) deallocate(mesh%adapt_children_item)
    end subroutine free_adapt


    subroutine free_sect(sect)
        type(hecmwST_section) :: sect

        if(associated(sect%sect_type)) deallocate(sect%sect_type)
        if(associated(sect%sect_opt)) deallocate(sect%sect_opt)
        if(associated(sect%sect_mat_ID_index)) deallocate(sect%sect_mat_ID_index)
        if(associated(sect%sect_mat_ID_item)) deallocate(sect%sect_mat_ID_item)
        if(associated(sect%sect_I_index)) deallocate(sect%sect_I_index)
        if(associated(sect%sect_I_item)) deallocate(sect%sect_I_item)
        if(associated(sect%sect_R_index)) deallocate(sect%sect_R_index)
        if(associated(sect%sect_R_item)) deallocate(sect%sect_R_item)
        if(associated(sect%sect_orien_ID)) deallocate(sect%sect_orien_ID)
    end subroutine free_sect


    subroutine free_mat(mat)
        type(hecmwST_material) :: mat

        if(associated(mat%mat_name)) deallocate(mat%mat_name)
        if(associated(mat%mat_item_index)) deallocate(mat%mat_item_index)
        if(associated(mat%mat_subitem_index)) deallocate(mat%mat_subitem_index)
        if(associated(mat%mat_table_index)) deallocate(mat%mat_table_index)
        if(associated(mat%mat_val)) deallocate(mat%mat_val)
        if(associated(mat%mat_temp)) deallocate(mat%mat_temp)
    end subroutine free_mat


    subroutine free_mpc(mpc)
        type(hecmwST_mpc) :: mpc

        if(associated(mpc%mpc_index)) deallocate(mpc%mpc_index)
        if(associated(mpc%mpc_item)) deallocate(mpc%mpc_item)
        if(associated(mpc%mpc_dof)) deallocate(mpc%mpc_dof)
        if(associated(mpc%mpc_val)) deallocate(mpc%mpc_val)
    end subroutine free_mpc


    subroutine free_amp(amp)
        type(hecmwST_amplitude) :: amp

        if(associated(amp%amp_name)) deallocate(amp%amp_name)
        if(associated(amp%amp_type_definition)) deallocate(amp%amp_type_definition)
        if(associated(amp%amp_type_time)) deallocate(amp%amp_type_time)
        if(associated(amp%amp_type_value)) deallocate(amp%amp_type_value)
        if(associated(amp%amp_index)) deallocate(amp%amp_index)
        if(associated(amp%amp_val)) deallocate(amp%amp_val)
        if(associated(amp%amp_table)) deallocate(amp%amp_table)
    end subroutine free_amp


    subroutine free_ngrp(grp)
        type(hecmwST_node_grp) :: grp

        if(associated(grp%grp_name)) deallocate(grp%grp_name)
        if(associated(grp%grp_index)) deallocate(grp%grp_index)
        if(associated(grp%grp_item)) deallocate(grp%grp_item)
        if(associated(grp%bc_grp_ID)) deallocate(grp%bc_grp_ID)
        if(associated(grp%bc_grp_type)) deallocate(grp%bc_grp_type)
        if(associated(grp%bc_grp_index)) deallocate(grp%bc_grp_index)
        if(associated(grp%bc_grp_dof)) deallocate(grp%bc_grp_dof)
        if(associated(grp%bc_grp_val)) deallocate(grp%bc_grp_val)
    end subroutine free_ngrp


    subroutine free_egrp(grp)
        type(hecmwST_elem_grp) :: grp

        if(associated(grp%grp_name)) deallocate(grp%grp_name)
        if(associated(grp%grp_index)) deallocate(grp%grp_index)
        if(associated(grp%grp_item)) deallocate(grp%grp_item)
        if(associated(grp%bc_grp_ID)) deallocate(grp%bc_grp_ID)
        if(associated(grp%bc_grp_type)) deallocate(grp%bc_grp_type)
        if(associated(grp%bc_grp_index)) deallocate(grp%bc_grp_index)
        if(associated(grp%bc_grp_val)) deallocate(grp%bc_grp_val)
    end subroutine free_egrp


    subroutine free_sgrp(grp)
        type(hecmwST_surf_grp) :: grp

        if(associated(grp%grp_name)) deallocate(grp%grp_name)
        if(associated(grp%grp_index)) deallocate(grp%grp_index)
        if(associated(grp%grp_item)) deallocate(grp%grp_item)
        if(associated(grp%bc_grp_ID)) deallocate(grp%bc_grp_ID)
        if(associated(grp%bc_grp_type)) deallocate(grp%bc_grp_type)
        if(associated(grp%bc_grp_index)) deallocate(grp%bc_grp_index)
        if(associated(grp%bc_grp_val)) deallocate(grp%bc_grp_val)
    end subroutine free_sgrp

    subroutine free_cpair(cpair)
        type(hecmwST_contact_pair) :: cpair

        if(associated(cpair%name)) deallocate(cpair%name)
        if(associated(cpair%type)) deallocate(cpair%type)
        if(associated(cpair%slave_grp_id)) deallocate(cpair%slave_grp_id)
        if(associated(cpair%master_grp_id)) deallocate(cpair%master_grp_id)
    end subroutine free_cpair

    subroutine free_reforg(reforg)
        type(hecmwST_refine_origin) :: reforg

        if(associated(reforg%index)) deallocate(reforg%index)
        if(associated(reforg%item_index)) deallocate(reforg%item_index)
        if(associated(reforg%item_item)) deallocate(reforg%item_item)
    end subroutine free_reforg

end module hecmw_dist_free_f

