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

module hecmw_dist_print_f
    use hecmw_util
    implicit none

    contains

    subroutine hecmw_dist_print(mesh)
        type(hecmwST_local_mesh) :: mesh 

        call hecmw_dist_print_etc(mesh)
        write(*,*) 
        call hecmw_dist_print_node(mesh)
        write(*,*) 
        call hecmw_dist_print_elem(mesh)
        write(*,*) 
        call hecmw_dist_print_pe(mesh)
        write(*,*)
        call hecmw_dist_print_adapt(mesh)
        write(*,*)
        call hecmw_dist_print_sect(mesh%section)
        write(*,*)
        call hecmw_dist_print_mat(mesh%material)
        write(*,*)
        call hecmw_dist_print_mpc(mesh%mpc)
        write(*,*)
        call hecmw_dist_print_amp(mesh%amp)
        write(*,*)
        call hecmw_dist_print_ngrp(mesh%node_group)
        write(*,*)
        call hecmw_dist_print_egrp(mesh%elem_group)
        write(*,*)
        call hecmw_dist_print_sgrp(mesh%surf_group)
    end subroutine hecmw_dist_print


    subroutine hecmw_dist_print_etc(mesh)
        type(hecmwST_local_mesh) :: mesh 

        write(*,*) 'FLAGS:'
        write(*,*) 'hecmw_flag_adapt:'
        write(*,*) mesh%hecmw_flag_adapt
        write(*,*) 'hecmw_flag_initcon:'
        write(*,*) mesh%hecmw_flag_initcon
        write(*,*) 'hecmw_flag_prattype:'
        write(*,*) mesh%hecmw_flag_parttype
        write(*,*) 'hecmw_flag_version:'
        write(*,*) mesh%hecmw_flag_version
        write(*,*) 'END of FLAGS'
        write(*,*) 
        write(*,*) 'HEADER:'
        write(*,*) trim(mesh%header)
        write(*,*) 'END of HEADER'
        write(*,*) 
        write(*,*) 'GRIDFILE:'
        write(*,*) trim(mesh%gridfile)
        write(*,*) 'END of GRIDFILE'
        write(*,*) 
        write(*,*) 'FILES:'
        write(*,*) 'mesh%hecmw_n_file:'
        write(*,*) mesh%hecmw_n_file
        if(associated(mesh%files)) then
            write(*,*) mesh%files
        endif
        write(*,*) 'END of FILES'
        write(*,*) 
        write(*,*) 'ZERO:'
        write(*,*) mesh%zero_temp
        write(*,*) 'END of ZERO'
    end subroutine hecmw_dist_print_etc


    subroutine hecmw_dist_print_node(mesh)
        type(hecmwST_local_mesh) :: mesh 

        write(*,*) 'NODE:'
        write(*,*) 'n_node:'
        write(*,*) mesh%n_node
        write(*,*) 'nn_internal:'
        write(*,*) mesh%nn_internal
        write(*,*) 'node_internal_list:'
        if(associated(mesh%node_internal_list)) then
            write(*,*) mesh%node_internal_list
        endif
        write(*,*) 'node_ID:'
        write(*,*) mesh%node_ID
        write(*,*) 'global_node_ID:'
        write(*,*) mesh%global_node_ID
        write(*,*) 'node:'
        write(*,*) mesh%node
        write(*,*) 'n_dof:'
        write(*,*) mesh%n_dof
        write(*,*) 'n_dof_grp:'
        write(*,*) mesh%n_dof_grp
        write(*,*) 'n_dof_tot:'
        write(*,*) mesh%n_dof_tot
        write(*,*) 'node_dof_index:'
        if(associated(mesh%node_dof_index)) then
            write(*,*) mesh%node_dof_index
        endif 
        write(*,*) 'node_dof_item:'
        if(associated(mesh%node_dof_item)) then
            write(*,*) mesh%node_dof_item
        endif 
        write(*,*) 'node_val_index:'
        if(associated(mesh%node_val_index)) then
            write(*,*) mesh%node_val_index
        endif
        write(*,*) 'node_val_item:'
        if(associated(mesh%node_val_item)) then
            write(*,*) mesh%node_val_item
        endif
        write(*,*) 'node_init_val_index:'
        if(associated(mesh%node_init_val_index)) then
            write(*,*) mesh%node_init_val_index
        endif
        write(*,*) 'node_init_val_item:'
        if(associated(mesh%node_init_val_item)) then
            write(*,*) mesh%node_init_val_item
        endif
        write(*,*) 'END of NODE'
    end subroutine hecmw_dist_print_node


    subroutine hecmw_dist_print_elem(mesh)
        type(hecmwST_local_mesh) :: mesh 

        write(*,*) 'ELEMENT:'
        write(*,*) 'n_elem:'
        write(*,*) mesh%n_elem
        write(*,*) 'ne_internal:'
        write(*,*) mesh%ne_internal
        write(*,*) 'elem_internal_list:'
        if(associated(mesh%elem_internal_list)) then
            write(*,*) mesh%elem_internal_list
        endif
        write(*,*) 'elem_ID:'
        write(*,*) mesh%elem_ID
        write(*,*) 'global_elem_ID:'
        write(*,*) mesh%global_elem_ID
        write(*,*) 'elem_type:'
        write(*,*) mesh%elem_type
        write(*,*) 'n_elem_type:'
        write(*,*) mesh%n_elem_type
        write(*,*) 'elem_type_index:'
        write(*,*) mesh%elem_type_index
        write(*,*) 'elem_type_item:'
        write(*,*) mesh%elem_type_item
        write(*,*) 'elem_node_index:'
        write(*,*) mesh%elem_node_index
        write(*,*) 'elem_node_item:'
        write(*,*) mesh%elem_node_item
        write(*,*) 'section_ID:'
        write(*,*) mesh%section_ID
        write(*,*) 'n_elem_mat_ID:'
        write(*,*) mesh%n_elem_mat_ID
        write(*,*) 'elem_mat_ID_index:'
        write(*,*) mesh%elem_mat_ID_index
        write(*,*) 'elem_mat_ID_item:'
        write(*,*) mesh%elem_mat_ID_item
        write(*,*) 'elem_mat_int_index:'
        if(associated(mesh%elem_mat_int_index)) then
            write(*,*) mesh%elem_mat_int_index
        endif
        write(*,*) 'elem_mat_int_val:'
        if(associated(mesh%elem_mat_int_val)) then
            write(*,*) mesh%elem_mat_int_val
        endif
        write(*,*) 'elem_val_index:'
        if(associated(mesh%elem_val_index)) then
            write(*,*) mesh%elem_val_index
        endif
        write(*,*) 'elem_val_item:'
        if(associated(mesh%elem_val_item)) then
            write(*,*) mesh%elem_val_item
        endif
        write(*,*) 'END of ELEMENT'
    end subroutine hecmw_dist_print_elem


    subroutine hecmw_dist_print_pe(mesh)
        type(hecmwST_local_mesh) :: mesh 

        write(*,*) 'PE:'
        write(*,*) 'zero:'
        write(*,*) mesh%zero
        write(*,*) 'MPI_COMM:'
        write(*,*) mesh%MPI_COMM
        write(*,*) 'PETOT:'
        write(*,*) mesh%PETOT
        write(*,*) 'PEsmpTOT:'
        write(*,*) mesh%PEsmpTOT
        write(*,*) 'my_rank:'
        write(*,*) mesh%my_rank
        write(*,*) 'errnof:'
        write(*,*) mesh%errnof
        write(*,*) 'n_subdomain:'
        write(*,*) mesh%n_subdomain
        write(*,*) 'n_neighbor_pe:'
        write(*,*) mesh%n_neighbor_pe
        write(*,*) 'neighbor_pe:'
        if(associated(mesh%neighbor_pe)) then
            write(*,*) mesh%neighbor_pe
        endif
        write(*,*) 'import_index:'
        if(associated(mesh%import_index)) then
            write(*,*) mesh%import_index
        endif
        write(*,*) 'import_item:'
        if(associated(mesh%import_item)) then
            write(*,*) mesh%import_item
        endif
        write(*,*) 'export_index:'
        if(associated(mesh%export_index)) then
            write(*,*) mesh%export_index
        endif
        write(*,*) 'export_item:'
        if(associated(mesh%export_item)) then
            write(*,*) mesh%export_item
        endif
        write(*,*) 'shared_index:'
        if(associated(mesh%shared_index)) then
            write(*,*) mesh%shared_index
        endif
        write(*,*) 'shared_item:'
        if(associated(mesh%shared_item)) then
            write(*,*) mesh%shared_item
        endif
        write(*,*) 'END of PE'
    end subroutine hecmw_dist_print_pe


    subroutine hecmw_dist_print_adapt(mesh)
        type(hecmwST_local_mesh) :: mesh 

        write(*,*) 'ADAPTATION:'
        write(*,*) 'coarse_grid_level:'
        write(*,*) mesh%coarse_grid_level
        write(*,*) 'when_i_was_refined_node:'
        if(associated(mesh%when_i_was_refined_node)) then
            write(*,*) mesh%when_i_was_refined_node
        endif
        write(*,*) 'when_i_was_refined_elem:'
        if(associated(mesh%when_i_was_refined_elem)) then
            write(*,*) mesh%when_i_was_refined_elem
        endif
        write(*,*) 'adapt_parent_type:'
        if(associated(mesh%adapt_parent_type)) then
            write(*,*) mesh%adapt_parent_type
        endif
        write(*,*) 'adapt_type:'
        if(associated(mesh%adapt_type)) then
            write(*,*) mesh%adapt_type
        endif
        write(*,*) 'adapt_level:'
        if(associated(mesh%adapt_level)) then
            write(*,*) mesh%adapt_level
        endif
        write(*,*) 'adapt_parent:'
        if(associated(mesh%adapt_parent)) then
            write(*,*) mesh%adapt_parent
        endif
        write(*,*) 'adapt_children_index:'
        if(associated(mesh%adapt_children_index)) then
            write(*,*) mesh%adapt_children_index
        endif
        write(*,*) 'adapt_children_item:'
        if(associated(mesh%adapt_children_item)) then
            write(*,*) mesh%adapt_children_item
        endif
        write(*,*) 'END of ADAPTATION'
    end subroutine hecmw_dist_print_adapt


    subroutine hecmw_dist_print_sect(sect)
        type(hecmwST_section) :: sect

        write(*,*) 'SECTION:'
        write(*,*) 'n_sect'
        write(*,*) sect%n_sect
        write(*,*) 'sect_type'
        if(associated(sect%sect_type)) then
            write(*,*) sect%sect_type
        endif
        write(*,*) 'sect_opt'
        if(associated(sect%sect_opt)) then
            write(*,*) sect%sect_opt
        endif
        write(*,*) 'sect_mat_ID_index'
        if(associated(sect%sect_mat_ID_index)) then
            write(*,*) sect%sect_mat_ID_index
        endif
        write(*,*) 'sect_mat_ID_item'
        if(associated(sect%sect_mat_ID_item)) then
            write(*,*) sect%sect_mat_ID_item
        endif
        write(*,*) 'sect_I_index'
        if(associated(sect%sect_I_index)) then
            write(*,*) sect%sect_I_index
        endif
        write(*,*) 'sect_I_item'
        if(associated(sect%sect_I_item)) then
            write(*,*) sect%sect_I_item
        endif
        write(*,*) 'sect_R_index'
        if(associated(sect%sect_R_index)) then
            write(*,*) sect%sect_R_index
        endif
        write(*,*) 'sect_R_item'
        if(associated(sect%sect_R_item)) then
            write(*,*) sect%sect_R_item
        endif
        write(*,*) 'END of SECTION'
    end subroutine hecmw_dist_print_sect


    subroutine hecmw_dist_print_mat(mat)
        type(hecmwST_material) :: mat

        write(*,*) 'MATERIAL:'
        write(*,*) 'n_mat'
        write(*,*) mat%n_mat
        write(*,*) 'n_mat_item'
        write(*,*) mat%n_mat_item
        write(*,*) 'n_mat_subitem'
        write(*,*) mat%n_mat_subitem
        write(*,*) 'n_mat_table'
        write(*,*) mat%n_mat_table
        write(*,*) 'mat_name'
        if(associated(mat%mat_name)) then
            write(*,*) mat%mat_name
        endif 
        write(*,*) 'mat_item_index'
        if(associated(mat%mat_item_index)) then
            write(*,*) mat%mat_item_index
        endif 
        write(*,*) 'mat_subitem_index'
        if(associated(mat%mat_subitem_index)) then
            write(*,*) mat%mat_subitem_index
        endif 
        write(*,*) 'mat_table_index'
        if(associated(mat%mat_table_index)) then
            write(*,*) mat%mat_table_index
        endif 
        write(*,*) 'mat_val'
        if(associated(mat%mat_val)) then
            write(*,*) mat%mat_val
        endif 
        write(*,*) 'mat_temp'
        if(associated(mat%mat_temp)) then
            write(*,*) mat%mat_temp
        endif 
        write(*,*) 'END of MATERIAL'
    end subroutine hecmw_dist_print_mat


    subroutine hecmw_dist_print_mpc(mpc)
        type(hecmwST_mpc) :: mpc

        write(*,*) 'MPC:'
        write(*,*) 'n_mpc'
        write(*,*) mpc%n_mpc
        write(*,*) 'mpc_index'
        if(associated(mpc%mpc_index)) then
            write(*,*) mpc%mpc_index
        endif
        write(*,*) 'mpc_item'
        if(associated(mpc%mpc_item)) then
            write(*,*) mpc%mpc_item
        endif
        write(*,*) 'mpc_dof'
        if(associated(mpc%mpc_dof)) then
            write(*,*) mpc%mpc_dof
        endif
        write(*,*) 'mpc_val'
        if(associated(mpc%mpc_val)) then
            write(*,*) mpc%mpc_val
        endif
        write(*,*) 'END of MPC'
    end subroutine hecmw_dist_print_mpc


    subroutine hecmw_dist_print_amp(amp)
        type(hecmwST_amplitude) :: amp

        write(*,*) 'AMPLITUDE:'
        write(*,*) 'n_amp'
        write(*,*) amp%n_amp
        write(*,*) 'amp_name'
        if(associated(amp%amp_name)) then
            write(*,*) amp%amp_name
        endif
        write(*,*) 'amp_type_definition'
        if(associated(amp%amp_type_definition)) then
            write(*,*) amp%amp_type_definition
        endif
        write(*,*) 'amp_type_time'
        if(associated(amp%amp_type_time)) then
            write(*,*) amp%amp_type_time
        endif
        write(*,*) 'amp_type_value'
        if(associated(amp%amp_type_value)) then
            write(*,*) amp%amp_type_value
        endif
        write(*,*) 'amp_index'
        if(associated(amp%amp_index)) then
            write(*,*) amp%amp_index
        endif
        write(*,*) 'amp_table'
        if(associated(amp%amp_table)) then
            write(*,*) amp%amp_table
        endif
        write(*,*) 'END of AMPLITUDE'
    end subroutine hecmw_dist_print_amp 


    subroutine hecmw_dist_print_ngrp(grp)
        type(hecmwST_node_grp) :: grp

        write(*,*) 'NGROUP:'
        write(*,*) 'n_grp'
        write(*,*) grp%n_grp
        write(*,*) 'grp_name'
        if(associated(grp%grp_name)) then
            write(*,*) grp%grp_name
        endif
        write(*,*) 'grp_index'
        if(associated(grp%grp_index)) then
            write(*,*) grp%grp_index
        endif
        write(*,*) 'grp_item'
        if(associated(grp%grp_item)) then
            write(*,*) grp%grp_item
        endif
        write(*,*) 'n_bc'
        write(*,*) grp%n_bc
        write(*,*) 'bc_grp_ID'
        if(associated(grp%bc_grp_ID)) then
            write(*,*) grp%bc_grp_ID
        endif
        write(*,*) 'bc_grp_type'
        if(associated(grp%bc_grp_type)) then
            write(*,*) grp%bc_grp_type
        endif
        write(*,*) 'bc_grp_index'
        if(associated(grp%bc_grp_index)) then
            write(*,*) grp%bc_grp_index
        endif
        write(*,*) 'bc_grp_dof'
        if(associated(grp%bc_grp_dof)) then
            write(*,*) grp%bc_grp_dof
        endif
        write(*,*) 'bc_grp_val'
        if(associated(grp%bc_grp_val)) then
            write(*,*) grp%bc_grp_val
        endif
        write(*,*) 'END of NGROUP'
    end subroutine hecmw_dist_print_ngrp


    subroutine hecmw_dist_print_egrp(grp)
        type(hecmwST_elem_grp) :: grp

        write(*,*) 'EGROUP:'
        write(*,*) 'n_grp'
        write(*,*) grp%n_grp
        write(*,*) 'grp_name'
        if(associated(grp%grp_name)) then
            write(*,*) grp%grp_name
        endif
        write(*,*) 'grp_index'
        if(associated(grp%grp_index)) then
            write(*,*) grp%grp_index
        endif
        write(*,*) 'grp_item'
        if(associated(grp%grp_item)) then
            write(*,*) grp%grp_item
        endif
        write(*,*) 'n_bc'
        write(*,*) grp%n_bc
        write(*,*) 'bc_grp_ID'
        if(associated(grp%bc_grp_ID)) then
            write(*,*) grp%bc_grp_ID
        endif
        write(*,*) 'bc_grp_type'
        if(associated(grp%bc_grp_type)) then
            write(*,*) grp%bc_grp_type
        endif
        write(*,*) 'bc_grp_index'
        if(associated(grp%bc_grp_index)) then
            write(*,*) grp%bc_grp_index
        endif
        write(*,*) 'bc_grp_val'
        if(associated(grp%bc_grp_val)) then
            write(*,*) grp%bc_grp_val
        endif
        write(*,*) 'END of EGROUP'
    end subroutine hecmw_dist_print_egrp


    subroutine hecmw_dist_print_sgrp(grp)
        type(hecmwST_surf_grp) :: grp

        write(*,*) 'SGROUP:'
        write(*,*) 'n_grp'
        write(*,*) grp%n_grp
        write(*,*) 'grp_name'
        if(associated(grp%grp_name)) then
            write(*,*) grp%grp_name
        endif
        write(*,*) 'grp_index'
        if(associated(grp%grp_index)) then
            write(*,*) grp%grp_index
        endif
        write(*,*) 'grp_item'
        if(associated(grp%grp_item)) then
            write(*,*) grp%grp_item
        endif
        write(*,*) 'n_bc'
        write(*,*) grp%n_bc
        write(*,*) 'bc_grp_ID'
        if(associated(grp%bc_grp_ID)) then
            write(*,*) grp%bc_grp_ID
        endif
        write(*,*) 'bc_grp_type'
        if(associated(grp%bc_grp_type)) then
            write(*,*) grp%bc_grp_type
        endif
        write(*,*) 'bc_grp_index'
        if(associated(grp%bc_grp_index)) then
            write(*,*) grp%bc_grp_index
        endif
        write(*,*) 'bc_grp_val'
        if(associated(grp%bc_grp_val)) then
            write(*,*) grp%bc_grp_val
        endif
        write(*,*) 'END of SGROUP'
    end subroutine hecmw_dist_print_sgrp

end module hecmw_dist_print_f
