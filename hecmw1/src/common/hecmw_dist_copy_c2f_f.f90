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

module hecmw_dist_copy_c2f_f
    use hecmw_util
    implicit none

    private
    character(len=100) :: sname,vname

    public :: hecmw_dist_copy_c2f

    contains

    subroutine hecmw_dist_copy_c2f(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        call get_flags(mesh, ierr)
        if(ierr /= 0) return

        call get_etc(mesh, ierr)
        if(ierr /= 0) return

        call get_node(mesh, ierr)
        if(ierr /= 0) return

        call get_elem(mesh, ierr)
        if(ierr /= 0) return

        call get_comm(mesh, ierr)
        if(ierr /= 0) return

        call get_adapt(mesh, ierr)
        if(ierr /= 0) return

        call get_refine(mesh, ierr)
        if(ierr /= 0) return

        call get_sect(mesh%section, ierr)
        if(ierr /= 0) return

        call get_mat(mesh%material, ierr)
        if(ierr /= 0) return

        call get_mpc(mesh%mpc, ierr)
        if(ierr /= 0) return

        call get_amp(mesh%amp, ierr)
        if(ierr /= 0) return

        call get_ngrp(mesh%node_group, ierr)
        if(ierr /= 0) return

        call get_egrp(mesh%elem_group, ierr)
        if(ierr /= 0) return

        call get_sgrp(mesh%surf_group, ierr)
        if(ierr /= 0) return

        call get_contact_pair(mesh%contact_pair, ierr)
        if(ierr /= 0) return
    end subroutine hecmw_dist_copy_c2f


    subroutine get_flags(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'hecmw_flag_adapt'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_flag_adapt, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_initcon'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_flag_initcon, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_parttype'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_flag_parttype, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_partdepth'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_flag_partdepth, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_version'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_flag_version, ierr)
        if(ierr /= 0) return
    end subroutine get_flags


    subroutine get_etc(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'gridfile'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%gridfile, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_n_file'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%hecmw_n_file, ierr)
        if(ierr /= 0) return

        if(mesh%hecmw_n_file > 0) then
            vname = 'files'
            allocate(mesh%files(mesh%hecmw_n_file))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%files, ierr)
            if(ierr /= 0) return
        endif

        vname = 'header'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%header, ierr)
        if(ierr /= 0) return

        vname = 'zero_temp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%zero_temp, ierr)
        if(ierr /= 0) return
    end subroutine get_etc


    subroutine get_node(mesh, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'n_node'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_node, ierr)
        if(ierr /= 0) return

        vname = 'n_node_gross'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_node_gross, ierr)
        if(ierr /= 0) return

        vname = 'nn_internal'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%nn_internal, ierr)
        if(ierr /= 0) return

        if(mesh%nn_internal > 0) then
            vname = 'node_internal_list'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%node_internal_list(mesh%nn_internal))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_internal_list, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_node_gross > 0) then
            vname = 'node_ID'
            allocate(mesh%node_ID(mesh%n_node_gross*2))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_ID, ierr)
            if(ierr /= 0) return

            vname = 'global_node_ID'
            allocate(mesh%global_node_ID(mesh%n_node_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%global_node_ID, ierr)
            if(ierr /= 0) return

            vname = 'node'
            allocate(mesh%node(mesh%n_node_gross*3))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node, ierr)
            if(ierr /= 0) return
        endif

        vname = 'n_dof'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_dof, ierr)
        if(ierr /= 0) return

        vname = 'n_dof_grp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_dof_grp, ierr)
        if(ierr /= 0) return

        vname = 'n_dof_tot'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_dof_tot, ierr)
        if(ierr /= 0) return

        if(mesh%n_dof_grp > 0) then
            vname = 'node_dof_index'
            allocate(mesh%node_dof_index(0:mesh%n_dof_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_dof_index, ierr)
            if(ierr /= 0) return

            vname = 'node_dof_item'
            allocate(mesh%node_dof_item(mesh%n_dof_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_dof_item, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_node_gross > 0) then
            vname = 'node_val_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%node_val_index(0:mesh%n_node_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_val_item'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                if(mesh%node_val_index(mesh%n_node_gross) > 0) then
                    allocate(mesh%node_val_item(mesh%node_val_index(mesh%n_node_gross)))
                    call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif

            vname = 'node_init_val_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%node_init_val_index(0:mesh%n_node_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_init_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_init_val_item'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                if(mesh%node_init_val_index(mesh%n_node_gross) > 0) then
                    allocate(mesh%node_init_val_item(mesh%node_init_val_index(mesh%n_node_gross)))
                    call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_init_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif
    end subroutine get_node


    subroutine get_elem(mesh, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'n_elem'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_elem, ierr)
        if(ierr /= 0) return

        vname = 'n_elem_gross'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_elem_gross, ierr)
        if(ierr /= 0) return

        vname = 'ne_internal'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%ne_internal, ierr)
        if(ierr /= 0) return

        if(mesh%ne_internal > 0) then
            vname = 'elem_internal_list'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%elem_internal_list(mesh%ne_internal))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_internal_list, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_ID'
            allocate(mesh%elem_ID(mesh%n_elem_gross*2))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_ID, ierr)
            if(ierr /= 0) return

            vname = 'global_elem_ID'
            allocate(mesh%global_elem_ID(mesh%n_elem_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%global_elem_ID, ierr)
            if(ierr /= 0) return

            vname = 'elem_type'
            allocate(mesh%elem_type(mesh%n_elem_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_type, ierr)
            if(ierr /= 0) return
        endif

        vname = 'n_elem_type'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_elem_type, ierr)
        if(ierr /= 0) return

        if(mesh%n_elem_type > 0) then
            vname = 'elem_type_index'
            allocate(mesh%elem_type_index(0:mesh%n_elem_type))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_type_index, ierr)
            if(ierr /= 0) return

            vname = 'elem_type_item'
            allocate(mesh%elem_type_item(mesh%n_elem_type))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_type_item, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_node_index'
            allocate(mesh%elem_node_index(0:mesh%n_elem_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_node_index, ierr)
            if(ierr /= 0) return

            vname = 'elem_node_item'
            if(mesh%elem_node_index(mesh%n_elem_gross) > 0) then
                allocate(mesh%elem_node_item(mesh%elem_node_index(mesh%n_elem_gross)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_node_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'section_ID'
            allocate(mesh%section_ID(mesh%n_elem_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%section_ID, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_mat_ID_index'
            allocate(mesh%elem_mat_ID_index(0:mesh%n_elem_gross))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_mat_ID_index, ierr)
            if(ierr /= 0) return

            if(mesh%elem_mat_ID_index(mesh%n_elem_gross) > 0) then
                vname = 'elem_mat_ID_item'
                allocate(mesh%elem_mat_ID_item(mesh%elem_mat_ID_index(mesh%n_elem_gross)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_mat_ID_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_elem_mat_ID'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_elem_mat_ID, ierr)
        if(ierr /= 0) return

        if(mesh%n_elem_mat_ID > 0) then
            vname = 'elem_mat_int_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%elem_mat_int_index(0:mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_mat_int_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_mat_int_val'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                if(mesh%elem_mat_int_index(mesh%n_elem_gross) > 0) then
                    allocate(mesh%elem_mat_int_val(mesh%elem_mat_int_index(mesh%n_elem_mat_ID)))
                    call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_mat_int_val, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_val_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%elem_val_index(0:mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_val_item'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                if(mesh%elem_val_index(mesh%n_elem_gross) > 0) then
                    allocate(mesh%elem_val_item(mesh%elem_val_index(mesh%n_elem_gross)))
                    call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif
    end subroutine get_elem


    subroutine get_comm(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh


        sname = 'hecmwST_local_mesh'

        vname = 'zero'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%zero, ierr)
        if(ierr /= 0) return

        vname = 'HECMW_COMM'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%MPI_COMM, ierr)
        if(ierr /= 0) return

        vname = 'PETOT'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%PETOT, ierr)
        if(ierr /= 0) return

        vname = 'PEsmpTOT'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%PEsmpTOT, ierr)
        if(ierr /= 0) return

        vname = 'my_rank'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%my_rank, ierr)
        if(ierr /= 0) return

        vname = 'errnof'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%errnof, ierr)
        if(ierr /= 0) return

        vname = 'n_subdomain'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_subdomain, ierr)
        if(ierr /= 0) return

        vname = 'n_neighbor_pe'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_neighbor_pe, ierr)
        if(ierr /= 0) return

        if(mesh%n_neighbor_pe > 0) then
            vname = 'neighbor_pe'
            allocate(mesh%neighbor_pe(mesh%n_neighbor_pe))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%neighbor_pe, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_neighbor_pe > 0) then
            vname = 'import_index'
            allocate(mesh%import_index(0:mesh%n_neighbor_pe))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%import_index, ierr)
            if(ierr /= 0) return

            if(mesh%import_index(mesh%n_neighbor_pe) > 0) then
                vname = 'import_item'
                allocate(mesh%import_item(mesh%import_index(mesh%n_neighbor_pe)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%import_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'export_index'
            allocate(mesh%export_index(0:mesh%n_neighbor_pe))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%export_index, ierr)
            if(ierr /= 0) return

            if(mesh%export_index(mesh%n_neighbor_pe) > 0) then
                vname = 'export_item'
                allocate(mesh%export_item(mesh%export_index(mesh%n_neighbor_pe)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%export_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'shared_index'
            allocate(mesh%shared_index(0:mesh%n_neighbor_pe))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%shared_index, ierr)
            if(ierr /= 0) return

            if(mesh%shared_index(mesh%n_neighbor_pe) > 0) then
                vname = 'shared_item'
                allocate(mesh%shared_item(mesh%shared_index(mesh%n_neighbor_pe)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%shared_item, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_comm


    subroutine get_adapt(mesh, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_local_mesh) :: mesh


        sname = 'hecmwST_local_mesh'

        vname = 'coarse_grid_level'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%coarse_grid_level, ierr)
        if(ierr /= 0) return

        vname = 'n_adapt'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_adapt, ierr)
        if(ierr /= 0) return

        if(mesh%n_node_gross > 0) then
            vname = 'when_i_was_refined_node'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%when_i_was_refined_node(mesh%n_node_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%when_i_was_refined_node, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'when_i_was_refined_elem'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%when_i_was_refined_elem(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%when_i_was_refined_elem, ierr)
                if(ierr /= 0) return
            endif

            vname = 'adapt_parent_type'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%adapt_parent_type(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_parent_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'adapt_type'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%adapt_type(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_type, ierr)
                if(ierr /= 0) return
            endif


            vname = 'adapt_level'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%adapt_level(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_level, ierr)
                if(ierr /= 0) return
            endif

            vname = 'adapt_parent'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%adapt_parent(mesh%n_elem_gross*2))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_parent, ierr)
                if(ierr /= 0) return
            endif

            vname = 'adapt_children_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%adapt_children_index(0:mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_children_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'adapt_children_item'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                if(mesh%adapt_children_index(mesh%n_elem_gross) > 0) then
                    allocate(mesh%adapt_children_item(mesh%adapt_children_index(mesh%n_elem_gross)*2))
                    call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%adapt_children_item, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif
    end subroutine get_adapt


    subroutine get_refine(mesh, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_local_mesh), target :: mesh
        type(hecmwST_refine_origin), pointer :: reforg

        sname = 'hecmwST_local_mesh'

        vname = 'n_refine'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%n_refine, ierr)
        if(ierr /= 0) return

        if(mesh%n_node_gross > 0) then
            vname = 'node_old2new'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%node_old2new(mesh%n_node_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_old2new, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_new2old'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%node_new2old(mesh%n_node_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%node_new2old, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_old2new'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%elem_old2new(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_old2new, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_new2old'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(mesh%elem_new2old(mesh%n_elem_gross))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mesh%elem_new2old, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_refine > 0) then
            sname = 'hecmwST_refine_origin'
            reforg => mesh%refine_origin

            vname = 'index'
            allocate(reforg%index(0:mesh%n_refine))
            call hecmw_dist_copy_c2f_set_if(sname, vname, reforg%index, ierr)
            if(ierr /= 0) return

            vname = 'item_index'
            allocate(reforg%item_index(0:reforg%index(mesh%n_refine)))
            call hecmw_dist_copy_c2f_set_if(sname, vname, reforg%item_index, ierr)
            if(ierr /= 0) return

            vname = 'item_item'
            allocate(reforg%item_item(reforg%item_index(reforg%index(mesh%n_refine))))
            call hecmw_dist_copy_c2f_set_if(sname, vname, reforg%item_item, ierr)
            if(ierr /= 0) return
        endif
    end subroutine get_refine


    subroutine get_sect(sect, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_section) :: sect

        sname = 'hecmwST_section'

        vname = 'n_sect'
        call hecmw_dist_copy_c2f_set_if(sname, vname, sect%n_sect, ierr)
        if(ierr /= 0) return

        if(sect%n_sect > 0) then
            vname = 'sect_type'
            allocate(sect%sect_type(sect%n_sect))
            call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_type, ierr)
            if(ierr /= 0) return

            vname = 'sect_opt'
            allocate(sect%sect_opt(sect%n_sect))
            call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_opt, ierr)
            if(ierr /= 0) return

            vname = 'sect_mat_ID_index'
            allocate(sect%sect_mat_ID_index(0:sect%n_sect))
            call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_mat_ID_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_mat_ID_index(sect%n_sect) > 0) then
                vname = 'sect_mat_ID_item'
                allocate(sect%sect_mat_ID_item(sect%sect_mat_ID_index(sect%n_sect)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_mat_ID_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'sect_I_index'
            allocate(sect%sect_I_index(0:sect%n_sect))
            call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_I_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_I_index(sect%n_sect) > 0) then
                vname = 'sect_I_item'
                allocate(sect%sect_I_item(sect%sect_I_index(sect%n_sect)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_I_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'sect_R_index'
            allocate(sect%sect_R_index(0:sect%n_sect))
            call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_R_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_R_index(sect%n_sect) > 0) then
                vname = 'sect_R_item'
                allocate(sect%sect_R_item(sect%sect_R_index(sect%n_sect)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, sect%sect_R_item, ierr)
                if(ierr /= 0) return
            endif
			
            allocate(sect%sect_orien_ID(sect%n_sect))
            sect%sect_orien_ID(:) = -1
        endif
    end subroutine get_sect


    subroutine get_mat(mat, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_material) :: mat

        sname = 'hecmwST_material'

        vname = 'n_mat'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mat%n_mat, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_item'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mat%n_mat_item, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_subitem'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mat%n_mat_subitem, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_table'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mat%n_mat_table, ierr)
        if(ierr /= 0) return

        if(mat%n_mat > 0) then
            vname = 'mat_name'
            allocate(mat%mat_name(mat%n_mat))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_name, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat > 0) then
            vname = 'mat_item_index'
            allocate(mat%mat_item_index(0:mat%n_mat))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_item_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_item > 0) then
            vname = 'mat_subitem_index'
            allocate(mat%mat_subitem_index(0:mat%n_mat_item))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_subitem_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_subitem > 0) then
            vname = 'mat_table_index'
            allocate(mat%mat_table_index(0:mat%n_mat_subitem))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_table_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_table > 0) then
            vname = 'mat_val'
            allocate(mat%mat_val(mat%n_mat_table))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_val, ierr)
            if(ierr /= 0) return

            vname = 'mat_temp'
            allocate(mat%mat_temp(mat%n_mat_table))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mat%mat_temp, ierr)
            if(ierr /= 0) return
        endif
    end subroutine get_mat


    subroutine get_mpc(mpc, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_mpc) :: mpc

        sname = 'hecmwST_mpc'

        vname = 'n_mpc'
        call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%n_mpc, ierr)
        if(ierr /= 0) return

        if(mpc%n_mpc > 0) then
            vname = 'mpc_index'
            allocate(mpc%mpc_index(0:mpc%n_mpc))
            call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%mpc_index, ierr)
            if(ierr /= 0) return

            if(mpc%mpc_index(mpc%n_mpc) > 0) then
                vname = 'mpc_item'
                allocate(mpc%mpc_item(mpc%mpc_index(mpc%n_mpc)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%mpc_item, ierr)
                if(ierr /= 0) return

                vname = 'mpc_dof'
                allocate(mpc%mpc_dof(mpc%mpc_index(mpc%n_mpc)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%mpc_dof, ierr)
                if(ierr /= 0) return

                vname = 'mpc_val'
                allocate(mpc%mpc_val(mpc%mpc_index(mpc%n_mpc)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%mpc_val, ierr)
                if(ierr /= 0) return

                vname = 'mpc_const'
                allocate(mpc%mpc_const(mpc%n_mpc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, mpc%mpc_const, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_mpc


    subroutine get_amp(amp, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_amplitude) :: amp

        sname = 'hecmwST_amplitude'

        vname = 'n_amp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, amp%n_amp, ierr)
        if(ierr /= 0) return

        if(amp%n_amp > 0) then
            vname = 'amp_name'
            allocate(amp%amp_name(amp%n_amp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_name, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_definition'
            allocate(amp%amp_type_definition(amp%n_amp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_type_definition, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_time'
            allocate(amp%amp_type_time(amp%n_amp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_type_time, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_value'
            allocate(amp%amp_type_value(amp%n_amp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_type_value, ierr)
            if(ierr /= 0) return

            vname = 'amp_index'
            allocate(amp%amp_index(0:amp%n_amp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_index, ierr)
            if(ierr /= 0) return

            if(amp%amp_index(amp%n_amp) > 0) then
                vname = 'amp_val'
                allocate(amp%amp_val(amp%amp_index(amp%n_amp)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_val, ierr)
                if(ierr /= 0) return

                vname = 'amp_table'
                allocate(amp%amp_table(amp%amp_index(amp%n_amp)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, amp%amp_table, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_amp


    subroutine get_ngrp(grp, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_node_grp) :: grp

        sname = 'hecmwST_node_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            allocate(grp%grp_name(grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_name, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            allocate(grp%grp_index(0:grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                allocate(grp%grp_item(grp%grp_index(grp%n_grp)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_ID(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_type(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_index(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_dof'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_dof(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_dof, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_val(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_ngrp


    subroutine get_egrp(grp, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_elem_grp) :: grp

        sname = 'hecmwST_elem_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            allocate(grp%grp_name(grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_name, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            allocate(grp%grp_index(0:grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                allocate(grp%grp_item(grp%grp_index(grp%n_grp)))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_ID(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_type(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_index(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_val(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_egrp


    subroutine get_sgrp(grp, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_surf_grp) :: grp

        sname = 'hecmwST_surf_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            allocate(grp%grp_name(grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_name, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            allocate(grp%grp_index(0:grp%n_grp))
            call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                allocate(grp%grp_item(grp%grp_index(grp%n_grp)*2))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_c2f_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_ID(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_type(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_index(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            call hecmw_dist_copy_c2f_isalloc_if(sname, vname, is_allocated, ierr)
            if(is_allocated == 1) then
                allocate(grp%bc_grp_val(grp%n_bc))
                call hecmw_dist_copy_c2f_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine get_sgrp


    subroutine get_contact_pair(cpair, ierr)
        integer(kind=kint) :: ierr,is_allocated
        type(hecmwST_contact_pair) :: cpair

        sname = 'hecmwST_contact_pair'

        vname = 'n_pair'
        call hecmw_dist_copy_c2f_set_if(sname, vname, cpair%n_pair, ierr)
        if(ierr /= 0) return

        if(cpair%n_pair > 0) then
            vname = 'name'
            allocate(cpair%name(cpair%n_pair))
            call hecmw_dist_copy_c2f_set_if(sname, vname, cpair%name, ierr)
            if(ierr /= 0) return

            vname = 'type'
            allocate(cpair%type(cpair%n_pair))
            call hecmw_dist_copy_c2f_set_if(sname, vname, cpair%type, ierr)
            if(ierr /= 0) return

            vname = 'slave_grp_id'
            allocate(cpair%slave_grp_id(cpair%n_pair))
            call hecmw_dist_copy_c2f_set_if(sname, vname, cpair%slave_grp_id, ierr)
            if(ierr /= 0) return

            vname = 'master_grp_id'
            allocate(cpair%master_grp_id(cpair%n_pair))
            call hecmw_dist_copy_c2f_set_if(sname, vname, cpair%master_grp_id, ierr)
            if(ierr /= 0) return
        endif
    end subroutine get_contact_pair

end module hecmw_dist_copy_c2f_f
