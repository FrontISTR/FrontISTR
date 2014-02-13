!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuaki Sakane (RIST)                          !
!                       Noboru Imai (Univ. of Tokyo)                   !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

! memo)
!   Intel 9 compiler generates codes to wast stack memory
!   when an array of string is passed to external subroutines defined with C.
!   Then the pointer of the head of the array is passed
!   to avoid consumptions of stack memory.


module hecmw_dist_copy_f2c_f
    use hecmw_util
    implicit none

    private
    character(len=100) :: sname,vname

    public :: hecmw_dist_copy_f2c

    contains

    subroutine hecmw_dist_copy_f2c(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        call put_flags(mesh, ierr)
        if(ierr /= 0) return

        call put_etc(mesh, ierr)
        if(ierr /= 0) return

        call put_node(mesh, ierr)
        if(ierr /= 0) return

        call put_elem(mesh, ierr)
        if(ierr /= 0) return

        call put_comm(mesh, ierr)
        if(ierr /= 0) return

        call put_adapt(mesh, ierr)
        if(ierr /= 0) return

        call put_refine(mesh, ierr)
        if(ierr /= 0) return

        call put_sect(mesh%section, ierr)
        if(ierr /= 0) return

        call put_mat(mesh%material, ierr)
        if(ierr /= 0) return

        call put_mpc(mesh%mpc, ierr)
        if(ierr /= 0) return

        call put_amp(mesh%amp, ierr)
        if(ierr /= 0) return

        call put_ngrp(mesh%node_group, ierr)
        if(ierr /= 0) return

        call put_egrp(mesh%elem_group, ierr)
        if(ierr /= 0) return

        call put_sgrp(mesh%surf_group, ierr)
        if(ierr /= 0) return

        call put_contact_pair(mesh%contact_pair, ierr)
        if(ierr /= 0) return
    end subroutine hecmw_dist_copy_f2c


    subroutine put_flags(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'hecmw_flag_adapt'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_flag_adapt, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_initcon'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_flag_initcon, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_parttype'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_flag_parttype, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_partdepth'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_flag_partdepth, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_flag_version'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_flag_version, ierr)
        if(ierr /= 0) return
    end subroutine put_flags


    subroutine put_etc(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'gridfile'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%gridfile, ierr)
        if(ierr /= 0) return

        vname = 'hecmw_n_file'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%hecmw_n_file, ierr)
        if(ierr /= 0) return

        if(mesh%hecmw_n_file > 0) then
            vname = 'files'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%files, ierr)
            if(ierr /= 0) return
        endif

        vname = 'header'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%header, ierr)
        if(ierr /= 0) return

        vname = 'zero_temp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%zero_temp, ierr)
        if(ierr /= 0) return
    end subroutine put_etc


    subroutine put_node(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'
        vname = 'n_node'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_node, ierr)
        if(ierr /= 0) return

        vname = 'n_node_gross'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_node_gross, ierr)
        if(ierr /= 0) return

        vname = 'nn_internal'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%nn_internal, ierr)
        if(ierr /= 0) return

        if((mesh%hecmw_flag_parttype == 0 .OR. mesh%hecmw_flag_parttype == 2) .AND. mesh%nn_internal > 0) then
            vname = 'node_internal_list'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_internal_list, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_node_gross > 0) then
            vname = 'node_ID'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_ID, ierr)
            if(ierr /= 0) return

            vname = 'global_node_ID'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%global_node_ID, ierr)
            if(ierr /= 0) return

            vname = 'node'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node, ierr)
            if(ierr /= 0) return
        endif

        vname = 'n_dof'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_dof, ierr)
        if(ierr /= 0) return

        vname = 'n_dof_grp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_dof_grp, ierr)
        if(ierr /= 0) return

        vname = 'n_dof_tot'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_dof_tot, ierr)
        if(ierr /= 0) return

        if(mesh%n_dof_grp > 0) then
            vname = 'node_dof_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_dof_index, ierr)
            if(ierr /= 0) return

            vname = 'node_dof_item'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_dof_item, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_node_gross > 0) then
            vname = 'node_val_index'
            if(associated(mesh%node_val_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_val_item'
            if(associated(mesh%node_val_item)) then
                if(mesh%node_val_index(mesh%n_node_gross) > 0) then
                    call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif

            vname = 'node_init_val_index'
            if(associated(mesh%node_init_val_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_init_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_init_val_item'
            if(associated(mesh%node_init_val_item)) then
                if(mesh%node_init_val_index(mesh%n_node_gross) > 0) then
                    call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_init_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif
    end subroutine put_node


    subroutine put_elem(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'n_elem'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_elem, ierr)
        if(ierr /= 0) return

        vname = 'n_elem_gross'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_elem_gross, ierr)
        if(ierr /= 0) return

        vname = 'ne_internal'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%ne_internal, ierr)
        if(ierr /= 0) return

        if((mesh%hecmw_flag_parttype == 0 .OR. mesh%hecmw_flag_parttype == 1) .AND. mesh%ne_internal > 0) then
            vname = 'elem_internal_list'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_internal_list, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_ID'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_ID, ierr)
            if(ierr /= 0) return

            vname = 'global_elem_ID'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%global_elem_ID, ierr)
            if(ierr /= 0) return

            vname = 'elem_type'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_type, ierr)
            if(ierr /= 0) return
        endif

        vname = 'n_elem_type'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_elem_type, ierr)
        if(ierr /= 0) return

        if(mesh%n_elem_type > 0) then
            vname = 'elem_type_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_type_index, ierr)
            if(ierr /= 0) return

            vname = 'elem_type_item'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_type_item, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_node_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_node_index, ierr)
            if(ierr /= 0) return

            vname = 'elem_node_item'
            if(mesh%elem_node_index(mesh%n_elem_gross) > 0) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_node_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'section_ID'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%section_ID, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_mat_ID_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_mat_ID_index, ierr)
            if(ierr /= 0) return

            if(mesh%elem_mat_ID_index(mesh%n_elem_gross) > 0) then
                vname = 'elem_mat_ID_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_mat_ID_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_elem_mat_ID'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_elem_mat_ID, ierr)
        if(ierr /= 0) return

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_mat_int_index'
            if(associated(mesh%elem_mat_int_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_mat_int_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_mat_int_val'
            if(associated(mesh%elem_mat_int_val)) then
                if(mesh%elem_mat_int_index(mesh%n_elem_gross) > 0) then
                    call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_mat_int_val, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_val_index'
            if(associated(mesh%elem_val_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_val_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_val_item'
            if(associated(mesh%elem_val_item)) then
                if(mesh%elem_val_index(mesh%n_elem_gross) > 0) then
                    call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_val_item, ierr)
                    if(ierr /= 0) return
                endif
            endif
        endif
    end subroutine put_elem


    subroutine put_comm(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh


        sname = 'hecmwST_local_mesh'

        vname = 'zero'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%zero, ierr)
        if(ierr /= 0) return

        vname = 'HECMW_COMM'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%MPI_COMM, ierr)
        if(ierr /= 0) return

        vname = 'PETOT'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%PETOT, ierr)
        if(ierr /= 0) return

        vname = 'PEsmpTOT'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%PEsmpTOT, ierr)
        if(ierr /= 0) return

        vname = 'my_rank'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%my_rank, ierr)
        if(ierr /= 0) return

        vname = 'errnof'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%errnof, ierr)
        if(ierr /= 0) return

        vname = 'n_subdomain'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_subdomain, ierr)
        if(ierr /= 0) return

        vname = 'n_neighbor_pe'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_neighbor_pe, ierr)
        if(ierr /= 0) return

        if(mesh%n_neighbor_pe > 0) then
            vname = 'neighbor_pe'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%neighbor_pe, ierr)
            if(ierr /= 0) return

            vname = 'import_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%import_index, ierr)
            if(ierr /= 0) return

            if(mesh%import_index(mesh%n_neighbor_pe) > 0) then
                vname = 'import_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%import_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'export_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%export_index, ierr)
            if(ierr /= 0) return

            if(mesh%export_index(mesh%n_neighbor_pe) > 0) then
                vname = 'export_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%export_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'shared_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%shared_index, ierr)
            if(ierr /= 0) return

            if(mesh%shared_index(mesh%n_neighbor_pe) > 0) then
                vname = 'shared_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%shared_item, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_comm


    subroutine put_adapt(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        if(mesh%hecmw_flag_adapt == 0) return;

        sname = 'hecmwST_local_mesh'

        vname = 'coarse_grid_level'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%coarse_grid_level, ierr)
        if(ierr /= 0) return

        vname = 'n_adapt'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_adapt, ierr)
        if(ierr /= 0) return

        if(mesh%n_node_gross > 0) then
            vname = 'when_i_was_refined_node'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%when_i_was_refined_node, ierr)
            if(ierr /= 0) return
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'when_i_was_refined_elem'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%when_i_was_refined_elem, ierr)
            if(ierr /= 0) return

            vname = 'adapt_parent_type'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_parent_type, ierr)
            if(ierr /= 0) return

            vname = 'adapt_type'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_type, ierr)
            if(ierr /= 0) return

            vname = 'adapt_level'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_level, ierr)
            if(ierr /= 0) return

            vname = 'adapt_parent'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_parent, ierr)
            if(ierr /= 0) return

            vname = 'adapt_children_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_children_index, ierr)
            if(ierr /= 0) return

            vname = 'adapt_children_item'
            if(mesh%adapt_children_index(mesh%n_elem_gross) > 0) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%adapt_children_item, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_adapt


    subroutine put_refine(mesh, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_local_mesh) :: mesh

        sname = 'hecmwST_local_mesh'

        vname = 'n_refine'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%n_refine, ierr)
        if(ierr /= 0) return

        if(mesh%n_refine == 0) return;

        if(mesh%n_node_gross > 0) then
            vname = 'node_old2new'
            if(associated(mesh%node_old2new)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_old2new, ierr)
                if(ierr /= 0) return
            endif

            vname = 'node_new2old'
            if(associated(mesh%node_new2old)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%node_new2old, ierr)
                if(ierr /= 0) return
            endif
        endif

        if(mesh%n_elem_gross > 0) then
            vname = 'elem_old2new'
            if(associated(mesh%elem_old2new)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_old2new, ierr)
                if(ierr /= 0) return
            endif

            vname = 'elem_new2old'
            if(associated(mesh%elem_new2old)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, mesh%elem_new2old, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_refine


    subroutine put_sect(sect, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_section) :: sect

        sname = 'hecmwST_section'

        vname = 'n_sect'
        call hecmw_dist_copy_f2c_set_if(sname, vname, sect%n_sect, ierr)
        if(ierr /= 0) return

        if(sect%n_sect > 0) then
            vname = 'sect_type'
            call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_type, ierr)
            if(ierr /= 0) return

            vname = 'sect_opt'
            call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_opt, ierr)
            if(ierr /= 0) return

            vname = 'sect_mat_ID_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_mat_ID_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_mat_ID_index(sect%n_sect) > 0) then
                vname = 'sect_mat_ID_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_mat_ID_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'sect_I_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_I_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_I_index(sect%n_sect) > 0) then
                vname = 'sect_I_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_I_item, ierr)
                if(ierr /= 0) return
            endif

            vname = 'sect_R_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_R_index, ierr)
            if(ierr /= 0) return

            if(sect%sect_R_index(sect%n_sect) > 0) then
                vname = 'sect_R_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, sect%sect_R_item, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_sect


    subroutine put_mat(mat, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_material),target :: mat
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_material'

        vname = 'n_mat'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mat%n_mat, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_item'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mat%n_mat_item, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_subitem'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mat%n_mat_subitem, ierr)
        if(ierr /= 0) return

        vname = 'n_mat_table'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mat%n_mat_table, ierr)
        if(ierr /= 0) return

        if(mat%n_mat > 0) then
            vname = 'mat_name'
            name_p => mat%mat_name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat > 0) then
            vname = 'mat_item_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mat%mat_item_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_item > 0) then
            vname = 'mat_subitem_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mat%mat_subitem_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_subitem > 0) then
            vname = 'mat_table_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mat%mat_table_index, ierr)
            if(ierr /= 0) return
        endif

        if(mat%n_mat_table > 0) then
            vname = 'mat_val'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mat%mat_val, ierr)
            if(ierr /= 0) return

            vname = 'mat_temp'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mat%mat_temp, ierr)
            if(ierr /= 0) return
        endif
    end subroutine put_mat


    subroutine put_mpc(mpc, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_mpc) :: mpc

        sname = 'hecmwST_mpc'

        vname = 'n_mpc'
        call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%n_mpc, ierr)
        if(ierr /= 0) return

        if(mpc%n_mpc > 0) then
            vname = 'mpc_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%mpc_index, ierr)
            if(ierr /= 0) return

            if(mpc%mpc_index(mpc%n_mpc) > 0) then
                vname = 'mpc_item'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%mpc_item, ierr)
                if(ierr /= 0) return

                vname = 'mpc_dof'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%mpc_dof, ierr)
                if(ierr /= 0) return

                vname = 'mpc_val'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%mpc_val, ierr)
                if(ierr /= 0) return

                vname = 'mpc_const'
                call hecmw_dist_copy_f2c_set_if(sname, vname, mpc%mpc_const, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_mpc


    subroutine put_amp(amp, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_amplitude) :: amp
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_amplitude'

        vname = 'n_amp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, amp%n_amp, ierr)
        if(ierr /= 0) return

        if(amp%n_amp > 0) then
            vname = 'amp_name'
            name_p => amp%amp_name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_definition'
            call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_type_definition, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_time'
            call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_type_time, ierr)
            if(ierr /= 0) return

            vname = 'amp_type_value'
            call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_type_value, ierr)
            if(ierr /= 0) return

            vname = 'amp_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_index, ierr)
            if(ierr /= 0) return

            if(amp%amp_index(amp%n_amp) > 0) then
                vname = 'amp_val'
                call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_val, ierr)
                if(ierr /= 0) return

                vname = 'amp_table'
                call hecmw_dist_copy_f2c_set_if(sname, vname, amp%amp_table, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_amp


    subroutine put_ngrp(grp, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_node_grp) :: grp
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_node_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            name_p => grp%grp_name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            if(associated(grp%bc_grp_ID)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            if(associated(grp%bc_grp_type)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            if(associated(grp%bc_grp_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_dof'
            if(associated(grp%bc_grp_dof)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_dof, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            if(associated(grp%bc_grp_val)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_ngrp


    subroutine put_egrp(grp, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_elem_grp) :: grp
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_elem_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            name_p => grp%grp_name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            if(associated(grp%bc_grp_ID)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            if(associated(grp%bc_grp_type)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            if(associated(grp%bc_grp_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            if(associated(grp%bc_grp_val)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_egrp


    subroutine put_sgrp(grp, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_surf_grp) :: grp
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_surf_grp'

        vname = 'n_grp'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_grp, ierr)
        if(ierr /= 0) return

        if(grp%n_grp > 0) then
            vname = 'grp_name'
            name_p => grp%grp_name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return

            vname = 'grp_index'
            call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_index, ierr)
            if(ierr /= 0) return

            vname = 'grp_item'
            if(grp%grp_index(grp%n_grp) > 0) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%grp_item, ierr)
                if(ierr /= 0) return
            endif
        endif

        vname = 'n_bc'
        call hecmw_dist_copy_f2c_set_if(sname, vname, grp%n_bc, ierr)
        if(ierr /= 0) return

        if(grp%n_bc > 0) then
            vname = 'bc_grp_ID'
            if(associated(grp%bc_grp_ID)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_ID, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_type'
            if(associated(grp%bc_grp_type)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_type, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_index'
            if(associated(grp%bc_grp_index)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_index, ierr)
                if(ierr /= 0) return
            endif

            vname = 'bc_grp_val'
            if(associated(grp%bc_grp_val)) then
                call hecmw_dist_copy_f2c_set_if(sname, vname, grp%bc_grp_val, ierr)
                if(ierr /= 0) return
            endif
        endif
    end subroutine put_sgrp


    subroutine put_contact_pair(cpair, ierr)
        integer(kind=kint) :: ierr
        type(hecmwST_contact_pair) :: cpair
        character(len=HECMW_NAME_LEN),pointer :: name_p

        sname = 'hecmwST_contact_pair'

        vname = 'n_pair'
        call hecmw_dist_copy_f2c_set_if(sname, vname, cpair%n_pair, ierr)
        if(ierr /= 0) return

        if(cpair%n_pair > 0) then
            vname = 'name'
            name_p => cpair%name(1)
            call hecmw_dist_copy_f2c_set_if(sname, vname, name_p, ierr)
            if(ierr /= 0) return

            vname = 'type'
            call hecmw_dist_copy_f2c_set_if(sname, vname, cpair%type, ierr)
            if(ierr /= 0) return

            vname = 'slave_grp_id'
            call hecmw_dist_copy_f2c_set_if(sname, vname, cpair%slave_grp_id, ierr)
            if(ierr /= 0) return

            vname = 'master_grp_id'
            call hecmw_dist_copy_f2c_set_if(sname, vname, cpair%master_grp_id, ierr)
            if(ierr /= 0) return
        endif
    end subroutine put_contact_pair

end module hecmw_dist_copy_f2c_f

