!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

!C
!C***
!C*** hecmw_adapt_new_mesh
!C***
!C

      subroutine hecmw_adapt_new_mesh (hecMESH, hecMESHnew)

      use  hecmw_util
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_local_mesh) :: hecMESHnew

      integer(kind=kint), pointer :: WORKaI(:)
      real(kind=kreal),   pointer :: WORKaR(:)

!C
!C-- START.
      if (hecMESH%my_rank.eq.3) write (*,'(/,a)') '#final mesh (n_node, nn_int, n_elem, ne_int)'
      call MPI_BARRIER (hecMESH%MPI_COMM, ierr)

      write (*,'(a,i5,4i8)') 'PE#', hecMESH%my_rank,                      &
     &                              hecMESH%n_node, hecMESH%nn_internal,  &
     &                              hecMESH%n_elem, hecMESH%ne_internal

!C
!C-- 1. ENTIRE INFO.
      hecMESHnew%hecmw_n_file= 0
      hecMESHnew%header              = hecMESH%header
      hecMESHnew%hecmw_flag_adapt    = hecMESH%hecmw_flag_adapt
      hecMESHnew%hecmw_flag_initcon  = hecMESH%hecmw_flag_initcon
      hecMESHnew%hecmw_flag_parttype = hecMESH%hecmw_flag_parttype
      hecMESHnew%hecmw_flag_partdepth= hecMESH%hecmw_flag_partdepth
      hecMESHnew%hecmw_flag_version  = hecMESH%hecmw_flag_version
      hecMESHnew%zero_temp           = hecMESH%zero_temp
      hecMESHnew%gridfile            = hecMESH%gridfile
      hecMESHnew%hecmw_n_file        = hecMESH%hecmw_n_file

!C
!C-- 2. NODE
      hecMESHnew%n_node     = hecMESH%n_node
      hecMESHnew%nn_internal= hecMESH%nn_internal
      hecMESHnew%n_dof      = 3
      hecMESHnew%n_dof_grp  = 1

      allocate (hecMESHnew%node   (3*hecMESHnew%n_node))
      allocate (hecMESHnew%node_ID(2*hecMESHnew%n_node))
      allocate (hecMESHnew%global_node_ID(hecMESHnew%n_node))
      hecMESHnew%global_node_ID= 1

      do i= 1, hecMESHnew%n_node
        in= hecMESH%adapt_OLDtoNEW_node(i)
        hecMESHnew%node(3*in-2)= hecMESH%node(3*i-2)
        hecMESHnew%node(3*in-1)= hecMESH%node(3*i-1)
        hecMESHnew%node(3*in  )= hecMESH%node(3*i  )

        hecMESHnew%node_ID(2*in-1)= hecMESH%node_ID(2*i-1)
        hecMESHnew%node_ID(2*in  )= hecMESH%node_ID(2*i  )
      enddo

      deallocate (hecMESH%node_ID, hecMESH%node, hecMESH%global_node_ID)

      allocate (hecMESHnew%node_dof_index(0:hecMESHnew%n_dof_grp))
      hecMESHnew%node_dof_index(0)= 0
      hecMESHnew%node_dof_index(1)= hecMESHnew%n_node
      allocate (hecMESHnew%node_dof_item(1))
      hecMESHnew%node_dof_item(1) =  3

!C
!C-- 3. ELEMENT
      hecMESHnew%n_elem     = hecMESH%n_elem
      hecMESHnew%ne_internal= hecMESH%ne_internal
      hecMESHnew%n_elem_type= hecMESH%n_elem_type

      hecMESHnew%n_elem_mat_ID= hecMESH%n_elem

      allocate (hecMESHnew%elem_type_index(0:hecMESHnew%n_elem_type))
      hecMESHnew%elem_type_index= hecMESH%elem_type_index
      allocate (hecMESHnew%elem_type_item(hecMESHnew%n_elem_type))

      hecMESHnew%elem_type_item = hecMESH%elem_type_item
      deallocate (hecMESH%elem_type_index, hecMESH%elem_type_item)

      allocate (hecMESHnew%elem_type        (hecMESHnew%n_elem))
      allocate (hecMESHnew%section_ID       (hecMESHnew%n_elem))
      allocate (hecMESHnew%elem_mat_ID_item (hecMESHnew%n_elem))

      do i= 1, hecMESHnew%n_elem
        in= hecMESH%adapt_OLDtoNEW_elem(i)
        hecMESHnew%elem_type       (in)= hecMESH%elem_type       (i)
        hecMESHnew%section_ID      (in)= hecMESH%section_ID      (i)
        hecMESHnew%elem_mat_ID_item(in)= hecMESH%elem_mat_ID_item(i)
      enddo

      allocate (hecMESHnew%elem_mat_ID_index (0:hecMESHnew%n_elem))
      hecMESHnew%elem_mat_ID_index= 0
      do i= 1, hecMESHnew%n_elem
        hecMESHnew%elem_mat_ID_index(i)= i
      enddo

      deallocate (hecMESH%elem_type, hecMESH%section_ID)
      deallocate (hecMESH%elem_mat_ID_index)
      deallocate (hecMESH%elem_mat_ID_item )

      nnn= hecMESH%elem_node_index(hecMESH%n_elem)
      allocate (hecMESHnew%elem_node_index(0:hecMESHnew%n_elem))
      allocate (hecMESHnew%elem_node_item (nnn))
      hecMESHnew%elem_node_index= 0
      hecMESHnew%elem_node_item = 0
      do i= 1, hecMESHnew%n_elem
        in= hecMESH%adapt_OLDtoNEW_elem(i)
        hecMESHnew%elem_node_index(in)= hecMESH%elem_node_index(i) - hecMESH%elem_node_index(i-1)
      enddo

      do i= 1, hecMESHnew%n_elem
        hecMESHnew%elem_node_index(i)= hecMESHnew%elem_node_index(i-1) + hecMESHnew%elem_node_index(i)
      enddo

      do i= 1, hecMESHnew%n_elem
        in = hecMESH%adapt_OLDtoNEW_elem(i)
        iS = hecMESH%elem_node_index(i-1)
        iE = hecMESH%elem_node_index(i  )
        iS0= hecMESHnew%elem_node_index(in-1)
        do k= iS+1, iE
          kk= k - iS + iS0
          nodK= hecMESH%adapt_OLDtoNEW_node(hecMESH%elem_node_item(k))
          hecMESHnew%elem_node_item(kk)= nodK
        enddo
      enddo
      deallocate (hecMESH%elem_node_index, hecMESH%elem_node_item)

      allocate (hecMESHnew%elem_ID(2*hecMESHnew%n_elem))
      do i= 1, hecMESHnew%n_elem
        in = hecMESH%adapt_OLDtoNEW_elem(i)
        hecMESHnew%elem_ID(2*in-1)= hecMESH%elem_ID(2*i-1)
        hecMESHnew%elem_ID(2*in  )= hecMESH%elem_ID(2*i  )
      enddo
      deallocate (hecMESH%elem_ID)
      allocate (hecMESHnew%global_elem_ID(hecMESHnew%n_elem))
      hecMESHnew%global_elem_ID= 1
      deallocate (hecMESH%global_elem_ID)
      allocate (hecMESHnew%elem_internal_list(hecMESHnew%ne_internal))
      do i= 1, hecMESH%ne_internal
        icel= hecMESH%elem_internal_list (i)
        in  = hecMESH%adapt_OLDtoNEW_elem(icel)
        hecMESHnew%elem_internal_list(i)= in
      enddo
      deallocate (hecMESH%elem_internal_list)

!C
!C-- 5. COMMUNICATION
      hecMESHnew%my_rank = hecMESH%my_rank
      hecMESHnew%zero    = hecMESH%zero
      hecMESHnew%PETOT   = hecMESH%PETOT
      hecMESHnew%PEsmpTOT= hecMESH%PEsmpTOT

      hecMESHnew%errnof  = hecMESH%errnof

      call MPI_COMM_DUP (hecMESH%MPI_COMM, hecMESHnew%MPI_COMM, ierr)

      hecMESHnew%n_subdomain  = hecMESH%n_subdomain
      hecMESHnew%n_neighbor_pe= hecMESH%n_neighbor_pe

      allocate (hecMESHnew%neighbor_pe(hecMESHnew%n_neighbor_pe))
      hecMESHnew%neighbor_pe= hecMESH%neighbor_pe

      nn1= hecMESH%adapt_import_new_index(hecMESH%n_neighbor_pe)
      nn2= hecMESH%adapt_export_new_index(hecMESH%n_neighbor_pe)

      allocate (hecMESHnew%import_index(0:hecMESHnew%n_neighbor_pe), hecMESHnew%import_item(nn1))
      allocate (hecMESHnew%export_index(0:hecMESHnew%n_neighbor_pe), hecMESHnew%export_item(nn2))

      hecMESHnew%import_index= 0
      hecMESHnew%export_index= 0

      do i= 1, hecMESHnew%n_neighbor_pe
        hecMESHnew%import_index(i)= hecMESH%adapt_import_new_index(i)
        hecMESHnew%export_index(i)= hecMESH%adapt_export_new_index(i)
      enddo

      do i= 1, nn1
        hecMESHnew%import_item(i) = hecMESH%adapt_import_new_item(i)
      enddo

      do i= 1, nn2
        hecMESHnew%export_item(i) = hecMESH%adapt_export_new_item(i)
      enddo

      allocate (hecMESHnew%shared_index(0:hecMESHnew%n_neighbor_pe))
      hecMESHnew%shared_index= 0

!C
!C-- 6. GRID ADAPTATION
      hecMESHnew%coarse_grid_level= hecMESH%coarse_grid_level
      hecMESHnew%n_adapt          = hecMESH%n_adapt

      allocate (hecMESHnew%when_i_was_refined_node(hecMESHnew%n_node))
      allocate (hecMESHnew%when_i_was_refined_elem(hecMESHnew%n_elem))
      allocate (hecMESHnew%adapt_parent_type(hecMESHnew%n_elem))
      allocate (hecMESHnew%adapt_type       (hecMESHnew%n_elem))
      allocate (hecMESHnew%adapt_level      (hecMESHnew%n_elem))
      allocate (hecMESHnew%adapt_parent     (hecMESHnew%n_elem*2))
      allocate (hecMESHnew%adapt_children_index (0: hecMESHnew%n_elem))
      allocate (hecMESHnew%adapt_children_item  (16*hecMESHnew%n_elem))

      do i= 1, hecMESHnew%n_node
        in= hecMESH%adapt_OLDtoNEW_node(i)
        hecMESHnew%when_i_was_refined_node(in)= hecMESH%when_i_was_refined_node(i)
      enddo

      hecMESHnew%adapt_children_index= 0
      hecMESHnew%adapt_children_item = 0
      do i= 1, hecMESHnew%n_elem
        in= hecMESH%adapt_OLDtoNEW_elem(i)
        hecMESHnew%when_i_was_refined_elem(in)= hecMESH%when_i_was_refined_elem(i)
        hecMESHnew%adapt_parent_type(in)= hecMESH%adapt_parent_type(i)
        hecMESHnew%adapt_type       (in)= hecMESH%adapt_type       (i)
        hecMESHnew%adapt_level      (in)= hecMESH%adapt_level      (i)

        hecMESHnew%adapt_parent(2*in-1)= hecMESH%adapt_parent(2*i-1)
        hecMESHnew%adapt_parent(2*in  )= hecMESH%adapt_parent(2*i  )

        hecMESHnew%adapt_children_index(i)= hecMESHnew%adapt_children_index(i-1) + 8

        iSorg= (i -1)*16
        iSnew= (in-1)*16
        do kk= 1, 16
          hecMESHnew%adapt_children_item(iSnew+kk)= hecMESH%adapt_children_item(iSorg+kk)
        enddo
      enddo

      deallocate (hecMESH%when_i_was_refined_node)
      deallocate (hecMESH%when_i_was_refined_elem)
      deallocate (hecMESH%adapt_type, hecMESH%adapt_parent_type, hecMESH%adapt_parent, hecMESH%adapt_level)
      deallocate (hecMESH%adapt_children_index, hecMESH%adapt_children_item)

!C
!C-- 7. SECTION
      hecMESHnew%section%n_sect= hecMESH%section%n_sect
     
      allocate (hecMESHnew%section%sect_type(hecMESHnew%section%n_sect))
      allocate (hecMESHnew%section%sect_opt (hecMESHnew%section%n_sect))
      hecMESHnew%section%sect_type= hecMESH%section%sect_type
      hecMESHnew%section%sect_opt = hecMESH%section%sect_opt

      nnn= hecMESH%section%sect_mat_ID_index(hecMESH%section%n_sect)
      allocate (hecMESHnew%section%sect_mat_ID_index(0:hecMESHnew%section%n_sect))
      allocate (hecMESHnew%section%sect_mat_ID_item (nnn))
      hecMESHnew%section%sect_mat_ID_index= hecMESH%section%sect_mat_ID_index
      hecMESHnew%section%sect_mat_ID_item = hecMESH%section%sect_mat_ID_item

      allocate (hecMESHnew%section%sect_I_index(0:hecMESHnew%section%n_sect))
      allocate (hecMESHnew%section%sect_R_index(0:hecMESHnew%section%n_sect))
      hecMESHnew%section%sect_I_index= 0
      hecMESHnew%section%sect_R_index= 0

!C
!C-- 8. MATERIAL
      hecMESHnew%material%n_mat        = hecMESH%material%n_mat
      hecMESHnew%material%n_mat_item   = hecMESH%material%n_mat_item
      hecMESHnew%material%n_mat_subitem= hecMESH%material%n_mat_subitem
      hecMESHnew%material%n_mat_table  = hecMESH%material%n_mat_table

      n_mat    = hecMESHnew%material%n_mat
      n_item   = hecMESHnew%material%n_mat_item
      n_subitem= hecMESHnew%material%n_mat_subitem
      nnn      = hecMESH%material%mat_table_index(n_subitem)

      allocate (hecMESHnew%material%mat_name(n_mat))
      hecMESHnew%material%mat_name= hecMESH%material%mat_name

      allocate (hecMESHnew%material%mat_item_index   (0:n_mat))
      allocate (hecMESHnew%material%mat_subitem_index(0:n_item))
      allocate (hecMESHnew%material%mat_table_index  (0:n_subitem))
      hecMESHnew%material%mat_item_index   = hecMESH%material%mat_item_index
      hecMESHnew%material%mat_subitem_index= hecMESH%material%mat_subitem_index
      hecMESHnew%material%mat_table_index  = hecMESH%material%mat_table_index

      allocate (hecMESHnew%material%mat_val(nnn), hecMESHnew%material%mat_temp(nnn))
      hecMESHnew%material%mat_val = hecMESH%material%mat_val
      hecMESHnew%material%mat_temp= hecMESH%material%mat_temp

!C
!C-- 9. MPC/AMPLITUDE
      hecMESHnew%mpc%n_mpc= hecMESH%mpc%n_mpc
      if (hecMESH%mpc%n_mpc.ne.0) then
        nnn= hecMESH%mpc%mpc_index(hecMESH%mpc%n_mpc)
        allocate (hecMESHnew%mpc%mpc_index(0:hecMESHnew%mpc%n_mpc))
        allocate (hecMESHnew%mpc%mpc_item (nnn))
        allocate (hecMESHnew%mpc%mpc_dof  (nnn))
        allocate (hecMESHnew%mpc%mpc_val  (nnn))
        hecMESHnew%mpc%mpc_index= hecMESH%mpc%mpc_index
        hecMESHnew%mpc%mpc_item = hecMESH%mpc%mpc_item
        hecMESHnew%mpc%mpc_dof  = hecMESH%mpc%mpc_dof
        hecMESHnew%mpc%mpc_val  = hecMESH%mpc%mpc_val
      endif
 
      hecMESHnew%amp%n_amp= 0
      if (hecMESH%amp%n_amp.ne.0) then
        nnn= hecMESH%amp%amp_index(hecMESH%amp%n_amp)
        allocate (hecMESHnew%amp%amp_type_definition(hecMESHnew%amp%n_amp))
        allocate (hecMESHnew%amp%amp_type_time      (hecMESHnew%amp%n_amp))
        allocate (hecMESHnew%amp%amp_type_value     (hecMESHnew%amp%n_amp))
        allocate (hecMESHnew%amp%amp_index(0:hecMESHnew%amp%n_amp))
        allocate (hecMESHnew%amp%amp_val  (nnn))
        allocate (hecMESHnew%amp%amp_table(nnn))
        hecMESHnew%amp%amp_type_definition= hecMESH%amp%amp_type_definition
        hecMESHnew%amp%amp_type_time      = hecMESH%amp%amp_type_time
        hecMESHnew%amp%amp_type_value     = hecMESH%amp%amp_type_value
        hecMESHnew%amp%amp_index          = hecMESH%amp%amp_index
        hecMESHnew%amp%amp_val            = hecMESH%amp%amp_val
        hecMESHnew%amp%amp_table          = hecMESH%amp%amp_table
      endif

!C
!C-- 10. NODE-GROUP
      n_grp  = hecMESH%node_group%n_grp

      if (n_grp.ne.0) then
        n_array= hecMESH%node_group%grp_index(n_grp)
        hecMESHnew%node_group%n_grp= hecMESH%node_group%n_grp

        allocate (hecMESHnew%node_group%grp_name (  n_grp))
        allocate (hecMESHnew%node_group%grp_index(0:n_grp))
        allocate (hecMESHnew%node_group%grp_item (n_array))

        hecMESHnew%node_group%grp_name = hecMESH%node_group%grp_name
        hecMESHnew%node_group%grp_index= hecMESH%node_group%grp_index

        do i= 1, n_array
          in0= hecMESH%node_group%grp_item(i)
          hecMESHnew%node_group%grp_item(i)= hecMESH%adapt_OLDtoNEW_node(in0)
        enddo
      endif

!C
!C-- 11. ELEM-GROUP
      n_grp  = hecMESH%elem_group%n_grp

      if (n_grp.ne.0) then
        n_array= hecMESH%elem_group%grp_index(n_grp)
        hecMESHnew%elem_group%n_grp= hecMESH%elem_group%n_grp

        allocate (hecMESHnew%elem_group%grp_name (  n_grp))
        allocate (hecMESHnew%elem_group%grp_index(0:n_grp))
        allocate (hecMESHnew%elem_group%grp_item (n_array))

        hecMESHnew%elem_group%grp_name = hecMESH%elem_group%grp_name
        hecMESHnew%elem_group%grp_index= hecMESH%elem_group%grp_index

        do i= 1, n_array
          in0= hecMESH%elem_group%grp_item(i)
          hecMESHnew%elem_group%grp_item(i)= hecMESH%adapt_OLDtoNEW_elem(in0)
        enddo
      endif

!C
!C-- 11. SURF-GROUP
      n_grp  = hecMESH%surf_group%n_grp

      if (n_grp.ne.0) then
        n_array= hecMESH%surf_group%grp_index(n_grp)
        hecMESHnew%surf_group%n_grp= hecMESH%surf_group%n_grp

        allocate (hecMESHnew%surf_group%grp_name (  n_grp))
        allocate (hecMESHnew%surf_group%grp_index(0:n_grp))
        allocate (hecMESHnew%surf_group%grp_item (2*n_array))

        hecMESHnew%surf_group%grp_name = hecMESH%surf_group%grp_name
        hecMESHnew%surf_group%grp_index= hecMESH%surf_group%grp_index

        do i= 1, n_array
          in0 = hecMESH%surf_group%grp_item(2*i-1)
          isuf= hecMESH%surf_group%grp_item(2*i  )
          hecMESHnew%surf_group%grp_item(2*i-1)= hecMESH%adapt_OLDtoNEW_elem(in0)
          hecMESHnew%surf_group%grp_item(2*i  )= isuf
        enddo
      endif

      end subroutine hecmw_adapt_new_mesh
    

