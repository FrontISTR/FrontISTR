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
!C*** hecmw_adapt_allocate
!C***
!C

      subroutine hecmw_adapt_allocate (hecMESH)

      use  hecmw_util
      type      (hecmwST_local_mesh) :: hecMESH

      integer(kind=kint), pointer    :: WORKaI(:)
      real(kind=kreal), pointer      :: WORKaR(:)

!C
!C +-------+
!C | ERROR |
!C +-------+
!C===
      if (hecMESH%hecmw_flag_initcon .eq.1) call hecmw_adapt_ERROR_EXIT (hecMESH, 201)
      if (hecMESH%hecmw_flag_parttype.eq.2) call hecmw_adapt_ERROR_EXIT (hecMESH, 202)
      if (hecMESH%n_dof_grp.ne.1)           call hecmw_adapt_ERROR_EXIT (hecMESH, 203)
      if (hecMESH%n_elem_mat_ID.ne.hecMESH%n_elem)                                     &
     &                                      call hecmw_adapt_ERROR_EXIT (hecMESH, 204)
!C===

!C
!C +--------------+
!C | INITIAL MESH |
!C +--------------+
!C===
      call MPI_BARRIER (hecMESH%MPI_COMM, ierr)
      write (*,'(a,i5,4i8)') 'PE#', hecMESH%my_rank,                      &
     &                              hecMESH%n_node, hecMESH%nn_internal,  &
     &                              hecMESH%n_elem, hecMESH%ne_internal

      hecMESH%nn_array= 9 * hecMESH%n_node + 1
      hecMESH%ne_array= 9 * hecMESH%n_elem + 1
      hecMESH%nx_array= max (hecMESH%nn_array,hecMESH%ne_array)

      if (hecMESH%hecmw_flag_adapt.eq.0) then
        allocate (hecMESH%when_i_was_refined_node (hecMESH%nn_array))
        allocate (hecMESH%when_i_was_refined_elem (hecMESH%ne_array))

        allocate (hecMESH%adapt_parent      (2*hecMESH%ne_array))
        allocate (hecMESH%adapt_parent_type (  hecMESH%ne_array))
        allocate (hecMESH%adapt_type        (  hecMESH%ne_array))
        allocate (hecMESH%adapt_level       (  hecMESH%ne_array))

        allocate (hecMESH%adapt_children_item  (2*8*hecMESH%ne_array))
        allocate (hecMESH%adapt_children_index (0:hecMESH%ne_array))

        hecMESH%when_i_was_refined_node= 0
        hecMESH%when_i_was_refined_elem= 0

        hecMESH%adapt_parent           = 0
        hecMESH%adapt_parent_type      = 0
        hecMESH%adapt_type             = 0
        hecMESH%adapt_level            = 0

        hecMESH%adapt_children_item = 0
        hecMESH%adapt_children_index= 0

        do icel= 1, hecMESH%n_elem
          hecMESH%adapt_parent        ( 2*icel)= -1
          hecMESH%adapt_children_index(   icel)=  8 * icel

          iS= 16*(icel-1)
          hecMESH%adapt_children_item (iS+ 2)= -1 
          hecMESH%adapt_children_item (iS+ 4)= -1 
          hecMESH%adapt_children_item (iS+ 6)= -1 
          hecMESH%adapt_children_item (iS+ 8)= -1 
          hecMESH%adapt_children_item (iS+10)= -1 
          hecMESH%adapt_children_item (iS+12)= -1 
          hecMESH%adapt_children_item (iS+14)= -1 
          hecMESH%adapt_children_item (iS+16)= -1 
        enddo
      endif
!C===

!C
!C +-------------------+
!C | RE-ALLOCATE: node |
!C +-------------------+
!C===
 
!C
!C-- COORDINATEs
      allocate (WORKaR(3*hecMESH%n_node))
      do i= 1, hecMESH%n_node
        WORKaR(3*i-2)= hecMESH%node(3*i-2)
        WORKaR(3*i-1)= hecMESH%node(3*i-1)
        WORKaR(3*i  )= hecMESH%node(3*i  )
      enddo

      deallocate (hecMESH%node)
        allocate (hecMESH%node(3*hecMESH%nn_array))
                 hecMESH%node= 0.d0
      do i= 1, hecMESH%n_node
        hecMESH%node(3*i-2)= WORKaR(3*i-2)
        hecMESH%node(3*i-1)= WORKaR(3*i-1)
        hecMESH%node(3*i  )= WORKaR(3*i  )
      enddo
      deallocate (WORKaR)

!C
!C-- node_ID
      allocate (WORKaI(2*hecMESH%n_node))
      do i= 1, hecMESH%n_node
        WORKaI(2*i-1)= hecMESH%node_ID(2*i-1)
        WORKaI(2*i  )= hecMESH%node_ID(2*i  )
      enddo

      deallocate (hecMESH%node_ID)
        allocate (hecMESH%node_ID(2*hecMESH%nn_array))
                  hecMESH%node_ID= 0
      do i= 1, hecMESH%n_node
        hecMESH%node_ID(2*i-1)= WORKaI(2*i-1)
        hecMESH%node_ID(2*i  )= WORKaI(2*i  )
      enddo
      deallocate (WORKaI)

!C
!C-- global_node_ID
      allocate (WORKaI(hecMESH%n_node))
      do i= 1, hecMESH%n_node
        WORKaI(i)= hecMESH%global_node_ID(i)
      enddo

      deallocate (hecMESH%global_node_ID)
        allocate (hecMESH%global_node_ID(hecMESH%nn_array))
                  hecMESH%global_node_ID= 0
      do i= 1, hecMESH%n_node
        hecMESH%global_node_ID(i)= WORKaI(i)
      enddo
      deallocate (WORKaI)
!C===

!C
!C +----------------------+
!C | RE-ALLOCATE: element |
!C +----------------------+
!C===

!C
!C-- elem_type, section_ID, global_elem_ID
      allocate (WORKaI(2*hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(2*i-1)= hecMESH%elem_type (i)
        WORKaI(2*i  )= hecMESH%section_ID(i)
      enddo

      deallocate (hecMESH%elem_type, hecMESH%section_ID)
        allocate (hecMESH%elem_type (hecMESH%ne_array))
        allocate (hecMESH%section_ID(hecMESH%ne_array))
                  hecMESH%elem_type = 0
                  hecMESH%section_ID= 0
      do i= 1, hecMESH%n_elem
        hecMESH%elem_type (i)= WORKaI(2*i-1)
        hecMESH%section_ID(i)= WORKaI(2*i  )
      enddo
      deallocate (WORKaI)

      allocate (WORKaI(hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(i)= hecMESH%global_elem_ID(i)
      enddo

      deallocate (hecMESH%global_elem_ID)
        allocate (hecMESH%global_elem_ID(hecMESH%ne_array))
                  hecMESH%global_elem_ID= 0
      do i= 1, hecMESH%n_elem
        hecMESH%global_elem_ID(i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

!C
!C-- elem_ID
      allocate (WORKaI(2*hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(2*i-1)= hecMESH%elem_ID(2*i-1)
        WORKaI(2*i  )= hecMESH%elem_ID(2*i  )
      enddo

      deallocate (hecMESH%elem_ID)
        allocate (hecMESH%elem_ID(2*hecMESH%ne_array))
                  hecMESH%elem_ID= 0
      do i= 1, hecMESH%n_elem
        hecMESH%elem_ID(2*i-1)= WORKaI(2*i-1)
        hecMESH%elem_ID(2*i  )= WORKaI(2*i  )
      enddo
      deallocate (WORKaI)

!C
!C-- elem_node_index
      allocate (WORKaI(hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(i)= hecMESH%elem_node_index(i)
      enddo

      deallocate (hecMESH%elem_node_index)
        allocate (hecMESH%elem_node_index(0:hecMESH%ne_array))
                  hecMESH%elem_node_index= 0
      do i= 1, hecMESH%n_elem
        hecMESH%elem_node_index (i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

!C
!C-- elem_node_item
      nnn= hecMESH%elem_node_index(hecMESH%n_elem)
      allocate (WORKaI(nnn))
      do i= 1, nnn
        WORKaI(i)= hecMESH%elem_node_item(i)
      enddo

      deallocate (hecMESH%elem_node_item)
        allocate (hecMESH%elem_node_item(6*hecMESH%ne_array))
                  hecMESH%elem_node_item= 0
      do i= 1, nnn
        hecMESH%elem_node_item (i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

!C
!C-- elem_mat_ID_index
      allocate (WORKaI(hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(i)= hecMESH%elem_mat_ID_index(i)
      enddo

      deallocate (hecMESH%elem_mat_ID_index)
        allocate (hecMESH%elem_mat_ID_index(0:hecMESH%ne_array))
                  hecMESH%elem_mat_ID_index= 0
      do i= 1, ne_array
        hecMESH%elem_mat_ID_index (i)= i
      enddo
      deallocate (WORKaI)

!C
!C-- elem_mat_ID_item
      nnn= hecMESH%n_elem
      allocate (WORKaI(nnn))
      do i= 1, nnn
        WORKaI(i)= hecMESH%elem_mat_ID_item(i)
      enddo

      deallocate (hecMESH%elem_mat_ID_item)
        allocate (hecMESH%elem_mat_ID_item(hecMESH%ne_array*nnn/hecMESH%n_elem))
                  hecMESH%elem_mat_ID_item= 0
      do i= 1, nnn
        hecMESH%elem_mat_ID_item (i)= WORKaI(i)
      enddo
      deallocate (WORKaI)
!C===
      
      if (hecMESH%hecmw_flag_adapt.eq.1) then
!C
!C +--------------------+
!C | RE-ALLOCATE: adapt |
!C +--------------------+
!C===

!C
!C-- when_i_was
      allocate (WORKaI(hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(i)= hecMESH%when_i_was_refined_elem(i) 
      enddo

      deallocate (hecMESH%when_i_was_refined_elem)
        allocate (hecMESH%when_i_was_refined_elem(hecMESH%ne_array))
                  hecMESH%when_i_was_refined_elem= 0
      do i= 1, hecMESH%n_elem
        hecMESH%when_i_was_refined_elem(i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

      allocate (WORKaI(hecMESH%n_node))
      do i= 1, hecMESH%n_node
        WORKaI(i)= hecMESH%when_i_was_refined_node(i) 
      enddo

      deallocate (hecMESH%when_i_was_refined_node)
        allocate (hecMESH%when_i_was_refined_node(hecMESH%nn_array))
                  hecMESH%when_i_was_refined_node= 0
      do i= 1, hecMESH%n_node
        hecMESH%when_i_was_refined_node(i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

!C
!C-- adapt_parent, parent_type, type, level

      allocate (WORKaI(5*hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(5*i-4)= hecMESH%adapt_parent     (2*i-1)
        WORKaI(5*i-3)= hecMESH%adapt_parent     (2*i  )
        WORKaI(5*i-2)= hecMESH%adapt_parent_type(  i)
        WORKaI(5*i-1)= hecMESH%adapt_type       (  i)
        WORKaI(4*i  )= hecMESH%adapt_level      (  i)
      enddo

      deallocate (hecMESH%adapt_parent, hecMESH%adapt_parent_type) 
      deallocate (hecMESH%adapt_type  , hecMESH%adapt_level      ) 

        allocate (hecMESH%adapt_parent     (2*hecMESH%ne_array))
        allocate (hecMESH%adapt_parent_type(  hecMESH%ne_array))
        allocate (hecMESH%adapt_type       (  hecMESH%ne_array))
        allocate (hecMESH%adapt_level      (  hecMESH%ne_array))
                  hecMESH%adapt_parent     = 0
                  hecMESH%adapt_parent_type= 0
                  hecMESH%adapt_type       = 0
                  hecMESH%adapt_level      = 0

      do i= 1, hecMESH%n_elem
        hecMESH%adapt_parent     (2*i-1)= WORKaI(5*i-4)
        hecMESH%adapt_parent     (2*i  )= WORKaI(5*i-3)
        hecMESH%adapt_parent_type(  i  )= WORKaI(5*i-2)
        hecMESH%adapt_type       (  i  )= WORKaI(5*i-1)
        hecMESH%adapt_level      (  i  )= WORKaI(5*i  )       
      enddo
      deallocate (WORKaI)
     
!C
!C-- adapt_children_index
      allocate (WORKaI(hecMESH%n_elem))
      do i= 1, hecMESH%n_elem
        WORKaI(i)= hecMESH%adapt_children_index(i) 
      enddo

      deallocate (hecMESH%adapt_children_index)
        allocate (hecMESH%adapt_children_index(0:hecMESH%ne_array))
                 hecMESH%adapt_children_index= 0
      do i= 1, hecMESH%n_elem
        hecMESH%adapt_children_index(i)= WORKaI(i)
      enddo
      deallocate (WORKaI)

!C
!C-- adapt_children_item
      nnn= hecMESH%adapt_children_index(hecMESH%n_elem)
      allocate (WORKaI(2*nnn))
      do i= 1, nnn
        WORKaI(2*i-1)= hecMESH%adapt_children_item(2*i-1) 
        WORKaI(2*i  )= hecMESH%adapt_children_item(2*i  ) 
      enddo

      deallocate (hecMESH%adapt_children_item)
        allocate (hecMESH%adapt_children_item(2*8*hecMESH%ne_array))
                 hecMESH%adapt_children_item= 0
      do i= 1, nnn
        hecMESH%adapt_children_item(2*i-1)= WORKaI(2*i-1)
        hecMESH%adapt_children_item(2*i  )= WORKaI(2*i  )
      enddo
      deallocate (WORKaI)
     
      endif

      hecMESH%hecmw_flag_adapt= 1
!C===

      if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#RE-ALLOCATE'
!C
!C +-------------------+
!C | REALLOCATE: group |
!C +-------------------+
!C===

!C
!C-- NODE group
      if (hecMESH%node_group%n_grp.ne.0) then
      if (hecMESH%node_group%grp_index(hecMESH%node_group%n_grp).ne.0) then
        nnn= hecMESH%node_group%grp_index(hecMESH%node_group%n_grp)
        allocate (WORKaI(nnn))
        do i= 1, nnn
          WORKaI(i)= hecMESH%node_group%grp_item(i)
        enddo

        deallocate (hecMESH%node_group%grp_item)
          allocate (hecMESH%node_group%grp_item(hecMESH%nx_array))
                    hecMESH%node_group%grp_item= 0
        do i= 1, nnn
          hecMESH%node_group%grp_item(i)= WORKaI(i)
        enddo
        deallocate (WORKaI)
      endif
      endif

!C
!C-- ELEMENT group
      if (hecMESH%elem_group%n_grp.ne.0) then
      if (hecMESH%elem_group%grp_index(hecMESH%elem_group%n_grp).ne.0) then
        nnn= hecMESH%elem_group%grp_index(hecMESH%elem_group%n_grp)
        allocate (WORKaI(nnn))
        do i= 1, nnn
          WORKaI(i)= hecMESH%elem_group%grp_item(i)
        enddo

          deallocate (hecMESH%elem_group%grp_item)
            allocate (hecMESH%elem_group%grp_item(hecMESH%nx_array))
                      hecMESH%elem_group%grp_item= 0
          do i= 1, nnn
            hecMESH%elem_group%grp_item(i)= WORKaI(i)
          enddo

        deallocate (WORKaI)
      endif
      endif

!C
!C-- SURFACE group
      if (hecMESH%surf_group%n_grp.ne.0) then
      if (hecMESH%surf_group%grp_index(hecMESH%surf_group%n_grp).ne.0) then
        nnn= hecMESH%surf_group%grp_index(hecMESH%surf_group%n_grp)
        allocate (WORKaI(2*nnn))
        do i= 1, nnn
          WORKaI(2*i-1)= hecMESH%surf_group%grp_item(2*i-1)
          WORKaI(2*i  )= hecMESH%surf_group%grp_item(2*i  )
        enddo

        deallocate (hecMESH%surf_group%grp_item)
          allocate (hecMESH%surf_group%grp_item(2*hecMESH%nx_array))
                    hecMESH%surf_group%grp_item= 0
        do i= 1, nnn
          hecMESH%surf_group%grp_item(2*i-1)= WORKaI(2*i-1)
          hecMESH%surf_group%grp_item(2*i  )= WORKaI(2*i  )
        enddo
        deallocate (WORKaI)
      endif
      endif
!C===

      if (hecMESH%my_rank.eq.0) write (*,'(a)') '#EDGE-INFO'
!C
!C +------------------+
!C | EDGE information |
!C +------------------+
!C===
      NE= max(hecMESH%n_node, hecMESH%n_elem)
  100 continue

      allocate(hecMESH%adapt_edge_node(2*NE))
               hecMESH%adapt_edge_node= 0
               hecMESH%n_adapt_edge   = 0

      do icel= 1, hecMESH%n_elem
        if (hecMESH%adapt_type(icel).eq.0) then
        ityp= hecMESH%elem_type(icel)
!C
!C-- 3D : tetrahedron
        if (ityp.eq.341) then
          iS = hecMESH%elem_node_index(icel-1)
          in1= hecMESH%elem_node_item (iS+1)
          in2= hecMESH%elem_node_item (iS+2)
          in3= hecMESH%elem_node_item (iS+3)
          in4= hecMESH%elem_node_item (iS+4)

          call hecmw_adapt_EDGE_INFO (hecMESH, in1,in2, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in1,in3, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in1,in4, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in2,in3, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in3,in4, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in4,in2, iedge, 0)
        endif        

!C
!C-- 3D : prisms
        if (ityp.eq.351) then
          iS = hecMESH%elem_node_index(icel-1)
          in1= hecMESH%elem_node_item (iS+1)
          in2= hecMESH%elem_node_item (iS+2)
          in3= hecMESH%elem_node_item (iS+3)
          in4= hecMESH%elem_node_item (iS+4)
          in5= hecMESH%elem_node_item (iS+5)
          in6= hecMESH%elem_node_item (iS+6)

          call hecmw_adapt_EDGE_INFO (hecMESH, in1,in2, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in2,in3, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in3,in1, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in4,in5, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in5,in6, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in6,in4, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in1,in4, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in2,in5, iedge, 0)
          call hecmw_adapt_EDGE_INFO (hecMESH, in3,in6, iedge, 0)
        endif        

        if (hecMESH%n_adapt_edge.ge.NE-6 .and. icel.lt.hecMESH%n_elem) then
          iii= hecMESH%n_elem/icel + 1
          NE = iii * NE + 1
          deallocate(hecMESH%adapt_edge_node)

          goto 100
        endif

        endif
      enddo

      call MPI_BARRIER (hecMESH%MPI_COMM, ierr)

      allocate(hecMESH%adapt_IWK(hecMESH%n_adapt_edge))
      hecMESH%adapt_IWK=   0

!C===

      hecMESH%n_adapt_elem_341= 0
      hecMESH%n_adapt_elem_351= 0

      do icel= 1, hecMESH%n_elem
        ityp= hecMESH%elem_type(icel)
        if (ityp.eq.341) then
          hecMESH%n_adapt_elem_341= hecMESH%n_adapt_elem_341 + 1
        endif
        if (ityp.eq.351) then
          hecMESH%n_adapt_elem_351= hecMESH%n_adapt_elem_351 + 1
        endif
      enddo

      end subroutine hecmw_adapt_allocate

