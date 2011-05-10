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
!C*** hecmw_adapt_bc_pointer
!C***
!C     
!C    creates NEW BOUNDARY POINTERs
!C
      subroutine hecmw_adapt_bc_pointer (hecMESH)

      use  hecmw_util
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(:  ), allocatable :: IW1, IW2, IW3
      integer(kind=kint), dimension(:,:), allocatable :: IUPD
      type(hecmwST_local_mesh)             :: hecMESH

      if (hecMESH%node_group%n_grp.ne.0) then
      if (hecMESH%node_group%grp_index(hecMESH%node_group%n_grp).ne.0) then
!C
!C +------------+
!C | NODE-GROUP |
!C +------------+
!C===
      allocate (IUPD(hecMESH%n_adapt_node_cur,                          &
     &               hecMESH%node_group%n_grp))
      allocate (IW1 (hecMESH%n_adapt_node_cur))
      allocate (IW2 (hecMESH%node_group%n_grp))

      IUPD = 0
      icou = 0
      IW2  = 0

      do ig= 1, hecMESH%node_group%n_grp
        icou0= 0
          IW1= 0
        do k= hecMESH%node_group%grp_index(ig-1)+1,                     &
     &        hecMESH%node_group%grp_index(ig)
               nod      = hecMESH%node_group%grp_item(k)
          icou          = icou  + 1
          icou0         = icou0 + 1
          IW1 (nod     )= 1
          IUPD(icou0,ig)= nod
        enddo

        do ie= 1, hecMESH%n_adapt_edge
          nod1= hecMESH%adapt_edge_node(2*ie-1)
          nod2= hecMESH%adapt_edge_node(2*ie  )

          if (IW1(nod1).eq.1  .and. IW1(nod2).eq.1 .and.                &
     &        hecMESH%adapt_iemb(ie ).ne.0) then          
            nod3          = hecMESH%adapt_IWK(ie)
            icou          = icou + 1
            icou0         = icou0 + 1
            IUPD(icou0,ig)= nod3
          endif
        enddo
        IW2(ig)= icou0
      enddo

!C
!C-- NEW POINTERs
      nnn1= hecMESH%node_group%grp_index(hecMESH%node_group%n_grp)
      deallocate (hecMESH%node_group%grp_item)
        allocate (hecMESH%node_group%grp_item(icou))

      hecMESH%node_group%grp_index= 0
      do ig= 1, hecMESH%node_group%n_grp
        hecMESH%node_group%grp_index(ig)=                               &
     &  hecMESH%node_group%grp_index(ig-1) + IW2(ig)
      enddo
      do ig= 1, hecMESH%node_group%n_grp
        do k= 1, IW2(ig)
          kk= hecMESH%node_group%grp_index(ig-1) + k
          hecMESH%node_group%grp_item(kk)= IUPD(k,ig)
        enddo
      enddo

      nnn2= hecMESH%node_group%grp_index(hecMESH%node_group%n_grp)
!      write (*,'(a, i5, 2i8)') '    node group', hecMESH%my_rank, nnn1, nnn2

      deallocate (IUPD, IW1, IW2)
!C=== 
      endif
      endif

      if (hecMESH%elem_group%n_grp.ne.0) then
      if (hecMESH%elem_group%grp_index(hecMESH%elem_group%n_grp).ne.0) then
!C
!C +---------------+
!C | ELEMENT-GROUP |
!C +---------------+
!C===
      allocate (IUPD(hecMESH%n_elem,                                    &
     &               hecMESH%elem_group%n_grp))
      allocate (IW1 (hecMESH%n_elem))
      allocate (IW2 (hecMESH%elem_group%n_grp))

      IUPD = 0
      icou = 0
      IW2  = 0

      do ig= 1, hecMESH%elem_group%n_grp
        icou0= 0
          IW1= 0
        do k= hecMESH%elem_group%grp_index(ig-1)+1,                     &
     &        hecMESH%elem_group%grp_index(ig)
               icel     = hecMESH%elem_group%grp_item(k)
          icou          = icou  + 1
          icou0         = icou0 + 1
          IW1 (icel    )= 1
          IUPD(icou0,ig)= icel
        enddo
        
        do icel= 1, hecMESH%n_elem
          if (hecMESH%adapt_type(icel).ne.0.and.IW1(icel).eq.1) then
            iS= hecMESH%adapt_children_index(icel-1) + 1
            iE= hecMESH%adapt_children_index(icel)
            icS= hecMESH%adapt_children_local(iS)

            if (hecMESH%when_i_was_refined_elem(icS).eq.                &
     &          hecMESH%n_adapt) then
              do k= iS, iE
                if (hecMESH%adapt_children_item(2*k-1).ne.0) then
                  iclocal= hecMESH%adapt_children_local(k)
                  icou          = icou  + 1
                  icou0         = icou0 + 1
                  IW1 (iclocal )= 1
                  IUPD(icou0,ig)= iclocal
                endif
              enddo
            endif
          endif
        enddo
        IW2(ig)= icou0
      enddo

!C
!C-- NEW POINTERs
      nnn1= hecMESH%elem_group%grp_index(hecMESH%elem_group%n_grp)
      deallocate (hecMESH%elem_group%grp_item)
        allocate (hecMESH%elem_group%grp_item(icou))

      hecMESH%elem_group%grp_index= 0
      do ig= 1, hecMESH%elem_group%n_grp
        hecMESH%elem_group%grp_index(ig)=                               &
     &  hecMESH%elem_group%grp_index(ig-1) + IW2(ig)
      enddo
      do ig= 1, hecMESH%elem_group%n_grp
        do k= 1, IW2(ig)
          kk= hecMESH%elem_group%grp_index(ig-1) + k
          hecMESH%elem_group%grp_item(kk)= IUPD(k,ig)
        enddo
      enddo
      nnn2= hecMESH%elem_group%grp_index(hecMESH%elem_group%n_grp)
!      write (*,'(a, i5, 2i8)') '    elem group', hecMESH%my_rank, nnn1, nnn2

      deallocate (IUPD, IW1, IW2)
!C=== 
      endif
      endif

      if (hecMESH%surf_group%n_grp.ne.0) then
      if (hecMESH%surf_group%grp_index(hecMESH%surf_group%n_grp).ne.0) then
!C
!C +---------------+
!C | SURFACE-GROUP |
!C +---------------+
!C===
      allocate (IUPD(2*hecMESH%n_node,                                  &
     &                 hecMESH%surf_group%n_grp))
      allocate (IW1 (2*hecMESH%n_elem))
      allocate (IW2 (hecMESH%surf_group%n_grp))
      allocate (IW3 (hecMESH%n_node))

      IUPD = 0
      icou = 0
      IW2  = 0

      do ig= 1, hecMESH%surf_group%n_grp
        icou0= 0
          IW1= 0
          IW3= 0
        do k= hecMESH%surf_group%grp_index(ig-1)+1,                     &
     &        hecMESH%surf_group%grp_index(ig)
               icel     = hecMESH%surf_group%grp_item(2*k-1)
               isuf     = hecMESH%surf_group%grp_item(2*k  )
          icou          = icou  + 1
          icou0         = icou0 + 1
          IW1 (icel    )= isuf
          IUPD(2*icou0-1,ig)= icel
          IUPD(2*icou0  ,ig)= isuf
          ip_type= hecMESH%elem_type(icel)
          if (ip_type.eq.351) then
            iSP= hecMESH%elem_node_index(icel-1)
            nP1= hecMESH%elem_node_item (iSP+1)
            nP2= hecMESH%elem_node_item (iSP+2)
            nP3= hecMESH%elem_node_item (iSP+3)
            nP4= hecMESH%elem_node_item (iSP+4)
            nP5= hecMESH%elem_node_item (iSP+5)
            nP6= hecMESH%elem_node_item (iSP+6)

            if       (isuf.eq.1) then
              IW3(nP2)= 1
              IW3(nP3)= 1
              IW3(nP5)= 1
              IW3(nP6)= 1
             else if (isuf.eq.2) then
              IW3(nP3)= 1
              IW3(nP1)= 1
              IW3(nP6)= 1
              IW3(nP4)= 1
             else if (isuf.eq.3) then
              IW3(nP1)= 1
              IW3(nP2)= 1
              IW3(nP4)= 1
              IW3(nP5)= 1
             else if (isuf.eq.4) then
              IW3(nP1)= 1
              IW3(nP2)= 1
              IW3(nP3)= 1
             else if (isuf.eq.5) then
              IW3(nP4)= 1
              IW3(nP5)= 1
              IW3(nP6)= 1
            endif
          endif

          if (ip_type.eq.341) then
            iSP= hecMESH%elem_node_index(icel-1)
            nP1= hecMESH%elem_node_item (iSP+1)
            nP2= hecMESH%elem_node_item (iSP+2)
            nP3= hecMESH%elem_node_item (iSP+3)
            nP4= hecMESH%elem_node_item (iSP+4)
            if       (isuf.eq.1) then
              IW3(nP2)= 1
              IW3(nP3)= 1
              IW3(nP4)= 1
             else if (isuf.eq.2) then
              IW3(nP1)= 1
              IW3(nP3)= 1
              IW3(nP4)= 1
             else if (isuf.eq.3) then
              IW3(nP1)= 1
              IW3(nP2)= 1
              IW3(nP4)= 1
             else if (isuf.eq.4) then
              IW3(nP1)= 1
              IW3(nP2)= 1
              IW3(nP3)= 1
            endif
          endif

        enddo
        
        do ie= 1, hecMESH%n_adapt_edge
          nod1= hecMESH%adapt_edge_node(2*ie-1)
          nod2= hecMESH%adapt_edge_node(2*ie  )

          if (IW3(nod1).eq.1  .and. IW3(nod2).eq.1 .and.                &
     &        hecMESH%adapt_iemb(ie).ne.0) then          
            nod3     = hecMESH%adapt_IWK(ie)
            IW3(nod3)= 1
          endif
        enddo

        do kkk= hecMESH%surf_group%grp_index(ig-1)+1,                   &
     &        hecMESH%surf_group%grp_index(ig)
          icel= hecMESH%surf_group%grp_item(2*kkk-1)
          if (hecMESH%adapt_type(icel).ne.0.and.IW1(icel).ne.0) then
            ip_type= hecMESH%elem_type(icel)
            iS= hecMESH%adapt_children_index(icel-1) + 1
            iE= hecMESH%adapt_children_index(icel)
            icS= hecMESH%adapt_children_local(iS)

            if (hecMESH%when_i_was_refined_elem(icS).eq.                &
     &          hecMESH%n_adapt) then
              do k= iS, iE
                if (hecMESH%adapt_children_item(2*k-1).ne.0) then
                  iclocal= hecMESH%adapt_children_local(k)
                  if (ip_type.eq.351) then
                    iSC= hecMESH%elem_node_index(iclocal-1)
                    nC1= hecMESH%elem_node_item (iSC+1)
                    nC2= hecMESH%elem_node_item (iSC+2)
                    nC3= hecMESH%elem_node_item (iSC+3)
                    nC4= hecMESH%elem_node_item (iSC+4)
                    nC5= hecMESH%elem_node_item (iSC+5)
                    nC6= hecMESH%elem_node_item (iSC+6)

                    nsuf1= IW3(nC2) + IW3(nC3) + IW3(nC5) + IW3(nC6)
                    nsuf2= IW3(nC1) + IW3(nC3) + IW3(nC4) + IW3(nC6)
                    nsuf3= IW3(nC1) + IW3(nC2) + IW3(nC4) + IW3(nC5)
                    nsuf4= IW3(nC1) + IW3(nC2) + IW3(nC3)
                    nsuf5= IW3(nC4) + IW3(nC5) + IW3(nC6)

                    if (nsuf1.eq.4) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 1
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 1
                    endif

                    if (nsuf2.eq.4) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 2
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 2
                    endif

                    if (nsuf3.eq.4) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 3
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 3
                    endif

                    if (nsuf4.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 4
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 4
                    endif

                    if (nsuf5.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 5
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 5
                    endif
                  endif

                  if (ip_type.eq.341) then
                    iSC= hecMESH%elem_node_index(iclocal-1)
                    nC1= hecMESH%elem_node_item (iSC+1)
                    nC2= hecMESH%elem_node_item (iSC+2)
                    nC3= hecMESH%elem_node_item (iSC+3)
                    nC4= hecMESH%elem_node_item (iSC+4)

                    nsuf1= IW3(nC2) + IW3(nC3) + IW3(nC4)
                    nsuf2= IW3(nC1) + IW3(nC3) + IW3(nC4)
                    nsuf3= IW3(nC1) + IW3(nC2) + IW3(nC4)
                    nsuf4= IW3(nC1) + IW3(nC2) + IW3(nC3)

                    if (nsuf1.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 1
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 1
                    endif

                    if (nsuf2.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 2
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 2
                    endif

                    if (nsuf3.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 3
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 3
                    endif

                    if (nsuf4.eq.3) then
                      icou          = icou  + 1
                      icou0         = icou0 + 1
                      IW1 (iclocal )= 4
                      IUPD(2*icou0-1,ig)= iclocal
                      IUPD(2*icou0  ,ig)= 4
                    endif

                  endif
                endif
              enddo
            endif
          endif
        enddo
        IW2(ig)= icou0
      enddo

!C
!C-- NEW POINTERs
      nnn1= hecMESH%surf_group%grp_index(hecMESH%surf_group%n_grp)
      deallocate (hecMESH%surf_group%grp_item)
        allocate (hecMESH%surf_group%grp_item(2*icou))

      hecMESH%surf_group%grp_index= 0
      do ig= 1, hecMESH%surf_group%n_grp
        hecMESH%surf_group%grp_index(ig)=                               &
     &  hecMESH%surf_group%grp_index(ig-1) + IW2(ig)
      enddo
      do ig= 1, hecMESH%surf_group%n_grp
        do k= 1, IW2(ig)
          kk= hecMESH%surf_group%grp_index(ig-1) + k
          hecMESH%surf_group%grp_item(2*kk-1)= IUPD(2*k-1,ig)
          hecMESH%surf_group%grp_item(2*kk  )= IUPD(2*k  ,ig)
        enddo
      enddo
      nnn2= hecMESH%surf_group%grp_index(hecMESH%surf_group%n_grp)
!      write (*,'(a, i5, 2i8)') '    surf group', hecMESH%my_rank, nnn1, nnn2

      deallocate (IUPD, IW1, IW2, IW3)
!C=== 
      endif
      endif

      return
      end



