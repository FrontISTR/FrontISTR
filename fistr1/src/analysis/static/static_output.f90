!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!C***
!> OUTPUT for FSTR solver
!C***
module m_static_output

   use m_static_get_prop
   use m_static_LIB_1d
   use m_static_LIB_2d
   use m_static_LIB_3d
   use m_static_LIB_shell

   contains

!C
!C 3D SOLID
!C
   subroutine fstr_output_3d(hecMESH, hecMAT, fstrSOLID, fstrPARAM)

      use m_fstr
   
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)         :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID
!C
!C-- SET DISPLACEMENT
!C
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(3*i-2)= hecMAT%X(3*i-2)
        fstrSOLID%unode(3*i-1)= hecMAT%X(3*i-1)
        fstrSOLID%unode(3*i  )= hecMAT%X(3*i  )
      enddo
!C
!C-- DISPLACEMENT
!C
      write(ILOG,*) '#### DISPLACEMENT 3D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      Z-DISP'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,3e12.4)') jj,(fstrSOLID%unode(3*(ii-1)+k),k=1,3)
      enddo
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 3D'
      call fstr_reaction_force_3d( hecMESH, fstrSOLID)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 3D'
      call fstr_nodal_stress_3d( hecMESH, fstrSOLID,fstrPARAM)

      end subroutine fstr_output_3d
!C
!C NODAL STRESS 3D
!C
      subroutine fstr_nodal_stress_3d(hecMESH,fstrSOLID,fstrPARAM, ifwrite)
      use m_fstr
      use m_static_lib
      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      logical, optional      :: ifwrite
!** Local variables
      REAL(kind=kreal) xx(20),yy(20),zz(20)
      REAL(kind=kreal) edisp(60),estrain(6),estress(6)
      REAL(kind=kreal) tt(20),tt0(20)
      integer(kind=kint) nodLOCAL(20)
      REAL(kind=kreal) edstrain(20,6),edstress(20,6)
!** Array for nodal recovery
      integer (kind=KINT), dimension(:), allocatable :: nnumber,tnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp 
      real(kind=KREAL), allocatable :: ndstrain(:,:),ndstress(:,:),tstrain(:,:),tstress(:,:)

!*Allocate array
      allocate ( nnumber(hecMESH%n_node) )
      allocate( ndstrain(hecMESH%n_node,6), ndstress(hecMESH%n_node,7) )
      allocate ( tnumber(hecMESH%n_node) )
      allocate( tstrain(hecMESH%n_node,6), tstress(hecMESH%n_node,6) )
!*ZERO CLEAR
      ndstrain=0.d0; ndstress=0.d0
      tstrain=0.d0; tstress=0.d0
      arrayTotal=0
      nnumber=0
      tnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node) )
      temp=0
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
        do ig0= 1, fstrSOLID%TEMP_ngrp_tot
          ig= fstrSOLID%TEMP_ngrp_ID(ig0)
          val=fstrSOLID%TEMP_ngrp_val(ig0)
!C START & END
          iS0= hecMESH%node_group%grp_index(ig-1) + 1
          iE0= hecMESH%node_group%grp_index(ig  )
           do ik= iS0, iE0
             in   = hecMESH%node_group%grp_item(ik)
             temp( in ) = val
          enddo
        enddo
      endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        if( ic_type==3422 ) ic_type=342
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            zz(j)=hecMESH%node(3*nodLOCAL(j)  )
            tt(j)=temp( nodLOCAL(j) )
            tt0(j)=ref_temp
            edisp(3*j-2)=fstrSOLID%unode(3*nodLOCAL(j)-2)
            edisp(3*j-1)=fstrSOLID%unode(3*nodLOCAL(j)-1)
            edisp(3*j  )=fstrSOLID%unode(3*nodLOCAL(j)  )
          enddo
!C** Create local stiffness
          if( fstrSOLID%elements(icel)%gausses(1)%pMaterial%nlgeom_flag ==0 ) then
          if( ic_type==361 ) then
            call UpdateST_C3D8IC(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn),        &
                           edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          else if( ic_type==301 ) then
            call UpdateST_C1( ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn), area,        &
                           edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          else
            call UpdateST_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn),            &
                           edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          endif
          endif
          if( ic_type==301 ) then
            call NodalStress_C1(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:))
          else
            call NodalStress_C3(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:))
          endif
!C** elem ID
!!!          ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area==hecMESH%my_rank ) then
            if( ic_type==301 ) then
              call ElementStress_C1(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            else
              call ElementStress_C3(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            endif
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-1) = estress
            s11=fstrSOLID%ESTRESS(7*ielem-6)
            s22=fstrSOLID%ESTRESS(7*ielem-5)
            s33=fstrSOLID%ESTRESS(7*ielem-4)
            s12=fstrSOLID%ESTRESS(7*ielem-3)
            s23=fstrSOLID%ESTRESS(7*ielem-2)
            s13=fstrSOLID%ESTRESS(7*ielem-1)
            ps=(s11+s22+s33)/3.0
            smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                  +s12**2+s23**2+s13**2
            fstrSOLID%ESTRESS(7*ielem)=dsqrt(3.d0*smises)
          endif

          do j=1,nn
            if( ic_type==301 ) then
              tstrain( nodLocal(j),: ) = tstrain(nodLOCAL(j),:) + edstrain(j,:)
              tstress( nodLocal(j),1:6 ) = tstress(nodLOCAL(j),1:6) + edstress(j,:)
              tnumber( nodLOCAL(j) )=tnumber( nodLOCAL(j) )+1
            else
              ndstrain( nodLocal(j),: ) = ndstrain(nodLOCAL(j),:) + edstrain(j,:)
              ndstress( nodLocal(j),1:6 ) = ndstress(nodLOCAL(j),1:6) + edstress(j,:)
              nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1
            endif
          enddo
        enddo
      enddo
!** Average over nodes
      do i=1,hecMESH%n_node
        if( tnumber(i).gt.0 ) then
          ndstrain(i,:)=ndstrain(i,:)/nnumber(i) + tstrain(i,:)/tnumber(i)
          ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i) + tstress(i,1:6)/tnumber(i)
        else
          ndstrain(i,:)=ndstrain(i,:)/nnumber(i)
          ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i)
        endif
      enddo
!** CALCULATE von MISES stress
      do i=1,hecMESH%n_node
        s11=ndstress(i,1)
        s22=ndstress(i,2)
        s33=ndstress(i,3)
        s12=ndstress(i,4)
        s23=ndstress(i,5)
        s13=ndstress(i,6)
        ps=(s11+s22+s33)/3.0
        smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
              +s12**2+s23**2+s13**2
        ndstress(i,7)=sqrt(3.d0*smises)
      enddo
!** Set Array 
      do i=1,hecMESH%n_node
        do j=1,6
          fstrSOLID%STRAIN(6*(i-1)+j)=ndstrain(i,j)
        enddo
        do j=1,7
          fstrSOLID%STRESS(7*(i-1)+j)=ndstress(i,j)
        enddo
      enddo
!** Show
      if( .not. present(ifwrite) ) then
      write(ILOG,*) '#### STRAIN'
      write(ILOG,'(a,a)')                                              &
                    '     NODE     E11         E22         E33     '   &
                             ,'    E12         E23         E13     '
      write(ILOG,'(a,a)')                                              &
                    '  -------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p6e12.4)') jj, (ndSTRAIN(ii,k),k=1,6)
      enddo
      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                              &
                    '     NODE     S11         S22         S33     '   &
                 ,'    S12         S23         S13        MISES'
      write(ILOG,'(a,a)')                                              &
                    '  -------+-----------+-----------+-----------+'   &
                 ,'-----------+-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p7e12.4)') jj, (ndSTRESS(ii,k),k=1,7)
      enddo
      endif
!*Deallocate array
      deallocate( nnumber,temp,ndstrain,ndstress)
      deallocate( tnumber,tstrain,tstress)
      end subroutine fstr_nodal_stress_3d
!C
!C REACTION FORCE 3D
!C
      subroutine fstr_reaction_force_3d(hecMESH,fstrSOLID)
      use m_fstr
      use m_static_lib
      use m_static_LIB_1d
      use m_static_LIB_3dIC
      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
!** Local variables
      REAL(kind=kreal) stiff(120,120), thick
      REAL(kind=kreal) edisp(60),force(60),ecoord(3,20)
      integer(kind=kint) nodLOCAL(20), isect, ihead
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
!** Freedom
      data NDOF/3/
!C Allocate array
      allocate ( id_spc( hecMESH%n_node ) )
      do i=1,hecMESH%n_node
        id_spc(i)=0
      enddo
!C
!C  SPC Nodal Check Flag
!C 
      icount=0
      do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          if( id_spc(in) .eq. 0 ) then
            icount=icount+1
            id_spc(in)=icount
          endif
        enddo
      enddo

      icount=0
      do i=1,hecMESH%n_node
        if( id_spc(i).ge.1 ) icount=icount+1
      enddo
      nspc=icount
!C Allocate array
      if( nspc .gt. 0 ) then
         allocate ( spcForce( nspc * NDOF ) )
         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc )
        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            ecoord(1,j)=hecMESH%node(3*nodLOCAL(j)-2)
            ecoord(2,j)=hecMESH%node(3*nodLOCAL(j)-1)
            ecoord(3,j)=hecMESH%node(3*nodLOCAL(j)  )
            edisp(3*j-2)=fstrSOLID%unode(3*nodLOCAL(j)-2)
            edisp(3*j-1)=fstrSOLID%unode(3*nodLOCAL(j)-1)
            edisp(3*j  )=fstrSOLID%unode(3*nodLOCAL(j)  )
          enddo
 !C** check node id
          iflag=0
          do j=1,nn
            if( id_spc(  nodLOCAL(j) ) .ge. 1 ) iflag=1
          enddo
          if( iflag == 1 ) then
!C** Create local stiffness
            if (ic_type==361) then   
              call STF_C3D8IC( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3))
            else if( ic_type==301 ) then
              isect= hecMESH%section_ID(icel)
              ihead = hecMESH%section%sect_R_index(isect-1)
              thick = hecMESH%section%sect_R_item(ihead+1)
              call STF_C1( ic_type,nn,ecoord(:,1:nn),thick,fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3) )
            else
              call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3), 1.d0)
            endif
            force(1:nn*3)= matmul( stiff(1:nn*3,1:nn*3), edisp(1:nn*3) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc(  nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k)       &
                  +force(NDOF*(j-1)+k)      
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show
      write(ILOG,'(a)')'    NODE     X-REAC      Y-REAC      Z-REAC'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p3e12.4)') j,     & 
          spcForce(NDOF*(nid-1)+1),           &                     
          spcForce(NDOF*(nid-1)+2),           &              
          spcForce(NDOF*(nid-1)+3)
        endif
      enddo
!C Deallocate array
      deallocate ( spcForce )
      deallocate ( id_spc )
      end subroutine fstr_reaction_force_3d
!C***
!C*** OUTPUT for FSTR solver
!C***
!C
!C 2D SOLID
!C
      subroutine FSTR_OUTPUT_2D(hecMESH, hecMAT, fstrSOLID, fstrPARAM)
      use m_fstr
   
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)      :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID

!C
!C-- SET DISPLACEMENT
!C
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(2*i-1)= hecMAT%X(2*i-1)
        fstrSOLID%unode(2*i  )= hecMAT%X(2*i  )
      enddo
!C
!C-- DISPLACEMENT
!C
      write(ILOG,*) '#### DISPLACEMENT 2D'
      write(ILOG,'(a)')                                      &
                    '     NODE    X-DISP      Y-DISP   '
      write(ILOG,'(a)')                                      &
                    ' --------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,2e12.4)') jj,(fstrSOLID%unode(2*(ii-1)+k),k=1,2)
      enddo
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 2D'
      call fstr_reaction_force_2d( hecMESH, fstrSOLID)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 2D'
      call fstr_nodal_stress_2d( hecMESH, fstrSOLID, fstrPARAM)

      end subroutine FSTR_OUTPUT_2D

!C
!C NODAL STRESS 2D
!C
      subroutine fstr_nodal_stress_2d(hecMESH,fstrSOLID,fstrPARAM)
      use m_fstr
      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
!** Local variables
      REAL(kind=kreal) xx(8),yy(8)
      REAL(kind=kreal) edisp(16),estrain(4),estress(4)
      REAL(kind=kreal) tt(8),tt0(8)
      integer(kind=kint) nodLOCAL(8)
      REAL(kind=kreal) edstrain(9,4),edstress(9,4)
!** Array for nodal recovery
      real    (kind=KREAL), dimension(:,:), allocatable :: arrayTotal
      integer (kind=KINT), dimension(:), allocatable :: nnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp
!*Allocate array
      allocate ( arrayTotal( hecMESH%n_node, 9  ) ) 
      allocate ( nnumber(hecMESH%n_node) )
!*ZERO CLEAR 
      arrayTotal=0.d0
      nnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node) )
      temp=0
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
        do ig0= 1, fstrSOLID%TEMP_ngrp_tot
          ig= fstrSOLID%TEMP_ngrp_ID(ig0)
          val=fstrSOLID%TEMP_ngrp_val(ig0)
!C START & END
          iS0= hecMESH%node_group%grp_index(ig-1) + 1
          iE0= hecMESH%node_group%grp_index(ig  )
           do ik= iS0, iE0
             in   = hecMESH%node_group%grp_item(ik)
             temp( in ) = val
          enddo
        enddo
      endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            tt(j)=temp( nodLOCAL(j) )
            tt0(j)=ref_temp
            edisp(2*j-1)=fstrSOLID%unode(2*nodLOCAL(j)-1)
            edisp(2*j  )=fstrSOLID%unode(2*nodLOCAL(j)  )
          enddo
!C** section  ID
          pa1=1.d0     
!C** Create local stiffness
          call UpdateST_C2(ic_type,nn,xx(1:nn),yy(1:nn),tt(1:nn),tt0(1:nn),pa1,              &
                           fstrSOLID%elements(icel)%iset,edisp(1:nn*2),fstrSOLID%elements(icel)%gausses )
          call NodalStress_C2(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:) )
!C** elem ID
!!!          ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area==hecMESH%my_rank ) then
            call ElementStress_C2(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem-2) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-3) = estress
            s11=fstrSOLID%ESTRESS(7*ielem-6)
            s22=fstrSOLID%ESTRESS(7*ielem-5)
            s33=fstrSOLID%ESTRESS(7*ielem-4)
            s12=fstrSOLID%ESTRESS(7*ielem-3)
            s23=0.0D0
            s13=0.0D0
            ps=(s11+s22+s33)/3.0
            smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )         &
                +s12**2+s23**2+s13**2
            fstrSOLID%ESTRESS(7*ielem)=sqrt(3.0*smises)
          endif

          do j=1,nn
            do k=1,4
              arrayTotal(nodLOCAL(j),k)=                                 &
              arrayTotal(nodLOCAL(j),k)+edstrain(j,k)
            enddo
            do k=1,4
              arrayTotal(nodLOCAL(j),k+4)=                                 &
              arrayTotal(nodLOCAL(j),k+4)+edstress(j,k)
            enddo
            nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1

          enddo
        enddo
      enddo
!** Average over nodes
      do i=1,hecMESH%n_node
        do j=1,8
          arrayTotal(i,j)=arrayTotal(i,j)/FLOAT(nnumber(i))
        enddo
      enddo
!** RE-CALCULATE von MISES stress
      do i=1,hecMESH%n_node
        s11=arrayTotal(i,5)
        s22=arrayTotal(i,6)
        s33=arrayTotal(i,8)
        s12=arrayTotal(i,7)
        s23=0.0
        s13=0.0
        ps=(s11+s22+s33)/3.0
        smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )      &  
              +s12**2+s23**2+s13**2
        arrayTotal(i,9)=sqrt(3.0*smises)
      enddo
!** Set Array 
      do i=1,hecMESH%n_node
        fstrSOLID%STRAIN(6*(i-1)+1)=arrayTotal(i,1)
        fstrSOLID%STRAIN(6*(i-1)+2)=arrayTotal(i,2)
        fstrSOLID%STRAIN(6*(i-1)+3)=arrayTotal(i,4)
        fstrSOLID%STRAIN(6*(i-1)+4)=arrayTotal(i,3)
        fstrSOLID%STRAIN(6*(i-1)+5)=0.0
        fstrSOLID%STRAIN(6*(i-1)+6)=0.0
        fstrSOLID%STRESS(7*(i-1)+1)=arrayTotal(i,5)
        fstrSOLID%STRESS(7*(i-1)+2)=arrayTotal(i,6)
        fstrSOLID%STRESS(7*(i-1)+3)=arrayTotal(i,8)
        fstrSOLID%STRESS(7*(i-1)+4)=arrayTotal(i,7)
        fstrSOLID%STRESS(7*(i-1)+5)=0.0
        fstrSOLID%STRESS(7*(i-1)+6)=0.0
        fstrSOLID%STRESS(7*(i-1)+7)=arrayTotal(i,9)
      enddo
!** Show
      write(ILOG,*) '#### STRIN'
      write(ILOG,'(a,a)')                                              &
               '   NODE       E11         E22         E33     '        &
                             ,'    E12' 
      write(ILOG,'(a,a)')                                              &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p4e12.4)') jj,                               &
                                 (fstrSOLID%STRAIN(6*(ii-1)+k),k=1,4)
      enddo
      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                              &
               '   NODE       S11         S22         S33     '        &
                             ,'    S12        MISES' 
      write(ILOG,'(a,a)')                                              &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p5e12.4)') jj,                               &
                                 (fstrSOLID%STRESS(7*(ii-1)+k),k=1,4)  &
                                 ,fstrSOLID%STRESS(7*(ii-1)+7)
      enddo
!*Deallocate array
      deallocate( arrayTotal,nnumber,temp)
      end subroutine fstr_nodal_stress_2d

!C
!C REACTION FORCE 2D
!C
      subroutine fstr_reaction_force_2d(hecMESH,fstrSOLID)
      use m_fstr
      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
!** Local variables
      REAL(kind=kreal) ecoord(2,8),stiff(18,18)
      REAL(kind=kreal) edisp(16),force(16)
      integer(kind=kint) nodLOCAL(8)
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
!** Freedom
      data NDOF/2/
!C Allocate array
      allocate ( id_spc( hecMESH%n_node) )
      do i=1,hecMESH%n_node
        id_spc(i)=0
      enddo
!C
!C  SPC Nodal Check Flag
!C 
      icount=0
      do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          if( id_spc(in) .eq. 0 ) then
            icount=icount+1
            id_spc(in)=icount
          endif
        enddo
      enddo

      icount=0
      do i=1,hecMESH%n_node
        if( id_spc(i).ge.1 ) icount=icount+1
      enddo
      nspc=icount
!C Allocate array
      if( nspc .gt. 0 ) then
         allocate ( spcForce( nspc * NDOF ) )
         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc )
        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            ecoord(1,j)=hecMESH%node(3*nodLOCAL(j)-2)
            ecoord(2,j)=hecMESH%node(3*nodLOCAL(j)-1)
            edisp(2*j-1)=fstrSOLID%unode(2*nodLOCAL(j)-1)
            edisp(2*j  )=fstrSOLID%unode(2*nodLOCAL(j)  )
          enddo
!C** check node id
          iflag=0
          do j=1,nn
            if( id_spc(  nodLOCAL(j) ) .ge. 1 ) iflag=1
          enddo
          if( iflag == 1 ) then
            pa1=1.d0
!C** Create local stiffness
            call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses,pa1    &
                        ,stiff(1:nn*2,1:nn*2),fstrSOLID%elements(icel)%iset)
            force(1:nn*2) = matmul( stiff(1:nn*2,1:nn*2), edisp(1:nn*2) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc( nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k)   &
                  +force(NDOF*(j-1)+k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show
      write(ILOG,'(a)') '   NODE    X-REACTION  Y-REACTION'
      write(ILOG,'(a)') ' --------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p2e12.4)') j,spcForce(NDOF*(nid-1)+1),   & 
                                   spcForce(NDOF*(nid-1)+2) 
        endif
      enddo
!C Deallocate array
      deallocate ( spcForce )
      deallocate ( id_spc )
      end subroutine fstr_reaction_force_2d
!C***
!C*** OUTPUT for FSTR solver
!C***
      subroutine fstr_output_6d( hecMESH,hecMAT,fstrSOLID,fstrPARAM )
      use m_fstr

      type ( hecmwST_matrix     ) :: hecMAT
      type ( hecmwST_local_mesh ) :: hecMESH
      type ( fstr_solid         ) :: fstrSOLID
      type ( fstr_param         ) :: fstrPARAM
      INTEGER(kind=kint) NDOF

      NDOF = hecMESH%n_dof
!C
!C-- DISPLACEMENT
  
      do i = 1, hecMESH%n_node
        do j = 1, NDOF
          fstrSOLID%unode(NDOF*(i-1)+j)= hecMAT%X(NDOF*(i-1)+j)
        enddo
      enddo
!C
      write(ILOG,*) '#### DISPLACEMENT SHELL'
      write(ILOG,'(a,a)')                                            & 
                    '     NODE    X-DISP      Y-DISP      Z-DISP   ' &
                             ,'   X-ROT       Y-ROT       Z-ROT'
      write(ILOG,'(a,a)')                                            & 
                    ' --------+-----------+-----------+-----------+' &
                             ,'-----------+-----------+-------------'
      do i = 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,6e12.4)') jj,(fstrSOLID%unode(6*(ii-1)+k),k=1,6)
      enddo

      call flush(ILOG)

!C
!C-- REACTION FORCE
  
      write(ILOG,*) '#### REACTION FORCE SHELL'
      call fstr_reaction_force_6d( hecMESH, fstrSOLID)

      call flush(ILOG)
!C
!C-- NODAL STRESS
  
      write(ILOG,*) '#### NODAL STRAIN/STRESS SHELL'
      call fstr_nodal_stress_6d( hecMESH, fstrSOLID,fstrPARAM)

      call flush(ILOG)

      end subroutine fstr_output_6d



!C
!C NODAL STRESS SHELL
!C
      subroutine fstr_nodal_stress_6d(hecMESH,fstrSOLID,fstrPARAM)
      use m_fstr

      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (fstr_param)      :: fstrPARAM

!** Local variables
      REAL(kind=kreal) xx(8),yy(8),zz(8),ri(4),si(4)
      REAL(kind=kreal) edisp(48), tt(8), tt0(8)
      integer(kind=kint) nodLOCAL(8)
      REAL(kind=kreal) array(11)

!** Array for nodal recovery
      real    (kind=KREAL),dimension(:,:),allocatable::arrayTotal_U
      real    (kind=KREAL),dimension(:,:),allocatable::arrayTotal_D
      integer (kind=KINT), dimension(:),  allocatable::nnumber
      real    (kind=KREAL),dimension(:),  allocatable::temp
      real    (kind=kreal) thick

      data ri / -1.0,  1.0, 1.0, -1.0 /
      data si / -1.0, -1.0, 1.0,  1.0 /

!*Allocate array
      allocate ( arrayTotal_U(hecMESH%n_node,11) )
      allocate ( arrayTotal_D(hecMESH%n_node,11) )
      allocate ( nnumber(hecMESH%n_node) )
 
      arrayTotal_U = 0.0d0
      arrayTotal_D = 0.0d0
      nnumber = 0
!C
!C Set Temperature

      allocate ( temp(hecMESH%n_node) )
      temp = 0.0D0
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
        do ig0 = 1, fstrSOLID%TEMP_ngrp_tot
          ig  = fstrSOLID%TEMP_ngrp_ID(ig0)
          val = fstrSOLID%TEMP_ngrp_val(ig0)
          iS  = hecMESH%node_group%grp_index(ig-1) + 1
          iE  = hecMESH%node_group%grp_index(ig  )
          do ik = iS, iE
            in = hecMESH%node_group%grp_item(ik)
            temp(in) = val
          enddo
        enddo
      endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+

      do itype  = 1, hecMESH%n_elem_type
        iS      = hecMESH%elem_type_index(itype-1) + 1
        iE      = hecMESH%elem_type_index(itype  )
        ic_type = hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_shell(ic_type)) cycle
!C        
        nn = hecmw_get_max_node(ic_type)
!C
        do icel = iS, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
            xx(j)  = hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)  = hecMESH%node(3*nodLOCAL(j)-1)
            zz(j)  = hecMESH%node(3*nodLOCAL(j))
            tt(j)  = temp( nodLOCAL(j) )
            tt0(j) = ref_temp

            edisp(6*j-5) = fstrSOLID%unode(6*nodLOCAL(j)-5)
            edisp(6*j-4) = fstrSOLID%unode(6*nodLOCAL(j)-4)
            edisp(6*j-3) = fstrSOLID%unode(6*nodLOCAL(j)-3)
            edisp(6*j-2) = fstrSOLID%unode(6*nodLOCAL(j)-2)
            edisp(6*j-1) = fstrSOLID%unode(6*nodLOCAL(j)-1)
            edisp(6*j  ) = fstrSOLID%unode(6*nodLOCAL(j))
          enddo

!C** Get Properties
          isect = hecMESH%section_ID(icel)
          call FSTR_GET_PROP( hecMESH,isect,ee,pp,rho,alpha,thick )
!C** Create local stiffness             
!C** elem ID
!!!          ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            fstrSOLID%ESTRAIN(10*ielem-9) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-8) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-7) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-6) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-5) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-4) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-3) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-2) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem-1) = 0.0D0
            fstrSOLID%ESTRAIN(10*ielem  ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-11) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-10) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-9 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-8 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-7 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-6 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-5 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-4 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-3 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-2 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem-1 ) = 0.0D0
            fstrSOLID%ESTRESS(12*ielem   ) = 0.0D0
          endif

          if     (ic_type.eq.731) then
            ti = 1.0
            call RCV_S3( xx,yy,zz,ee,pp,thick,edisp,ti,array )
            do j = 1, 3
              jj = nodLOCAL(j)
              do k = 1, 11
                arrayTotal_U(jj,k) = arrayTotal_U(jj,k) + array(k)
              enddo
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              fstrSOLID%ESTRAIN(10*ielem-9 ) = array(1)/nn
              fstrSOLID%ESTRAIN(10*ielem-8 ) = array(2)/nn
              fstrSOLID%ESTRAIN(10*ielem-7 ) = array(3)/nn
              fstrSOLID%ESTRAIN(10*ielem-6 ) = array(4)/nn
              fstrSOLID%ESTRAIN(10*ielem-5 ) = array(5)/nn
              fstrSOLID%ESTRESS(12*ielem-11) = array(6)/nn
              fstrSOLID%ESTRESS(12*ielem-10) = array(7)/nn
              fstrSOLID%ESTRESS(12*ielem-9 ) = array(8)/nn
              fstrSOLID%ESTRESS(12*ielem-8 ) = array(9)/nn
              fstrSOLID%ESTRESS(12*ielem-7 ) = array(10)/nn
              s11=fstrSOLID%ESTRESS(12*ielem-11)
              s22=fstrSOLID%ESTRESS(12*ielem-10)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-9 )
              s23=fstrSOLID%ESTRESS(12*ielem-8 )
              s13=fstrSOLID%ESTRESS(12*ielem-7 )
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )   &   
                                      +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem-6)=sqrt(3.0*smises)
            endif

            ti = -1.0
            call RCV_S3( xx,yy,zz,ee,pp,thick,edisp,ti,array )
            do j = 1, 3
              jj = nodLOCAL(j)
              do k = 1, 11
                arrayTotal_D(jj,k) = arrayTotal_D(jj,k) + array(k)
              enddo
              nnumber(jj) = nnumber(jj) + 1
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              fstrSOLID%ESTRAIN(10*ielem-4) = array(1)/nn
              fstrSOLID%ESTRAIN(10*ielem-3) = array(2)/nn
              fstrSOLID%ESTRAIN(10*ielem-2) = array(3)/nn
              fstrSOLID%ESTRAIN(10*ielem-1) = array(4)/nn
              fstrSOLID%ESTRAIN(10*ielem  ) = array(5)/nn
              fstrSOLID%ESTRESS(12*ielem-5) = array(6)/nn
              fstrSOLID%ESTRESS(12*ielem-4) = array(7)/nn
              fstrSOLID%ESTRESS(12*ielem-3) = array(8)/nn
              fstrSOLID%ESTRESS(12*ielem-2) = array(9)/nn
              fstrSOLID%ESTRESS(12*ielem-1) = array(10)/nn
              s11=fstrSOLID%ESTRESS(12*ielem-5)
              s22=fstrSOLID%ESTRESS(12*ielem-4)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-3)
              s23=fstrSOLID%ESTRESS(12*ielem-2)
              s13=fstrSOLID%ESTRESS(12*ielem-1)
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )  &   
                                    +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem)=sqrt(3.0*smises)
            endif

          elseif (ic_type.eq.741) then
            do j = 1, 4
              jj = nodLOCAL(j)
              ti = 1.0
              call RCV_S4( xx,yy,zz,ee,pp,thick,                  &
                                  edisp,ri(j),si(j),ti,array )
              do k = 1, 11
                arrayTotal_U(jj,k) = arrayTotal_U(jj,k) + array(k)
              enddo

              if( ID_area.eq.hecMESH%my_rank ) then
                fstrSOLID%ESTRAIN(10*ielem-9 ) =                  & 
                fstrSOLID%ESTRAIN(10*ielem-9 ) + array(1)/nn
                fstrSOLID%ESTRAIN(10*ielem-8 ) =                  & 
                fstrSOLID%ESTRAIN(10*ielem-8 ) + array(2)/nn
                fstrSOLID%ESTRAIN(10*ielem-7 ) =                  & 
                fstrSOLID%ESTRAIN(10*ielem-7 ) + array(3)/nn
                fstrSOLID%ESTRAIN(10*ielem-6 ) =                  & 
                fstrSOLID%ESTRAIN(10*ielem-6 ) + array(4)/nn
                fstrSOLID%ESTRAIN(10*ielem-5 ) =                  & 
                fstrSOLID%ESTRAIN(10*ielem-5 ) + array(5)/nn
                fstrSOLID%ESTRESS(12*ielem-11) =                  & 
                fstrSOLID%ESTRESS(12*ielem-11) + array(6)/nn
                fstrSOLID%ESTRESS(12*ielem-10) =                  & 
                fstrSOLID%ESTRESS(12*ielem-10) + array(7)/nn
                fstrSOLID%ESTRESS(12*ielem-9 ) =                  & 
                fstrSOLID%ESTRESS(12*ielem-9 ) + array(8)/nn
                fstrSOLID%ESTRESS(12*ielem-8 ) =                  & 
                fstrSOLID%ESTRESS(12*ielem-8 ) + array(9)/nn
                fstrSOLID%ESTRESS(12*ielem-7 ) =                  & 
                fstrSOLID%ESTRESS(12*ielem-7 ) + array(10)/nn
              endif
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              s11=fstrSOLID%ESTRESS(12*ielem-11)
              s22=fstrSOLID%ESTRESS(12*ielem-10)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-9 )
              s23=fstrSOLID%ESTRESS(12*ielem-8 )
              s13=fstrSOLID%ESTRESS(12*ielem-7 )
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &   
                                           +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem-6)=sqrt(3.0*smises)
            endif

            do j = 1, 4
              jj = nodLOCAL(j)
              ti = -1.0
              call RCV_S4( xx,yy,zz,ee,pp,thick,                    &
                                  edisp,ri(j),si(j),ti,array )
              do k = 1, 11
                arrayTotal_D(jj,k) = arrayTotal_D(jj,k) + array(k)
              enddo
              nnumber(jj) = nnumber(jj) + 1

              if( ID_area.eq.hecMESH%my_rank ) then
                fstrSOLID%ESTRAIN(10*ielem-4) =                & 
                fstrSOLID%ESTRAIN(10*ielem-4) + array(1)/nn
                fstrSOLID%ESTRAIN(10*ielem-3) =                &
                fstrSOLID%ESTRAIN(10*ielem-3) + array(2)/nn
                fstrSOLID%ESTRAIN(10*ielem-2) =                &
                fstrSOLID%ESTRAIN(10*ielem-2) + array(3)/nn
                fstrSOLID%ESTRAIN(10*ielem-1) =                &
                fstrSOLID%ESTRAIN(10*ielem-1) + array(4)/nn
                fstrSOLID%ESTRAIN(10*ielem  ) =                &
                fstrSOLID%ESTRAIN(10*ielem  ) + array(5)/nn
                fstrSOLID%ESTRESS(12*ielem-5) =                &
                fstrSOLID%ESTRESS(12*ielem-5) + array(6)/nn
                fstrSOLID%ESTRESS(12*ielem-4) =                &
                fstrSOLID%ESTRESS(12*ielem-4) + array(7)/nn
                fstrSOLID%ESTRESS(12*ielem-3) =                &
                fstrSOLID%ESTRESS(12*ielem-3) + array(8)/nn
                fstrSOLID%ESTRESS(12*ielem-2) =                &
                fstrSOLID%ESTRESS(12*ielem-2) + array(9)/nn
                fstrSOLID%ESTRESS(12*ielem-1) =                &
                fstrSOLID%ESTRESS(12*ielem-1) + array(10)/nn
              endif
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              s11=fstrSOLID%ESTRESS(12*ielem-5)
              s22=fstrSOLID%ESTRESS(12*ielem-4)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-3)
              s23=fstrSOLID%ESTRESS(12*ielem-2)
              s13=fstrSOLID%ESTRESS(12*ielem-1)
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )  &   
                                    +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem)=sqrt(3.0*smises)
            endif

          endif
        enddo
      enddo

!** Average over nodes

      do i = 1, hecMESH%n_node

        do j = 1, 11
          arrayTotal_U(i,j) = arrayTotal_U(i,j) / FLOAT(nnumber(i))
          arrayTotal_D(i,j) = arrayTotal_D(i,j) / FLOAT(nnumber(i))

          if( j < 6 ) then 
            fstrSOLID%STRAIN(10*(i-1)+j  ) = arrayTotal_U(i,j)        
            fstrSOLID%STRAIN(10*(i-1)+j+5) = arrayTotal_D(i,j)        
          else
            fstrSOLID%STRESS(12*(i-1)+j-5) = arrayTotal_U(i,j)        
            fstrSOLID%STRESS(12*(i-1)+j+1) = arrayTotal_D(i,j)        
          endif

        enddo

!** RE-CALCULATE von MISES stress

        s1 = fstrSOLID%STRESS(12*(i-1)+1)
        s2 = fstrSOLID%STRESS(12*(i-1)+2)
        s3 = fstrSOLID%STRESS(12*(i-1)+3)
        s4 = fstrSOLID%STRESS(12*(i-1)+4)
        s5 = fstrSOLID%STRESS(12*(i-1)+5)
        fstrSOLID%STRESS(12*(i-1)+6)                                  & 
           = DSQRT( s1*s1 + s2*s2 - s1*s2 + 3.0d0*s3*s3 )

        s1 = fstrSOLID%STRESS(12*(i-1)+7)
        s2 = fstrSOLID%STRESS(12*(i-1)+8)
        s3 = fstrSOLID%STRESS(12*(i-1)+9)
        s4 = fstrSOLID%STRESS(12*(i-1)+10)
        s5 = fstrSOLID%STRESS(12*(i-1)+11)
        fstrSOLID%STRESS(12*i)                                        & 
           = DSQRT( s1*s1 + s2*s2 - s1*s2 + 3.0d0*s3*s3 )

      enddo

!** Show
      write(ILOG,*) '#### STRAIN'
      write(ILOG,'(a,a)')                                             & 
               '     NODE          E11         E22         E12     '  &
                             ,'    E23         E13'
      write(ILOG,'(a,a)')                                             & 
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,a5,1p5e12.4)') jj,' (+) ',                     &
                                 (fstrSOLID%STRAIN(10*(ii-1)+k),k=1,5)
        write(ILOG,'(10x,a5,1p5e12.4)')   ' (-) ',                      &
                                 (fstrSOLID%STRAIN(10*(ii-1)+k),k=6,10)
      enddo
      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                               &
               '     NODE          S11         S22         S12     '    &
                             ,'    S23         S13        MISES'
      write(ILOG,'(a,a)')                                               &
               '  ------------+-----------+-----------+-----------+'    &
                             ,'-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,a5,1p6e12.4)') jj,' (+) ',                     &
                                 (fstrSOLID%STRESS(12*(ii-1)+k),k=1,6)
        write(ILOG,'(10x,a5,1p6e12.4)')   ' (-) ',                      &
                                 (fstrSOLID%STRESS(12*(ii-1)+k),k=7,12)
      enddo

      deallocate( arrayTotal_U,arrayTotal_D,nnumber,temp)

      end subroutine fstr_nodal_stress_6d


!C
!C REACTION FORCE SHELL
!C
      subroutine fstr_reaction_force_6d(hecMESH,fstrSOLID)
      use m_fstr

      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID

!** Local variables
      REAL(kind=kreal)    xx(8),yy(8),zz(8)
      REAL(kind=kreal)    edisp(48),force(48),local_stf(1176)
      integer(kind=kint) nodLOCAL(8)

      real   (kind=KREAL),dimension(:), allocatable :: spcForce
      integer(kind=KINT), dimension(:), allocatable :: id_spc

!C Allocate array
      allocate ( id_spc( hecMESH%n_node ) )
      id_spc = 0

      NDOF = hecMESH%n_dof
!C
!C  SPC Nodal Check Flag
!C 
      icount = 0
      do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot

        ig = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS = hecMESH%node_group%grp_index(ig-1) + 1
        iE = hecMESH%node_group%grp_index(ig  )

        do ik = iS, iE
          in = hecMESH%node_group%grp_item(ik)
          if( id_spc(in) .eq. 0 ) then
            icount = icount + 1
            id_spc(in) = icount
          endif
        enddo

      enddo

      icount = 0
      do i = 1, hecMESH%n_node
        if( id_spc(i).ge.1 ) icount = icount + 1
      enddo
      nspc = icount

!C Allocate array
      if( nspc .gt. 0 ) then
         allocate ( spcForce(nspc*NDOF) )
         spcForce = 0.0
      else
        deallocate ( id_spc )
        return
      endif

!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+

      do itype = 1, hecMESH%n_elem_type
        iS = hecMESH%elem_type_index(itype-1) + 1
        iE = hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_shell(ic_type)) cycle

!C** Set number of nodes        
        nn = hecmw_get_max_node(ic_type)

        do icel = iS, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
            xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
            yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
            zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
            edisp( 6*j-5 ) = fstrSOLID%unode( 6*nodLOCAL(j)-5 )
            edisp( 6*j-4 ) = fstrSOLID%unode( 6*nodLOCAL(j)-4 )
            edisp( 6*j-3 ) = fstrSOLID%unode( 6*nodLOCAL(j)-3 )
            edisp( 6*j-2 ) = fstrSOLID%unode( 6*nodLOCAL(j)-2 )
            edisp( 6*j-1 ) = fstrSOLID%unode( 6*nodLOCAL(j)-1 )
            edisp( 6*j   ) = fstrSOLID%unode( 6*nodLOCAL(j)   )
          enddo

          ifspc = 0
          do j = 1, nn
            if( id_spc( nodLOCAL(j) ) .ge. 1 ) ifspc = 1
          enddo

          if( ifspc .eq. 1 ) then
!C** Get Properties
            isect = hecMESH%section_ID(icel)
            call FSTR_GET_PROP(hecMESH,isect,ee,pp,rho,alpha,thick)

!C** Create local stiffness             
            if     ( ic_type.EQ.731 ) then
              call STF_S3 ( xx,yy,zz,ee,pp,thick,local_stf )
            elseif ( ic_type.eq.741 ) then
              call STF_S4 ( xx,yy,zz,ee,pp,thick,local_stf )
            endif

!C*** [K]{U} ****
            force = 0.0
            do jj = 1, nn*NDOF
              do ii = 1, jj-1
                num = ( (jj-1)*jj + 2*ii ) / 2
                force(ii) = force(ii) + local_stf(num)*edisp(jj)
                force(jj) = force(jj) + local_stf(num)*edisp(ii)
              enddo
              num = ( (jj-1)*jj + 2*jj ) / 2
              force(jj) = force(jj) + local_stf(num)*edisp(jj)
            enddo

!C*** Add SPC FORCE ****
            do j = 1, nn
              isid = id_spc( nodLOCAL(j) )
              if( isid .gt. 0 ) then
                do k = 1, NDOF
                  kk = NDOF*(isid-1) + k
                  spcForce(kk) = spcForce(kk) + force( NDOF*(j-1)+k )
                enddo
              endif
            enddo

          endif

        enddo
      enddo

!C*** Show
      write(ILOG,'(a,a)')                                              &
                    '   NODE      X-REAC      Y-REAC      Z-REAC   '   &
                             ,'  RX-REAC     RY-REAC     RZ-REAC'
      write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i = 1, hecMESH%nn_internal
        in = hecMESH%global_node_ID(i)
        isid = id_spc(i)
        if( isid .gt. 0 ) then
          write(ILOG,'(i10,1p6e12.4)') in                              &
          , spcForce( NDOF*(isid-1)+1 ), spcForce( NDOF*(isid-1)+2 )   &
          , spcForce( NDOF*(isid-1)+3 ), spcForce( NDOF*(isid-1)+4 )   &
          , spcForce( NDOF*(isid-1)+5 ), spcForce( NDOF*(isid-1)+6 )
        endif
      enddo

      deallocate ( id_spc,spcForce )

      end subroutine fstr_reaction_force_6d
!C 
!C OUTPUT FOR ELEMENTS
!C
      subroutine fstr_output_elem(hecMESH,fstrSOLID)
      use m_fstr
      implicit REAL(kind=kreal) (A-H,O-Z)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID

      if( hecMESH%n_dof .eq. 3 ) then
        write(ILOG,*) '#### 3D STRAIN @Element'
        write(ILOG,'(a,a)')                                            &
                    '  ELEMENT     E11         E22         E33     '   &
                             ,'    E12         E23         E13     '
        write(ILOG,'(a,a)')                                            &
                    '  -------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p6e12.4)') jj,                           &
                                 (fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,6)
          endif
        enddo
        write(ILOG,*) '#### 3D STRESS @Element'
        write(ILOG,'(a,a)')                                            &
                    '  ELEMENT     S11         S22         S33     '   &
                 ,'    S12         S23         S13        MISES'
        write(ILOG,'(a,a)')                                            &
                    '  -------+-----------+-----------+-----------+'   &
                 ,'-----------+-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p7e12.4)') jj,                            &
                                 (fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,7)
          endif
        enddo

      elseif( hecMESH%n_dof .eq. 2 ) then
        write(ILOG,*) '#### 2D STRAIN @Element'
        write(ILOG,'(a,a)')                                             &
               'ELEMENT       E11         E22         E33     '         &
                             ,'    E12' 
        write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'    &
                             ,'-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p4e12.4)') jj,                            &
                                 (fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,4)
          endif
        enddo
        write(ILOG,*) '#### 2D STRESS @Element'
        write(ILOG,'(a,a)')                                             &
               'ELEMENT       S11         S22         S33     '         &
                             ,'    S12        MISES' 
        write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'    &
                             ,'-----------+-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p5e12.4)') jj,                              &
                                 (fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,4)    &
                                 ,fstrSOLID%ESTRESS(7*(ii-1)+7)
          endif
        enddo

      elseif( hecMESH%n_dof .eq. 6 ) then
        write(ILOG,*) '#### SHELL STRAIN @Element'
        write(ILOG,'(a,a)')                                               &
               '  ELEMENT          E11         E22         E12     '      &
                             ,'    E23         E13'
        write(ILOG,'(a,a)')                                               &
               '  ------------+-----------+-----------+-----------+'      &
                             ,'-----------+-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,a5,1p5e12.4)') jj,' (+) ',                   &
                                 (fstrSOLID%ESTRAIN(10*(ii-1)+k),k=1,5)
            write(ILOG,'(10x,a5,1p5e12.4)')   ' (-) ',                    &
                                 (fstrSOLID%ESTRAIN(10*(ii-1)+k),k=6,10)
          endif
        enddo
        write(ILOG,*) '#### SHELL STRESS @Element'
        write(ILOG,'(a,a)')                                               &
               '  ELEMENT          S11         S22         S12     '      &
                             ,'    S23         S13        MISES'
        write(ILOG,'(a,a)')                                               &
               '  ------------+-----------+-----------+-----------+'      &
                             ,'-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
          jj = hecMESH%global_elem_ID(i)
!!!          ii = hecMESH%elem_ID(i*2-1)
          ii = i
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,a5,1p6e12.4)') jj,' (+) ',                   &
                                 (fstrSOLID%ESTRESS(12*(ii-1)+k),k=1,6)
            write(ILOG,'(10x,a5,1p6e12.4)')   ' (-) ',                    &
                                 (fstrSOLID%ESTRESS(12*(ii-1)+k),k=7,12)
          endif
        enddo
      endif
   end subroutine fstr_output_elem
end module m_static_output
