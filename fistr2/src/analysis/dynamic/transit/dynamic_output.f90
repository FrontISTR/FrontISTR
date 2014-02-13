!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                       Tomotaka Ogasawara (Univ. of Tokyo)            !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!C================================================================C
!> This module contains functions for output in dynamic analysis
!C================================================================C

module m_dynamic_output
use m_fstr
use m_static_lib
use m_static_get_prop
use m_dynamic_post
use m_hecmw2fstr_mesh_conv
use m_static_make_result

implicit none


contains

!C***
!C*** OUTPUT for FSTR solver
!C***
!C
!C 3D SOLID
!C
      subroutine dynamic_output_3d(hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)

      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)         :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID
      type (fstr_dynamic)       :: fstrDYNAMIC
      integer(kind=kint) :: idx, i, ii, jj, k
!C
!C-- displacement, velocity and accelaration
!C
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!C
!C-- SET DISPLACEMENT
!C
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(3*i-2)= fstrDYNAMIC%DISP(3*i-2, idx)
        fstrSOLID%unode(3*i-1)= fstrDYNAMIC%DISP(3*i-1, idx)
        fstrSOLID%unode(3*i  )= fstrDYNAMIC%DISP(3*i  , idx)
      enddo

!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then

      write(ILOG,*) '#### DISPLACEMENT 3D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      Z-DISP'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%DISP(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!C
!C-- VELOCITY
!C
      if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### VELOCITY 3D'
      write(ILOG,'(a)')'    NODE     X-VELO      Y-VELO      Z-VELO'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%VEL(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### ACCELERATION 3D'
      write(ILOG,'(a)')'    NODE     X-ACCE      Y-ACCE      Z-ACCE'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%ACC(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!!
!C
!C-- end of displacement, velocity and accelaration
!C
!C
!C-- reaction force, strain and stress
!C
!      if(fstrDYNAMIC%i_step .gt. 0) then
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 3D'
      call dynamic_reaction_force_3d( hecMESH, fstrSOLID, fstrDYNAMIC)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 3D'
      call dynamic_nodal_stress_3d( hecMESH, fstrSOLID,fstrPARAM, fstrDYNAMIC)
!!
!      end if
!C
!C-- end of reaction force, strain and stress
!C
      end subroutine dynamic_output_3d
!C
!C NODAL STRESS 3D
!C
      subroutine dynamic_nodal_stress_3d(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      REAL(kind=kreal) xx(20),yy(20),zz(20)
      REAL(kind=kreal) edisp(60),estrain(6),estress(6)
      REAL(kind=kreal) val, tt(20),tt0(20)                 
      integer(kind=kint) nodLOCAL(20),ierror 
      REAL(kind=kreal) edstrain(20,6),edstress(20,6)
!** Array for nodal recovery
      integer (kind=KINT), dimension(:), allocatable :: nnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp
      real(kind=KREAL), allocatable :: ndstrain(:,:),ndstress(:,:)
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, ielem, ID_area, ii,jj
      real(kind=KREAL) :: s11, s22, s33, s12, s23, s13, ps, smises

!*Allocate array
      allocate ( ndstrain(hecMESH%n_node,6), ndstress(hecMESH%n_node,7) ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( nnumber(hecMESH%n_node)           ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
!*ZERO CLEAR 
      ndstrain=0.d0; ndstress=0.d0
      nnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node) ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
     !   if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        nn = getNumberOfNodes(ic_type)
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
          if (ic_type==361) then
            call UpdateST_C3D8IC(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn),        &
                           edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          else
            call UpdateST_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn)   &
                          ,edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          endif
          call NodalStress_C3(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:) )
!C** elem ID
!!        ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            call ElementStress_C3(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-1) = estress
            s11=fstrSOLID%ESTRESS(7*ielem-6)
            s22=fstrSOLID%ESTRESS(7*ielem-5)
            s33=fstrSOLID%ESTRESS(7*ielem-4)
            s12=fstrSOLID%ESTRESS(7*ielem-3)
            s23=fstrSOLID%ESTRESS(7*ielem-2)
            s13=fstrSOLID%ESTRESS(7*ielem-1)
            ps=(s11+s22+s33)/3.d0
            smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                  +s12**2+s23**2+s13**2
            fstrSOLID%ESTRESS(7*ielem)=sqrt(3.0*smises)
          end if

          do j=1,nn
            ndstrain( nodLocal(j),: ) = ndstrain(nodLOCAL(j),:) + edstrain(j,:)
            ndstress( nodLocal(j),1:6 ) = ndstress(nodLOCAL(j),1:6) + edstress(j,:)
            nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1  
          enddo
        enddo
      enddo
!** Average over nodes
      do i=1,hecMESH%n_node
        ndstrain(i,:)=ndstrain(i,:)/nnumber(i)
        ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i)
      enddo
!** RE-CALCULATE von MISES stress
      do i=1,hecMESH%n_node
        s11=ndstress(i,1)
        s22=ndstress(i,2)
        s33=ndstress(i,3)
        s12=ndstress(i,4)
        s23=ndstress(i,5)
        s13=ndstress(i,6)
        ps=(s11+s22+s33)/3.d0
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

      if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then

      write(ILOG,*) '#### STRAIN'
      write(ILOG,'(a,a)')                                                &
                    '     NODE     E11         E22         E33     '     &
                             ,'    E12         E23         E13     '
      write(ILOG,'(a,a)')                                                &
                    '  -------+-----------+-----------+-----------+'     &
                             ,'-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p6e12.4)') jj,(ndSTRAIN(ii,k),k=1,6)
      enddo

      end if
      if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then

      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                                &
                    '     NODE     S11         S22         S33     '     &
                 ,'    S12         S23         S13        MISES'
      write(ILOG,'(a,a)')                                                &
                    '  -------+-----------+-----------+-----------+'     &
                 ,'-----------+-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p7e12.4)') jj,(ndSTRESS(ii,k),k=1,7)
      enddo

      end if

!*Deallocate array
      deallocate( ndstrain,ndstress   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( nnumber      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( temp         ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      end subroutine dynamic_nodal_stress_3d
!C
!C REACTION FORCE 3D
!C
      subroutine dynamic_reaction_force_3d(hecMESH,fstrSOLID, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      integer, parameter :: NDOF=3
      REAL(kind=kreal) stiff(120,120)
      REAL(kind=kreal) edisp(60),force(60),ecoord(3,20)
      integer(kind=kint) nodLOCAL(20),ierror
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
	  integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, iflag, icount, nspc, nid


!C Allocate array
      allocate ( id_spc( hecMESH%n_node )     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
         allocate ( spcForce( nspc * NDOF )   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_3d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
    !    if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        nn = getNumberOfNodes(ic_type)
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
            else
              call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3),1.d0)
            endif
            force(1:nn*3)= matmul( stiff(1:nn*3,1:nn*3), edisp(1:nn*3) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc(  nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k) +force(NDOF*(j-1)+k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show

      if( fstrDYNAMIC%iout_list(4) .eq. 1 ) then

      write(ILOG,'(a)')'    NODE     X-REAC      Y-REAC      Z-REAC'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p3e12.4)') j,  &                          
          spcForce(NDOF*(nid-1)+1),        &                         
          spcForce(NDOF*(nid-1)+2),        &                        
          spcForce(NDOF*(nid-1)+3)
        endif
      enddo

      end if

!C Deallocate array
      deallocate ( spcForce     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate ( id_spc       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_reaction_force_3d
!C***
!C*** OUTPUT for FSTR solver
!C***
!C
!C 2D SOLID
!C
      subroutine DYNAMIC_OUTPUT_2D(hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)         :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC
      integer(kind=kint) :: idx
      integer :: i,k, ii,jj
!C
!C-- displacement, velocity and accelaration
!C
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!C
!C-- SET DISPLACEMENT
!C
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(2*i-1)= fstrDYNAMIC%DISP(2*i-1, idx)
        fstrSOLID%unode(2*i  )= fstrDYNAMIC%DISP(2*i  , idx)
      enddo

!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then

      write(ILOG,*) '#### DISPLACEMENT 2D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%DISP(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- VELOCITY
!C
      if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### VELOCITY 2D'
      write(ILOG,'(a)')'    NODE     X-VELO      Y-VELO      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%VEL(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### ACCELERATION 2D'
      write(ILOG,'(a)')'    NODE     X-ACCE      Y-ACCE      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%ACC(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- end of displacement, velocity and accelaration
!C
!C
!C-- reaction force, strain and stress
!C
!      if(fstrDYNAMIC%i_step .gt. 0) then
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 2D'
      call dynamic_reaction_force_2d( hecMESH, fstrSOLID, fstrDYNAMIC)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 2D'
      call dynamic_nodal_stress_2d( hecMESH, fstrSOLID, fstrPARAM, fstrDYNAMIC)
!!
!      end if
!C
!C-- end of reaction force, strain and stress
!C
      end subroutine DYNAMIC_OUTPUT_2D

!C
!C NODAL STRESS 2D
!C
      subroutine dynamic_nodal_stress_2d(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      REAL(kind=kreal) xx(8),yy(8)
      REAL(kind=kreal) edisp(16),estrain(4),estress(4)
      REAL(kind=kreal) val, tt(8),tt0(8)
      integer(kind=kint) nodLOCAL(8),ierror
      REAL(kind=kreal) edstrain(9,4),edstress(9,4)
!** Array for nodal recovery
      real    (kind=KREAL), dimension(:,:), allocatable :: arrayTotal
      integer (kind=KINT), dimension(:), allocatable :: nnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, ielem, ID_area, ii,jj
      real(kind=KREAL) :: s11, s22, s33, s12, s23, s13, ps, smises
	  
!*Allocate array
      allocate ( arrayTotal( hecMESH%n_node, 9  )    ,STAT=ierror) 
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( nnumber(hecMESH%n_node)             ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
!*ZERO CLEAR
      arrayTotal=0
      nnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node)      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
    !    if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = getNumberOfNodes(ic_type)
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
            edisp(2*j  )=fstrSOLID%unode(2*nodLOCAL(j))
          enddo  

!C** Create local stiffness
          call UpdateST_C2(ic_type,nn,xx(1:nn),yy(1:nn),tt(1:nn),tt0(1:nn),1.d0,              &
                     fstrSOLID%elements(icel)%iset,edisp(1:nn*2),fstrSOLID%elements(icel)%gausses )
          call NodalStress_C2(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:) )
!C** elem ID
!!        ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            call ElementStress_C2(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-1) = estress
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
          end if

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
        do j=1,9
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
        smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
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
      if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
      write(ILOG,*) '#### STRIN'
      write(ILOG,'(a,a)')                                         & 
               '   NODE       E11         E22         E33     '   &
                             ,'    E12' 
      write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p4e12.4)') jj,(fstrSOLID%STRAIN(6*(ii-1)+k),k=1,4)
      enddo
      end if

      if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                         &
               '   NODE       S11         S22         S33     '   &
                             ,'    S12        MISES' 
      write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p5e12.4)') jj,                               &
                                 (fstrSOLID%STRESS(7*(ii-1)+k),k=1,4)  &
                                 ,fstrSOLID%STRESS(7*(ii-1)+7)
      enddo
      end if
!*Deallocate array
      deallocate( arrayTotal      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( nnumber         ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( temp            ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_nodal_stress_2d

!C
!C REACTION FORCE 2D
!C
      subroutine dynamic_reaction_force_2d(hecMESH,fstrSOLID, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      integer, parameter :: NDOF=3
      REAL(kind=kreal) ecoord(2,8),stiff(18,18)
      REAL(kind=kreal) val, edisp(16),force(16)
      integer(kind=kint) nodLOCAL(8),ierror
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, jS, icel, iflag, nn
      integer :: i,j,k, in, ielem, ii,jj, icount, nspc, nid, itype

!C Allocate array
      allocate ( id_spc( hecMESH%n_node )    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
         allocate ( spcForce( nspc * NDOF )    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_2d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
    !    if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = getNumberOfNodes(ic_type)
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
!C** Create local stiffness
            call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses,1.d0        &
                        ,stiff(1:nn*2,1:nn*2),fstrSOLID%elements(icel)%iset)
            force = matmul( stiff(1:nn*2,1:nn*2), edisp(1:nn*2) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc( nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k) +force(NDOF*(j-1)+k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show
      if( fstrDYNAMIC%iout_list(4) .eq. 1 ) then

      write(ILOG,'(a)') '   NODE    X-REACTION  Y-REACTION'
      write(ILOG,'(a)') ' --------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p2e12.4)') j,spcForce(NDOF*(nid-1)+1), spcForce(NDOF*(nid-1)+2)
        endif
      enddo

      end if
!C Deallocate array
      deallocate ( spcForce    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate ( id_spc      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_reaction_force_2d
!C***
!C*** OUTPUT for FSTR solver
!C***
      subroutine dynamic_output_6d( hecMESH,hecMAT,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type ( hecmwST_matrix     ) :: hecMAT
      type ( hecmwST_local_mesh ) :: hecMESH
      type ( fstr_solid         ) :: fstrSOLID
      type ( fstr_param         ) :: fstrPARAM
      type (fstr_dynamic        ) :: fstrDYNAMIC

      INTEGER(kind=kint) NDOF
      integer(kind=kint) :: idx
      integer :: ig, ig0, iE0, ik, iE, icel
      integer :: i,j,k, ielem, ii,jj

      NDOF = hecMESH%n_dof
!C
!C-- displacement, velocity and accelaration
!C
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!C
!C-- SET DISPLACEMENT
!C
      do i = 1, hecMESH%n_node
        do j = 1, NDOF
          fstrSOLID%unode(NDOF*(i-1)+j)= fstrDYNAMIC%DISP(NDOF*(i-1)+j, idx)
        enddo
      enddo

!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then

      write(ILOG,*) '#### DISPLACEMENT SHELL'
      write(ILOG,'(a,a)')                                              &
                    '     NODE    X-DISP      Y-DISP      Z-DISP   '   &
                             ,'   X-ROT       Y-ROT       Z-ROT'
      write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,6e12.4)') jj,(fstrDYNAMIC%DISP(6*(ii-1)+k, idx),k=1,6)
      enddo

      end if
!C
!C-- VELOCITY
!C
      if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### VELOCITY SHELL'
      write(ILOG,'(a,a)')                                              &
                    '     NODE    X-VELO      Y-VELO      Z-VELO   '   &
                             ,'   X-ROT-VEL   Y-ROT-VEL   Z-ROT-VEL'
      write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,6e12.4)') jj,(fstrDYNAMIC%VEL(6*(ii-1)+k, idx),k=1,6)
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### ACCELERATION SHELL'
      write(ILOG,'(a,a)')                                              &
                    '     NODE    X-ACCE      Y-ACCE      Z-ACCE   '   &
                             ,'   X-ROT-ACC   Y-ROT-ACC   Z-ROT-ACC'
      write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,6e12.4)') jj,(fstrDYNAMIC%ACC(6*(ii-1)+k, idx),k=1,6)
      enddo

      end if
!C
!C-- end of displacement, velocity and accelaration
!C
!C
!C-- reaction force, strain and stress
!C
!      if(fstrDYNAMIC%i_step .gt. 0) then
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE SHELL'
      call dynamic_reaction_force_6d( hecMESH, fstrSOLID, fstrDYNAMIC)
      call flush(ILOG)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRAIN/STRESS SHELL'
      call dynamic_nodal_stress_6d( hecMESH, fstrSOLID,fstrPARAM, fstrDYNAMIC)
      call flush(ILOG)
!!
!      end if
!C
!C-- end of reaction force, strain and stress
!C
      end subroutine dynamic_output_6d



!C
!C NODAL STRESS SHELL
!C
      subroutine dynamic_nodal_stress_6d(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC

!** Local variables
      REAL(kind=kreal) xx(8),yy(8),zz(8),ri(4),si(4)
      REAL(kind=kreal) val, edisp(48), tt(8), tt0(8)
      integer(kind=kint) nodLOCAL(8),ierror
      REAL(kind=kreal) array(11)

!** Array for nodal recovery
      real    (kind=KREAL),dimension(:,:),allocatable::arrayTotal_U
      real    (kind=KREAL),dimension(:,:),allocatable::arrayTotal_D
      integer (kind=KINT), dimension(:),  allocatable::nnumber
      real    (kind=KREAL),dimension(:),  allocatable::temp
      real    (kind=kreal) ee,pp,rho,alpha,thick, ti
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, ielem, ID_area, ii,jj, isect
      real(kind=KREAL) :: s11, s22, s33, s12, s23, s13, ps, smises
      real(kind=KREAL) :: s1, s2, s3, s4,s5

      data ri / -1.0,  1.0, 1.0, -1.0 /
      data si / -1.0, -1.0, 1.0,  1.0 /

!*Allocate array
      allocate ( arrayTotal_U(hecMESH%n_node,11)   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_6d, arrayTotal_U>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( arrayTotal_D(hecMESH%n_node,11)   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_6d, arrayTotal_D>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( nnumber(hecMESH%n_node)           ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_6d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
 
      arrayTotal_U = 0.0d0
      arrayTotal_D = 0.0d0
      nnumber = 0
!C
!C Set Temperature

      allocate ( temp(hecMESH%n_node)               ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_6d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
    !    if (.not. hecmw_is_etype_shell(ic_type)) cycle
!C        
        nn = getNumberOfNodes(ic_type)
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
!!        ielem = hecMESH%elem_ID(icel*2-1)
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
          end if

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
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem-6)=sqrt(3.0*smises)
            end if

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
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem)=sqrt(3.0*smises)
            end if

          elseif (ic_type.eq.741) then
            do j = 1, 4
              jj = nodLOCAL(j)
              ti = 1.0
              call RCV_S4( xx,yy,zz,ee,pp,thick,edisp,ri(j),si(j),ti,array )
              do k = 1, 11
                arrayTotal_U(jj,k) = arrayTotal_U(jj,k) + array(k)
              enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              fstrSOLID%ESTRAIN(10*ielem-9 ) = fstrSOLID%ESTRAIN(10*ielem-9 ) + array(1)/nn
              fstrSOLID%ESTRAIN(10*ielem-8 ) = fstrSOLID%ESTRAIN(10*ielem-8 ) + array(2)/nn
              fstrSOLID%ESTRAIN(10*ielem-7 ) = fstrSOLID%ESTRAIN(10*ielem-7 ) + array(3)/nn
              fstrSOLID%ESTRAIN(10*ielem-6 ) = fstrSOLID%ESTRAIN(10*ielem-6 ) + array(4)/nn
              fstrSOLID%ESTRAIN(10*ielem-5 ) = fstrSOLID%ESTRAIN(10*ielem-5 ) + array(5)/nn
              fstrSOLID%ESTRESS(12*ielem-11) = fstrSOLID%ESTRESS(12*ielem-11) + array(6)/nn
              fstrSOLID%ESTRESS(12*ielem-10) = fstrSOLID%ESTRESS(12*ielem-10) + array(7)/nn
              fstrSOLID%ESTRESS(12*ielem-9 ) = fstrSOLID%ESTRESS(12*ielem-9 ) + array(8)/nn
              fstrSOLID%ESTRESS(12*ielem-8 ) = fstrSOLID%ESTRESS(12*ielem-8 ) + array(9)/nn
              fstrSOLID%ESTRESS(12*ielem-7 ) = fstrSOLID%ESTRESS(12*ielem-7 ) + array(10)/nn
            end if
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              s11=fstrSOLID%ESTRESS(12*ielem-11)
              s22=fstrSOLID%ESTRESS(12*ielem-10)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-9 )
              s23=fstrSOLID%ESTRESS(12*ielem-8 )
              s13=fstrSOLID%ESTRESS(12*ielem-7 )
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem-6)=sqrt(3.0*smises)
            end if

            do j = 1, 4
              jj = nodLOCAL(j)
              ti = -1.0
              call RCV_S4( xx,yy,zz,ee,pp,thick,edisp,ri(j),si(j),ti,array )
              do k = 1, 11
                arrayTotal_D(jj,k) = arrayTotal_D(jj,k) + array(k)
              enddo
              nnumber(jj) = nnumber(jj) + 1

              if( ID_area.eq.hecMESH%my_rank ) then
                fstrSOLID%ESTRAIN(10*ielem-4) = fstrSOLID%ESTRAIN(10*ielem-4) + array(1)/nn
                fstrSOLID%ESTRAIN(10*ielem-3) = fstrSOLID%ESTRAIN(10*ielem-3) + array(2)/nn
                fstrSOLID%ESTRAIN(10*ielem-2) = fstrSOLID%ESTRAIN(10*ielem-2) + array(3)/nn
                fstrSOLID%ESTRAIN(10*ielem-1) = fstrSOLID%ESTRAIN(10*ielem-1) + array(4)/nn
                fstrSOLID%ESTRAIN(10*ielem  ) = fstrSOLID%ESTRAIN(10*ielem  ) + array(5)/nn
                fstrSOLID%ESTRESS(12*ielem-5) = fstrSOLID%ESTRESS(12*ielem-5) + array(6)/nn
                fstrSOLID%ESTRESS(12*ielem-4) = fstrSOLID%ESTRESS(12*ielem-4) + array(7)/nn
                fstrSOLID%ESTRESS(12*ielem-3) = fstrSOLID%ESTRESS(12*ielem-3) + array(8)/nn
                fstrSOLID%ESTRESS(12*ielem-2) = fstrSOLID%ESTRESS(12*ielem-2) + array(9)/nn
                fstrSOLID%ESTRESS(12*ielem-1) = fstrSOLID%ESTRESS(12*ielem-1) + array(10)/nn
              end if
            enddo

            if( ID_area.eq.hecMESH%my_rank ) then
              s11=fstrSOLID%ESTRESS(12*ielem-5)
              s22=fstrSOLID%ESTRESS(12*ielem-4)
              s33=0.0D0
              s12=fstrSOLID%ESTRESS(12*ielem-3)
              s23=fstrSOLID%ESTRESS(12*ielem-2)
              s13=fstrSOLID%ESTRESS(12*ielem-1)
              ps=(s11+s22+s33)/3.0
              smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
              fstrSOLID%ESTRESS(12*ielem)=sqrt(3.0*smises)
            end if

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
        fstrSOLID%STRESS(12*(i-1)+6) = DSQRT( s1*s1 + s2*s2 - s1*s2 + 3.0d0*s3*s3 )

        s1 = fstrSOLID%STRESS(12*(i-1)+7)
        s2 = fstrSOLID%STRESS(12*(i-1)+8)
        s3 = fstrSOLID%STRESS(12*(i-1)+9)
        s4 = fstrSOLID%STRESS(12*(i-1)+10)
        s5 = fstrSOLID%STRESS(12*(i-1)+11)
        fstrSOLID%STRESS(12*i) = DSQRT( s1*s1 + s2*s2 - s1*s2 + 3.0d0*s3*s3 )

      enddo

!** Show
      if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
      write(ILOG,*) '#### STRAIN'
      write(ILOG,'(a,a)')                                              &
               '     NODE          E11         E22         E12     '   &
                             ,'    E23         E13'
      write(ILOG,'(a,a)')                                              &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,a5,1p5e12.4)') jj,' (+) ',(fstrSOLID%STRAIN(10*(ii-1)+k),k=1,5)
        write(ILOG,'(10x,a5,1p5e12.4)')   ' (-) ',(fstrSOLID%STRAIN(10*(ii-1)+k),k=6,10)
      enddo
      end if

      if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
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
        write(ILOG,'(i10,a5,1p6e12.4)') jj,' (+) ',(fstrSOLID%STRESS(12*(ii-1)+k),k=1,6)
        write(ILOG,'(10x,a5,1p6e12.4)')   ' (-) ',(fstrSOLID%STRESS(12*(ii-1)+k),k=7,12)
      enddo
      end if

      deallocate( arrayTotal_U   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_6d, arrayTotal_U>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( arrayTotal_D   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_6d, arrayTotal_D>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( nnumber        ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_6d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( temp           ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_6d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_nodal_stress_6d


!C
!C REACTION FORCE SHELL
!C
      subroutine dynamic_reaction_force_6d(hecMESH,fstrSOLID, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC

!** Local variables
      REAL(kind=kreal)    xx(8),yy(8),zz(8)
      REAL(kind=kreal)    edisp(48),force(48),local_stf(1176)
      integer(kind=kint) nodLOCAL(8),ierror

      real   (kind=KREAL),dimension(:), allocatable :: spcForce
      integer(kind=KINT), dimension(:), allocatable :: id_spc
	  real    (kind=kreal) ee,pp,rho,alpha,thick, ti
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel, isid, isect, num
      integer :: i,j,k,kk, in, itype, ielem, ID_area, ii,jj, NDOF, icount, nspc, ifspc

!C Allocate array
      allocate ( id_spc( hecMESH%n_node )      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_6d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

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
         allocate ( spcForce(nspc*NDOF)    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_6d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

         spcForce = 0.0
      else
        deallocate ( id_spc                ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_6d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

        return
      endif

!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+

      do itype = 1, hecMESH%n_elem_type
        iS = hecMESH%elem_type_index(itype-1) + 1
        iE = hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
    !    if (.not. hecmw_is_etype_shell(ic_type)) cycle

!C** Set number of nodes        
        nn = getNumberOfNodes(ic_type)

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
      if( fstrDYNAMIC%iout_list(4) .eq. 1 ) then

      write(ILOG,'(a,a)')                                             &
                    '   NODE      X-REAC      Y-REAC      Z-REAC   '  &
                             ,'  RX-REAC     RY-REAC     RZ-REAC'
      write(ILOG,'(a,a)')                                             &
                    ' --------+-----------+-----------+-----------+'  &
                             ,'-----------+-----------+-------------'
      do i = 1, hecMESH%nn_internal
        in = hecMESH%global_node_ID(i)
        isid = id_spc(i)
        if( isid .gt. 0 ) then
          write(ILOG,'(i10,1p6e12.4)') in                             &
          , spcForce( NDOF*(isid-1)+1 ), spcForce( NDOF*(isid-1)+2 )  &
          , spcForce( NDOF*(isid-1)+3 ), spcForce( NDOF*(isid-1)+4 )  &
          , spcForce( NDOF*(isid-1)+5 ), spcForce( NDOF*(isid-1)+6 ) 
        endif
      enddo

      end if

      deallocate ( id_spc          ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_6d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate ( spcForce        ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_6d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_reaction_force_6d
!C 
!C OUTPUT FOR ELEMENTS
!C
      subroutine dynamic_output_elem(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel, isid, isect, num
      integer :: i,j,k,kk, in, itype, ielem, ID_area, ii,jj, NDOF, icount, nspc, ifspc

      if( hecMESH%n_dof .eq. 3 ) then

        if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then

        write(ILOG,*) '#### 3D STRAIN @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     E11         E22         E33     '  &
                             ,'    E12         E23         E13     '
        write(ILOG,'(a,a)')                                            &
                    '  -------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p6e12.4)') jj,(fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,6)
          end if
        enddo
        end if

        if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then

        write(ILOG,*) '#### 3D STRESS @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     S11         S22         S33     '  &
                 ,'    S12         S23         S13        MISES'
        write(ILOG,'(a,a)')                                           &
                    '  -------+-----------+-----------+-----------+'  &
                 ,'-----------+-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p7e12.4)') jj,(fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,7)
          end if
        enddo
        end if

      elseif( hecMESH%n_dof .eq. 2 ) then

        if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
        write(ILOG,*) '#### 2D STRAIN @Element'
        write(ILOG,'(a,a)')                                      &
               'ELEMENT       E11         E22         E33     '  &
                             ,'    E12' 
        write(ILOG,'(a,a)')                                           &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p4e12.4)') jj,(fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,4)
          end if
        enddo
        end if

        if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
        write(ILOG,*) '#### 2D STRESS @Element'
        write(ILOG,'(a,a)')                                     &
               'ELEMENT       S11         S22         S33     ' &
                             ,'    S12        MISES' 
        write(ILOG,'(a,a)')                                            &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p5e12.4)') jj,(fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,4)  &
                                   ,fstrSOLID%ESTRESS(7*(ii-1)+7)
          end if
        enddo
        end if

      elseif( hecMESH%n_dof .eq. 6 ) then

        if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
        write(ILOG,*) '#### SHELL STRAIN @Element'
        write(ILOG,'(a,a)')                                            &
               '  ELEMENT          E11         E22         E12     '   &
                             ,'    E23         E13'
        write(ILOG,'(a,a)')                                           &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,a5,1p5e12.4)') jj,' (+) ',(fstrSOLID%ESTRAIN(10*(ii-1)+k),k=1,5)
            write(ILOG,'(10x,a5,1p5e12.4)')    ' (-) ',(fstrSOLID%ESTRAIN(10*(ii-1)+k),k=6,10)
          end if
        enddo
        end if

        if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
        write(ILOG,*) '#### SHELL STRESS @Element'
        write(ILOG,'(a,a)')                                             &
               '  ELEMENT          S11         S22         S12     '    &
                             ,'    S23         S13        MISES'
        write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'    &
                             ,'-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,a5,1p6e12.4)') jj,' (+) ',(fstrSOLID%ESTRESS(12*(ii-1)+k),k=1,6)
            write(ILOG,'(10x,a5,1p6e12.4)')    ' (-) ',(fstrSOLID%ESTRESS(12*(ii-1)+k),k=7,12)
          end if
        enddo
        end if

      endif

      end subroutine dynamic_output_elem
	  
!C================================================================C
!C-- subroutine  matvec
!C================================================================C
  subroutine matvec(y,x,hecMAT,ndof,D,AU,AL)
    use hecmw
    type ( hecmwST_matrix      ) :: hecMAT

    integer(kind=kint) :: ndof,i,is,ie,j,icol,ki,kj,ix,iy,ip,nn
    real(kind=kreal) :: D(ndof*ndof*hecMAT%NP)
    real(kind=kreal) :: AU(ndof*ndof*hecMAT%NPU)
    real(kind=kreal) :: AL(ndof*ndof*hecMAT%NPL)
    real(kind=kreal) :: x(ndof*hecMAT%NP)
    real(kind=kreal) :: y(ndof*hecMAT%NP)

    nn=ndof*ndof

    y=0.

    do i=1,hecMAT%NP
      is=hecMAT%indexU(i-1)+1
      ie=hecMAT%indexU(i)
      do j=is,ie
        icol=hecMAT%itemU(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AU(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      is=hecMAT%indexL(i-1)+1
      ie=hecMAT%indexL(i)
      do j=is,ie
        icol=hecMAT%itemL(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AL(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      do ki=1,ndof
        iy=ndof*(i-1)+ki
        do kj=1,ndof
          ix=ndof*(i-1)+kj
          ip=nn*(i-1)+ndof*(ki-1)+kj
          y(iy)=y(iy)+D(ip)*x(ix)
        enddo
      enddo
    enddo

  end subroutine matvec

  subroutine dynamic_output_and_post(istep,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      use hecmw

      integer (kind = kint) :: istep
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      integer (kind = kint) :: idummy

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if
!
!C-- output new displacement, velocity and accelaration
      if( mod(istep,fstrDYNAMIC%nout) /= 0 ) return
!
      write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP =',istep,'TIME =',fstrDYNAMIC%t_delta*istep
!C-- output

      call flush(ILOG)
      IF(myrank .EQ. 0) THEN
        IF(hecMAT%Iarray(99) .EQ. 1) THEN
          write(IMSG,*) '*----------------------------*'
          write(IMSG,*) '## No.of ITER:',hecMAT%Iarray(1)
          write(IMSG,*) '*----------------------------*'
        ENDIF
        call flush(IDBG)
        call flush(IMSG)
      ENDIF

!C-- output
      if( hecMESH%n_dof .eq. 3 ) then
        call dynamic_output_3d(hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrDYNAMIC)
      else if( hecMESH%n_dof .eq. 2 ) then
        call dynamic_output_2d(hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrDYNAMIC)
      else if( hecMESH%n_dof.EQ.6) THEN
        call dynamic_output_6d(hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrDYNAMIC)
      endif

!C-- OUTPUT FOR ELEMENTS
      call dynamic_output_elem(hecMESH,fstrSOLID,fstrPARAM,fstrDYNAMIC)
!
      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- POST PROCESSING
!C
      call dynamic_post(hecMESH,hecMAT,fstrSOLID,istep,fstrDYNAMIC)

!C
!C-- POST PROCESSING VIA MEMORY
!C
      if( IVISUAL.eq. 1 ) then
        call fstr_make_result(fstrSOLID,fstrRESULT)
        call fstr2hecmw_mesh_conv(hecMESH)
     !   call hecmw_visualize_init
        idummy=1
     !   call hecmw_visualize(hecMESH,fstrRESULT,istep,fstrDYNAMIC%n_step,idummy)
     !   call hecmw_visualize_finalize
     !   call hecmw2fstr_mesh_conv(hecMESH)
      !  call hecmw_result_free(fstrRESULT)
      endif

  end subroutine dynamic_output_and_post


!C================================================================C
!C-- subroutine dynamic_output_monit
!C================================================================C
  subroutine dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, dynamic_IW4, dynamic_IW5,  &
                                    dynamic_IW6, my_rank_monit_1)
      use hecmw
!C
!C-- global variable
!C
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      type ( fstr_param          ) :: fstrPARAM
!C
!C-- local variable
!C
    integer(kind=kint) :: dynamic_IW4
    integer(kind=kint) :: dynamic_IW5
    integer(kind=kint) :: dynamic_IW6
    integer(kind=kint) :: ii,jj
    integer(kind=kint) :: my_rank_monit_1

    if( mod(fstrDYNAMIC%i_step,fstrDYNAMIC%nout_monit) /= 0 ) return
!!
    if( hecMESH%my_rank .eq. my_rank_monit_1) then
!
      jj=fstrPARAM%global_local_id(1,fstrDYNAMIC%node_monit_1)
      ii=fstrPARAM%global_local_id(2,fstrDYNAMIC%node_monit_1)
!!
!!    step = 0
!!
      if(fstrDYNAMIC%i_step .eq. 0) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!!
!!    step = 1,2,.......
!!
      else if(fstrDYNAMIC%i_step .gt. 0) then
!
!C-- implicit dynamic analysis
        if(fstrDYNAMIC%idx_eqa .eq. 1) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )
        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )
        end if
!
!C-- explicit dynamic analysis
        else if(fstrDYNAMIC%idx_eqa .eq. 11) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &
                 ,fstrDYNAMIC%t_delta*fstrDYNAMIC%i_step &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
!                ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )
        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step -1 &
           write(dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step -1 &
                 ,fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1) &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step -1 &
           write(dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step -1 &
                 ,fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1) &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!
        end if
!
      end if
    end if

  end subroutine dynamic_output_monit

end module m_dynamic_output
