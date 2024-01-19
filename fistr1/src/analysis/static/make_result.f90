!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to prepare output of static analysis
module m_make_result
  private

  public:: fstr_write_result
  public:: fstr_make_result
  public:: fstr_reorder_node_shell
  public:: fstr_reorder_rot_shell
  public:: fstr_reorder_node_beam
  public:: setup_contact_output_variables


contains

  !C***
  !>  OUTPUT result file for static and dynamic analysis
  !C***
  subroutine fstr_write_result( hecMESH, fstrSOLID, fstrPARAM, istep, time, flag, fstrDYNAMIC)
    use m_fstr
    use m_out
    use m_static_lib
    use mMaterial
    use hecmw_util

    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type (fstr_param       )  :: fstrPARAM    !< analysis control parameters
    integer(kind=kint)        :: istep, flag
    type (fstr_dynamic), intent(in), optional  :: fstrDYNAMIC
    real(kind=kreal)          :: time        !< current time
    integer(kind=kint) :: n_lyr, ntot_lyr, tmp, is_33shell, is_33beam, cid
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, mm, ngauss, it
    real(kind=kreal), pointer :: tnstrain(:), testrain(:), yield_ratio(:)
    integer(kind=kint) :: idx
    real(kind=kreal), allocatable   :: work(:), unode(:), rnode(:)
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN)    :: comment
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname, cnum
    character(len=6), allocatable   :: clyr(:)
    logical :: is_dynamic

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN
    yield_ratio => fstrSOLID%YIELD_RATIO

    is_dynamic = present(fstrDYNAMIC)

    if( is_dynamic ) then
      idx = 1
      if( fstrDYNAMIC%idx_eqa==1 .and. istep>0 ) idx = 2
    endif

    ndof = hecMESH%n_dof
    mm   = hecMESH%n_node
    if( hecMESH%n_elem > hecMESH%n_node ) mm = hecMESH%n_elem
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==4 ) mdof = 6
    if( ndof==6 ) mdof = 6

    ntot_lyr   = fstrSOLID%max_lyr
    is_33shell = fstrSOLID%is_33shell
    is_33beam  = fstrSOLID%is_33beam

    nn = mm * mdof
    allocate( work(nn) )

    ! --- INITIALIZE
    header = '*fstrresult'
    if( present(fstrDYNAMIC) ) then
      comment = 'dynamic_result'
    else
    comment = 'static_result'
    endif
    call hecmw_result_init( hecMESH, istep, header, comment )

    ! --- TIME
    id = HECMW_RESULT_DTYPE_GLOBAL
    label = 'TOTALTIME'
    work(1) = time
    call hecmw_result_add( id, 1, label, work )

    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(3)%outinfo%on(1) ) then
      if(ndof == 4) then
        id = HECMW_RESULT_DTYPE_NODE
        ! for VELOCITY
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), 3 )
        allocate( unode(3*hecMESH%n_node) )
        unode = 0.0d0
        do i=1, hecMESH%n_node
          do j = 1, 3
            unode((i-1)*3 + j) = fstrDYNAMIC%DISP((i-1)*4 + j, idx)
          enddo
        enddo
        label = 'VELOCITY'
        call hecmw_result_add( id, nitem, label, unode )
        deallocate( unode )
        ! for PRESSURE
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), 1 )
        allocate( unode(hecMESH%n_node) )
        unode = 0.0d0
        do i=1, hecMESH%n_node
          unode(i) = fstrDYNAMIC%DISP(i*4, idx)
        enddo
        label = 'PRESSURE'
        call hecmw_result_add( id, nitem, label, unode )
        deallocate( unode )
      else if(ndof == 6) then
        id = HECMW_RESULT_DTYPE_NODE
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), 3 )
        allocate( unode(hecMESH%n_node*3) )
        unode = 0.0d0
        if( is_dynamic ) then
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrDYNAMIC%DISP((i-1)*ndof+1:(i-1)*ndof+3, idx)
          enddo
        else
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrSOLID%unode((i-1)*ndof+1:(i-1)*ndof+3)
          enddo
        endif
        label = 'DISPLACEMENT'
        call hecmw_result_add( id, nitem, label, unode )
        deallocate( unode )
      else
        id = HECMW_RESULT_DTYPE_NODE
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), ndof )
        allocate( unode(hecMESH%n_node*ndof) )
        unode = 0.0d0
        if( is_dynamic ) then
          unode(:) = fstrDYNAMIC%DISP(:,idx)
        else
          unode(:) = fstrSOLID%unode
        endif
        label = 'DISPLACEMENT'
        if(is_33beam == 1)then
          call fstr_reorder_node_beam(fstrSOLID, hecMESH, unode)
        endif
        if(is_33shell == 1)then
          call fstr_reorder_node_shell(fstrSOLID, hecMESH, unode)
        endif
        call hecmw_result_add( id, nitem, label, unode )
        deallocate( unode )
      endif
    endif

    ! --- ROTATION
    if (fstrSOLID%output_ctrl(3)%outinfo%on(18)) then
      if(ndof == 6) then
        id = HECMW_RESULT_DTYPE_NODE
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), 3 )
        label = 'ROTATION'
        allocate( rnode(hecMESH%n_node*3) )
        rnode = 0.0d0
        if( is_dynamic ) then
          do i=1, hecMESH%n_node
            rnode((i-1)*3+1:(i-1)*3+3) = fstrDYNAMIC%DISP((i-1)*ndof+4:(i-1)*ndof+6, idx)
          enddo
        else
          do i=1, hecMESH%n_node
            rnode((i-1)*3+1:(i-1)*3+3) = fstrSOLID%unode((i-1)*ndof+4:(i-1)*ndof+6)
          enddo
        endif
        call hecmw_result_add( id, nitem, label, rnode )
        deallocate( rnode )
      else
        if ( is_33shell == 1) then
          id = HECMW_RESULT_DTYPE_NODE
          nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), ndof )
          label = 'ROTATION'
          allocate( rnode(hecMESH%n_node*ndof) )
          rnode = 0.0d0
          call fstr_reorder_rot_shell(fstrSOLID, hecMESH, rnode)
          call hecmw_result_add( id, nitem, label, rnode )
          deallocate( rnode )
        end if
      end if
    endif

    ! --- VELOCITY
    if( is_dynamic .and. fstrSOLID%output_ctrl(3)%outinfo%on(15) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(15), ndof )
      label = 'VELOCITY'
      call hecmw_result_add( id, nitem, label, fstrDYNAMIC%VEL(:,idx) )
    endif

    ! --- ACCELERATION
    if( is_dynamic .and. fstrSOLID%output_ctrl(3)%outinfo%on(16) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(16), ndof )
      label = 'ACCELERATION'
      call hecmw_result_add( id, nitem, label, fstrDYNAMIC%ACC(:,idx) )
    endif

    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(3)%outinfo%on(2) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(2), ndof )
      label = 'REACTION_FORCE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%REACTION )
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(is_33shell == 1 .or. ndof == 6)then
      call fstr_write_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL, "      " )
    else
      call fstr_write_result_main( hecMESH, fstrSOLID, fstrSOLID%SOLID, "      " )
    endif

    !laminated shell
    if( associated(fstrSOLID%SHELL) .and. fstrSOLID%output_ctrl(3)%outinfo%on(27) ) then
      allocate(clyr(2*ntot_lyr))
      do i=1,ntot_lyr
        write(cnum,"(i0)")i
        clyr(2*i-1)="_L"//trim(cnum)//"+"
        clyr(2*i  )="_L"//trim(cnum)//"-"
      enddo
      do i=1,ntot_lyr
        call fstr_write_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL%LAYER(i)%PLUS,  clyr(2*i-1) )
        call fstr_write_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL%LAYER(i)%MINUS, clyr(2*i  ) )
      enddo
      deallocate(clyr)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! --- STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(9) .and. ndof/=6 ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(9), ndof )
      ngauss = fstrSOLID%maxn_gauss
      work(:) = 0.d0
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( associated(fstrSOLID%elements(i)%gausses) ) then
            if( k <= size(fstrSOLID%elements(i)%gausses) ) then
              do j = 1, nitem
                work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%strain_out(j)
              enddo
            endif
          end if
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif

    ! --- STRESS @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(10) .and. ndof/=6 ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(10), ndof )
      ngauss = fstrSOLID%maxn_gauss
      work(:) = 0.d0
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRESS',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( associated(fstrSOLID%elements(i)%gausses) ) then
            if( k <= size(fstrSOLID%elements(i)%gausses) ) then
              do j = 1, nitem
                work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%stress_out(j)
              enddo
            endif
          end if
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif

    ! --- PLASTIC STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(11) .and. fstrSOLID%StaticType/=3 ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(11), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'PLASTIC_GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( k > size(fstrSOLID%elements(i)%gausses) ) then
            work(i) = 0.d0
          else
            work(i) = fstrSOLID%elements(i)%gausses(k)%plstrain
          endif
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif

    ! --- THERMAL STRAIN @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(12) .and. associated(tnstrain) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(12), ndof )
      label = 'THERMAL_NodalSTRAIN'
      call hecmw_result_add( id, nitem, label, tnstrain )
    endif

    ! --- THERMAL STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(13) .and. associated(testrain) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(13), ndof )
      label = 'THERMAL_ElementalSTRAIN'
      call hecmw_result_add( id, nitem, label, testrain )
    endif

    ! --- THERMAL STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(14) .and. associated(testrain) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(14), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'THERMAL_GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( k > ngauss ) then
            do j = 1, nitem
              work(nitem*(i-1)+j) = 0.d0
            enddo
          else
            do j = 1, nitem
              !                work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%tstrain(j)
            enddo
          end if
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif

    ! --- YIELD RATIO
    if( fstrSOLID%output_ctrl(3)%outinfo%on(29) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(29), ndof )
      label = "YIELD_RATIO"
      call hecmw_result_add( id, nitem, label, yield_ratio )
    endif

    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(30), ndof )
      label = 'CONTACT_NFORCE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_NFORCE )
    endif

    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(31), ndof )
      label = 'CONTACT_FRICTION'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_FRIC )
    endif

    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(32), ndof )
      label = 'CONTACT_RELVEL'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_RELVEL )
    endif

    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(33), ndof )
      label = 'CONTACT_STATE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_STATE )
    endif

    ! --- CONTACT NORMAL TRACTION @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(36) .and. associated(fstrSOLID%CONT_NTRAC) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(36), ndof )
      label = 'CONTACT_NTRACTION'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_NTRAC )
    endif

    ! --- CONTACT FRICTION TRACTION @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(37) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(37), ndof )
      label = 'CONTACT_FTRACTION'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_FTRAC )
    endif

    ! --- TEMPERATURE @node
    !if(present(fstrPARAM))then
    !  if( fstrSOLID%output_ctrl(4)%outinfo%on(38) .and. fstrPARAM%solution_type == kstHEATSTATIC)then
    !    ncomp = ncomp + 1
    !    nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(38), ndof )
    !    fstrRESULT%nn_dof(ncomp) = nn
    !    fstrRESULT%node_label(ncomp) = 'TEMPERATURE'
    !    do i = 1, hecMESH%n_node
    !      do j = 1, nn
    !        fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%temperature(nn*(i-1)+j)
    !      enddo
    !    enddo
    !    iitem = iitem + nn
    !  endif
    !endif

    ! --- WRITE
    nameID = 'fstrRES'
    if( flag==0 ) then
      call hecmw_result_write_by_name( nameID )
    else
      addfname = '_dif'
      call hecmw_result_write_by_addfname( nameID, addfname )
    endif

    ! --- FINALIZE
    call hecmw_result_finalize

    deallocate( work )
  end subroutine fstr_write_result

  subroutine fstr_write_result_main( hecMESH, fstrSOLID, RES, clyr )
    use m_fstr
    use m_out
    use m_static_lib
    use mMaterial
    use hecmw_util

    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type (fstr_solid_physic_val) :: RES
    integer(kind=kint)        :: istep, flag
    integer(kind=kint)        :: n_lyr, cid

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname
    character(len=6)                :: clyr
    character(len=4)                :: cnum
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, mm, ngauss, it
    real(kind=kreal), allocatable   :: work(:)

    ndof = hecMESH%n_dof
    allocate( work(hecMESH%n_elem) )

    ! --- STRAIN @node
    if (fstrSOLID%output_ctrl(3)%outinfo%on(3)) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(3), ndof )
      label = 'NodalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%STRAIN )
    endif

    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(4) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(4), ndof )
      label = 'NodalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%STRESS )
    endif

    ! --- MISES @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(5) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(5), ndof )
      label = 'NodalMISES'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%MISES )
    endif

    ! --- NODAL PRINC STRESS
    if( fstrSOLID%output_ctrl(3)%outinfo%on(19) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(19), ndof )
      label = 'NodalPrincipalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%PSTRESS )
    endif

    ! --- NODAL PRINC STRAIN
    if( fstrSOLID%output_ctrl(3)%outinfo%on(21) ) then
      id = HECMW_RESULT_DTYPE_NODE
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(21), ndof )
      label = 'NodalPrincipalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%PSTRAIN )
    endif

    ! --- NODAL PRINC STRESS VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(23) ) then
      id = HECMW_RESULT_DTYPE_NODE
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(23), ndof )
        label = 'NodalPrincipalSTRESSVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%PSTRESS_VECT(:,k) )
      end do
    endif

    ! --- NODAL PRINC STRAIN VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(25) ) then
      id = HECMW_RESULT_DTYPE_NODE
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(25), ndof )
        label = 'NodalPrincipalSTRAINVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%PSTRAIN_VECT(:,k) )
      end do
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! --- STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(6) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(6), ndof )
      label = 'ElementalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%ESTRAIN )
    endif

    ! --- STRESS @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(7) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(7), ndof )
      label = 'ElementalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%ESTRESS )
    endif

    ! --- NQM @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(35) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(35), ndof )
      label = 'ElementalNQM'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%ENQM )
    endif

    ! --- MISES @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(8)) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(8), ndof )
      label = 'ElementalMISES'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EMISES )
    endif

    ! --- Principal_STRESS @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(20) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(20), ndof )
      label = 'ElementalPrincipalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EPSTRESS )
    endif

    ! --- Principal_STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(22) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(22), ndof )
      label = 'ElementalPrincipalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EPSTRAIN )
    endif

    ! --- ELEM PRINC STRESS VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(24) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(24), ndof )
        label = 'ElementalPrincipalSTRESSVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%EPSTRESS_VECT(:,k) )
      end do
    endif

    !ELEM PRINC STRAIN VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(26) ) then
      id = HECMW_RESULT_DTYPE_ELEM
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(26), ndof )
        label = 'ElementalPrincipalSTRAINVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%EPSTRAIN_VECT(:,k) )
      end do
    endif

    ! --- PLASTIC STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(43) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(43), ndof )
      ngauss = fstrSOLID%maxn_gauss
      label = 'ElementalPLSTRAIN'//trim(clyr)
      do i = 1, hecMESH%n_elem
        work(i) = 0.d0
        do j = 1, size(fstrSOLID%elements(i)%gausses) 
          work(i) = work(i) + fstrSOLID%elements(i)%gausses(j)%plstrain
        enddo
        work(i) = work(i) / size(fstrSOLID%elements(i)%gausses) 
        RES%EPLSTRAIN(i) = work(i)
      enddo
      call hecmw_result_add( id, nitem, label, work )
    endif
    deallocate( work )
  
  end subroutine fstr_write_result_main

  !C***
  !>  MAKE RESULT for static and dynamic analysis (WITHOUT ELEMENTAL RESULTS) --------------------------------------------------------------
  !C***
  subroutine fstr_make_result( hecMESH, fstrSOLID, fstrRESULT, istep, time, fstrDYNAMIC )
    use m_fstr
    use hecmw_util

    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type(hecmwST_result_data) :: fstrRESULT
    integer(kind=kint)        :: istep
    real(kind=kreal) :: time
    type(fstr_dynamic), intent(in), optional  :: fstrDYNAMIC
    integer(kind=kint) :: n_lyr, ntot_lyr, it, coef33, is_33shell, is_33beam
    integer(kind=kint) :: i, j, k, ndof, mdof, gcomp, gitem, ncomp, nitem, iitem, ecomp, eitem, jitem, nn, mm
    integer(kind=kint) :: idx
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)
    real(kind=kreal), allocatable   ::unode(:)
    character(len=4) :: cnum
    character(len=6), allocatable   :: clyr(:)
    logical :: is_dynamic

    is_dynamic = present(fstrDYNAMIC)

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN

    ntot_lyr   = fstrSOLID%max_lyr
    is_33shell = fstrSOLID%is_33shell
    is_33beam  = fstrSOLID%is_33beam

    mm = hecMESH%n_node
    if( hecMESH%n_elem>hecMESH%n_node ) mm = hecMESH%n_elem

    if( is_dynamic ) then
      idx = 1
      if( fstrDYNAMIC%idx_eqa==1 .and. istep>0 ) idx = 2
    endif

    ndof = hecMESH%n_dof
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==4 ) mdof = 6
    if( ndof==6 ) mdof = 6

    if(is_33shell == 1 .and. fstrSOLID%output_ctrl(4)%outinfo%on(27) )then
      coef33 = 1 + 2*ntot_lyr
    else
      coef33 = 1
    endif

    call hecmw_nullify_result_data( fstrRESULT )
    gcomp = 0
    gitem = 0
    ncomp = 0
    nitem = 0
    ecomp = 0
    eitem = 0

    ! --- COUNT SUM OF ALL NITEM
    ! --- TIME
    gcomp = gcomp + 1
    gitem = gitem + 1
    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
      if(ndof == 4) then
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 1 )
      else if(ndof == 6) then
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
      else
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
      endif
    endif
    ! --- VELOCITY
    if( is_dynamic .and. fstrSOLID%output_ctrl(4)%outinfo%on(15) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(15), ndof )
    endif
    ! --- ACCELERATION
    if( is_dynamic .and. fstrSOLID%output_ctrl(4)%outinfo%on(16) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(16), ndof )
    endif
    ! --- ROTATION (Only for 781 shell)
    if( fstrSOLID%output_ctrl(4)%outinfo%on(18) ) then
      if(ndof == 6) then
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(18), 3 )
      else
        if( is_33shell == 1 ) then
          ncomp = ncomp + 1
          nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(18), ndof )
        endif
      endif
    endif
    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(4)%outinfo%on(2) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(2), ndof )
    endif
    ! --- STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(3) ) then
      ncomp = ncomp + 1*coef33
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(3), ndof )*coef33
    endif
    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(4) ) then
      ncomp = ncomp + 1*coef33
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(4), ndof )*coef33
    endif
    ! --- MISES @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(5) ) then
      ncomp = ncomp + 1*coef33
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(5), ndof )*coef33
    endif
    ! --- Principal Stress @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(19) ) then
      ncomp = ncomp + 1*coef33
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(19), ndof )*coef33
    endif
    ! --- Principal Strain @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(21) ) then
      ncomp = ncomp + 1*coef33
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(21), ndof )*coef33
    endif
    ! --- Principal Stress Vector @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(23) ) then
      ncomp = ncomp + 3*coef33
      nitem = nitem + 3*n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(23), ndof )*coef33
    endif
    ! --- Principal Strain Vector @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(25) ) then
      ncomp = ncomp + 3*coef33
      nitem = nitem + 3*n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(25), ndof )*coef33
    endif
    ! --- THERMAL STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(12) .and. associated(tnstrain) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(12), ndof )
    endif
    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(30), ndof )
    endif
    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(31), ndof )
    endif
    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(32), ndof )
    endif
    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(33), ndof )
    endif
    ! --- CONTACT NORMAL TRACTION @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(36) .and. associated(fstrSOLID%CONT_NTRAC) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(36), ndof )
    endif
    ! --- CONTACT FRICTION TRACTION @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(37) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(37), ndof )
    endif
    ! --- NODE ID @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(38) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(38), ndof )
    endif
    ! --- TEMPERATURE @node
    !if( fstrSOLID%output_ctrl(4)%outinfo%on(41) .and. associated(fstrSOLID%CONT_FTRAC) ) then
    !  ncomp = ncomp + 1
    !  nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(41), ndof )
    !endif

    ! --- STRAIN @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(6) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(6), ndof )
    endif
    ! --- STRESS @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(7) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(7), ndof )
    endif
    ! --- MISES @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(8) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(8), ndof )
    endif
    ! --- Principal Stress @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(20) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(20), ndof )
    endif
    ! --- Principal Strain @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(22) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(22), ndof )
    endif
    ! --- Principal Stress Vector @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(24) ) then
      ecomp = ecomp + 3
      eitem = eitem + 3*n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(24), ndof )
    endif
    ! --- Principal Strain Vector @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(26) ) then
      ecomp = ecomp + 3
      eitem = eitem + 3*n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(26), ndof )
    endif
    ! --- PLASTIC STRAIN @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(43) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(43), ndof )
    endif
    ! --- MATERIAL @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(34) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(34), ndof )
    endif
    ! --- ELEM ID @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(39) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(39), ndof )
    endif
    ! --- SECTION ID @element
    if( fstrSOLID%output_ctrl(4)%outinfo%on(40) ) then
      ecomp = ecomp + 1
      eitem = eitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(40), ndof )
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fstrRESULT%ng_component = gcomp
    fstrRESULT%nn_component = ncomp
    fstrRESULT%ne_component = ecomp
    allocate( fstrRESULT%ng_dof(gcomp) )
    allocate( fstrRESULT%global_label(gcomp) )
    allocate( fstrRESULT%global_val_item(gitem) )
    allocate( fstrRESULT%nn_dof(ncomp) )
    allocate( fstrRESULT%node_label(ncomp) )
    allocate( fstrRESULT%node_val_item(nitem*hecMESH%n_node) )
    allocate( fstrRESULT%ne_dof(ecomp) )
    allocate( fstrRESULT%elem_label(ecomp) )
    allocate( fstrRESULT%elem_val_item(eitem*hecMESH%n_elem) )
    ncomp = 0
    iitem = 0
    ecomp = 0
    jitem = 0

    ! --- TIME
    fstrRESULT%ng_dof(1) = 1
    fstrRESULT%global_label(1) = "TOTALTIME"
    fstrRESULT%global_val_item(1) = time

    ! --- DISPLACEMENT
    if (fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
      if(ndof == 4) then
        ! DIPLACEMENT
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'VELOCITY'
        do i = 1, hecMESH%n_node
          do j = 1, 3
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrDYNAMIC%DISP(4*(i-1)+j,idx)
          enddo
        enddo
        iitem = iitem + nn
        ! PRESSURE
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 1 )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'PRESSURE'
        do i = 1, hecMESH%n_node
          fstrRESULT%node_val_item(nitem*(i-1)+1+iitem) = fstrDYNAMIC%DISP(4*i,idx)
        enddo
        iitem = iitem + nn
      else if(ndof == 6) then
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'DISPLACEMENT'
        allocate( unode(3*hecMESH%n_node) )
        unode = 0.0d0
        if( is_dynamic ) then
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrDYNAMIC%DISP((i-1)*ndof+1:(i-1)*ndof+3, idx)
          enddo
        else
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrSOLID%unode((i-1)*ndof+1:(i-1)*ndof+3)
          enddo
        endif
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = unode(nn*(i-1)+j)
          enddo
        enddo
        deallocate( unode )
        iitem = iitem + nn

      else
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'DISPLACEMENT'
        allocate( unode(ndof*hecMESH%n_node) )
        unode = 0.0d0
        if( is_dynamic ) then
          unode(:) = fstrDYNAMIC%DISP(:,idx)
        else
          unode(:) = fstrSOLID%unode(:)
        endif
        if(is_33beam == 1)then
          call fstr_reorder_node_beam(fstrSOLID, hecMESH, unode)
        endif
        if(is_33shell == 1)then
          call fstr_reorder_node_shell(fstrSOLID, hecMESH, unode)
        endif
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = unode(nn*(i-1)+j)
          enddo
        enddo
        deallocate( unode )
        iitem = iitem + nn
      endif
    endif

    ! --- VELOCITY
    if( is_dynamic .and. fstrSOLID%output_ctrl(4)%outinfo%on(15) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(15), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'VELOCITY'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrDYNAMIC%VEL(nn*(i-1)+j,idx)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- ACCELERATION
    if( is_dynamic .and. fstrSOLID%output_ctrl(4)%outinfo%on(16) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(16), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'ACCELERATION'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrDYNAMIC%ACC(nn*(i-1)+j,idx)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- ROTATION
    if( fstrSOLID%output_ctrl(4)%outinfo%on(18) ) then

      if(ndof == 6) then
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'ROTATION'
        allocate( unode(3*hecMESH%n_node) )
        unode = 0.0d0
        if( is_dynamic ) then
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrDYNAMIC%DISP((i-1)*ndof+4:(i-1)*ndof+6, idx)
          enddo
        else
          do i=1, hecMESH%n_node
            unode((i-1)*3+1:(i-1)*3+3) = fstrSOLID%unode((i-1)*ndof+4:(i-1)*ndof+6)
          enddo
        endif
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = unode(nn*(i-1)+j)
          enddo
        enddo
        deallocate( unode )
        iitem = iitem + nn
      else
        if ( is_33shell == 1) then
          ncomp = ncomp + 1
          nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
          fstrRESULT%nn_dof(ncomp) = nn
          fstrRESULT%node_label(ncomp) = 'ROTATION'
          allocate( unode(ndof*hecMESH%n_node) )
          unode = 0.0d0
          call fstr_reorder_rot_shell(fstrSOLID, hecMESH, unode)
          do i = 1, hecMESH%n_node
            do j = 1, nn
              fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = unode(nn*(i-1)+j)
            enddo
          enddo
          deallocate( unode )
          iitem = iitem + nn
        end if
      end if
    endif

    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(4)%outinfo%on(2) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(2), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'REACTION_FORCE'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%REACTION(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(is_33shell == 1 .or. ndof == 6)then
      call fstr_make_result_main( hecMESH, fstrSOLID, fstrRESULT, &
        & fstrSOLID%SHELL, nitem, iitem, ncomp, eitem, jitem, ecomp, 1, "      " )
    else
      call fstr_make_result_main( hecMESH, fstrSOLID, fstrRESULT, &
        & fstrSOLID%SOLID, nitem, iitem, ncomp, eitem, jitem, ecomp, 1, "      " )
    endif

    !laminated shell
    if( associated(fstrSOLID%SHELL) .and. fstrSOLID%output_ctrl(4)%outinfo%on(27) .and. is_33shell == 1 ) then
      allocate(clyr(2*ntot_lyr))
      do i=1,ntot_lyr
        write(cnum,"(i0)")i
        clyr(2*i-1)="_L"//trim(cnum)//"+"
        clyr(2*i  )="_L"//trim(cnum)//"-"
      enddo
      do i=1,ntot_lyr
        call fstr_make_result_main( hecMESH, fstrSOLID, fstrRESULT, &
          & fstrSOLID%SHELL%LAYER(i)%PLUS,  nitem, iitem, ncomp, eitem, jitem, ecomp, i+1, clyr(2*i-1) )
        call fstr_make_result_main( hecMESH, fstrSOLID, fstrRESULT, &
          & fstrSOLID%SHELL%LAYER(i)%MINUS, nitem, iitem, ncomp, eitem, jitem, ecomp, i+1, clyr(2*i  ) )
      enddo
      deallocate(clyr)
    endif

    ! --- THERMAL STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(12) .and. associated(tnstrain) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(12), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'THERMAL_NodalSTRAIN'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = tnstrain(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(30), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_NFORCE'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_NFORCE(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(31), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_FRICTION'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_FRIC(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(32), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_RELVEL'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_RELVEL(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(33), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_STATE'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_STATE(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT NORMAL TRACTION @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(36) .and. associated(fstrSOLID%CONT_NTRAC) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(36), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_NTRACTION'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_NTRAC(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- CONTACT FRICTION TRACTION @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(37) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(37), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'CONTACT_FTRACTION'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%CONT_FTRAC(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- NODE ID @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(38) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(38), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NODE_ID'
      do i = 1, hecMESH%n_node
        fstrRESULT%node_val_item(nitem*(i-1)+1+iitem) = hecMESH%global_node_ID(i)
      enddo
      iitem = iitem + nn
    endif

    ! --- TEMPERATURE @node
    !if( fstrSOLID%output_ctrl(4)%outinfo%on(41) .and. associated(fstrSOLID%CONT_FTRAC) ) then
    !  ncomp = ncomp + 1
    !  nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(41), ndof )
    !  fstrRESULT%nn_dof(ncomp) = nn
    !  fstrRESULT%node_label(ncomp) = 'TEMPERATURE'
    !  do i = 1, hecMESH%n_node
    !    fstrRESULT%node_val_item(nitem*(i-1)+1+iitem) = fstrSOLID%global_node_ID(i)
    !  enddo
    !  iitem = iitem + nn
    !endif


    ! --- MATERIAL @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(34)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(34), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'MATERIAL_ID'
      do i = 1, hecMESH%n_elem
        j = hecMESH%section_ID(i)
        fstrRESULT%elem_val_item(eitem*(i-1)+1+jitem) = hecMESH%section%sect_mat_ID_item(j)
      enddo
      jitem = jitem + nn
    endif

    ! --- ELEM ID @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(39)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(39), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ELEM_ID'
      do i = 1, hecMESH%n_elem
        fstrRESULT%elem_val_item(eitem*(i-1)+1+jitem) = hecMESH%global_elem_ID(i)
      enddo
      jitem = jitem + nn
    endif

    ! --- SECTION ID @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(40)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(40), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'SECTION_ID'
      do i = 1, hecMESH%n_elem
        fstrRESULT%elem_val_item(eitem*(i-1)+1+jitem) = hecMESH%section_ID(i)
      enddo
      jitem = jitem + nn
    endif
  end subroutine fstr_make_result

  subroutine fstr_make_result_main( hecMESH, fstrSOLID, fstrRESULT, RES, nitem, &
     &                              iitem, ncomp, eitem, jitem, ecomp, nlyr, clyr )
    use m_fstr
    use m_out
    use m_static_lib
    use mMaterial
    use hecmw_util

    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_result_data):: fstrRESULT
    type (fstr_solid_physic_val) :: RES
    integer(kind=kint)        :: istep, flag
    integer(kind=kint)        :: n_lyr, cid

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname
    character(len=6)                :: clyr
    character(len=4)                :: cnum
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, eitem, nn, mm, ngauss, it
    integer(kind=kint) :: iitem, ncomp, jitem, ecomp, nlyr

    ndof = hecMESH%n_dof

    ! --- STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(3)) then
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(3), ndof )
      ncomp = ncomp + 1
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalSTRAIN'//trim(clyr)
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%STRAIN(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- STRESS @node
    if(fstrSOLID%output_ctrl(4)%outinfo%on(4)) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(4), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalSTRESS'//trim(clyr)
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%STRESS((nn)*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- MISES @node
    if(fstrSOLID%output_ctrl(4)%outinfo%on(5)) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(5), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalMISES'//trim(clyr)
      do i = 1, hecMESH%n_node
        fstrRESULT%node_val_item(nitem*(i-1)+1+iitem) = RES%MISES(i)
      enddo
      iitem = iitem + nn
    endif

    ! --- Princ STRESS @node
    if(fstrSOLID%output_ctrl(4)%outinfo%on(19)) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(19), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalPrincipalSTRESS'//trim(clyr)
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%PSTRESS((nn)*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- Princ STRESS Vector @node
    if(fstrSOLID%output_ctrl(4)%outinfo%on(23)) then
      do k=1,3
        write(cnum, '(i0)') k
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(23), ndof )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'NodalPrincipalSTRESSVector'//trim(cnum)//trim(clyr)
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%PSTRESS_VECT((nn)*(i-1)+j,k)
          enddo
        enddo
        iitem = iitem + nn
      end do
    endif

    ! --- Princ STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(21)) then
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(21), ndof )
      ncomp = ncomp + 1
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalPrincipalSTRAIN'//trim(clyr)
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%PSTRAIN(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif

    ! --- Princ STRAIN Vector @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(25)) then
      do k=1,3
        write(cnum, '(i0)') k
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(25), ndof )
        ncomp = ncomp + 1
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'NodalPrincipalSTRAINVector'//trim(cnum)//trim(clyr)
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = RES%PSTRAIN_VECT(nn*(i-1)+j,k)
          enddo
        enddo
        iitem = iitem + nn
      enddo
    endif

    ! --- STRAIN @elem
    if( fstrSOLID%output_ctrl(4)%outinfo%on(6)) then
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(6), ndof )
      ecomp = ecomp + 1
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalSTRAIN'
      do i = 1, hecMESH%n_elem
        do j = 1, nn
          fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%ESTRAIN(nn*(i-1)+j)
        enddo
      enddo
      jitem = jitem + nn
    endif

    ! --- STRESS @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(7)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(7), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalSTRESS'
      do i = 1, hecMESH%n_elem
        do j = 1, nn
          fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%ESTRESS((nn)*(i-1)+j)
        enddo
      enddo
      jitem = jitem + nn
    endif

    ! --- MISES @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(8)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(8), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalMISES'
      do i = 1, hecMESH%n_elem
        fstrRESULT%elem_val_item(eitem*(i-1)+1+jitem) = RES%EMISES(i)
      enddo
      jitem = jitem + nn
    endif

    ! --- Principal_STRESS @element
    if(fstrSOLID%output_ctrl(4)%outinfo%on(20)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(20), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalPrincipalSTRESS'
      do i = 1, hecMESH%n_elem
        do j = 1, nn
          fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%EPSTRESS((nn)*(i-1)+j)
        enddo
      enddo
      jitem = jitem + nn
    endif

    ! --- Principal_STRAIN @element
    if(fstrSOLID%output_ctrl(4)%outinfo%on(22)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(22), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalPrincipalSTRAIN'
      do i = 1, hecMESH%n_elem
        do j = 1, nn
          fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%EPSTRAIN((nn)*(i-1)+j)
        enddo
      enddo
      jitem = jitem + nn
    endif

    ! --- ELEM PRINC STRESS VECTOR
    if(fstrSOLID%output_ctrl(4)%outinfo%on(24)) then
      do k = 1, 3
        write(cnum,'(i0)')k
        ecomp = ecomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(24), ndof )
        fstrRESULT%ne_dof(ecomp) = nn
        fstrRESULT%elem_label(ecomp) = 'ElementalPrincipalSTRESSVector'//trim(cnum)
        do i = 1, hecMESH%n_elem
          do j = 1, nn
            fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%EPSTRESS_VECT((nn)*(i-1)+j,k)
          enddo
        enddo
        jitem = jitem + nn
      enddo
    endif

    ! --- ELEM PRINC STRAIN VECTOR
    if(fstrSOLID%output_ctrl(4)%outinfo%on(26)) then
      do k = 1, 3
        write(cnum,'(i0)')k
        ecomp = ecomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(26), ndof )
        fstrRESULT%ne_dof(ecomp) = nn
        fstrRESULT%elem_label(ecomp) = 'ElementalPrincipalSTRAINVector'//trim(cnum)
        do i = 1, hecMESH%n_elem
          do j = 1, nn
            fstrRESULT%elem_val_item(eitem*(i-1)+j+jitem) = RES%EPSTRAIN_VECT((nn)*(i-1)+j,k)
          enddo
        enddo
        jitem = jitem + nn
      enddo
    endif

    ! --- PLSTRAIN @elem
    if(fstrSOLID%output_ctrl(4)%outinfo%on(43)) then
      ecomp = ecomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(43), ndof )
      fstrRESULT%ne_dof(ecomp) = nn
      fstrRESULT%elem_label(ecomp) = 'ElementalPLSTRAIN'
      do i = 1, hecMESH%n_elem
        RES%EPLSTRAIN(i) = 0.d0
        do j = 1, size(fstrSOLID%elements(i)%gausses) 
          RES%EPLSTRAIN(i) = RES%EPLSTRAIN(i) + fstrSOLID%elements(i)%gausses(j)%plstrain
        enddo
        RES%EPLSTRAIN(i) = RES%EPLSTRAIN(i) / size(fstrSOLID%elements(i)%gausses)
        fstrRESULT%elem_val_item(eitem*(i-1)+1+jitem) = RES%EPLSTRAIN(i)
      enddo
      jitem = jitem + nn
    endif

  end subroutine fstr_make_result_main

  subroutine fstr_reorder_node_shell(fstrSOLID, hecMESH, unode)
    use m_fstr
    use m_out
    use m_static_lib

    implicit none
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j, k, itype, is, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, n1, n2
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 781)then
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, 4
            n1 = hecMESH%elem_node_item(jS+j  )
            n2 = hecMESH%elem_node_item(jS+j+4)
            unode(3*n2-2) = unode(3*n1-2)
            unode(3*n2-1) = unode(3*n1-1)
            unode(3*n2  ) = unode(3*n1  )
          enddo
        enddo
      elseif(ic_type == 761)then
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, 3
            n1 = hecMESH%elem_node_item(jS+j  )
            n2 = hecMESH%elem_node_item(jS+j+3)
            unode(3*n2-2) = unode(3*n1-2)
            unode(3*n2-1) = unode(3*n1-1)
            unode(3*n2  ) = unode(3*n1  )
          enddo
        enddo
      endif
    enddo

  end subroutine fstr_reorder_node_shell

  subroutine fstr_reorder_rot_shell(fstrSOLID, hecMESH, unode)
    use m_fstr
    use m_out
    use m_static_lib

    implicit none
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j, k, itype, is, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, n1, n2
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 781)then
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, 4
            n1 = hecMESH%elem_node_item(jS+j)
            n2 = hecMESH%elem_node_item(jS+j+4)
            unode(3*n1-2) = fstrSOLID%unode(3*n2-2)
            unode(3*n1-1) = fstrSOLID%unode(3*n2-1)
            unode(3*n1  ) = fstrSOLID%unode(3*n2  )
            unode(3*n2-2) = fstrSOLID%unode(3*n2-2)
            unode(3*n2-1) = fstrSOLID%unode(3*n2-1)
            unode(3*n2  ) = fstrSOLID%unode(3*n2  )
          enddo
        enddo
      elseif(ic_type == 761)then
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, 3
            n1 = hecMESH%elem_node_item(jS+j)
            n2 = hecMESH%elem_node_item(jS+j+3)

            unode(3*n1-2) = fstrSOLID%unode(3*n2-2)
            unode(3*n1-1) = fstrSOLID%unode(3*n2-1)
            unode(3*n1  ) = fstrSOLID%unode(3*n2  )
            unode(3*n2-2) = fstrSOLID%unode(3*n2-2)
            unode(3*n2-1) = fstrSOLID%unode(3*n2-1)
            unode(3*n2  ) = fstrSOLID%unode(3*n2  )
          enddo
        enddo
      endif
    enddo

  end subroutine fstr_reorder_rot_shell

  subroutine fstr_reorder_node_beam(fstrSOLID, hecMESH, unode)
    use m_fstr
    use m_out
    use m_static_lib

    implicit none
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j, k, itype, is, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, a, b
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 641)then
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, 2
            a = hecMESH%elem_node_item(jS+j)
            b = hecMESH%elem_node_item(jS+j+2)
            unode(3*b-2) = unode(3*a-2)
            unode(3*b-1) = unode(3*a-1)
            unode(3*b  ) = unode(3*a  )
          enddo
        enddo
      endif
    enddo

  end subroutine fstr_reorder_node_beam

  subroutine setup_contact_output_variables( hecMESH, fstrSOLID, phase )
    use m_fstr
    use hecmw_util
    use mContact
    implicit none
    type(hecmwST_local_mesh), intent(in)  :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    integer(kind=kint), intent(in)        :: phase !< -1:clear,3:result,4:vis

    integer(kind=kint), parameter :: nval = 10
    logical, save :: updated(nval) = .false.
    integer(kind=kint) :: ndof, i
    real(kind=kreal) :: area

    ndof = hecMESH%n_dof

    if( phase == -1 ) then
      updated(1:nval) = .false.
      return
    else
      if( phase /= 3 .and. phase /= 4 ) return !irregular case
    end if

    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      if( paraContactFlag .and. .not. updated(1)) then
        call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_NFORCE,1)
      end if
      updated(1) = .true.
    endif

    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      if( paraContactFlag .and. .not. updated(2)) then
        call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_FRIC,1)
      end if
      updated(2) = .true.
    endif

    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      if( paraContactFlag .and. .not. updated(3)) then
        call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_RELVEL,1)
      end if
      updated(3) = .true.
    endif

    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      if( paraContactFlag .and. .not. updated(4)) then
        call fstr_setup_parancon_contactvalue(hecMESH,1,fstrSOLID%CONT_STATE,2)
      end if
      updated(4) = .true.
    endif

    ! --- CONTACT AREA for CONTACT TRACTION
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(36) .or. fstrSOLID%output_ctrl(phase)%outinfo%on(37) ) then
      if( .not. updated(5)) call calc_contact_area( hecMESH, fstrSOLID, 0 )
      ! fstr_setup_parancon_contactvalue is not necessary because
      ! contact area is calculated from original surface group
    end if

    ! --- CONTACT NORMAL TRACTION @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(36) .and. associated(fstrSOLID%CONT_NTRAC) ) then
      if( paraContactFlag .and. .not. updated(6)) then
        if( .not. updated(1)) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_NFORCE,1)
      end if
      fstrSOLID%CONT_NTRAC(:) = 0.d0
      do i=1,hecMESH%nn_internal
        area = fstrSOLID%CONT_AREA(i)
        if( area < 1.d-16 ) cycle
        fstrSOLID%CONT_NTRAC(3*i-2:3*i) = fstrSOLID%CONT_NFORCE(3*i-2:3*i)/area
      end do
      updated(6) = .true.
    endif

    ! --- CONTACT FRICTION TRACTION @node
    if( fstrSOLID%output_ctrl(phase)%outinfo%on(37) .and. associated(fstrSOLID%CONT_FTRAC) ) then
      if( paraContactFlag .and. .not. updated(7)) then
        if( .not. updated(1)) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_FRIC,1)
      end if
      fstrSOLID%CONT_FTRAC(:) = 0.d0
      do i=1,hecMESH%nn_internal
        area = fstrSOLID%CONT_AREA(i)
        if( area < 1.d-16 ) cycle
        fstrSOLID%CONT_FTRAC(3*i-2:3*i) = fstrSOLID%CONT_FRIC(3*i-2:3*i)/area
      end do
      updated(7) = .true.
    endif

  end subroutine

  subroutine fstr_setup_parancon_contactvalue(hecMESH,ndof,vec,vtype)
  use m_fstr
  implicit none
  type(hecmwST_local_mesh), intent(in)      :: hecMESH
  integer(kind=kint), intent(in)            :: ndof
  real(kind=kreal), pointer, intent(inout)  :: vec(:)
  integer(kind=kint), intent(in)            :: vtype !1:value, 2:state
  !
  real(kind=kreal) ::  rhsB
  integer(kind=kint) ::  i,j,N,i0,N_loc,nndof
  integer(kind=kint) :: offset, pid, lid
  integer(kind=kint), allocatable :: displs(:)
  real(kind=kreal), allocatable   :: vec_all(:)

  !
  N_loc = hecMESH%nn_internal
  allocate(displs(0:nprocs))
  displs(:) = 0
  displs(myrank+1) = N_loc
  call hecmw_allreduce_I(hecMESH, displs, nprocs+1, hecmw_sum)
  do i=1,nprocs
    displs(i) = displs(i-1) + displs(i)
  end do
  offset = displs(myrank)
  N = displs(nprocs)

  allocate(vec_all(ndof*N))

  if( vtype == 1 ) then
    vec_all(:) = 0.d0
    do i= hecMESH%nn_internal+1,hecMESH%n_node
      pid = hecMESH%node_ID(i*2)
      lid = hecMESH%node_ID(i*2-1)
      i0 = (displs(pid) + (lid-1))*ndof
      vec_all(i0+1:i0+ndof) = vec((i-1)*ndof+1:i*ndof)
      vec((i-1)*ndof+1:i*ndof) = 0.d0
    enddo

    call hecmw_allreduce_R(hecMESH, vec_all, N*ndof, hecmw_sum)

    do i=1,ndof*N_loc
      vec(i) = vec(i) + vec_all(offset*ndof+i)
    end do
  else if( vtype == 2 ) then
    vec_all(:) = -1000.d0
    do i= hecMESH%nn_internal+1,hecMESH%n_node
      if( vec(i) == 0.d0 ) cycle
      pid = hecMESH%node_ID(i*2)
      lid = hecMESH%node_ID(i*2-1)
      i0 = displs(pid) + lid
      vec_all(i0) = vec(i)
    enddo

    call hecmw_allreduce_R(hecMESH, vec_all, N, hecmw_max)

    do i=1,N_loc
      if( vec_all(offset+i) == -1000.d0 ) cycle
      if( vec(i) < vec_all(offset+i) ) vec(i) = vec_all(offset+i)
    end do
  end if

  deallocate(displs,vec_all)
  end subroutine


end module m_make_result
