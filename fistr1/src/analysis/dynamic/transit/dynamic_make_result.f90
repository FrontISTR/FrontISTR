!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to prepare output of dynamic analysis
module m_dynamic_make_result
  use m_static_make_result
  implicit none
contains

  !C***
  !>  OUTPUT result file for dynamic analysis
  !C***
  subroutine fstr_write_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, maxstep, istep, time )
    use m_fstr
    use m_out
    use m_static_lib
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type (fstr_dynamic)       :: fstrDYNAMIC
    integer(kind=kint)        :: maxstep, istep
    real(kind=kreal)          :: time        !< current time
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN)    :: comment
    character(len=HECMW_NAME_LEN)   :: label, s, nameID
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, idx, ngauss
    integer(kind=kint) :: n_lyr, ntot_lyr, is_33shell, is_33beam
    real(kind=kreal), allocatable   :: work(:), unode(:)

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN

    if( fstrDYNAMIC%idx_eqa==1 .and. istep>0 ) then
      idx = 2
    else
      idx = 1
    endif

    ndof = hecMESH%n_dof
    nn = hecMESH%n_node
    if( hecMESH%n_elem>hecMESH%n_node ) nn = hecMESH%n_elem
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==4 ) mdof = 6
    if( ndof==6 ) mdof = 6

    ntot_lyr   = fstrSOLID%max_lyr
    is_33shell = fstrSOLID%is_33shell
    is_33beam  = fstrSOLID%is_33beam

    nn = nn * mdof
    allocate( work(nn) )

    ! --- INITIALIZE
    header = '*fstrresult'
    comment = 'dynamic_result'
    call hecmw_result_init( hecMESH, istep, header, comment )

    ! --- TIME
    id = 3 !global data
    label = 'TOTALTIME'
    work(1) = time
    call hecmw_result_add( id, 1, label, work )

    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(3)%outinfo%on(1) ) then
      if(ndof /= 4) then
        id = 1
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), ndof )
        allocate( unode(hecMESH%n_dof*hecMESH%n_node) )
        unode = 0.0d0
        unode(:) = fstrDYNAMIC%DISP(:,idx)
        label = 'DISPLACEMENT'
        if(is_33beam == 1)then
          call fstr_reorder_node_beam(fstrSOLID, hecMESH, unode)
        endif
        if(is_33shell == 1)then
          call fstr_reorder_node_shell(fstrSOLID, hecMESH, unode)
        endif
        call hecmw_result_add( id, nitem, label, unode )
        deallocate( unode )
      else
        id = 1
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
      endif
    endif
    ! --- VELOCITY
    if( fstrSOLID%output_ctrl(3)%outinfo%on(15) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(15), ndof )
      label = 'VELOCITY'
      call hecmw_result_add( id, nitem, label, fstrDYNAMIC%VEL(:,idx) )
    endif
    ! --- ACCELERATION
    if( fstrSOLID%output_ctrl(3)%outinfo%on(16) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(16), ndof )
      label = 'ACCELERATION'
      call hecmw_result_add( id, nitem, label, fstrDYNAMIC%ACC(:,idx) )
    endif
    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(3)%outinfo%on(2) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(2), ndof )
      label = 'REACTION_FORCE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%REACTION )
    endif
    ! --- STRAIN @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(3) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(3), ndof )
      label = 'NodalSTRAIN'
      call hecmw_result_add( id, nitem, label, fstrSOLID%STRAIN )
    endif
    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(4) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(4), ndof )
      label = 'NodalSTRESS'
      call hecmw_result_add( id, nitem, label, fstrSOLID%STRESS )
    endif
    ! --- MISES @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(5) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(5), ndof )
      label = 'NodalMISES'
      call hecmw_result_add( id, nitem, label, fstrSOLID%MISES )
    endif
    ! --- STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(6) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(6), ndof )
      label = 'ElementalSTRAIN'
      call hecmw_result_add( id, nitem, label, fstrSOLID%ESTRAIN )
    endif
    ! --- STRESS @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(7) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(7), ndof )
      label = 'ElementalSTRESS'
      call hecmw_result_add( id, nitem, label, fstrSOLID%ESTRESS )
    endif
    ! --- MISES @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(8) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(8), ndof )
      label = 'ElementalMISES'
      call hecmw_result_add( id, nitem, label, fstrSOLID%EMISES )
    endif
    ! --- STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(9) .and. ndof/=6 ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(9), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRESS',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( k > size(fstrSOLID%elements(i)%gausses) ) then
            do j = 1, nitem
              work(nitem*(i-1)+j) = 0.0D0
            enddo
          else
            do j = 1, nitem
              work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%strain_out(j)
            enddo
          endif
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif
    ! --- STRESS @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(10) .and. ndof/=6 ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(10), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRESS',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( k > size(fstrSOLID%elements(i)%gausses) ) then
            do j = 1, nitem
              work(nitem*(i-1)+j) = 0.0D0
            enddo
          else
            do j = 1, nitem
              work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%stress_out(j)
            enddo
          endif
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif
    ! --- PLASTIC STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(11) .and. fstrSOLID%StaticType/=3 ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(11), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRESS',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          if( k > size(fstrSOLID%elements(i)%gausses) ) then
            do j = 1, nitem
              work(nitem*(i-1)+j) = 0.0D0
            enddo
          else
            do j = 1, nitem
              work(i) = fstrSOLID%elements(i)%gausses(k)%plstrain
            enddo
          endif
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif
    ! --- THERMAL STRAIN @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(12) .and. associated(tnstrain) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(12), ndof )
      label = 'THERMAL_NodalSTRAIN'
      call hecmw_result_add( id, nitem, label, tnstrain )
    endif
    ! --- THERMAL STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(13) .and. associated(testrain) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(13), ndof )
      label = 'THERMAL_ElementalSTRAIN'
      call hecmw_result_add( id, nitem, label, testrain )
    endif
    ! --- THERMAL STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(14) .and. associated(testrain) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(14), ndof )
      ngauss = fstrSOLID%maxn_gauss
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'THERMAL_GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          do j = 1, nitem
            !                work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%tstrain(j)
          enddo
        enddo
        call hecmw_result_add( id, nitem, label, work )
      enddo
    endif
    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_NFORCE,1)
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(30), ndof )
      label = 'CONTACT_NFORCE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_NFORCE )
    endif

    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_FRIC,1)
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(31), ndof )
      label = 'CONTACT_FRICTION'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_FRIC )
    endif

    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_RELVEL,1)
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(32), ndof )
      label = 'CONTACT_RELVEL'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_RELVEL )
    endif

    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,1,fstrSOLID%CONT_STATE,2)
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(33), ndof )
      label = 'CONTACT_STATE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%CONT_STATE )
    endif

    ! --- WRITE
    nameID = 'fstrRES'
    call hecmw_result_write_by_name( nameID )
    ! --- FINALIZE
    call hecmw_result_finalize

    deallocate( work )
  end subroutine fstr_write_dynamic_result

  !C***
  !>  MAKE RESULT for dynamic analysis (WITHOUT ELEMENTAL RESULTS)
  !C***
  subroutine fstr_make_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, fstrRESULT )
    use m_fstr
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_solid)          :: fstrSOLID
    type(fstr_dynamic)        :: fstrDYNAMIC
    type(hecmwST_result_data) :: fstrRESULT
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)

    integer(kind=kint) :: i, j, ndof, mdof, gcomp, gitem, ncomp, nitem, iitem, ecomp, eitem, jitem, nn, idx

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN

    if( fstrDYNAMIC%idx_eqa==1 .and. fstrDYNAMIC%i_step>0 ) then
      idx = 2
    else
      idx = 1
    endif

    ndof = hecMESH%n_dof
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==4 ) mdof = 6
    if( ndof==6 ) mdof = 6

    call hecmw_nullify_result_data( fstrRESULT )
    gcomp = 0
    gitem = 0
    ncomp = 0
    nitem = 0
    ecomp = 0
    eitem = 0

    ! --- TIME
    gcomp = gcomp + 1
    gitem = gitem + 1

    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
      if(ndof /= 4) then
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
      else
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 3 )
        ncomp = ncomp + 1
        nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), 1 )
      endif
    endif
    ! --- VELOCITY
    if( fstrSOLID%output_ctrl(4)%outinfo%on(15) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(15), ndof )
    endif
    ! --- ACCELERATION
    if( fstrSOLID%output_ctrl(4)%outinfo%on(16) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(16), ndof )
    endif
    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(4)%outinfo%on(2) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(2), ndof )
    endif
    ! --- STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(3) .and. ndof==6 ) then
      ncomp = ncomp + 2
      nitem = nitem + mdof
    else if( fstrSOLID%output_ctrl(4)%outinfo%on(3) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(3), ndof )
    endif
    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(4) .and. ndof==6 ) then
      ncomp = ncomp + 2
      nitem = nitem + mdof
    else if( fstrSOLID%output_ctrl(4)%outinfo%on(4) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(4), ndof )
    endif
    ! --- MISES @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(5) .and. ndof==6 ) then
      ncomp = ncomp + 2
      nitem = nitem + 2
    else if( fstrSOLID%output_ctrl(4)%outinfo%on(5) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(5), ndof )
    endif
    ! --- THERMAL STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(12) .and. associated(tnstrain) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(12), ndof )
    endif
    ! --- CONTACT NORMAL FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(30), ndof )
    endif
    ! --- CONTACT FRICTION FORCE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(31), ndof )
    endif
    ! --- CONTACT RELATIVE VELOCITY @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(32), ndof )
    endif
    ! --- CONTACT STATE @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      ncomp = ncomp + 1
      nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(33), ndof )
    endif

    fstrRESULT%ng_component = gcomp
    fstrRESULT%nn_component = ncomp
    fstrRESULT%ne_component = ecomp
    allocate( fstrRESULT%ng_dof(gcomp) )
    allocate( fstrRESULT%global_label(gcomp) )
    allocate( fstrRESULT%global_val_item(gitem) )
    allocate( fstrRESULT%nn_dof(ncomp) )
    allocate( fstrRESULT%node_label(ncomp) )
    allocate( fstrRESULT%node_val_item(nitem*hecMESH%n_node) )
    !allocate( fstrRESULT%ne_dof(ecomp) )
    !allocate( fstrRESULT%elem_label(ecomp) )
    !allocate( fstrRESULT%elem_val_item(eitem*hecMESH%n_elem) )
    ncomp = 0
    iitem = 0
    ecomp = 0
    jitem = 0

    ! --- TIME
    fstrRESULT%ng_dof(1) = 1
    fstrRESULT%global_label(1) = "TOTALTIME"
    fstrRESULT%global_val_item(1) = fstrDYNAMIC%t_curr

    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
      if(ndof /= 4) then
        ncomp = ncomp + 1
        nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
        fstrRESULT%nn_dof(ncomp) = nn
        fstrRESULT%node_label(ncomp) = 'DISPLACEMENT'
        do i = 1, hecMESH%n_node
          do j = 1, nn
            fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrDYNAMIC%DISP(nn*(i-1)+j,idx)
          enddo
        enddo
        iitem = iitem + nn
      else
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
      endif
    endif
    ! --- VELOCITY
    if( fstrSOLID%output_ctrl(4)%outinfo%on(15) ) then
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
    if( fstrSOLID%output_ctrl(4)%outinfo%on(16) ) then
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
    ! --- STRAIN @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(3) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(3), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalSTRAIN'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%STRAIN(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif
    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(4) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(4), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalSTRESS'
      do i = 1, hecMESH%n_node
        do j = 1, nn
          fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%STRESS(nn*(i-1)+j)
        enddo
      enddo
      iitem = iitem + nn
    endif
    ! --- MISES @node
    if( fstrSOLID%output_ctrl(4)%outinfo%on(5) ) then
      ncomp = ncomp + 1
      nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(5), ndof )
      fstrRESULT%nn_dof(ncomp) = nn
      fstrRESULT%node_label(ncomp) = 'NodalMISES'
      do i = 1, hecMESH%n_node
        fstrRESULT%node_val_item(nitem*(i-1)+1+iitem) = fstrSOLID%MISES(i)
      enddo
      iitem = iitem + nn
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
    if( fstrSOLID%output_ctrl(3)%outinfo%on(30) .and. associated(fstrSOLID%CONT_NFORCE) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_NFORCE,1)
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
    if( fstrSOLID%output_ctrl(3)%outinfo%on(31) .and. associated(fstrSOLID%CONT_FRIC) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_FRIC,1)
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
    if( fstrSOLID%output_ctrl(3)%outinfo%on(32) .and. associated(fstrSOLID%CONT_RELVEL) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,ndof,fstrSOLID%CONT_RELVEL,1)
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
    if( fstrSOLID%output_ctrl(3)%outinfo%on(33) .and. associated(fstrSOLID%CONT_STATE) ) then
      if( paraContactFlag ) call fstr_setup_parancon_contactvalue(hecMESH,1,fstrSOLID%CONT_STATE,2)
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

  end subroutine fstr_make_dynamic_result

end module m_dynamic_make_result
