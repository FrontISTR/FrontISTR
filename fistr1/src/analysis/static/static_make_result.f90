!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.7                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Suemitsu(Advancesoft)                       !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provide a function to prepare output of static analysis
module m_static_make_result
  private

  public:: fstr_write_static_result
  public:: fstr_make_static_result
  public:: fstr_reorder_node_shell
  public:: fstr_reorder_rot_shell
  public:: fstr_reorder_node_beam
  
  
  contains

!C***
!>  OUTPUT result file for static analysis
!C***
  subroutine fstr_write_static_result( hecMESH, fstrSOLID, maxstep, istep, flag )
    use m_fstr
    use m_out
    use m_static_lib
    use mMaterial
    use hecmw_util
    
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    integer(kind=kint)        :: maxstep, istep, flag
    integer(kind=kint) :: n_lyr, ntot_lyr, tmp, is_33shell, is_33beam, cid
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, mm, ngauss, it
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)
    real(kind=kreal), allocatable   :: work(:), unode(:), rnode(:)
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname, cnum
    character(len=6), allocatable   :: clyr(:)

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN

    ndof = hecMESH%n_dof
    mm   = hecMESH%n_node
    if( hecMESH%n_elem > hecMESH%n_node ) mm = hecMESH%n_elem
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==6 ) mdof = 6

    ntot_lyr   = fstrSOLID%max_lyr
    is_33shell = fstrSOLID%is_33shell
    is_33beam  = fstrSOLID%is_33beam

    nn = mm * mdof
    
    allocate( work(nn) )

    ! --- INITIALIZE
    header = '*fstrresult'
    call hecmw_result_init( hecMESH, maxstep, istep, header )

    ! --- DISPLACEMENT
    if( fstrSOLID%output_ctrl(3)%outinfo%on(1)) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), ndof )
      allocate( unode(hecMESH%n_node*ndof) )
      unode = 0.0d0
      unode = fstrSOLID%unode
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

    ! --- ROTATION
    if (fstrSOLID%output_ctrl(3)%outinfo%on(18)) then
      if ( is_33shell == 1) then
        id = 1
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(1), ndof )
        label = 'ROTATION'
        allocate( rnode(hecMESH%n_node*ndof) )
        rnode = 0.0d0
        call fstr_reorder_rot_shell(fstrSOLID, hecMESH, rnode)
        call hecmw_result_add( id, nitem, label, rnode )
        deallocate( rnode )
      end if
    endif

    ! --- REACTION FORCE
    if( fstrSOLID%output_ctrl(3)%outinfo%on(2) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(2), ndof )
      label = 'REACTION_FORCE'
      call hecmw_result_add( id, nitem, label, fstrSOLID%QFORCE )
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(is_33shell == 1 .or. ndof == 6)then
      call fstr_write_static_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL, "      " )
    else
      call fstr_write_static_result_main( hecMESH, fstrSOLID, fstrSOLID%SOLID, "      " )
    endif

    !laminated shell
    if( fstrSOLID%output_ctrl(3)%outinfo%on(27) ) then
      allocate(clyr(2*ntot_lyr))
      do i=1,ntot_lyr
        write(cnum,"(i0)")i
        clyr(2*i-1)="_L"//trim(cnum)//"+"
        clyr(2*i  )="_L"//trim(cnum)//"-"
      enddo
      do i=1,ntot_lyr
        call fstr_write_static_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL%LAYER(i)%PLUS,  clyr(2*i-1) )
        call fstr_write_static_result_main( hecMESH, fstrSOLID, fstrSOLID%SHELL%LAYER(i)%MINUS, clyr(2*i  ) )
      enddo
      deallocate(clyr)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! --- STRAIN @gauss
    if( fstrSOLID%output_ctrl(3)%outinfo%on(9) .and. ndof/=6 ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(9), ndof )
      ngauss = NumOfQuadPoints( hecMESH%elem_type_item(1) )
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
         if( size(fstrSOLID%elements(i)%gausses) .lt. ngauss ) then
           do j = 1, nitem
            work(nitem*(i-1)+j) = 0.0D0
           enddo
         else
           do j = 1, nitem
            work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%strain(j)
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
      ngauss = NumOfQuadPoints( hecMESH%elem_type_item(1) )
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'GaussSTRESS',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
         if( size(fstrSOLID%elements(i)%gausses) .lt. ngauss ) then
           do j = 1, nitem
            work(nitem*(i-1)+j) = 0.0D0
           enddo
         else
           do j = 1, nitem
            work(nitem*(i-1)+j) = fstrSOLID%elements(i)%gausses(k)%stress(j)
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
      ngauss = NumOfQuadPoints( hecMESH%elem_type_item(1) )
      do k = 1, ngauss
        write(s,*) k
        write(label,'(a,a)') 'PLASTIC_GaussSTRAIN',trim(adjustl(s))
        label = adjustl(label)
        do i = 1, hecMESH%n_elem
          work(i) = fstrSOLID%elements(i)%gausses(k)%plstrain
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
      ngauss = NumOfQuadPoints( hecMESH%elem_type_item(1) )
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
  end subroutine fstr_write_static_result

  subroutine fstr_write_static_result_main( hecMESH, fstrSOLID, RES, clyr )
    use m_fstr
    use m_out
    use m_static_lib
    use mMaterial
    use hecmw_util
    
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type (fstr_solid_physic_val) :: RES
    integer(kind=kint)        :: maxstep, istep, flag
    integer(kind=kint)        :: n_lyr, cid

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname
    character(len=6)                :: clyr
    character(len=4)                :: cnum
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, mm, ngauss, it

    ndof = hecMESH%n_dof

    ! --- STRAIN @node
    if (fstrSOLID%output_ctrl(3)%outinfo%on(3)) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(3), ndof )
      label = 'NodalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%STRAIN )
    endif

    ! --- STRESS @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(4) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(4), ndof )
      label = 'NodalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%STRESS )
    endif

    ! --- MISES @node
    if( fstrSOLID%output_ctrl(3)%outinfo%on(5) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(5), ndof )
      label = 'NodalMISES'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%MISES )
    endif
    
    ! --- NODAL PRINC STRESS
    if( fstrSOLID%output_ctrl(3)%outinfo%on(19) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(19), ndof )
      label = 'NodalPrincipalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%PSTRESS )
    endif

    ! --- NODAL PRINC STRAIN
    if( fstrSOLID%output_ctrl(3)%outinfo%on(21) ) then
      id = 1
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(21), ndof )
      label = 'NodalPrincipalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%PSTRAIN )
    endif

    ! --- NODAL PRINC STRESS VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(23) ) then
      id = 1
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(23), ndof )
        label = 'NodalPrincipalSTRESSVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%PSTRESS_VECT(:,k) )
      end do
    endif

    ! --- NODAL PRINC STRAIN VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(25) ) then
      id = 1
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
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(6), ndof )
      label = 'ElementalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%ESTRAIN )
    endif

    ! --- STRESS @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(7) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(7), ndof )
      label = 'ElementalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%ESTRESS )
    endif

    ! --- MISES @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(8)) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(8), ndof )
      label = 'ElementalMISES'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EMISES )
    endif

    ! --- Principal_STRESS @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(20) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(20), ndof )
      label = 'ElementalPrincipalSTRESS'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EPSTRESS )
    endif

    ! --- Principal_STRAIN @element
    if( fstrSOLID%output_ctrl(3)%outinfo%on(22) ) then
      id = 2
      nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(22), ndof )
      label = 'ElementalPrincipalSTRAIN'//trim(clyr)
      call hecmw_result_add( id, nitem, label, RES%EPSTRAIN )
    endif
    
    ! --- ELEM PRINC STRESS VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(24) ) then
      id = 2
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(24), ndof )
        label = 'ElementalPrincipalSTRESSVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%EPSTRESS_VECT(:,k) )
      end do
    endif

    !ELEM PRINC STRAIN VECTOR
    if( fstrSOLID%output_ctrl(3)%outinfo%on(26) ) then
      id = 2
      do k=1,3
        write(cnum,'(i0)')k
        nitem = n_comp_valtype( fstrSOLID%output_ctrl(3)%outinfo%vtype(26), ndof )
        label = 'ElementalPrincipalSTRAINVector'//trim(cnum)//trim(clyr)
        call hecmw_result_add( id, nitem, label, RES%EPSTRAIN_VECT(:,k) )
      end do
    endif

  end subroutine fstr_write_static_result_main

!C***
!>  MAKE RESULT for static analysis (WITHOUT ELEMENTAL RESULTS) --------------------------------------------------------------
!C***
  subroutine fstr_make_static_result( hecMESH, fstrSOLID, fstrRESULT )
    use m_fstr
    use hecmw_util
    
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_result_data):: fstrRESULT
    integer(kind=kint) :: n_lyr, ntot_lyr, it, coef33, is_33shell, is_33beam
    integer(kind=kint) :: i, j, k, ndof, mdof, ncomp, nitem, iitem, ecomp, eitem, jitem, nn, mm
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)
    real(kind=kreal), allocatable   ::unode(:)
    character(len=4) :: cnum
    character(len=6), allocatable   :: clyr(:)

    tnstrain => fstrSOLID%TNSTRAIN
    testrain => fstrSOLID%TESTRAIN

    ntot_lyr   = fstrSOLID%max_lyr
    is_33shell = fstrSOLID%is_33shell
    is_33beam  = fstrSOLID%is_33beam

    mm = hecMESH%n_node
    if( hecMESH%n_elem>hecMESH%n_node ) mm = hecMESH%n_elem

    ndof = hecMESH%n_dof
    if( ndof==2 ) mdof = 3
    if( ndof==3 ) mdof = 6
    if( ndof==6 ) mdof = 6

    if(is_33shell == 1 .and. fstrSOLID%output_ctrl(4)%outinfo%on(27) )then
      coef33 = 1 + 2*ntot_lyr
    else
      coef33 = 1
    endif

    call hecmw_nullify_result_data( fstrRESULT )
    ncomp = 0
    nitem = 0
    ecomp = 0
    eitem = 0

! --- COUNT SUM OF ALL NITEM
! --- DISPLACEMENT
        if( fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
          ncomp = ncomp + 1
          nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
        endif
! --- ROTATION (Only for 781 shell)
        if( fstrSOLID%output_ctrl(4)%outinfo%on(18) .and. is_33shell == 1 ) then
          ncomp = ncomp + 1
          nitem = nitem + n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(18), ndof )
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fstrRESULT%nn_component = ncomp
        fstrRESULT%ne_component = ecomp
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

! --- DISPLACEMENT
        if (fstrSOLID%output_ctrl(4)%outinfo%on(1) ) then
          ncomp = ncomp + 1
          nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(1), ndof )
          fstrRESULT%nn_dof(ncomp) = nn
          fstrRESULT%node_label(ncomp) = 'DISPLACEMENT'
          allocate( unode(ndof*hecMESH%n_node) )
          unode = 0.0d0
          unode(:) = fstrSOLID%unode(:)
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
! --- ROTATION
        if( fstrSOLID%output_ctrl(4)%outinfo%on(18) .and. is_33shell == 1 ) then
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
        endif
! --- REACTION FORCE
        if( fstrSOLID%output_ctrl(4)%outinfo%on(2) ) then
          ncomp = ncomp + 1
          nn = n_comp_valtype( fstrSOLID%output_ctrl(4)%outinfo%vtype(2), ndof )
          fstrRESULT%nn_dof(ncomp) = nn
          fstrRESULT%node_label(ncomp) = 'REACTION_FORCE'
          do i = 1, hecMESH%n_node
            do j = 1, nn
              fstrRESULT%node_val_item(nitem*(i-1)+j+iitem) = fstrSOLID%QFORCE(nn*(i-1)+j)
            enddo
          enddo
          iitem = iitem + nn
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(is_33shell == 1 .or. ndof == 6)then
          call fstr_make_static_result_main( hecMESH, fstrSOLID, fstrRESULT, &
          & fstrSOLID%SHELL, nitem, iitem, ncomp, 1, "      " )
        else
          call fstr_make_static_result_main( hecMESH, fstrSOLID, fstrRESULT, &
          & fstrSOLID%SOLID, nitem, iitem, ncomp, 1, "      " )
        endif

        !laminated shell
        if( fstrSOLID%output_ctrl(3)%outinfo%on(27) .and. is_33shell == 1 ) then
          allocate(clyr(2*ntot_lyr))
          do i=1,ntot_lyr
            write(cnum,"(i0)")i
            clyr(2*i-1)="_L"//trim(cnum)//"+"
            clyr(2*i  )="_L"//trim(cnum)//"-"
          enddo
          do i=1,ntot_lyr
            call fstr_make_static_result_main( hecMESH, fstrSOLID, fstrRESULT, &
            & fstrSOLID%SHELL%LAYER(i)%PLUS,  nitem, iitem, ncomp, i+1, clyr(2*i-1) )
            call fstr_make_static_result_main( hecMESH, fstrSOLID, fstrRESULT, &
            & fstrSOLID%SHELL%LAYER(i)%MINUS, nitem, iitem, ncomp, i+1, clyr(2*i  ) )
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

    end subroutine fstr_make_static_result

  subroutine fstr_make_static_result_main( hecMESH, fstrSOLID, fstrRESULT, RES, nitem, iitem, ncomp, nlyr, clyr )
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
    integer(kind=kint)        :: maxstep, istep, flag
    integer(kind=kint)        :: n_lyr, cid

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: s, label, nameID, addfname
    character(len=6)                :: clyr
    character(len=4)                :: cnum
    integer(kind=kint) :: i, j, k, ndof, mdof, id, nitem, nn, mm, ngauss, it
    integer(kind=kint) :: iitem, ncomp, nlyr

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

  end subroutine fstr_make_static_result_main

  subroutine fstr_reorder_node_shell(fstrSOLID, hecMESH, unode)
    use m_fstr
    use m_out
    use m_static_lib
    
    implicit none
    type (fstr_solid)         :: fstrSOLID
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j, k, itype, iS, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, n1, n2
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 781)then
        do icel = iS, iE
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
        do icel = iS, iE
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
    integer(kind=kint) :: i, j, k, itype, iS, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, n1, n2
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 781)then
        do icel = iS, iE
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
        do icel = iS, iE
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
    integer(kind=kint) :: i, j, k, itype, iS, iE, ic_type, jS, icel
    integer(kind=kint) :: mm, a, b
    real(kind=kreal), allocatable   :: unode(:)

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 641)then
        do icel = iS, iE
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
  
end module m_static_make_result
