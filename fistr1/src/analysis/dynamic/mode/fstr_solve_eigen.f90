!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function to control eigen analysis
module m_fstr_solve_eigen
  contains

!> SOLVE EIGENVALUE PROBLEM
  subroutine fstr_solve_eigen( hecMESH, hecMAT, fstrEIG, fstrSOLID, &
                             & fstrRESULT, fstrPARAM, fstrMAT)
    use hecmw_util
    use m_fstr
    use m_fstr_StiffMatrix
    use m_fstr_AddBC
    use m_fstr_EIG_getamat
    use m_fstr_EIG_lanczos
    use m_fstr_EIG_lanczos_util
    use m_fstr_EIG_matmult
    use m_fstr_EIG_mgs1
    use m_fstr_EIG_output
    use m_fstr_EIG_setMASS
    use m_fstr_EIG_tridiag
    use m_static_lib
    use m_static_make_result
    use m_eigen_lib
    use m_hecmw2fstr_mesh_conv
    use fstr_matrix_con_contact
    use lczeigen

    implicit none

    type (hecmwST_local_mesh ) :: hecMESH
    type (hecmwST_matrix     ) :: hecMAT
    type (fstr_solid         ) :: fstrSOLID
    type (hecmwST_result_data) :: fstrRESULT
    type (fstr_param         ) :: fstrPARAM
    type (fstr_eigen         ) :: fstrEIG
    type (fstrST_matrix_contact_lagrange) :: fstrMAT

    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    integer(kind=kint) :: i, j, k , ii, iii, ik, in, in1, in2, in3, nstep, istep
    integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0, jiter, iiter, kiter
    integer(kind=kint) :: kk, jjiter, ppc
    real(kind=kreal)   :: t1, t2, aalf, tmp, tmp2, gm, gm2, r1, r2, r3, r4, r5, r6

    integer(kind=kint) :: IOUT,IREOR,eITMAX,itype,iS,iE,ic_type,icel,jS,nn, maxItr
    integer(kind=kint), allocatable :: isnode33(:)
    real(kind=kreal),   allocatable :: gmass(:)

    t1 = hecmw_Wtime()

    numnp  = hecMAT%NP
    numn   = hecMAT%N
    NDOF   = hecMESH%n_dof
    ntotal = numnp*NDOF

    fstrSOLID%dunode = 0.d0
    call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, 0.d0, 0.d0 )
    call fstr_AddBC(1,  hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, 2)

    allocate(fstrEIG%mass(ntotal), STAT=ierror)
    call setMASS(fstrSOLID, hecMESH, hecMAT, fstrEIG)

    if( myrank == 0 ) then
      write(IMSG,*) 'fstr_mat_ass: OK'
    endif

    allocate(isnode33(numnp))
    isnode33 = 0

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(hecmw_is_etype_33struct(ic_type))then
        nn = HECMW_get_max_node(ic_type)/2
        do icel = iS, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            ii = hecMESH%elem_node_item(jS+j+nn)
            isnode33(ii)=1
          enddo
        enddo
      endif
    enddo

    allocate( EFILT( ntotal) )

    efilt  = 1.0
    kcount = 0

    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0  = hecMESH%node_group%grp_index(ig-1) + 1
      iE0  = hecMESH%node_group%grp_index(ig  )
      it   = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      itS0 = (it - mod(it,10))/10
      itE0 = mod(it,10)

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        do i = itS0,itE0
          EFILT((in-1)*NDOF+i)=0.0
        enddo
      enddo
    enddo

    call hecmw_update_m_R(hecMESH,EFILT,numnp,NDOF)

    Gtotal  = numnp*NDOF
    Gntotal = numn*NDOF
    novl    = ( numnp - numn )*NDOF

    kcount = 0
    DO I = 1, Gntotal
      IF(EFILT(I).EQ.0) kcount = kcount + 1
    enddo

    call hecmw_allreduce_I1(hecMESH,Gtotal,hecmw_sum)
    call hecmw_allreduce_I1(hecMESH,kcount,hecmw_sum)

    eITMAX  = fstrEIG%lczmax
    ITLIMIT = Gtotal - kcount
    IF(eITMAX.GT.ITLIMIT) THEN
      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*) '*-------------------------------------------*'
        WRITE(IMSG,*) '  WARNING: LCZMAX exceeds system matrix size.'
        WRITE(IMSG,*) '  Resetting LCZMAX to system matrix size.'
        WRITE(IMSG,*) '*-------------------------------------------*'
      endif
      eITMAX = ITLIMIT
    endif

    CTOL = fstrEIG%lcztol
    NGET = fstrEIG%nget
    neig = eITMAX + NGET
    !maxItr = ITLIMIT
    maxItr = neig - 1

    write(*,*)"fstr_solve_lanczos"
    write(*,*)"neig", neig
    write(*,*)"maxItr", maxItr
    call fstr_solve_lanczos(hecMESH, hecMAT, fstrEIG, maxItr)

    call fstr_eigen_output(hecMESH,hecMAT,ILOG,fstrEIG)

    nstep = fstrEIG%nget

    do istep = 1, nstep
      do i=1,ntotal
        hecMAT%X(i) = ewk(i,istep)
      enddo

      IF(NDOF.EQ.3) THEN
        call hecmw_update_3_R(hecMESH,hecMAT%X,hecMAT%NP)
      ELSE IF( NDOF .EQ. 2 ) THEN
        call hecmw_update_2_R(hecMESH,hecMAT%X,hecMAT%NP)
      ELSE IF(ndof.EQ.6) THEN
        call hecmw_update_m_R(hecMESH,hecMAT%X,hecMAT%NP,ndof)
      endif

      if( IRESULT.eq.1 ) then
        header = "*fstrresult"
        call hecmw_result_init(hecMESH,nstep,istep,header)
        label = "DISPLACEMENT"
        call hecmw_result_add(1,NDOF,label,hecMAT%X)
        nameID = "fstrRES"
        call hecmw_result_write_by_name(nameID)
        call hecmw_result_finalize
      endif

      if( IVISUAL.eq.1 ) then
        call hecmw_nullify_result_data(fstrRESULT)
        fstrRESULT%nn_component = 1
        fstrRESULT%ne_component = 0
        allocate(fstrRESULT%nn_dof(1))
        allocate(fstrRESULT%node_label(1))
        allocate(fstrRESULT%node_val_item(NDOF*hecMAT%NP))
        fstrRESULT%nn_dof(1) = NDOF
        fstrRESULT%node_label(1) = 'DISPLACEMENT'
        fstrRESULT%node_val_item = hecMAT%X
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize(hecMESH,fstrRESULT,istep,nstep,1)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif
    enddo

    IF(myrank == 0)THEN
      WRITE(IMSG,'("### FSTR_SOLVE_EIGEN FINISHED!")')
      WRITE(*,'("### FSTR_SOLVE_EIGEN FINISHED!")')
    endif

  end subroutine fstr_solve_eigen

end module m_fstr_solve_eigen

