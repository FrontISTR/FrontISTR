!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains steady state frequency analysis
module fstr_frequency_analysis

  use m_fstr
  use m_fstr_CreateMatrix_and_DampingForce
  use m_fstr_AddBC
  use m_fstr_EIG_setMASS
  use m_fstr_NodalKinematics
  use m_fstr_Update
  use m_fstr_Residual, only: fstr_Update_REACTION_SPC
  use m_fstr_NodalStress
  use m_make_result
  use m_hecmw2fstr_mesh_conv
  use hecmw_setup_util, only: node_global_to_local

  implicit none

contains

  subroutine fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, fstrEIG,  &
      fstrDYNAMIC, fstrRESULT, fstrPARAM, &
      fstrCPL, fstrFREQ, hecLagMAT, restart_step_num)
    !C
    !C-- global variable
    !C
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_result_data)            :: fstrRESULT
    type(fstr_param)                     :: fstrPARAM
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_couple)                    :: fstrCPL
    type(fstr_freqanalysis)              :: fstrFREQ
    type(hecmwST_matrix_lagrange)        :: hecLagMAT
    integer(kind=kint)                   :: restart_step_num

    !C
    !C-- local variable
    !C
    integer(kind=kint)            :: numnode, numelm, startmode, endmode, nummode, ndof, im, in, ntotal, vistype
    integer(kind=kint)            :: numfreq, idnode, numdisp, imonit, nmonit, monitnode(1)
    integer(kind=kint)            :: freqiout(3)
    integer(kind=kint)            :: ierr
    integer(kind=kint), parameter :: ilogin = 9056
    real(kind=kreal), allocatable :: eigenvalue(:), loadvecRe(:), loadvecIm(:)
    real(kind=kreal), allocatable :: bjre(:), bjim(:)
    real(kind=kreal), allocatable :: dispRe(:), dispIm(:), velRe(:), velIm(:), accRe(:), accIm(:)
    real(kind=kreal)              :: freq, omega, val, dx, dy, dz, f_start, f_end
    real(kind=kreal)              :: t_start, t_end, time, dxi, dyi, dzi
    type(fstr_freqanalysis_data)  :: freqData

    numnode   = hecMESH%n_node
    numelm    = hecMESH%n_elem
    ndof      = hecMESH%n_dof
    startmode = fstrFREQ%start_mode
    endmode   = fstrFREQ%end_mode
    nummode   = endmode - startmode +1

    freqData%numMode    = nummode
    freqData%numNodeDOF = numnode*ndof
    allocate(freqData%eigOmega(nummode))
    allocate(freqData%eigVector(numnode*ndof, nummode))

    call setupFREQParam(fstrDYNAMIC, f_start, f_end, numfreq, freqData%rayAlpha, freqData%rayBeta, idnode, vistype, freqiout)
    call setupOutputItem(fstrSOLID, freqiout)

    monitnode(1) = idnode
    imonit       = 0
    if( node_global_to_local(hecMESH, monitnode, 1) == 1 ) then
      if( monitnode(1) <= hecMESH%nn_internal ) imonit = monitnode(1)  !monitor is reported by the rank owning it
    end if
    nmonit = 0
    if( imonit > 0 ) nmonit = 1
    call hecmw_allreduce_I1(hecMESH, nmonit, hecmw_sum)
    if( nmonit == 0 ) then
      if( myrank == 0 ) then
        write(*,   *) "Error: monitor node not found:", idnode
        write(ilog,*) "Error: monitor node not found:", idnode
      end if
      call hecmw_abort( hecmw_comm_get_comm() )
    end if

    if( myrank == 0 ) then
      write(*,*) "Rayleigh alpha:", freqData%rayAlpha
      write(*,*) "Rayleigh beta:", freqData%rayBeta
      write(ilog,*) "Rayleigh alpha:", freqData%rayAlpha
      write(ilog,*) "Rayleigh beta:", freqData%rayBeta
    end if


    allocate(eigenvalue(nummode))
    allocate(loadvecRe(numnode*ndof))
    allocate(loadvecIm(numnode*ndof))
    allocate(bjre(nummode))
    allocate(bjim(nummode))
    allocate(dispRe(numnode*ndof))
    allocate(dispIm(numnode*ndof))
    allocate(velRe(numnode*ndof))
    allocate(velIm(numnode*ndof))
    allocate(accRe(numnode*ndof))
    allocate(accIm(numnode*ndof))

    loadvecRe(:) = 0.0D0
    loadvecIm(:) = 0.0D0

    if( myrank == 0 ) then
      write(*,*) "--frequency analysis--"
      write(*,   *) "read from=", trim(fstrFREQ%eigenlog_filename)
      write(ilog,*) "read from=", trim(fstrFREQ%eigenlog_filename)
      write(*,   *) "start mode=", startmode
      write(ilog,*) "start mode=", startmode
      write(*,   *) "end mode=", endmode
      write(ilog,*) "end mode=", endmode
    end if
    open(unit=ilogin, file=trim(fstrFREQ%eigenlog_filename), status="OLD", action="READ", iostat=ierr)
    if( ierr /= 0 ) then
      if( myrank == 0 ) write(*,*) "Error: cannot open eigenlog file: ", trim(fstrFREQ%eigenlog_filename)
      call hecmw_abort( hecmw_comm_get_comm() )
    endif
    call read_eigen_values(ilogin, startmode, endmode, eigenvalue, freqData%eigOmega)
    !call read_eigen_vector(ilogin, startmode, endmode, ndof, numnode, freqData%eigVector)
    close(ilogin)

    call read_eigen_vector_res(hecMESH, startmode, endmode, ndof, numnode, freqData%eigVector)

    call extract_surf2node(hecMESH, fstrFREQ, ndof, loadvecRe, loadvecIm)
    call assemble_nodeload(hecMESH, fstrFREQ, ndof, loadvecRe, loadvecIm)

    if( myrank == 0 ) write(*,*) "calc mass matrix"
    call calcMassMatrix(fstrPARAM, hecMESH, hecMAT, fstrSOLID, fstrEIG, hecLagMAT)
    if( myrank == 0 ) write(*,*) "scale eigenvector"
    call scaleEigenVector(hecMESH, fstrEIG, ndof*numnode, nummode, freqData%eigVector)

    if( myrank == 0 ) then
      write(*,   *) "start frequency:", f_start
      write(ilog,*) "start frequency:", f_start
      write(*,   *) "end frequency:", f_end
      write(ilog,*) "end frequency:", f_end
      write(*   ,*) "number of the sampling points", numfreq
      write(ilog,*) "number of the sampling points", numfreq
      write(*   ,*) "monitor nodeid=", idnode
      write(ilog,*) "monitor nodeid=", idnode
    end if

    do im=1, numfreq
      freq = (f_end-f_start)/dble(numfreq)*dble(im) + f_start
      omega = 2.0D0 * 3.14159265358979D0 * freq

      call calcFreqCoeff(hecMESH, freqData, loadvecRe, loadvecIm, omega, bjRe, bjIm)
      call calcDispVector(freqData, bjRe, bjIm, dispRe, dispIm)
      call calcVelVector(freqData, omega, bjRe, bjIm, velRe, velIm)
      call calcAccVector(freqData, omega, bjRe, bjIm, accRe, accIm)

      if( imonit > 0 ) then
        dx  = sqrt(dispRe(3*(imonit-1)+1)**2 + dispIm(3*(imonit-1)+1)**2)
        dy  = sqrt(dispRe(3*(imonit-1)+2)**2 + dispIm(3*(imonit-1)+2)**2)
        dz  = sqrt(dispRe(3*(imonit-1)+3)**2 + dispIm(3*(imonit-1)+3)**2)
        val = sqrt(dx**2 + dy**2 + dz**2)
        write(*,    *) freq, "[Hz] : ", val
        write(ilog, *) freq, "[Hz] : ", val
      end if

      if( myrank == 0 ) then
        if(IRESULT==1) then
          write(*,   *) freq, "[Hz] : ", im, ".res"
          write(ilog,*) freq, "[Hz] : ", im, ".res"
        end if
        if(IVISUAL==1 .and. vistype==1) then
          write(*,   *) freq, "[Hz] : ", im, ".vis"
          write(ilog,*) freq, "[Hz] : ", im, ".vis"
        end if
      end if
      call output_result(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, 'fstrRES', 'frequency', freq, im, &
        IRESULT==1, IVISUAL==1 .and. vistype==1, dispRe, dispIm, velRe, velIm, accRe, accIm)
    end do

    call setupDYNAParam(fstrDYNAMIC, t_start, t_end, freq, numdisp)
    if( myrank == 0 ) then
      write(*,   *) "start time:", t_start
      write(ilog,*) "start time:", t_start
      write(*,   *) "end time:", t_end
      write(ilog,*) "end time:", t_end
      write(*,   *) "frequency:", freq
      write(ilog,*) "frequency:", freq
      write(*,   *) "node id:", idnode
      write(ilog,*) "node id:", idnode
      write(*,   *) "num disp:", numdisp
      write(ilog,*) "num disp:", numdisp
    end if

    omega = 2.0D0 * 3.14159265358979D0 * freq
    call calcFreqCoeff(hecMESH, freqData, loadvecRe, loadvecIm, omega, bjRe, bjIm)

    do im=1, numdisp
      time = (t_end-t_start)/dble(numdisp)*dble(im-1) + t_start
      call calcDispVectorTime(freqData, time, omega, bjRe, bjIm, dispRe, dispIm)
      call calcVelVectorTime(freqData, time, omega, bjRe, bjIm, velRe, velIm)
      call calcAccVectorTime(freqData, time, omega, bjRe, bjIm, accRe, accIm)
      if( imonit > 0 ) then
        dx  = dispRe(3*(imonit-1)+1)
        dy  = dispRe(3*(imonit-1)+2)
        dz  = dispRe(3*(imonit-1)+3)
        dxi = dispIm(3*(imonit-1)+1)
        dyi = dispIm(3*(imonit-1)+2)
        dzi = dispIm(3*(imonit-1)+3)
      end if

      if( myrank == 0 ) then
        if(IRESULT==1) then
          write(*,   *) "time=", time, " : ", im, ".res"
          write(ilog,*) "time=", time, " : ", im, ".res"
        end if
        if(IVISUAL==1 .and. vistype==2) then
          write(*,   *) "time=", time, " : ", im, ".vis"
          write(ilog,*) "time=", time, " : ", im, ".vis"
        end if
      end if
      call output_result(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, 'fstrDYNA', 'TOTALTIME', time, im, &
        IRESULT==1, IVISUAL==1 .and. vistype==2, dispRe, dispIm, velRe, velIm, accRe, accIm)
    end do

    deallocate(freqData%eigOmega)
    deallocate(freqData%eigVector)
    deallocate(eigenvalue)
    deallocate(loadvecRe)
    deallocate(loadvecIm)
    deallocate(bjre)
    deallocate(bjim)
    deallocate(dispRe)
    deallocate(dispIm)
    deallocate(velRe)
    deallocate(velIm)
    deallocate(accRe)
    deallocate(accIm)

  end subroutine

  subroutine read_eigen_values(logfile, startmode, endmode, eigenvalue, anglfreq)
    !---- args
    integer(kind=kint), intent(in)  :: logfile
    integer(kind=kint), intent(in)  :: startmode
    integer(kind=kint), intent(in)  :: endmode
    real(kind=kreal), intent(inout) :: eigenvalue(:) !intend(endmode-startmode+1)
    real(kind=kreal), intent(inout) :: anglfreq(:)
    !---- vals
    integer(kind=kint)           :: im, endflag, id
    character(len=HECMW_MSG_LEN) :: line
    real(kind=kreal)             :: freq
    !---- body

    rewind(logfile)
    endflag = 0
    !Find eigenvalue header.
    do
      read(logfile, '(A80)', err=119) line
      if(trim(adjustl(line)) == "NO.  EIGENVALUE  FREQUENCY   (HZ)        X           Y           Z           X") then
        endflag = 1
        exit
      end if
    end do
    read(logfile, '(A80)') line
    !read eigenvalue
    do im=1, startmode-1
      read(logfile, '(A80)') line
    end do
    do im=1, (endmode-startmode+1)
      read(logfile, '(i5,3e12.4,a)', err=119) id, eigenvalue(im), anglfreq(im), freq, line
    end do
    return

    !error handling
    119 write(*,*) "Error to find eigenvalue information from logfile"
    write(ilog,*) "Error to find eigenvalue information from logfile"
    stop
  end subroutine

  subroutine read_eigen_vector(logfile, startmode, endmode, numdof, numnode, eigenvector)
    !---- args
    integer(kind=kint), intent(in)  :: logfile
    integer(kind=kint), intent(in)  :: startmode
    integer(kind=kint), intent(in)  :: endmode
    integer(kind=kint), intent(in)  :: numdof
    integer(kind=kint), intent(in)  :: numnode
    real(kind=kreal), intent(inout) :: eigenvector(:, :) !intend (numdof*NN,nmode)
    !---- vals
    integer(kind=kint)           :: im, in, gblid, j, idx
    real(kind=kreal)             :: vec(6)
    character(len=HECMW_MSG_LEN) :: line
    !---- body

    rewind(logfile)
    !Find first eigenvector header
    do im=1, startmode-1
      do
        read(logfile, '(a80)', err=119, end=119) line
        if(line(1:9) == " Mode No.") then
          exit ! goto next mode
        end if
      end do
    end do

    !read eigenvector
    do im=1, (endmode-startmode+1)
      !find header
      do
        read(logfile, '(a80)', err=119, end=119) line
        if(line(1:9) == " Mode No.") then
          exit !find eigenmode
        end if
      end do

      read(logfile, '(a80)', err=119) line
      read(logfile, '(a80)', err=119) line
      read(logfile, '(a80)', err=119) line
      read(logfile, '(a80)', err=119) line
      !read eigenvector component
      do in=1, numnode
        select case(numdof)
          case(2)
            read(logfile, '(i10,2e12.4)', err=119) gblid, (vec(j), j=1,2)

          case(3)
            !read(logfile, '(i10,3e12.4)', ERR=119) gblid, (vec(j), j=1,3)
            read(logfile, '(i10,3e16.8)', err=119) gblid, (vec(j), j=1,3)

          case(6)
            read(logfile, '(i10,6e12.4)', err=119) gblid, (vec(j), j=1,6)
          case default
            !error
            goto 119
        end select

        do j=1, numdof
          idx = (in-1)*numdof + j
          eigenvector(idx,im) = vec(j)
        end do
      end do
    end do
    return

    !error handling
    119 write(*,*) "Error to find eigenvector from logfile"
    write(ilog,*) "Error to find eigenvector from logfile"
    stop

  end subroutine

  subroutine read_eigen_vector_res(hecMESH, startmode, endmode, numdof, numnode, eigenvector)
    !---- args
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in)       :: startmode
    integer(kind=kint), intent(in)       :: endmode
    integer(kind=kint), intent(in)       :: numdof
    integer(kind=kint), intent(in)       :: numnode
    real(kind=kreal), intent(inout)      :: eigenvector(:, :) !intend (numdof*NN,nmode)
    !---- vals
    integer(kind=kint), parameter :: compidx = 1 !Component index of displacement
    integer(kind=kint)            :: imode, idx, ind, a, b, nallcomp, j
    type(hecmwST_result_data)     :: eigenres
    character(len=HECMW_NAME_LEN) :: name
    !---- body

    name = 'result-in'
    do imode=startmode, endmode
      call nullify_result_data(eigenres)
      call hecmw_result_read_by_name(hecMESH, name, imode, eigenres)

      nallcomp = 0
      do ind=1,eigenres%nn_component
        nallcomp = nallcomp + eigenres%nn_dof(ind)
      end do

      idx = imode - startmode + 1
      do ind=1, numnode
        do j=1, numdof
          a = (ind-1)*nallcomp + j  !src vector index
          b = (ind-1)*numdof   + j
          eigenvector(b,idx) = eigenres%node_val_item(a)
        end do
      end do
      call free_result_data(eigenres)
    end do

  contains

    subroutine free_result_data(res)
      !---- args
      type(hecmwST_result_data) :: res
      !---- vals
      !---- body
      if(associated(res%nn_dof))        deallocate(res%nn_dof)
      if(associated(res%ne_dof))        deallocate(res%ne_dof)
      if(associated(res%node_label))    deallocate(res%node_label)
      if(associated(res%elem_label))    deallocate(res%elem_label)
      if(associated(res%node_val_item)) deallocate(res%node_val_item)
      if(associated(res%elem_val_item))  deallocate(res%elem_val_item)
    end subroutine

    subroutine nullify_result_data(res)
      !---- args
      type(hecmwST_result_data) :: res
      !---- vals
      !---- body
      nullify(res%nn_dof)
      nullify(res%ne_dof)
      nullify(res%node_label)
      nullify(res%elem_label)
      nullify(res%node_val_item)
      nullify(res%elem_val_item)
    end subroutine


  end subroutine

  !
  ! The displacement, the velocity and the acceleration of the frequency response
  ! analysis are selected by iout_list of !DYNAMIC, and the other items follow
  ! !OUTPUT_RES as in the static and the transient dynamic analysis
  !
  subroutine setupOutputItem(fstrSOLID, iout)
    !---- args
    type(fstr_solid), intent(inout) :: fstrSOLID
    integer(kind=kint), intent(in)  :: iout(3)
    !---- vals
    integer(kind=kint) :: i
    !---- body

    do i=3, 4
      fstrSOLID%output_ctrl(i)%outinfo%on(1)  = (iout(1) == 1)
      fstrSOLID%output_ctrl(i)%outinfo%on(15) = (iout(2) == 1)
      fstrSOLID%output_ctrl(i)%outinfo%on(16) = (iout(3) == 1)
    end do
  end subroutine

  !
  ! Set one phase of the response to fstrSOLID and fstrDYNAMIC, and recover the
  ! strain and the stress of it with the routines of the static analysis.
  ! The incremental form keeps the incompatible mode of the 361 element consistent
  ! with the displacement, as the material is linear elastic here.
  !
  subroutine updateSolidResult(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, istp, disp, vel, acc)
    !---- args
    type(hecmwST_local_mesh)       :: hecMESH
    type(hecmwST_matrix)           :: hecMAT
    type(fstr_solid)               :: fstrSOLID
    type(fstr_dynamic)             :: fstrDYNAMIC
    integer(kind=kint), intent(in) :: istp
    real(kind=kreal), intent(in)   :: disp(:) !intend (numnodeDOF)
    real(kind=kreal), intent(in)   :: vel(:) !intend (numnodeDOF)
    real(kind=kreal), intent(in)   :: acc(:) !intend (numnodeDOF)
    !---- vals
    integer(kind=kint) :: ndof, ntotal, idx
    !---- body

    ndof   = hecMESH%n_dof
    ntotal = ndof*hecMESH%n_node

    idx = 1
    if( fstrDYNAMIC%idx_eqa == 1 .and. istp > 0 ) idx = 2
    fstrDYNAMIC%DISP(1:ntotal,idx) = disp(1:ntotal)
    fstrDYNAMIC%VEL (1:ntotal,idx) = vel(1:ntotal)
    fstrDYNAMIC%ACC (1:ntotal,idx) = acc(1:ntotal)

    hecMAT%X(1:ntotal) = disp(1:ntotal) - fstrSOLID%unode(1:ntotal) - fstrSOLID%dunode(1:ntotal)
    call fstr_apply_solution_increment( hecMESH, fstrSOLID, ndof, hecMAT%X )
    call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID, 0.0D0, 0.0D0, 1 )
    call fstr_Update_REACTION_SPC( 1, hecMESH, fstrSOLID )

    if( ndof == 2 ) then
      call fstr_NodalStress2D( hecMESH, fstrSOLID )
    else if( ndof == 6 ) then
      call fstr_NodalStress6D( hecMESH, fstrSOLID )
    else
      call fstr_NodalStress3D( hecMESH, fstrSOLID )
    end if
  end subroutine

  !
  ! The response is complex, so the real part and the imaginary part are written as
  ! two sets of the same items, distinguished by the label suffix. Both are recovered
  ! once and shared by the result file and the visualization data.
  !
  subroutine output_result(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, nameid, glabel, gvalue, istp, &
      wres, wvis, dispRe, dispIm, velRe, velIm, accRe, accIm)
    !---- args
    type(hecmwST_local_mesh)       :: hecMESH
    type(hecmwST_matrix)           :: hecMAT
    type(fstr_solid)               :: fstrSOLID
    type(fstr_dynamic)             :: fstrDYNAMIC
    character(len=*), intent(in)   :: nameid
    character(len=*), intent(in)   :: glabel
    real(kind=kreal), intent(in)   :: gvalue
    integer(kind=kint), intent(in) :: istp
    logical, intent(in)            :: wres, wvis
    real(kind=kreal), intent(in)   :: dispRe(:), dispIm(:) !intend (numnodeDOF)
    real(kind=kreal), intent(in)   :: velRe(:), velIm(:) !intend (numnodeDOF)
    real(kind=kreal), intent(in)   :: accRe(:), accIm(:) !intend (numnodeDOF)
    !---- vals
    type(hecmwST_result_data)       :: visRe, visIm
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN)    :: comment
    character(len=HECMW_NAME_LEN)   :: label, resid
    real(kind=kreal)                :: gval(1)
    !---- body

    if( .not. (wres .or. wvis) ) return

    if( wres ) then
      header='*fstrresult'
      comment='frequency_result'
      call hecmw_result_init(hecMESH, istp, header, comment)

      label   = glabel
      gval(1) = gvalue
      call hecmw_result_add(HECMW_RESULT_DTYPE_GLOBAL, 1, label, gval)
    end if

    call updateSolidResult(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, istp, dispRe, velRe, accRe)
    if( wres ) call fstr_write_result_add(hecMESH, fstrSOLID, istp, '_real', fstrDYNAMIC)
    if( wvis ) call fstr_make_result(hecMESH, fstrSOLID, visRe, istp, gvalue, fstrDYNAMIC, '_real')

    call updateSolidResult(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, istp, dispIm, velIm, accIm)
    if( wres ) call fstr_write_result_add(hecMESH, fstrSOLID, istp, '_imag', fstrDYNAMIC)
    if( wvis ) call fstr_make_result(hecMESH, fstrSOLID, visIm, istp, gvalue, fstrDYNAMIC, '_imag')

    if( wres ) then
      resid = nameid
      call hecmw_result_write_by_name(resid)
      call hecmw_result_finalize()
    end if

    if( wvis ) then
      call fstr_merge_result(hecMESH, visRe, visIm)
      call hecmw_result_free(visIm)
      call fstr2hecmw_mesh_conv(hecMESH)
      call hecmw_visualize_init
      call hecmw_visualize( hecMESH, visRe, istp )
      call hecmw2fstr_mesh_conv(hecMESH)
      call hecmw_result_free(visRe)
    end if
    return
  end subroutine

  subroutine extract_surf2node(hecMESH, freqData, numdof, loadvecRe, loadvecIm)
    !---- args
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_freqanalysis), intent(in)  :: freqData
    integer(kind=kint), intent(in)       :: numdof
    real(kind=kreal), intent(inout)      :: loadvecRe(:) !intend(numnode*ndof)
    real(kind=kreal), intent(inout)      :: loadvecIm(:) !intend(numnode*ndof)
    !---- vals
    integer(kind=kint), parameter   :: MAXNODE = 100
    integer(kind=kint)              :: sgrpID, is, ie, ic, nsurf, ic_type, outtype, node_index(MAXNODE)
    integer(kind=kint)              :: nn, iss, nodeid, dof_index, ndof
    integer(kind=kint)              :: i, j, k, l, m, isn, nsize
    integer(kind=kint)              :: iwk(60), nodLOCAL(20)
    real(kind=kreal)                :: vect(60), xx(20), yy(20), zz(20), forcere(3), forceim(3)
    !---- body

    ndof = 3
    do i=1,freqData%FLOAD_ngrp_tot
      if(freqData%FLOAD_ngrp_TYPE(i) == kFLOADTYPE_SURF) then  !FLOAD type=surface
        sgrpID     = freqData%FLOAD_ngrp_ID(i)
        dof_index  = freqData%FLOAD_ngrp_DOF(i)
        forcere(:) = 0.0D0
        forceim(:) = 0.0D0
        forcere(dof_index) = freqData%FLOAD_ngrp_valre(i)
        forceim(dof_index) = freqData%FLOAD_ngrp_valim(i)

        is = hecMESH%surf_group%grp_index(sgrpID-1) + 1
        ie = hecMESH%surf_group%grp_index(sgrpID)
        do j=is, ie
          ic      = hecMESH%surf_group%grp_item(2*j-1)
          nsurf   = hecMESH%surf_group%grp_item(2*j)
          ic_type = hecMESH%elem_type(ic)
          nn  = hecmw_get_max_node(ic_type)
          isn = hecMESH%elem_node_index(ic-1)
          do k=1, nn
            nodLOCAL(k) = hecMESH%elem_node_item(isn+k)
            xx(k) = hecMESH%node(3*nodLOCAL(k)-2)
            yy(k) = hecMESH%node(3*nodLOCAL(k)-1)
            zz(k) = hecMESH%node(3*nodLOCAL(k)  )
            do l=1, ndof
              iwk(ndof*(k-1)+l) = ndof*(nodLOCAL(k)-1)+l
            end do
          end do

          call DL_C3_freq(ic_type, nn, xx, yy, zz, nsurf, forcere, vect, nsize)
          do k=1,nsize
            loadvecRe(iwk(k)) = loadvecRe(iwk(k)) + vect(k)
          end do

          call DL_C3_freq(ic_type, nn, xx, yy, zz, nsurf, forceim, vect, nsize)
          do k=1,nsize
            loadvecIm(iwk(k)) = loadvecIm(iwk(k)) + vect(k)
          end do
        end do
      end if
    end do

    return
  end subroutine

  subroutine DL_C3_freq(ETYPE, NN, XX, YY, ZZ, LTYPE, force, VECT, nsize)
    !---- args
    integer(kind=kint), intent(in)    :: ETYPE !--solid element type
    integer(kind=kint), intent(in)    :: NN !--node num
    integer(kind=kint), intent(in)    :: LTYPE !--solid element face
    real(kind=kreal), intent(in)      :: XX(:) !--node x pos
    real(kind=kreal), intent(in)      :: YY(:) !--node y pos
    real(kind=kreal), intent(in)      :: ZZ(:) !--node z pos
    real(kind=kreal), intent(in)      :: force(3) !--node surfforce
    real(kind=kreal), intent(inout)   :: VECT(:)
    integer(kind=kint), intent(inout) :: nsize
    !---- vals
    integer(kind=kint), parameter :: NDOF = 3
    real(kind=kreal)              :: WG
    integer(kind=kint)            :: NOD(NN)
    real(kind=kreal)              :: elecoord(3, NN), localcoord(3)
    real(kind=kreal)              :: H(NN)
    integer(kind=kint)            :: I, IG2, NSUR, SURTYPE
    !---- body

    call getSubFace( ETYPE, LTYPE, SURTYPE, NOD )
    NSUR = getNumberOfNodes( SURTYPE )

    do I=1,NSUR
      elecoord(1,i)=XX(NOD(I))
      elecoord(2,i)=YY(NOD(i))
      elecoord(3,i)=ZZ(NOD(i))
    end do
    nsize         = NN*NDOF
    VECT(1:NSIZE) = 0.0D0
    do IG2=1,NumOfQuadPoints( SURTYPE )
      call getQuadPoint( SURTYPE, IG2, localcoord(1:2) )
      call getShapeFunc( SURTYPE, localcoord(1:2), H(1:NSUR) )

      WG=getWeight( SURTYPE, IG2 )
      do I=1,NSUR
        VECT(3*NOD(I)-2)=VECT(3*NOD(I)-2)+WG*H(I)*force(1)
        VECT(3*NOD(I)-1)=VECT(3*NOD(I)-1)+WG*H(I)*force(2)
        VECT(3*NOD(I)  )=VECT(3*NOD(I)  )+WG*H(I)*force(3)
      end do
    end do
  end subroutine

  subroutine assemble_nodeload(hecMESH, freqData, numdof, loadvecRe, loadvecIm)
    !---- args
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_freqanalysis), intent(in)  :: freqData
    integer(kind=kint), intent(in)       :: numdof
    real(kind=kreal), intent(inout)      :: loadvecRe(:)
    real(kind=kreal), intent(inout)      :: loadvecIm(:)
    !---- vals
    integer(kind=kint) :: i, vecsize, ig, is, ie, in, nodeid, dof_index

    !---- body

    do i=1, freqData%FLOAD_ngrp_tot
      if(freqData%FLOAD_ngrp_TYPE(i) == kFLOADTYPE_NODE) then
        ig = freqData%FLOAD_ngrp_ID(i)
        is = hecMESH%node_group%grp_index(ig-1) + 1
        ie = hecMESH%node_group%grp_index(ig)
        do in=is, ie
          nodeid    = hecMESH%node_group%grp_item(in)
          dof_index = freqData%FLOAD_ngrp_DOF(i)
          loadvecRe((nodeid-1)*numdof + dof_index) = loadvecRe((nodeid-1)*numdof + dof_index) + freqData%FLOAD_ngrp_valre(i)
          loadvecIm((nodeid-1)*numdof + dof_index) = loadvecIm((nodeid-1)*numdof + dof_index) + freqData%FLOAD_ngrp_valim(i)
        end do
      end if
    end do
    return
  end subroutine

  subroutine calcMassMatrix(fstrPARAM, hecMESH, hecMAT, fstrSOLID, fstrEIG, hecLagMAT)
    !---- args
    type(fstr_param), intent(in)         :: fstrPARAM
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(inout)  :: hecMAT
    type(fstr_solid), intent(inout)      :: fstrSOLID
    type(fstr_eigen), intent(inout)      :: fstrEIG
    type(hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT
    !---- vals
    integer(kind=kint) :: ntotal
    !---- body


    fstrSOLID%dunode = 0.d0
    call fstr_CreateMatrix_and_DampingForce( hecMESH, hecMAT, fstrSOLID, 0.d0, 0.d0 )
    call fstr_AddBC(1, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 2)

    call setMASS(fstrSOLID, hecMESH, hecMAT, fstrEIG)

  end subroutine

  subroutine scaleEigenVector(hecMESH, fstrEIG, ntotaldof, nmode, eigenvector)
    !---- args
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_eigen), intent(in)    :: fstrEIG
    integer(kind=kint), intent(in)  :: ntotaldof
    integer(kind=kint), intent(in)  :: nmode
    real(kind=kreal), intent(inout) :: eigenvector(:, :)
    !---- vals
    integer(kind=kint) :: imode, idof, nintdof
    real(kind=kreal)   :: mas
    !---- body

    nintdof = hecMESH%nn_internal*hecMESH%n_dof

    do imode=1,nmode
      mas = 0.0D0
      do idof=1,nintdof
        mas = mas + fstrEIG%mass(idof)*eigenvector(idof,imode)**2
      end do
      call hecmw_allreduce_R1(hecMESH, mas, hecmw_sum)
      do idof=1,ntotaldof
        eigenvector(idof,imode) =  eigenvector(idof,imode) / sqrt(mas)
      end do
    end do
  end subroutine

  subroutine checkOrthVector(fstrEIG, eigenvector, imode, jmode, prod)
    !---- args
    type(fstr_eigen), intent(in)    :: fstrEIG
    real(kind=kreal), intent(in)    :: eigenvector(:, :)
    integer(kind=kint), intent(in)  :: imode
    integer(kind=kint), intent(in)  :: jmode
    real(kind=kreal), intent(inout) :: prod
    !---- vals
    integer(kind=kint) :: idof, s
    !---- body
    s = size(eigenvector(:,1))
    prod = 0.0D0

    do idof=1, s
      prod = prod + eigenvector(idof,imode)*fstrEIG%mass(idof)*eigenvector(idof,jmode)
    end do
    return
  end subroutine

  subroutine writeOutvector(im, vector)
    !---- args
    integer(kind=kint), intent(in) :: im
    real(kind=kreal), intent(in)   :: vector(:)
    !---- vals
    integer(kind=kint) :: i, s
    !---- body
    s = size(vector)
    do i=1, s
      if(i == 1) then
        write(*,'("eigenvec",i2.2,":[",e12.5,", ")') im, vector(i)
      else if(i /= s) then
        write(*,'(e12.5,", ")') vector(i)
      else
        write(*,'(e12.5,"];")') vector(i)
      end if
    end do
    write(*,*)
    return
  end subroutine

  subroutine calcDotProduct(hecMESH, a, b, c)
    !---- args
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    real(kind=kreal), intent(in)    :: a(:)
    real(kind=kreal), intent(in)    :: b(:)
    real(kind=kreal), intent(inout) :: c
    !---- vals
    integer(kind=kint) :: nintdof
    !---- body
    nintdof = hecMESH%nn_internal*hecMESH%n_dof  !external nodes are owned by another rank
    c = dot_product(a(1:nintdof), b(1:nintdof))
    call hecmw_allreduce_R1(hecMESH, c, hecmw_sum)
    return
  end subroutine

  subroutine calcFreqCoeff(hecMESH, freqData, loadRe, loadIm, inpOmega, bjRe, bjIm)
    !---- args
    type(hecmwST_local_mesh), intent(in)     :: hecMESH
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: loadRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(in)             :: loadIm(:) !intend (numNodeDOF)
    real(kind=kreal), intent(in)             :: inpOmega
    real(kind=kreal), intent(inout)          :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: bjIm(:) !intend (numMode)
    !---- vals
    integer(kind=kint) :: imode
    real(kind=kreal)   :: ujfr, ujfi, a, b, alp, beta
    !---- body

    alp  = freqData%rayAlpha
    beta = freqData%rayBeta

    do imode=1, freqData%numMode
      call calcDotProduct(hecMESH, freqData%eigVector(:,imode), loadRe, ujfr)
      call calcDotProduct(hecMESH, freqData%eigVector(:,imode), loadIm, ujfi)

      a = ujfr*(freqData%eigOmega(imode)**2 - inpOmega**2) + ujfi*(alp + beta*freqData%eigOmega(imode)**2)*inpOmega
      b = (freqData%eigOmega(imode)**2 - inpOmega**2)**2 + ((alp + beta*freqData%eigOmega(imode)**2)*inpOmega)**2
      bjRe(imode) = a / b

      a = ujfi*(freqData%eigOmega(imode)**2 -inpOmega**2) - ujfr*(alp + beta*freqData%eigOmega(imode)**2)*inpOmega
      b = (freqData%eigOmega(imode)**2 - inpOmega**2)**2 + ((alp + beta*freqData%eigOmega(imode)**2)*inpOmega)**2
      bjIm(imode) = a / b
    end do
    return
  end subroutine

  subroutine calcDispVector(freqData, bjRe, bjIm, dispRe, dispIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: dispRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: dispIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint) :: imode
    !---- body

    dispRe(:) = 0.0D0
    dispIm(:) = 0.0D0

    do imode=1, freqData%numMode
      dispRe(:) = dispRe(:) + bjRe(imode)*freqData%eigVector(:,imode)
      dispIm(:) = dispIm(:) + bjIm(imode)*freqData%eigVector(:,imode)
    end do
    return
  end subroutine

  subroutine calcVelVector(freqData, omega, bjRe, bjIm, velRe, velIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: omega
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: velRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: velIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint) :: imode
    !---- body

    velRe(:) = 0.0D0
    velIm(:) = 0.0D0

    do imode=1, freqData%numMode
      velRe(:) = velRe(:) - omega * bjIm(imode) * freqData%eigVector(:,imode)
      velIm(:) = velIm(:) + omega * bjRe(imode) * freqData%eigVector(:,imode)
    end do
  end subroutine

  subroutine calcAccVector(freqData, omega, bjRe, bjIm, accRe, accIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: omega
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: accRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: accIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint) :: imode
    !---- body

    accRe(:) = 0.0D0
    accIm(:) = 0.0D0

    do imode=1, freqData%numMode
      accRe(:) = accRe(:) - omega**2 * bjRe(imode) * freqData%eigVector(:,imode)
      accIm(:) = accIm(:) - omega**2 * bjIm(imode) * freqData%eigVector(:,imode)
    end do

  end subroutine

  subroutine setupFREQParam(fstrDYNAMIC, f_start, f_end, numfreq, raym, rayk, idnode, vistype, ioutl)
    !---- args
    type(fstr_dynamic), intent(in)    :: fstrDYNAMIC
    real(kind=kreal), intent(inout)   :: f_start
    real(kind=kreal), intent(inout)   :: f_end
    integer(kind=kint), intent(inout) :: numfreq
    real(kind=kreal), intent(inout)   :: raym
    real(kind=kreal), intent(inout)   :: rayk
    integer(kind=kint), intent(inout) :: idnode
    integer(kind=kint), intent(inout) :: vistype
    integer(kind=kint), intent(inout) :: ioutl(3)
    !---- vals

    !---- body
    f_start    = fstrDYNAMIC%t_start
    f_end      = fstrDYNAMIC%t_end
    numfreq    = fstrDYNAMIC%n_step
    raym       = fstrDYNAMIC%ray_m
    rayk       = fstrDYNAMIC%ray_k
    idnode     = fstrDYNAMIC%nout_monit
    vistype    = fstrDYNAMIC%ngrp_monit
    ioutl(1:3) = fstrDYNAMIC%iout_list(1:3)
    return
  end subroutine

  subroutine calcDispVectorTime(freqData, time, omega, bjRe, bjIm, dispRe, dispIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: time
    real(kind=kreal), intent(in)             :: omega
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: dispRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: dispIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint)  :: imode, idf, s
    complex(kind=kreal) :: a, b, c
    !---- body

    dispRe(:) = 0.0D0
    dispIm(:) = 0.0D0
    a = exp(dcmplx(0.0D0, omega*time))

    do imode=1, freqData%numMode
      s = size(freqData%eigvector(:,imode))
      b = dcmplx(bjRe(imode), bjIm(imode)) * a
      do idf=1, s
        c = b*dcmplx(freqData%eigVector(idf,imode), 0.0D0)
        dispRe(idf) = dispRe(idf) + dble(c)
        dispIm(idf) = dispIm(idf) + imag(c)
      end do
    end do
    return
  end subroutine

  subroutine calcVelVectorTime(freqData, time, omega, bjRe, bjIm, velRe, velIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: time
    real(kind=kreal), intent(in)             :: omega
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: velRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: velIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint)  :: imode, idf, s
    complex(kind=kreal) :: a, b, c
    !---- body

    velRe(:) = 0.0D0
    velIm(:) = 0.0D0
    a = dcmplx(0.0D0, 1.0D0)*dcmplx(omega, 0.0D0)*exp(dcmplx(0.0D0, omega*time))

    do imode=1, freqData%numMode
      s = size(freqData%eigvector(:,imode))
      b = dcmplx(bjRe(imode), bjIm(imode)) * a
      do idf=1, s
        c = b*dcmplx(freqData%eigVector(idf,imode), 0.0D0)
        velRe(idf) = velRe(idf) + dble(c)
        velIm(idf) = velIm(idf) + imag(c)
      end do
    end do
    return
  end subroutine

  subroutine calcAccVectorTime(freqData, time, omega, bjRe, bjIm, accRe, accIm)
    !---- args
    type(fstr_freqanalysis_data), intent(in) :: freqData
    real(kind=kreal), intent(in)             :: time
    real(kind=kreal), intent(in)             :: omega
    real(kind=kreal), intent(in)             :: bjRe(:) !intend (numMode)
    real(kind=kreal), intent(in)             :: bjIm(:) !intend (numMode)
    real(kind=kreal), intent(inout)          :: accRe(:) !intend (numNodeDOF)
    real(kind=kreal), intent(inout)          :: accIm(:) !intend (numNodeDOF)
    !---- vals
    integer(kind=kint)  :: imode, idf, s
    complex(kind=kreal) :: a, b, c
    !---- body

    accRe(:) = 0.0D0
    accIm(:) = 0.0D0
    a = dcmplx(-1.0D0, 0.0D0)*dcmplx(omega**2, 0.0D0)*exp(dcmplx(0.0D0, omega*time))

    do imode=1, freqData%numMode
      s = size(freqData%eigvector(:,imode))
      b = dcmplx(bjRe(imode), bjIm(imode)) * a
      do idf=1, s
        c = b*dcmplx(freqData%eigVector(idf,imode), 0.0D0)
        accRe(idf) = accRe(idf) + dble(c)
        accIm(idf) = accIm(idf) + imag(c)
      end do
    end do
    return
  end subroutine

  subroutine setupDYNAParam(fstrDYNAMIC, t_start, t_end, dynafreq, numdisp)
    !---- args
    type(fstr_dynamic), intent(in)    :: fstrDYNAMIC
    real(kind=kreal), intent(inout)   :: t_start
    real(kind=kreal), intent(inout)   :: t_end
    real(kind=kreal), intent(inout)   :: dynafreq
    integer(kind=kint), intent(inout) :: numdisp
    !---- vals
    !---- body
    t_start  = fstrDYNAMIC%gamma
    t_end    = fstrDYNAMIC%beta
    dynafreq = fstrDYNAMIC%t_delta
    numdisp  = fstrDYNAMIC%nout
    return
  end subroutine

end module
