!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides functions to initialize heat analysis
module m_heat_init
contains

  subroutine heat_init(hecMESH, fstrHEAT)
    use m_fstr
    implicit none
    type(fstr_heat) :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: i, j

    allocate(fstrHEAT%TEMP0(hecMESH%n_node))
    allocate(fstrHEAT%TEMPC(hecMESH%n_node))
    allocate(fstrHEAT%TEMP (hecMESH%n_node))
    fstrHEAT%TEMP0 = 0.0d0
    fstrHEAT%TEMPC = 0.0d0
    fstrHEAT%TEMP  = 0.0d0

    if(hecMESH%hecmw_flag_initcon == 1)then
      do i = 1, hecMESH%n_node
        j = hecMESH%node_init_val_index(i)
        fstrHEAT%TEMP0(i) = hecMESH%node_init_val_item(j)
        fstrHEAT%TEMPC(i) = fstrHEAT%TEMP0(i)
        fstrHEAT%TEMP (i) = fstrHEAT%TEMP0(i)
      enddo
      write(ILOG,*) ' Initial condition of temperatures: OK'
    endif
    if( associated(g_InitialCnd) ) then
        do j=1,size(g_InitialCnd)
          if( g_InitialCnd(j)%cond_name=="temperature" ) then
            do i= 1, hecMESH%n_node
              fstrHEAT%TEMP0(i)= g_InitialCnd(j)%realval(i)
              fstrHEAT%TEMPC(i)= fstrHEAT%TEMP0(i)
              fstrHEAT%TEMP (i)= fstrHEAT%TEMP0(i)
            enddo
            exit
          endif
        enddo
        write(ILOG,*) ' Initial condition of temperatures: OK'
    endif

    call heat_fix_dummy_nodes(hecMESH, fstrHEAT)
  end subroutine heat_init

  !> Mark the rotational dummy nodes of 641/761/781 elements and fix them to zero
  !! temperature. The elements are assembled on their translational nodes only, so
  !! the rows of the dummy nodes would otherwise stay empty.
  subroutine heat_fix_dummy_nodes(hecMESH, fstrHEAT)
    use m_fstr
    use hecmw_etype, only : hecmw_is_etype_33struct, hecmw_get_max_node
    use hecmw_setup_util, only : hecmw_expand_integer_array, hecmw_expand_real_array
    implicit none
    type(fstr_heat) :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: icel, i, j, in0, nn, n_dummy, old_size, new_size

    allocate(fstrHEAT%is_dummy_node(hecMESH%n_node))
    fstrHEAT%is_dummy_node = .false.

    do icel = 1, hecMESH%n_elem
      if( .not. hecmw_is_etype_33struct(hecMESH%elem_type(icel)) ) cycle
      nn = hecmw_get_max_node(hecMESH%elem_type(icel))/2
      in0 = hecMESH%elem_node_index(icel-1)
      do i = 1, nn
        fstrHEAT%is_dummy_node(hecMESH%elem_node_item(in0+nn+i)) = .true.
      enddo
    enddo

    n_dummy = count(fstrHEAT%is_dummy_node)
    if( n_dummy == 0 ) return

    old_size = fstrHEAT%T_FIX_tot
    new_size = old_size + n_dummy
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_node, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_ampl, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%T_FIX_val,  old_size, new_size )
    fstrHEAT%T_FIX_tot = new_size

    j = old_size
    do i = 1, hecMESH%n_node
      if( .not. fstrHEAT%is_dummy_node(i) ) cycle
      fstrHEAT%TEMP0(i) = 0.0d0
      fstrHEAT%TEMPC(i) = 0.0d0
      fstrHEAT%TEMP (i) = 0.0d0
      j = j + 1
      fstrHEAT%T_FIX_node(j) = i
      fstrHEAT%T_FIX_ampl(j) = 0
      fstrHEAT%T_FIX_val(j)  = 0.0d0
    enddo
  end subroutine heat_fix_dummy_nodes

  subroutine heat_init_log(hecMESH, fstrHEAT)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_heat) :: fstrHEAT
    integer(kind=kint) :: nradiate

    if(hecMESH%my_rank == 0)then
      write(IMSG,*) '============================='
      write(IMSG,*) '  H E A T   T R A N S F E R  '
      write(IMSG,*) '============================='
      write(ISTA,*)
      write(ISTA,*)'  ISTEP    INCR    ITER     RESIDUAL     IITER   '
      write(ISTA,*)'-------------------------------------------------'
    endif

    nradiate = fstrHEAT%R_SUF_tot
    call hecmw_allreduce_I1(hecMESH, nradiate, hecmw_sum)
    if(nradiate > 0 .and. hecMESH%my_rank == 0)then
      write(IMSG,"(a,f12.4)")"* Absolute zero temperature for radiation: ", hecMESH%zero_temp
    endif
  end subroutine heat_init_log

  subroutine heat_finalize(fstrHEAT)
    use m_fstr
    implicit none
    type(fstr_heat) :: fstrHEAT
    if(associated(fstrHEAT%STEP_DLTIME)) deallocate(fstrHEAT%STEP_DLTIME)
    if(associated(fstrHEAT%STEP_EETIME)) deallocate(fstrHEAT%STEP_EETIME)
    if(associated(fstrHEAT%STEP_DELMIN)) deallocate(fstrHEAT%STEP_DELMIN)
    if(associated(fstrHEAT%STEP_DELMAX)) deallocate(fstrHEAT%STEP_DELMAX)
    if(associated(fstrHEAT%mat_name)) deallocate(fstrHEAT%mat_name)
    if(associated(fstrHEAT%RHOtab)) deallocate(fstrHEAT%RHOtab)
    if(associated(fstrHEAT%CPtab)) deallocate(fstrHEAT%CPtab)
    if(associated(fstrHEAT%CONDtab)) deallocate(fstrHEAT%CONDtab)
    if(associated(fstrHEAT%RHO)) deallocate(fstrHEAT%RHO)
    if(associated(fstrHEAT%RHOtemp)) deallocate(fstrHEAT%RHOtemp)
    if(associated(fstrHEAT%CP)) deallocate(fstrHEAT%CP)
    if(associated(fstrHEAT%CPtemp)) deallocate(fstrHEAT%CPtemp)
    if(associated(fstrHEAT%COND)) deallocate(fstrHEAT%COND)
    if(associated(fstrHEAT%CONDtemp)) deallocate(fstrHEAT%CONDtemp)
    if(associated(fstrHEAT%RHOfuncA)) deallocate(fstrHEAT%RHOfuncA)
    if(associated(fstrHEAT%RHOfuncB)) deallocate(fstrHEAT%RHOfuncB)
    if(associated(fstrHEAT%CPfuncA)) deallocate(fstrHEAT%CPfuncA)
    if(associated(fstrHEAT%CPfuncB)) deallocate(fstrHEAT%CPfuncB)
    if(associated(fstrHEAT%CONDfuncA)) deallocate(fstrHEAT%CONDfuncA)
    if(associated(fstrHEAT%CONDfuncB)) deallocate(fstrHEAT%CONDfuncB)
    if(associated(fstrHEAT%AMPLtab)) deallocate(fstrHEAT%AMPLtab)
    if(associated(fstrHEAT%AMPL)) deallocate(fstrHEAT%AMPL)
    if(associated(fstrHEAT%AMPLtime)) deallocate(fstrHEAT%AMPLtime)
    if(associated(fstrHEAT%AMPLfuncA)) deallocate(fstrHEAT%AMPLfuncA)
    if(associated(fstrHEAT%AMPLfuncB)) deallocate(fstrHEAT%AMPLfuncB)
    if(associated(fstrHEAT%TEMP0)) deallocate(fstrHEAT%TEMP0)
    if(associated(fstrHEAT%TEMPC)) deallocate(fstrHEAT%TEMPC)
    if(associated(fstrHEAT%TEMP)) deallocate(fstrHEAT%TEMP)
    if(associated(fstrHEAT%is_dummy_node)) deallocate(fstrHEAT%is_dummy_node)
    if(associated(fstrHEAT%T_FIX_node)) deallocate(fstrHEAT%T_FIX_node)
    if(associated(fstrHEAT%T_FIX_ampl)) deallocate(fstrHEAT%T_FIX_ampl)
    if(associated(fstrHEAT%T_FIX_val)) deallocate(fstrHEAT%T_FIX_val)
    if(associated(fstrHEAT%Q_NOD_node)) deallocate(fstrHEAT%Q_NOD_node)
    if(associated(fstrHEAT%Q_NOD_ampl)) deallocate(fstrHEAT%Q_NOD_ampl)
    if(associated(fstrHEAT%Q_NOD_val)) deallocate(fstrHEAT%Q_NOD_val)
    if(associated(fstrHEAT%Q_VOL_elem)) deallocate(fstrHEAT%Q_VOL_elem)
    if(associated(fstrHEAT%Q_VOL_ampl)) deallocate(fstrHEAT%Q_VOL_ampl)
    if(associated(fstrHEAT%Q_VOL_val)) deallocate(fstrHEAT%Q_VOL_val)
    if(associated(fstrHEAT%Q_SUF_elem)) deallocate(fstrHEAT%Q_SUF_elem)
    if(associated(fstrHEAT%Q_SUF_ampl)) deallocate(fstrHEAT%Q_SUF_ampl)
    if(associated(fstrHEAT%Q_SUF_surf)) deallocate(fstrHEAT%Q_SUF_surf)
    if(associated(fstrHEAT%Q_SUF_val)) deallocate(fstrHEAT%Q_SUF_val)
    if(associated(fstrHEAT%R_SUF_elem)) deallocate(fstrHEAT%R_SUF_elem)
    if(associated(fstrHEAT%R_SUF_ampl)) deallocate(fstrHEAT%R_SUF_ampl)
    if(associated(fstrHEAT%R_SUF_surf)) deallocate(fstrHEAT%R_SUF_surf)
    if(associated(fstrHEAT%R_SUF_val)) deallocate(fstrHEAT%R_SUF_val)
    if(associated(fstrHEAT%H_SUF_elem)) deallocate(fstrHEAT%H_SUF_elem)
    if(associated(fstrHEAT%H_SUF_ampl)) deallocate(fstrHEAT%H_SUF_ampl)
    if(associated(fstrHEAT%H_SUF_surf)) deallocate(fstrHEAT%H_SUF_surf)
    if(associated(fstrHEAT%H_SUF_val)) deallocate(fstrHEAT%H_SUF_val)
    if(associated(fstrHEAT%weldline)) deallocate(fstrHEAT%weldline)
  end subroutine heat_finalize

end module m_heat_init
