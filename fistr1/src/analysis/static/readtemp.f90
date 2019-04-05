!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module mReadTemp
contains

  !> Read in temperature distribution from external file
  subroutine read_temperature_result(hecMESH, nstep, sstep, interval, factor, temp, temp_bak)
    use hecmw_result
    use m_fstr
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: nstep, sstep, interval
    real(kind=kreal) :: factor
    real(kind=kreal),pointer :: temp(:), temp_bak(:)

    type(hecmwST_result_data) :: result
    character(len=HECMW_NAME_LEN) :: name_ID
    integer(kind=kint) :: i, k0, kt, nstep_active, tstep
    real(kind=kreal) :: w

    real(kind=kreal), pointer :: temp0(:,:)

    allocate(temp0(2,size(temp)))

    nstep_active = (nstep-sstep)/interval + 1
    kt = floor(factor*dble(nstep_active)-1.d-10)
    w = factor*dble(nstep_active)-dble(kt)

    do k0=1,2
      tstep = sstep+(kt+k0-2)*interval
      if( tstep == 0 ) then
        do i=1, hecMESH%n_node
          temp0(k0,i) = temp_bak(i)
        enddo
        cycle
      end if

      call hecmw_nullify_result_data( result )
      name_ID = 'fstrTEMP'
      call hecmw_result_read_by_name(hecMESH, name_ID, tstep, result)

      if (result%nn_component /= 1 .or. result%nn_dof(1) /= 1) then
        write(*,*) ' Read temperature result failed; not heat analysis result'
      endif

      do i=1, hecMESH%n_node
        temp0(k0,i) = result%node_val_item(i)
      enddo

      call hecmw_result_free(result)
    enddo

    do i=1, hecMESH%n_node
      temp(i) = (1.d0-w)*temp0(1,i)+w*temp0(2,i)
    enddo

    deallocate(temp0)
    write(IDBG,*) ' Read temperature from result file : OK'
  end subroutine read_temperature_result

end module
