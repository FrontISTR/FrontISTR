MODULE mReadTemp

  contains

!> Read in temperature distribution from external file
    subroutine read_temperature_result(hecMESH, nstep, tstep, temp)
      use hecmw_result
      use m_fstr
      type(hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: nstep, tstep
      real(kind=kreal),pointer :: temp(:)
      type(hecmwST_result_data) :: result
      character(len=HECMW_NAME_LEN) :: name_ID
      integer(kind=kint) :: i

      call hecmw_nullify_result_data( result )
      name_ID = 'fstrTEMP'
      call hecmw_result_read_by_name(hecMESH, name_ID, nstep, tstep, result)

      if (result%nn_component /= 1 .or. result%nn_dof(1) /= 1) then
        write(*,*) ' Read temperature result failed; not heat analysis result'
      endif

      do i= 1, hecMESH%n_node
        temp(i) = result%node_val_item(i)
      enddo

      call hecmw_result_free(result)

      write(IDBG,*) ' Read temperature from result file : OK'
    end subroutine read_temperature_result
	
END MODULE
