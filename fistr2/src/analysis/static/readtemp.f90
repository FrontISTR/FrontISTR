MODULE mReadTemp

  contains

!> Read in temperature distribution from external file
    subroutine read_temperature_result(hecMESH, tstep, temp)
      use hecmw_result
	  use hecmw
      type(hecmwST_local_mesh) :: hecMESH
      integer(kind=kint) :: tstep
      real(kind=kreal),pointer :: temp(:)
      type(hecmwST_result_data) :: result
      integer(kind=kint) :: i

      call hecmw_result_read(tstep, result)

      if (result%nn_component /= 1 .or. result%nn_dof(1) /= 1) then
        write(*,*) ' Read temperature result failed; not heat analysis result'
      endif

      do i= 1, hecMESH%n_node
        temp(i) = result%node_val_item(i)
      enddo

      call hecmw_result_free(result)

      write(IDBG,*) ' Read temperature from result file: OK'
    end subroutine read_temperature_result
	
END MODULE
