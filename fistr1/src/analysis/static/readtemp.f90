!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module mReadTemp

  private
  public :: read_temperature_result

contains

  !> Read in temperature distribution from external file
  subroutine read_temperature_result(hecMESH, nstep, sstep, rtype, interval, factor, ctime, temp, temp_bak)
    use m_fstr
    type(hecmwST_local_mesh), intent(in) :: hecMESH  !< hecmw mesh
    integer(kind=kint), intent(in) :: nstep          !< num of total steps of heat analysis
    integer(kind=kint), intent(in) :: sstep          !< step to start reading
    integer(kind=kint), intent(in) :: rtype          !< read type; 1: step-based, 2: time-based
    integer(kind=kint), intent(in) :: interval       !< interval of reading next result
    real(kind=kreal), intent(in) :: factor           !< readtemp factor
    real(kind=kreal), intent(in) :: ctime            !< current time
    real(kind=kreal), pointer :: temp(:)             !< temperature result
    real(kind=kreal), pointer :: temp_bak(:)         !< temperature at previous step
    if( rtype == 1 ) then
      call read_temperature_result_by_step(hecMESH, nstep, sstep, interval, factor, temp, temp_bak)
    elseif( rtype == 2 ) then
      call read_temperature_result_by_time(hecMESH, nstep, sstep, ctime, temp)
    endif
    write(IDBG,*) ' Read temperature from result file : OK'
  end subroutine read_temperature_result

  subroutine read_temperature_result_by_step(hecMESH, nstep, sstep, interval, factor, temp, temp_bak)
    use hecmw_result
    use m_fstr
    type(hecmwST_local_mesh), intent(in) :: hecMESH  !< hecmw mesh
    integer(kind=kint), intent(in) :: nstep          !< num of total steps of heat analysis
    integer(kind=kint), intent(in) :: sstep          !< step to start reading
    integer(kind=kint), intent(in) :: interval       !< interval of reading next result
    real(kind=kreal), intent(in) :: factor           !< readtemp factor
    real(kind=kreal), pointer :: temp(:)             !< temperature result
    real(kind=kreal), pointer :: temp_bak(:)         !< temperature at previous step
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
  end subroutine read_temperature_result_by_step

  !> Read in temperature distribution from external file
  subroutine read_temperature_result_by_time(hecMESH, nstep, sstep, ctime, temp)
    use hecmw_result
    use m_fstr
    type(hecmwST_local_mesh), intent(in) :: hecMESH  !< hecmw mesh
    integer(kind=kint), intent(in) :: nstep          !< num of total steps of heat analysis
    integer(kind=kint), intent(in) :: sstep          !< step to start reading
    real(kind=kreal), intent(in) :: ctime            !< current time
    real(kind=kreal), pointer :: temp(:)             !< temperature result

    logical, save :: flg_first_call = .true.
    logical, save :: flg_constant = .false.
    logical, save :: flg_finished = .false.
    integer(kind=kint), save ::  tstep1, tstep2
    real(kind=kreal), save :: ttime1, ttime2

    integer(kind=kint) :: ierr
    real(kind=kreal) :: w
    logical :: flg_read1, flg_read2
    real(kind=kreal), allocatable :: temp0(:,:)

    allocate(temp0(2,size(temp)))
    flg_read1 = .false.
    flg_read2 = .false.
    ! read first 2 results at first call
    if (flg_first_call) then
      call read_next_result(hecMESH, nstep, sstep, tstep1, ttime1, temp0(1,:), ierr)
      !  if no more result; error: no result file
      if (ierr /= 0) then
        write(*,*) ' Read temperature result failed; no result file'
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      flg_read1 = .true.
      call read_next_result(hecMESH, nstep, tstep1+1, tstep2, ttime2, temp0(2,:), ierr)
      !  if no more result (only 1 result exists); assume temperature is constant
      if (ierr /= 0) flg_constant = .true.
      flg_read2 = .true.
      flg_first_call = .false.
    endif
    ! when temperature is assumed constant
    if (flg_constant) then
      if (.not. flg_read1) call read_result_at(hecMESH, tstep1, ttime1, temp0(1,:))
      temp(:) = temp0(1,:)
      return
    endif
    ! when ctime is out-of-range (left side)
    if (ctime <= ttime1) then
      if (.not. flg_read1) call read_result_at(hecMESH, tstep1, ttime1, temp0(1,:))
      temp(:) = temp0(1,:)
      return
    endif
    ! proceed until ctime is within [ttime1, ttime2]
    do while (ttime2 < ctime)
      tstep1 = tstep2
      ttime1 = ttime2
      if (flg_read2) then
        temp0(1,:) = temp0(2,:)
        flg_read1 = .true.
      endif
      flg_read2 = .false.
      call read_next_result(hecMESH, nstep, tstep1+1, tstep2, ttime2, temp0(2,:), ierr)
      !  if no more result; finished reading (ctime is out-of-range)
      if (ierr /= 0) then
        flg_finished = .true.
        exit
      endif
      flg_read2 = .true.
    enddo
    ! when ctime is out-of-range (right side)
    if (flg_finished) then
      if (.not. flg_read1) call read_result_at(hecMESH, tstep1, ttime1, temp0(1,:))
      temp(:) = temp0(1,:)
      return
    endif
    ! now, ttime1 <= ctime <= ttime2
    if (.not. flg_read1) call read_result_at(hecMESH, tstep1, ttime1, temp0(1,:))
    if (.not. flg_read2) call read_result_at(hecMESH, tstep2, ttime2, temp0(2,:))
    w = (ctime - ttime1) / (ttime2 - ttime1)
    temp(:) = (1-w)*temp0(1,:) + w*temp0(2,:)
    !write(IDBG,'(a,f10.2,a,i0,a,f10.2,a,i0,a,f10.2,a)') &
    !     ' Interpolated temp at',ctime,' from step ',tstep1,' (time',ttime1,') and step ',tstep2,' (time',ttime2,')'
    return
  end subroutine read_temperature_result_by_time

  subroutine read_result_at(hecMESH, tstep, ttime, temp)
    use m_fstr
    use hecmw_result
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: tstep
    real(kind=kreal), intent(out) :: ttime
    real(kind=kreal), intent(out) :: temp(:)
    type(hecmwST_result_data) :: result
    character(len=HECMW_NAME_LEN) :: name_ID
    integer(kind=kint) :: i
    call hecmw_nullify_result_data( result )
    name_ID = 'fstrTEMP'
    call hecmw_result_read_by_name(hecMESH, name_ID, tstep, result)
    if (result%nn_component /= 1 .or. result%nn_dof(1) /= 1) then
      write(*,*) ' Read temperature result failed; not heat analysis result'
    endif
    ttime = result%global_val_item(1)
    do i=1, hecMESH%n_node
      temp(i) = result%node_val_item(i)
    enddo
    call hecmw_result_free(result)
  end subroutine read_result_at

  subroutine read_next_result(hecMESH, nstep, sstep, tstep, ttime, temp, ierr)
    use m_fstr
    use hecmw_result
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: nstep
    integer(kind=kint), intent(in) :: sstep
    integer(kind=kint), intent(out) :: tstep
    real(kind=kreal), intent(out) :: ttime
    real(kind=kreal), intent(out) :: temp(:)
    integer(kind=kint), intent(out) :: ierr
    character(len=HECMW_NAME_LEN) :: name_ID
    integer(kind=kint) :: istep
    name_ID = 'fstrTEMP'
    tstep = -1
    do istep = sstep, nstep
      call hecmw_result_checkfile_by_name(name_ID, istep, ierr)
      if (ierr == 0) then
        tstep = istep
        exit
      endif
    enddo
    if (tstep < 0) then
      ierr = 1
      return
    endif
    call read_result_at(hecMESH, tstep, ttime, temp)
  end subroutine read_next_result

end module mReadTemp
