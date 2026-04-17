!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides the entry point for ELEMCHECK (pre-analysis input validation)
module m_fstr_precheck

  use hecmw
  use m_fstr
  use m_inputcheck_mesh_quality
  use m_inputcheck_make_result
  use m_hecmw2fstr_mesh_conv

  implicit none

contains

  !> Main entry point (called from fistr_main)
  subroutine fstr_input_precheck(hecMESH, hecMAT, fstrSOLID)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(in)        :: hecMAT
    type(fstr_solid), intent(in)            :: fstrSOLID

    integer(kind=kint) :: i, cid
    real(kind=kreal), allocatable :: elem_vol(:), elem_asp(:)
    type(hecmwST_result_data) :: fstrRESULT

    ! Density check for eigen/dynamic analysis
    if(fstrPR%solution_type == kstEIGEN .or. &
       fstrPR%solution_type == kstSTATICEIGEN .or. &
       fstrPR%solution_type == kstDYNAMIC) then
      do i = 1, hecMESH%section%n_sect
        if(hecMESH%section%sect_type(i) == 4) cycle
        cid = hecMESH%section%sect_mat_ID_item(i)
        if(fstrSOLID%materials(cid)%variables(M_DENSITY) == 0.0d0) then
          write(*,*) "*** error: density is not assigned or set to zero"
          call hecmw_abort(hecmw_comm_get_comm())
        endif
      enddo
    endif

    if(fstrPR%solution_type == kstPRECHECK) then
      if(myrank == 0) then
        write(IMSG,*)
        write(IMSG,*) ' ****   STAGE Precheck  **'
      endif

      ! Allocate element quality arrays
      allocate(elem_vol(hecMESH%n_elem))
      allocate(elem_asp(hecMESH%n_elem))

      ! Mesh quality check
      call inputcheck_mesh_quality(hecMESH, hecMAT, elem_vol, elem_asp)
      write(IDBG,*) 'inputcheck_mesh_quality: OK'

      ! Build visualization result data and output
      call inputcheck_make_result(hecMESH, fstrRESULT, elem_vol, elem_asp)
      if(IVISUAL == 1) then
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize(hecMESH, fstrRESULT, 1)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        write(IDBG,*) 'elemcheck visualization: OK'
      endif
      call hecmw_result_free(fstrRESULT)

      ! Write result file
      if(IRESULT == 1) then
        call inputcheck_write_result(hecMESH, elem_vol, elem_asp)
        write(IDBG,*) 'elemcheck result file: OK'
      endif

      deallocate(elem_vol)
      deallocate(elem_asp)
    endif

    if(fstrPR%solution_type == kstNZPROF) then
      call hecmw_nonzero_profile(hecMESH, hecMAT)
    endif
  end subroutine fstr_input_precheck

  subroutine fstr_get_thickness(hecMESH, mid, thick)
    use hecmw
    use m_fstr
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: mid, ihead
    real(kind=kreal)   :: thick

    ihead = hecMESH%section%sect_R_index(mid-1)
    thick = hecMESH%section%sect_R_item(ihead+1)
  end subroutine fstr_get_thickness

  !> Non-zero profile output for gnuplot visualization
  subroutine hecmw_nonzero_profile(hecMESH, hecMAT)
    use hecmw_util
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT

    integer(kind=kint) :: i, j, in, jS, jE, ftype, n, ndof, nnz, fio
    real(kind=kreal) :: rnum, dens, cond
    character :: fileid*3

    fio = 70 + hecMESH%my_rank
    write(fileid,"(i3.3)")hecMESH%my_rank

    ftype = 4

    n = hecMAT%N
    ndof = 3*n
    nnz = 9*n + 9*2*hecMAT%indexL(n)
    dens = 100*dble(nnz)/dble(9*n*n)
    rnum = (7.21d0+0.01*dlog10(dble(hecMAT%N)))*10.0d0/dble(hecMAT%N)
    cond = 1.0d0

    open(fio,file='nonzero.dat.'//fileid, status='replace')
    do i= 1, n
      jS= hecMAT%indexL(i-1) + 1
      jE= hecMAT%indexL(i  )
      write(fio,"(i0,a,i0)")i,"  ",i
      do j= jS, jE
        in = hecMAT%itemL(j)
        write(fio,"(i0,a,i0)")i, "  ",in
        write(fio,"(i0,a,i0)")in,"  ",i
      enddo
    enddo
    close(fio)

    open(fio,file='nonzero.plt.'//fileid, status='replace')
    if(ftype == 4)then
      write(fio,"(a)")'set terminal png size 1500,1500'
    else
      write(fio,"(a)")'set terminal postscript eps enhanced color solid "TimesNewRomanPSMT" 20'
    endif

    write(fio,"(a)")'unset key'
    write(fio,"(a)")'unset xtics'
    write(fio,"(a)")'unset ytics'
    write(fio,"(a)")'set size ratio 1.0'
    write(fio,"(a)")'set border lw 1.0'
    write(fio,"(a,i0,a)")'set xrange[0.5:',n,'.5]'
    write(fio,"(a,i0,a)")'set yrange[0.5:',n,'.5] reverse '

    if(ftype == 4)then
      write(fio,"(a)")'set out "image.'//fileid//'.png"'
    else
      write(fio,"(a)")'set out "image.'//fileid//'.eps"'
      write(fio,"(a)"     )'set label 1 "Name" at graph 1.1,0.9'
      write(fio,"(a)")'set label 2 "N" at graph 1.1,0.85'
      write(fio,"(a)")'set label 3 "Non-Zero Elem." at graph 1.1,0.8'
      write(fio,"(a)")'set label 4 "Density [%]" at graph 1.1,0.75'
      write(fio,"(a)")'set label 9 "Condition Num." at graph 1.1,0.7'
      write(fio,"(a)"     )'set label 5 ":  matrix" at graph 1.4,0.9'
      write(fio,"(a,i0,a)")'set label 6 ":  ',ndof,'" at graph 1.4,0.85'
      write(fio,"(a,i0,a)")'set label 7 ":  ',nnz,'" at graph 1.4,0.8'
      write(fio,"(a,1pe9.2,a)")'set label 8 ": ',dens,'" at graph 1.4,0.75'
      write(fio,"(a,1pe9.2,a)")'set label 10 ": ',cond,'" at graph 1.4,0.7'
    endif

    write(fio,"(a,f12.5,a)")'plot "nonzero.dat.'//fileid//'" pointtype 5 pointsize ',rnum,' linecolor rgb "#F96566"'
    close(fio)

    write(*,*)''
    write(*,*)' ### Command recommendation'
    write(*,*)' gnuplot -persist "nonzero.plt"'
  end subroutine hecmw_nonzero_profile

end module m_fstr_precheck
