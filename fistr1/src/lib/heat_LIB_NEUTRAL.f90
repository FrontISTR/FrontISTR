!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide functions to generate Femap neutral output
module m_heat_lib_neutral
contains

  !----------------------------------------------------------------------
  !     heat_put_neutral_601 ( hecMESH )                                l
  !     heat_put_neutral_402 ( hecMESH )                                l
  !     heat_put_neutral_403 ( hecMESH )                                l
  !     heat_put_neutral_404 ( hecMESH )                                l
  !     heat_put_neutral_409 ( hecMESH )                                l
  !     heat_put_neutral_450 ( hecMESH )                                l
  !     heat_put_neutral_451 ( hecMESH )                                l
  !----------------------------------------------------------------------
  subroutine heat_put_neutral_601 ( hecMESH )

    use m_fstr

    implicit none
    integer(kind=kint) :: icou, i, im1, im2, is
    real(kind=kreal) :: CD, CP, RHO, val
    real(kind=kreal) :: rdum(0:9)
    character(len=80) :: line

    type(hecmwST_local_mesh) :: hecMESH

    CP = 0.0d0
    CD = 0.0d0
    RHO = 0.0d0

    !C
    !C==put MATERIAL : ( BLOCK NO. = 601 )
    !C
    write(INEU,'(a)') '   -1'
    write(INEU,*) '  601'

    icou = 0
    do im1 = 1, hecMESH%material%n_mat
      line = ' '
      line(7:23) = ',-601,55,0,0,1,0,'
      write(line(1:6),'(i6)') im1
      write(INEU,'(a23)') line(1:23)
      write(INEU,*) '<NULL>'
      write(INEU,*) '10,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '25,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,'
      write(INEU,*) '200,'

      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'

      do im2= 1, 3
        icou = icou  + 1
        is= hecMESH%material%mat_TABLE_index(icou -1) + 1
        val  = hecMESH%material%mat_VAL (is)

        if (im2.eq.3) CD  = val
        if (im2.eq.2) CP  = val
        if (im2.eq.1) RHO = val
      enddo

      rdum( 0) = 0.d0
      rdum( 1) = 0.d0
      rdum( 2) = CD
      rdum( 3) = CD
      rdum( 4) = CD
      rdum( 5) = CD
      rdum( 6) = CD
      rdum( 7) = CD
      rdum( 8) = CP
      rdum( 9) = RHO

      write(INEU,'(1p,10(E9.2,'',''))') ( rdum(i),i=0,9 )

      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,'
      write(INEU,*) '50,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '70,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'

    enddo
    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_601

  !----------------------------------------------------------------------

  !C
  subroutine heat_put_neutral_402 ( hecMESH )

    use m_fstr

    implicit none
    integer(kind=kint) :: i, im
    integer(kind=kint) :: idum(0:9)

    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C==put Section : ( BLOCK NO. = 402 )
    !C
    write(INEU,'(a)') '   -1'
    write(INEU,*) '  402'

    do im = 1, hecMESH%section%n_sect
      idum(0) = im
      idum(1) = 110
      idum(2) = hecMESH%section%sect_mat_ID_item(im)
      if    ( hecMESH%section%sect_type(im) == 1 ) then
        idum(3) = 25
      elseif( hecMESH%section%sect_type(im) == 2 ) then
        idum(3) = 17
      elseif( hecMESH%section%sect_type(im) == 3 ) then
        idum(3) = 5
      elseif( hecMESH%section%sect_type(im) == 4 ) then
        idum(3) = 9
      endif
      idum(4) = 1
      idum(5) = 0

      write(INEU,'(6(i6,'',''))') ( idum(i),i=0,5 )
      write(INEU,*) '<NULL>'
      write(INEU,*) '0,0,0,0,'
      write(INEU,*) '90,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,'
      write(INEU,*) '0,0,'
      write(INEU,*) '190,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0.,0.,0.,0.,0.,'
      write(INEU,*) '0,'
      write(INEU,*) '0,'

    enddo

    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_402

  !----------------------------------------------------------------------
  !C
  subroutine heat_put_neutral_403 ( hecMESH )

    use m_fstr

    implicit none
    integer(kind=kint) :: i, inod
    real(kind=kreal)   :: xx, yy, zz 

    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C==put NODE : ( BLOCK NO. = 403 )
    !C
    write(INEU,'(a)') '   -1'
    write(INEU,*) '  403'

    do i = 1, hecMESH%n_node

      inod = hecMESH%global_node_ID(i)
      xx = hecMESH%node( 3*i-2 )
      yy = hecMESH%node( 3*i-1 )
      zz = hecMESH%node( 3*i   )

      write(INEU,'(i8,a,3(e15.7,'',''))')  &
        &          inod,',0,0,1,46,0,0,0,0,0,0,',xx,yy,zz
    enddo

    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_403

  !----------------------------------------------------------------------
  !C
  subroutine heat_put_neutral_404 ( hecMESH )

    use m_fstr

    implicit none
    integer(kind=kint) :: i, j, k, ij, jj, is, ie, ielm, icol, isid, isop, ietyp, istyp, itopo
    integer(kind=kint) :: nn(20), nna(10), nnb(10)

    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C==put ELEMENT : ( BLOCK NO. = 404 )
    !C

    write(INEU,'(a)') '   -1'
    write(INEU,*) '  404'

    do i = 1, hecMESH%n_elem

      ielm = hecMESH%global_elem_ID(i)
      icol = 124
      isid = hecMESH%section_ID(i)
      isop = hecMESH%section%sect_opt(isid)

      is = hecMESH%elem_node_index(i-1) + 1
      ie = hecMESH%elem_node_index(i)

      k = 0
      do j = is, ie
        k = k + 1
        jj    = hecMESH%elem_node_item(j)
        nn(k) = hecMESH%global_node_ID(jj)
      enddo

      nna = 0  ;  nnb = 0
      ietyp = hecMESH%elem_type(i)
      if( ietyp == 231 ) then
        istyp = 25
        if( isop == 1 ) istyp = 19
        if( isop == 2 ) istyp = 35
        itopo = 2
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
      elseif( ietyp == 731 ) then
        istyp = 17
        itopo = 2
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
      elseif( ietyp == 232 ) then
        istyp = 26
        if( isop == 1 ) istyp = 20
        if( isop == 2 ) istyp = 36
        itopo = 3
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(5) = nn(4)
        nna(6) = nn(5)
        nna(7) = nn(6)
      elseif( ietyp == 732 ) then
        istyp = 18
        itopo = 3
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(5) = nn(4)
        nna(6) = nn(5)
        nna(7) = nn(6)
      elseif( ietyp == 241 ) then
        istyp = 25
        if( isop == 1 ) istyp = 19
        if( isop == 2 ) istyp = 35
        itopo = 4
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(4) = nn(4)
      elseif( ietyp == 741 ) then
        istyp = 17
        itopo = 4
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(4) = nn(4)
      elseif( ietyp == 242 ) then
        istyp = 26
        if( isop == 1 ) istyp = 20
        if( isop == 2 ) istyp = 36
        itopo = 5
        do ij = 1, 8
          nna(ij) = nn(ij)
        enddo
      elseif( ietyp == 742 ) then
        istyp = 18
        itopo = 5
        do ij = 1, 8
          nna(ij) = nn(ij)
        enddo
      elseif( ietyp == 341 ) then
        istyp = 25
        itopo = 6
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(5) = nn(4)
      elseif( ietyp == 351 ) then
        istyp = 25
        itopo = 7
        nna(1) = nn(1)
        nna(2) = nn(2)
        nna(3) = nn(3)
        nna(5) = nn(4)
        nna(6) = nn(5)
        nna(7) = nn(6)
      elseif( ietyp == 361 ) then
        istyp = 25
        itopo = 8
        do ij = 1, 8
          nna(ij) = nn(ij)
        enddo
      elseif( ietyp == 342 ) then
        istyp = 26
        itopo = 10
        nna( 1) = nn( 1)
        nna( 2) = nn( 2)
        nna( 3) = nn( 3)
        nna( 5) = nn( 4)
        nna( 9) = nn( 5)
        nna(10) = nn( 6)
        nnb( 1) = nn( 7)
        nnb( 3) = nn( 8)
        nnb( 4) = nn( 9)
        nnb( 5) = nn(10)
      elseif( ietyp == 352 ) then
        istyp = 26
        itopo = 11
        nna( 1) = nn( 1)
        nna( 2) = nn( 2)
        nna( 3) = nn( 3)
        nna( 5) = nn( 4)
        nna( 6) = nn( 5)
        nna( 7) = nn( 6)
        nna( 9) = nn( 7)
        nna(10) = nn( 8)
        nnb( 1) = nn( 9)
        nnb( 3) = nn(13)
        nnb( 4) = nn(14)
        nnb( 5) = nn(15)
        nnb( 7) = nn(10)
        nnb( 8) = nn(11)
        nnb( 9) = nn(12)
      elseif( ietyp == 362 ) then
        istyp = 26
        itopo = 12
        do ij = 1, 10
          nna(ij) = nn(ij)
        enddo
        nnb( 1) = nn(11)
        nnb( 2) = nn(12)
        nnb( 3) = nn(17)
        nnb( 4) = nn(18)
        nnb( 5) = nn(19)
        nnb( 6) = nn(20)
        nnb( 7) = nn(13)
        nnb( 8) = nn(14)
        nnb( 9) = nn(15)
        nnb(10) = nn(16)
      endif

      write(INEU,'(5(i8,'',''),a)') &
        &          ielm,icol,isid,istyp,itopo,'1,0,0,0,0,0,0,0,'
      write(INEU,'(10(i8,'',''))') (nna(j),j=1,10)
      write(INEU,'(10(i8,'',''))') (nnb(j),j=1,10)
      write(INEU,*) '0,0,0,'
      write(INEU,*) '0,0,0,'
      write(INEU,*) '0,0,0,'
      write(INEU,*) '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'

    enddo

    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_404

  !----------------------------------------------------------------------
  !C
  subroutine heat_put_neutral_409

    use m_fstr
    !C
    !C==put View : BLOCK NO. = 409 )
    !C
    write(INEU,'(a)') '   -1'
    write(INEU,*) ' 409'
    write(INEU,*) '1,'
    write(INEU,*) 'Default XY View'
    write(INEU,*) '2,0,1,'
    !C      write(INEU,*) '35.2644,-45.,0.,'
    write(INEU,*) '0.,0.,0.,'
    write(INEU,*) '2.5,1.25,1.5,'
    write(INEU,*) '1.,1.,0,0.,0.,0.,0.,0.,0.,'
    write(INEU,*) '1.03572,0.51035,0.,'
    write(INEU,*) '1.2574,'
    write(INEU,*) '0.,0.,1.,1.,'
    write(INEU,*) '2,0,1,1,0,'
    write(INEU,*) '-1,-1,0,1,0,1,1,60031,0,4000000,'
    write(INEU,*) '9,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '0,0,0,'
    write(INEU,*) '100.,100.,1,7,'
    write(INEU,*) '0,1,1.,'
    write(INEU,*) '0.,0.,0.,'
    write(INEU,*) '0.,0.,1.,'
    write(INEU,*) '0,1,0,0,'
    write(INEU,*) '0.,0.,0.,'
    write(INEU,*) '1.,0.,0.,'
    write(INEU,*) '0.,0.,0.,'
    write(INEU,*) '0.,1.,0.,'
    write(INEU,*) '2.,1.,70.,0.5,'
    write(INEU,*) '0.,0.,0.,0.,0.,0,100.,1000.,0.,0.,0.,100.,-100.,1,'
    write(INEU,*) '5.,90.,10.,10.,1.,'
    write(INEU,*) '4,176,0,0,0,0,0,0,0,0,0.,0.,0.,'
    write(INEU,*) '0,0,0,0,0,0,14,110,'
    write(INEU,*) '0,1,1,1,0,1,0,1,1,0,1,1,0,1,1,1,0,0,1,'
    write(INEU,*) '0,0,0.00000001,25.,100.,0.,0.,0.,20,'
    write(INEU,*) '0,1,1,0,0,1,20.,0,'
    write(INEU,*) '12,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,0.,'
    write(INEU,*) '0,5,0,0,0,0.,25.,'
    write(INEU,*) '4,16408,20,16504,100,16488,'
    write(INEU,*) '0.,0.,'
    write(INEU,*) '0.,0.,0.,0.,'
    write(INEU,*) '9,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '1.,'
    write(INEU,*) '2,'
    write(INEU,*) '<NULL>'
    write(INEU,*) '<NULL>'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '0.,0.,0.,0.,'
    write(INEU,*) '0.,0.,0.,0.,'
    write(INEU,*) '90,1,124,1,0,'
    write(INEU,*) '0,60,0,0,'
    write(INEU,*) '0,24,0,0,'
    write(INEU,*) '0,100,0,0,'
    write(INEU,*) '0,2,0,0,'
    write(INEU,*) '0,24580,0,0,'
    write(INEU,*) '0,124,0,0,'
    write(INEU,*) '0,46,0,0,'
    write(INEU,*) '0,120,0,0,'
    write(INEU,*) '0,124,0,1,'
    write(INEU,*) '0,124,0,0,'
    write(INEU,*) '0,12,0,1,'
    write(INEU,*) '0,62,0,0,'
    write(INEU,*) '0,62,0,0,'
    write(INEU,*) '0,10,0,0,'
    write(INEU,*) '0,52,0,0,'
    write(INEU,*) '0,4,0,0,'
    write(INEU,*) '0,120,0,0,'
    write(INEU,*) '0,12,0,0,'
    write(INEU,*) '0,2,0,0,'
    write(INEU,*) '0,120,0,0,'
    write(INEU,*) '0,8312,0,0,'
    write(INEU,*) '0,24600,0,0,'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '1,74,0,1,'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '3,124,0,1,'
    write(INEU,*) '0,24636,0,0,'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '0,4,0,0,'
    write(INEU,*) '0,100,0,0,'
    write(INEU,*) '0,124,0,1,'
    write(INEU,*) '0,60,0,1,'
    write(INEU,*) '0,56,0,1,'
    write(INEU,*) '0,24,0,0,'
    write(INEU,*) '0,8216,0,1,'
    write(INEU,*) '0,4,0,0,'
    write(INEU,*) '0,124,2,0,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,0,1,'
    write(INEU,*) '1,124,5,1,'
    write(INEU,*) '0,0,0,1,'
    write(INEU,*) '0,24,0,1,'
    write(INEU,*) '0,124,0,0,'
    write(INEU,*) '0,100,0,1,'
    write(INEU,*) '1,100,0,1,'
    write(INEU,*) '0,0,0,1,'
    write(INEU,*) '0,16,0,0,'
    write(INEU,*) '0,124,4,1,'
    write(INEU,*) '0,62,0,0,'
    write(INEU,*) '2,124,1,1,'
    write(INEU,*) '1,8254,0,0,'
    write(INEU,*) '0,124,1,1,'
    write(INEU,*) '1,0,5,1,'
    write(INEU,*) '0,124,0,1,'
    write(INEU,*) '0,100,0,1,'
    write(INEU,*) '0,100,0,1,'
    write(INEU,*) '1,46,0,1,'
    write(INEU,*) '1,120,0,1,'
    write(INEU,*) '1,4,0,1,'
    write(INEU,*) '1,52,0,1,'
    write(INEU,*) '1,24,0,1,'
    write(INEU,*) '1,93,0,1,'
    write(INEU,*) '1,12,0,1,'
    write(INEU,*) '1,10,0,1,'
    write(INEU,*) '1,104,0,1,'
    write(INEU,*) '0,100,0,0,'
    write(INEU,*) '0,24,0,0,'
    write(INEU,*) '0,60,0,0,'
    write(INEU,*) '0,104,0,0,'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,4,0,0,'
    write(INEU,*) '0,0,1,0,'
    write(INEU,*) '0,0,0,0,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,62,1,1,'
    write(INEU,*) '0,60,4,0,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,*) '0,0,1,1,'
    write(INEU,'(a)') '-1,'
    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_409

  !----------------------------------------------------------------------
  !C
  subroutine heat_put_neutral_450

    use m_fstr
    !C
    !C==put RESULT : BLOCK NO. = 450 )
    !C
    write(INEU,'(a)') '   -1'
    write(INEU,*) '  450'
    write(INEU,*) '1,'
    write(INEU,*) 'hecmw_FSTR_heat_result'
    write(INEU,*) '0,0,'
    write(INEU,*) '0.,'
    write(INEU,*) '1,'
    write(INEU,*) '<NULL>'
    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_450

  !----------------------------------------------------------------------
  !C
  subroutine heat_put_neutral_451 ( hecMESH, hecHEAT )

    use m_fstr

    implicit none
    integer(kind=kint) :: i, inod
    real(kind=kreal)   :: tmin, tmax, tt, absmin, absmax

    type(fstr_heat)          :: hecHEAT
    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C==put DISPLACEMENT : ( BLOCK NO. = 451 )
    !C    & STRESS/STRAIN
    !C

    tmin = 1.0e9
    tmax =-1.0e9

    !C

    do i = 1, hecMESH%n_node

      tt = hecHEAT%TEMP(i)
      if( tt > tmax ) tmax = tt
      if( tt < tmin ) tmin = tt

    enddo

    absmax = dabs(tmax)
    absmin = dabs(tmin)
    if( absmin > absmax ) absmax = absmin

    write(INEU,'(a)') '   -1'
    write(INEU,*) '  451'

    write(INEU,*) '1,1,1,'
    write(INEU,*) 'Temperature'
    write(INEU,'(3(e15.7,'',''))') tmin, tmax, absmax
    write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
    write(INEU,*) '0,0,0,0,0,0,0,0,0,0,'
    write(INEU,*) '0,0,6,7,'
    write(INEU,*) '0,0,1,'

    do i = 1, hecMESH%n_node
      inod = hecMESH%global_node_ID(i)
      tt = hecHEAT%TEMP(i)
      write(INEU,'(i8,'','',E15.7,'','')') inod,tt
    enddo
    write(INEU,*) '-1,0.'
    !C
    !C==TERMINATION
    !C
    write(INEU,'(a)') '   -1'

    return
  end subroutine heat_put_neutral_451
end module m_heat_lib_neutral

