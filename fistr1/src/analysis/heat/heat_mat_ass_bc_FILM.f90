!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine for setting heat transfer
!! boundary conditions
module m_heat_mat_ass_bc_FILM
contains
  !C
  !C***
  !C*** MAT_ASS_FILM
  !C***
  !C
  subroutine heat_mat_ass_bc_FILM( hecMESH, hecMAT, fstrHEAT, CTIME, DTIME, beta )

    use m_fstr
    use m_heat_get_amplitude
    use m_heat_LIB_FILM

    implicit none
    integer(kind=kint) :: k, icel, isuf, iam1, iam2, ic_type, isect, nn, is, j, mm, m, ic, ip
    integer(kind=kint) :: inod, jp, jnod, isU, ieU, ik, isL, ieL
    real(kind=kreal)   :: CTIME, DTIME, QQ, HH, SINK, thick, beta
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH

    real(kind=kreal)   :: xx(20), yy(20), zz(20)
    real(kind=kreal)   :: term1(64), term2(20), stiff(8,8)
    integer(kind=kint) :: nodLocal(20), nsuf(8), nodSurf(8)

!C
    !$omp parallel default(none), &
      !$omp&  private(k,icel,isuf,iam1,iam2,QQ,HH,SINK,ic_type,isect,nn,is,j,nodLocal, &
      !$omp&  xx,yy,zz,thick,mm,term1,term2,stiff,nsuf,nodSurf,ip,inod,jnod,ic,isU,ieU,ik,jp,isL,ieL), &
      !$omp&  shared(fstrHEAT,CTIME,hecMAT,hecMESH)
    !$omp do
    do k = 1, fstrHEAT%H_SUF_tot
      icel    = fstrHEAT%H_SUF_elem(k)
      isuf    = fstrHEAT%H_SUF_surf(k)
      iam1    = fstrHEAT%H_SUF_ampl(k,1)
      iam2    = fstrHEAT%H_SUF_ampl(k,2)
      call heat_get_amplitude ( fstrHEAT, iam1, CTIME, QQ )
      HH      = fstrHEAT%H_SUF_val (k,1) * QQ
      call heat_get_amplitude ( fstrHEAT, iam2, CTIME, QQ )
      SINK    = fstrHEAT%H_SUF_val (k,2) * QQ
      ic_type = hecMESH%elem_type(icel)
      isect   = hecMESH%section_ID(icel)
      !C**
      nn = hecmw_get_max_node(ic_type)
      !C**
      is = hecMESH%elem_node_index(icel-1)
      do j = 1, nn
        nodLOCAL(j) = hecMESH%elem_node_item(is+j)
        xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
        yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
        zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
      enddo
      !C**
      if    ( ic_type.eq.231 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm = 2
        call heat_FILM_231(nn,xx,yy,zz,thick,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.232 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm = 3
        call heat_FILM_232(nn,xx,yy,zz,thick,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.241 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm = 2
        call heat_FILM_241(nn,xx,yy,zz,thick,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.242 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        mm = 3
        call heat_FILM_242(nn,xx,yy,zz,thick,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.341 ) then
        mm = 3
        call heat_FILM_341(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.342 ) then
        mm = 6
        call heat_FILM_342(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.351 ) then
        mm = 4
        if( isuf.eq.1 .or. isuf.eq.2 ) mm = 3
        call heat_FILM_351(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.352 ) then
        mm = 8
        if( isuf.eq.1 .or. isuf.eq.2 ) mm = 6
        call heat_FILM_352(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.361 ) then
        mm = 4
        call heat_FILM_361(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.362 ) then
        mm = 8
        call heat_FILM_362(nn,xx,yy,zz,isuf,HH,SINK,mm,term1,term2,nsuf)
      elseif( ic_type.eq.731 ) then
        mm = 3
        call heat_FILM_731(nn,xx,yy,zz,HH,SINK,term1,term2)
        do m = 1, mm
          nsuf(m) = m
        enddo

      elseif( ic_type.eq.741 ) then
        mm = 4
        call heat_FILM_741(nn,xx,yy,zz,HH,SINK,term1,term2)
        do m = 1, mm
          nsuf(m) = m
        enddo

      endif
      !C

      do ip = 1, mm
        nodSurf(ip) = nodLOCAL(nsuf(ip))
      end do

      ic = 0
      stiff = 0.d0
      do ip = 1, mm
        do jp = 1, mm
          ic = ic + 1
          stiff(ip,jp) = -term1(ic)
        enddo
      enddo

      call hecmw_mat_ass_elem(hecMAT, mm, nodSurf, stiff)

      do ip = 1, mm
        !$omp atomic
        hecMAT%B(nodSurf(ip)) = hecMAT%B(nodSurf(ip)) - term2(ip)
      end do

      !C
      !C
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine  heat_mat_ass_bc_FILM
end module m_heat_mat_ass_bc_FILM
