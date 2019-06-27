!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine for setting distributed
!!  heat flux boundary conditions
module m_heat_mat_ass_bc_DFLUX
  use m_fstr

  implicit none
contains
  !C
  !C***
  !C*** MAT_ASS_DFLUX
  !C***
  !C
  subroutine heat_mat_ass_bc_DFLUX( hecMESH, hecMAT, fstrHEAT, CTIME, DTIME, beta )
    use m_heat_get_amplitude
    use m_heat_LIB_DFLUX
    use m_static_LIB_3d
    integer(kind=kint) :: k, icel, ic_type, isect, isuf, iamp, nn, is, j
    real(kind=kreal)   :: CTIME, DTIME, QQ, val, asect, thick, vol, beta
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    real(kind=kreal)   :: xx(20), yy(20), zz(20)
    real(kind=kreal)   :: vect(20), ss(2000)
    integer(kind=kint) :: ig0, ig, iS0, iE0, nodLocal(20)
    real(kind=kreal), allocatable :: Bbak(:)

    !C
    !$omp parallel default(none), &
      !$omp&  private(k,icel,ic_type,isect,isuf,iamp,QQ,val,nn,is,j,nodLOCAL,xx,yy,zz,asect,vect,thick), &
      !$omp&  shared(fstrHEAT,CTIME,hecMAT,hecMESH)
    !$omp do
    do k = 1, fstrHEAT%Q_SUF_tot

      icel    = fstrHEAT%Q_SUF_elem(k)
      ic_type = hecMESH%elem_type(icel)
      isect   = hecMESH%section_ID(icel)
      isuf    = fstrHEAT%Q_SUF_surf(k)
      iamp    = fstrHEAT%Q_SUF_ampl(k)

      call heat_get_amplitude( fstrHEAT,iamp,CTIME,QQ )
      val     = fstrHEAT%Q_SUF_val (k) * QQ
      if( dabs(val) < 1.d-16 ) cycle
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
      if    ( ic_type.eq.111 ) then
        is = hecMesh%section%sect_R_index(isect)
        asect = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_111(nn,xx,yy,zz,asect,isuf,val,vect)

      elseif( ic_type.eq.231 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_231(nn,xx,yy,zz,thick,isuf,val,vect)

      elseif( ic_type.eq.232 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_232(nn,xx,yy,zz,thick,isuf,val,vect)

      elseif( ic_type.eq.241 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_241(nn,xx,yy,zz,thick,isuf,val,vect)

      elseif( ic_type.eq.242 ) then
        is = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_242(nn,xx,yy,zz,thick,isuf,val,vect)

      elseif( ic_type.eq.341 ) then
        call heat_DFLUX_341(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.342 ) then
        call heat_DFLUX_342(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.351 ) then
        call heat_DFLUX_351(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.352 ) then
        call heat_DFLUX_352(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.361 ) then
        call heat_DFLUX_361(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.362 ) then
        call heat_DFLUX_362(nn,xx,yy,zz,isuf,val,vect)

      elseif( ic_type.eq.731 ) then
        is    = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_731(nn,xx,yy,zz,thick,isuf,val,vect)

      elseif( ic_type.eq.741 ) then
        is    = hecMesh%section%sect_R_index(isect)
        thick = hecMESH%section%sect_R_item(is)
        call heat_DFLUX_741(nn,xx,yy,zz,thick,isuf,val,vect)

      endif

      do j = 1, nn
        !$omp atomic
        hecMAT%B( nodLOCAL(j) ) = hecMAT%B( nodLOCAL(j) ) - vect(j)
      enddo
    enddo
    !$omp end do
    !$omp end parallel

    if( fstrHEAT%WL_tot>0 ) then
      allocate( Bbak(size(hecMAT%B)) )

      do ig0 = 1, fstrHEAT%WL_tot
        ig = fstrHEAT%weldline(ig0)%egrpid
        iS0= hecMESH%elem_group%grp_index(ig-1) + 1
        iE0= hecMESH%elem_group%grp_index(ig  )
        vol =0.d0; Bbak=0.d0
        do k=iS0, iE0
          icel   = hecMESH%elem_group%grp_item(k)
          ic_type = hecMESH%elem_type(icel)
          isect   = hecMESH%section_ID(icel)
          isuf    = 0

          nn = hecmw_get_max_node(ic_type)
          is = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item(is+j)
            xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
            yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
            zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
          enddo

          if( .not. isWeldActive(nn, xx,yy,zz, fstrHEAT%weldline(ig0), ctime-0.5d0*dtime) ) cycle

          vol = vol + VOLUME_C3(ic_type,NN,XX,YY,ZZ)
          val = fstrHEAT%weldline(ig0)%I * fstrHEAT%weldline(ig0)%U *   &
            fstrHEAT%weldline(ig0)%coe
          write(IDBG,*) "Element", hecMESH%global_elem_id(icel),"with dflux", val, vol

          if( ic_type.eq.341 ) then
            call heat_DFLUX_341(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.342 ) then
            call heat_DFLUX_342(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.351 ) then
            call heat_DFLUX_351(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.352 ) then
            call heat_DFLUX_352(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.361 ) then
            call heat_DFLUX_361(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.362 ) then
            call heat_DFLUX_362(nn,xx,yy,zz,isuf,val,vect)


          else
            write(*,*) '###ERROR### : Element type not supported for heat analysis'
            write(*,*) ' ic_type = ', ic_type
            call hecmw_abort(hecmw_comm_get_comm())

          endif

          do j = 1, nn
            Bbak( nodLOCAL(j) ) = Bbak( nodLOCAL(j) ) - vect(j)
          enddo
        enddo

        if( vol>0 ) then
          Bbak = Bbak/vol
          hecMAT%B = hecMAT%B + Bbak
        endif
      enddo

      deallocate(Bbak)
    endif

  end subroutine heat_mat_ass_bc_DFLUX

  logical function isWeldActive(nn, xx,yy,zz, weldline, ctime)
    integer, intent(in)          ::  nn
    real(kind=kreal), intent(in) :: xx(:), yy(:), zz(:)
    type(tWeldLine), intent(in)  :: weldline
    real(kind=kreal), intent(in) :: ctime

    integer :: i
    real(kind=kreal) :: cpos, wpos, tend

    isWeldActive = .false.
    tend = weldline%tstart+ (weldline%n2-weldline%n1)/weldline%v
    if( ctime<weldline%tstart .or. ctime>tend ) return
    cpos = 0
    do i=1,nn
      if( weldline%xyz==1 ) then
        cpos = cpos+xx(i)
      elseif( weldline%xyz==2 ) then
        cpos = cpos+yy(i)
      else
        cpos = cpos+zz(i)
      endif
    enddo
    cpos = cpos/nn
    wpos = weldline%n1 + weldline%v*(ctime-weldline%tstart)
    if( dabs(cpos-wpos)<weldline%distol ) isWeldActive = .true.
  end function
end module m_heat_mat_ass_bc_DFLUX
