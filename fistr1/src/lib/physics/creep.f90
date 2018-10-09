!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions for creep calculation
module mCreep
  use hecmw_util
  use mMaterial
  use m_ElasticLinear

  implicit none

contains

  !> This subroutine calculates stiffness for elastically isotropic
  !>     materials with isotropic creep
  subroutine iso_creep(matl, sectType, stress, strain, extval,plstrain,            &
      dtime,ttime,stiffness, temp)
    type( tMaterial ), intent(in)    :: matl      !< material properties
    integer, intent(in)              :: sectType  !< not used currently
    real(kind=kreal), intent(in)     :: stress(6) !< Piola-Kirchhoff stress
    real(kind=kreal), intent(in)     :: strain(6) !< strain
    real(kind=kreal), intent(in)     :: extval(:) !< plastic strain
    real(kind=kreal), intent(in)     :: plstrain  !< plastic strain increment
    real(kind=kreal), intent(in)     :: ttime     !< total time at the start of the current increment
    real(kind=kreal), intent(in)     :: dtime     !< time length of the increment
    real(kind=kreal), intent(out)    :: stiffness(6,6) !< stiffness
    real(kind=kreal), optional       :: temp      !> temprature

    integer :: i, j
    logical :: ierr
    real(kind=kreal) :: ina(1), outa(3)
    real(kind=kreal) :: xxn, aa

    real(kind=kreal) :: c3,e,un,G,ddg,stri(6),p,dstri,c4,c5,f,df, eqvs

    aa = 0.0d0
    xxn = 0.0d0
    !
    !     elastic
    !
    call calElasticMatrix( matl, sectTYPE, stiffness )
    if( dtime==0.d0 .or. all(stress==0.d0) ) return
    !
    !     elastic constants
    !
    if( present(temp) ) then
      ina(1) = temp
      call fetch_TableData( MC_ISOELASTIC, matl%dict, outa(1:2), ierr, ina )
    else
      call fetch_TableData( MC_ISOELASTIC, matl%dict, outa(1:2), ierr )
    endif
    if( ierr ) then
      stop "error in isotropic elasticity definition"
    else
      e=outa(1)
      un=outa(2)
    endif

    !      Norton
    if( matl%mtype==NORTON ) then         ! those with no yield surface
      if( present( temp ) ) then
        ina(1) = temp
        call fetch_TableData( MC_NORTON, matl%dict, outa, ierr, ina )
      else
        call fetch_TableData( MC_NORTON, matl%dict, outa, ierr )
      endif
      xxn=outa(2)
      aa=outa(1)*((ttime+dtime)**(outa(3)+1.d0)-ttime**(outa(3)+1.d0))/(outa(3)+1.d0)
    endif

    G=0.5d0*e/(1.d0+un)
    !
    !     creep
    !
    stri(:)=stress(:)
    p=-(stri(1)+stri(2)+stri(3))/3.d0
    do i=1,3
      stri(i)=stri(i)+p
    enddo
    !
    dstri=dsqrt(1.5d0*(stri(1)*stri(1)+stri(2)*stri(2)+stri(3)*stri(3)+    &
      2.d0*(stri(4)*stri(4)+stri(5)*stri(5)+stri(6)*stri(6))) )
    !
    !     unit trial vector
    !
    stri(:)=stri(:)/dstri

    eqvs = dstri
    if( eqvs<1.d-10 ) eqvs=1.d-10
    f=aa*eqvs**xxn
    df=xxn*f/eqvs
    !
    !     stiffness matrix
    !
    c3=6.d0*G*G
    c4=c3*plstrain/(dstri+3.d0*G*plstrain)
    c3=c4-c3*df/(3.d0*G*df+1.d0)
    c5=c4/3.d0

    do i=1,6
      do j=1,6
        stiffness(i,j) = stiffness(i,j) +c3*stri(i)*stri(j)
      enddo
    enddo
    do i=1,3
      stiffness(i,i) = stiffness(i,i) - c4
      do j=1,3
        stiffness(i,j) = stiffness(i,j) +c5
      enddo
    enddo
    do i=4,6
      stiffness(i,i) = stiffness(i,i) - c4/2.d0
    enddo

  end subroutine

  !> This subroutine calculates stresses and creep status for an elastically isotropic
  !>     material with isotropic creep
  subroutine update_iso_creep(matl, sectType, strain, stress, extval,plstrain,                &
      dtime,ttime,temp)
    type( tMaterial ), intent(in)    :: matl      !< material properties
    integer, intent(in)              :: sectType  !< not used currently
    real(kind=kreal), intent(in)     :: strain(6) !< strain
    real(kind=kreal), intent(inout)  :: stress(6) !< Piola-Kirchhoff stress
    real(kind=kreal), intent(inout)  :: extval(:) !< plastic strain
    real(kind=kreal), intent(out)    :: plstrain  !< plastic strain increment
    real(kind=kreal), intent(in)     :: ttime     !< total time at the start of the current increment
    real(kind=kreal), intent(in)     :: dtime     !< time length of the increment
    real(kind=kreal), optional       :: temp      !> temprature

    integer :: i
    logical :: ierr
    real(kind=kreal) :: ina(1), outa(3)
    real(kind=kreal) :: xxn, aa

    real(kind=kreal) :: c3,e,un,G,dg,ddg,stri(6),p,dstri,c4,c5,f,df, eqvs

    aa = 0.0d0
    xxn = 0.0d0

    if( dtime==0.d0 ) return
    !
    !     elastic constants
    !
    if( present(temp) ) then
      ina(1) = temp
      call fetch_TableData( MC_ISOELASTIC, matl%dict, outa(1:2), ierr, ina )
    else
      call fetch_TableData( MC_ISOELASTIC, matl%dict, outa(1:2), ierr )
    endif
    if( ierr ) then
      stop "error in isotropic elasticity definition"
    else
      e=outa(1)
      un=outa(2)
    endif

    !      Norton
    if( matl%mtype==NORTON ) then         ! those with no yield surface
      if( present( temp ) ) then
        ina(1) = temp
        call fetch_TableData( MC_NORTON, matl%dict, outa, ierr, ina )
      else
        call fetch_TableData( MC_NORTON, matl%dict, outa, ierr )
      endif
      if( ierr ) then
        stop "error in isotropic elasticity definition"
      else
        xxn=outa(2)
        aa=outa(1)*((ttime+dtime)**(outa(3)+1.d0)-ttime**(outa(3)+1.d0))/(outa(3)+1.d0)
      endif
    endif

    G=0.5d0*e/(1.d0+un)

    !
    !     creep
    !
    stri(:)=stress(:)
    p=-(stri(1)+stri(2)+stri(3))/3.d0
    do i=1,3
      stri(i)=stri(i)+p
    enddo
    !
    dstri=dsqrt(1.5d0*(stri(1)*stri(1)+stri(2)*stri(2)+stri(3)*stri(3)+    &
      2.d0*(stri(4)*stri(4)+stri(5)*stri(5)+stri(6)*stri(6))) )
    !
    !     determination of the consistency parameter
    !
    dg=0.d0
    do
      if( matl%mtype==NORTON ) then
        eqvs = dstri-3.d0*G*dg
        f=aa*eqvs**xxn
        df=xxn*f/eqvs
        ddg = (f-dg)/(3.d0*G*df+1.d0)
        dg = dg+ddg
        if((ddg<dg*1.d-6).or.(ddg<1.d-12)) exit
      endif
    enddo

    stri(:) = stri(:)-3.d0*G*dg*stri(:)/dstri
    stress(1:3) = stri(1:3)-p
    stress(4:6) = stri(4:6)

    !
    !     state variables
    !
    plstrain= dg
    extval(1)=eqvs

  end subroutine

  !> Update viscoplastic state
  subroutine updateViscoState( gauss )
    use mMechGauss
    type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point

    gauss%fstatus(2) = gauss%fstatus(2)+gauss%plstrain
  end subroutine

end module
