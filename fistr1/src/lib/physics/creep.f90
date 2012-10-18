!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions for creep calculation
!
!>  \author     X.Yuan(Advancesoft)
!>  \date       2010/10/06
!>  \version    0.00
!
!======================================================================!
module mCreep

  use mMaterial
  use m_ElasticLinear
  
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
  
    !> This subrooutine calculates stiffness for elastically isotropic
    !>     materials with isotropic creep
    subroutine iso_creep(matl, sectType, stress, strain, extval,plstrain,            &
             dtime,ttime,stiffness, temp)
      TYPE( tMaterial ), INTENT(IN)    :: matl      !< material properties
      INTEGER, INTENT(IN)              :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)     :: stress(6) !< Piola-Kirchhoff stress
      REAL(KIND=kreal), INTENT(IN)     :: strain(6) !< strain
      REAL(KIND=kreal), INTENT(IN)     :: extval(:) !< plastic strain
      REAL(KIND=kreal), INTENT(IN)     :: plstrain  !< plastic strain increment
      REAL(KIND=kreal), INTENT(IN)     :: ttime     !< total time at the start of the current increment
      REAL(KIND=kreal), INTENT(IN)     :: dtime     !< time length of the increment
      REAL(KIND=kreal), INTENT(out)    :: stiffness(6,6) !< stiffness
      REAL(KIND=kreal), OPTIONAL       :: temp      !> temprature

      integer :: i, j
      logical :: ierr
      real(kind=kreal) :: ina(1), outa(3)
      real(kind=kreal) :: xxa, xxn, aa

      real(kind=kreal) :: c1,c2,c3,e,un,um2,dg,    &
       ddg,stri(6),p,eeq,dstri,c4,c5,f,df, eqvs

      c1 =dsqrt(2.d0/3.d0)
      c2 = 2.d0/3.d0
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
        xxa=outa(1)*(ttime+dtime)**outa(3)
        if(xxa<1.d-20) xxa=1.d-20
        xxn=outa(2)
        aa=xxa*dtime
      endif

      um2=e/(1.d0+un)
!
!     creep
!
      stri(:)=stress(:)
      p=-(stri(1)+stri(2)+stri(3))/3.d0
      do i=1,3
         stri(i)=stri(i)+p
      enddo
!
      dstri=dsqrt(stri(1)*stri(1)+stri(2)*stri(2)+stri(3)*stri(3)+    &
           2.d0*(stri(4)*stri(4)+stri(5)*stri(5)+stri(6)*stri(6)))
!
!     unit trial vector
!
      stri(:)=stri(:)/dstri

      eqvs = extval(1)
      if( eqvs<1.d-10 ) eqvs=1.d-10
      f=aa*eqvs**xxn
      df=xxn*f/eqvs
      dg=(f-plstrain)/c1
!
!     stiffness matrix
!
      c3=um2*um2
      c4=c3*dg/dstri
      c3=c4-c3*df/(um2*df+c2)
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
   
   !> This subrooutine calculates stresses and creep status for an elastically isotropic
   !>     material with isotropic creep
   subroutine update_iso_creep(matl, sectType, strain, stress, extval,plstrain,                &
             dtime,ttime,temp)
      TYPE( tMaterial ), INTENT(IN)    :: matl      !< material properties
      INTEGER, INTENT(IN)              :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)     :: strain(6) !< strain
      REAL(KIND=kreal), INTENT(INOUT)  :: stress(6) !< Piola-Kirchhoff stress
      REAL(KIND=kreal), INTENT(INOUT)  :: extval(:) !< plastic strain
      REAL(KIND=kreal), INTENT(OUT)    :: plstrain  !< plastic strain increment
      REAL(KIND=kreal), INTENT(IN)     :: ttime     !< total time at the start of the current increment
      REAL(KIND=kreal), INTENT(IN)     :: dtime     !< time length of the increment
      REAL(KIND=kreal), OPTIONAL       :: temp      !> temprature

      integer :: i
      logical :: ierr
      real(kind=kreal) :: ina(1), outa(3)
      real(kind=kreal) :: xxa, xxn, aa

      real(kind=kreal) :: c1,c2,c3,e,un,um2,dg,    &
       ddg,stri(6),p,dstri,c4,c5,f,df, eqvs
	   
      if( dtime==0.d0 ) return
!
!     state variables
!
      c1 = dsqrt(2.d0/3.d0)
	  c2 = 2.d0/3.d0
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
          xxa=outa(1)*(ttime+dtime)**outa(3)
          if(xxa<1.d-20) xxa=1.d-20
          xxn=outa(2)
          aa=xxa*dtime
        endif
      endif

      um2=e/(1.d0+un)

!
!     creep
!
      stri(:)=stress(:)
      p=-(stri(1)+stri(2)+stri(3))/3.d0
      do i=1,3
         stri(i)=stri(i)+p
      enddo
!
      dstri=dsqrt(stri(1)*stri(1)+stri(2)*stri(2)+stri(3)*stri(3)+    &
           2.d0*(stri(4)*stri(4)+stri(5)*stri(5)+stri(6)*stri(6)))    
!
!     determination of the consistency parameter
!
      dg=0.d0
      do
        if( matl%mtype==NORTON ) then
          eqvs = (dstri-um2*dg)/c1
          if( eqvs<1.d-10 ) eqvs=1.d-10
          f=aa*eqvs**xxn
          df=xxn*f/eqvs
          ddg=(c1*f-c2*dg)/(um2*df+c2)
          dg=dg+ddg
          if((ddg<dg*1.d-4).or.(ddg<1.d-10)) exit
        endif
      enddo
	  
	  stri(:) = stri(:)-um2*dg*stri(:)/dstri
      stress(1:3) = stri(1:3)-p
      stress(4:6) = stri(4:6)
		
!
!     state variables
!
      plstrain= c1*dg
      extval(1)=eqvs
	  
   end subroutine
   
   !> Update viscoplastic state 
   subroutine updateViscoState( gauss )
      use mMechGauss
      type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point

      gauss%fstatus(2) = gauss%fstatus(2)+gauss%plstrain
   end subroutine
   
end module
