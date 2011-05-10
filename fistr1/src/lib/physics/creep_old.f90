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

      real(kind=kreal) :: c1,c2,c3,ep0(6),eei(6),e,un,al,am1,um2,ep(6),dg,    &
       ddg,stri(6),p,eeq0,eeq,um,dstri,c4,c5,f,df, eqvs, stiff(21)

      c1 =dsqrt(2.d0/3.d0)
      c2 = 2.d0/3.d0
	  
!
!     elastic
!
      if( all(stress==0.d0) ) then
        call calElasticMatrix( matl, sectTYPE, stiffness )
        return
      endif
!
!     state variables
!
      eeq0=extval(1)
      ep0(1:6)=extval(2:7)
!
!     elastic strains
!
      eei(:)=strain(:)-ep0(:)
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
      al=un*um2/(1.d0-2.d0*un)
      am1=al+um2
      um=um2/2.d0
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

      dg=eeq0/c1
      eqvs = (dstri-um2*dg)/c1
      if( eqvs<1.d-10 ) eqvs=1.d-10
      f=aa*eqvs**xxn
      df=xxn*f/eqvs
!
!     stiffness matrix
!
         c3=um2*um2
         c4=c3*dg/dstri
         c3=c4-c3*df/(um2*df+c2)
         c5=c4/3.d0
!print *, dg		 
!
         stiff(1)=am1+c3*stri(1)*stri(1)+c5-c4
         stiff(2)=al+c3*stri(1)*stri(2)+c5
         stiff(3)=am1+c3*stri(2)*stri(2)+c5-c4
         stiff(4)=al+c3*stri(1)*stri(3)+c5
         stiff(5)=al+c3*stri(2)*stri(3)+c5
         stiff(6)=am1+c3*stri(3)*stri(3)+c5-c4
         stiff(7)=0.d0+c3*stri(1)*stri(4)
         stiff(8)=0.d0+c3*stri(2)*stri(4)
         stiff(9)=0.d0+c3*stri(3)*stri(4)
         stiff(10)=um+c3*stri(4)*stri(4)-c4/2.d0
         stiff(11)=0.d0+c3*stri(1)*stri(5)
         stiff(12)=0.d0+c3*stri(2)*stri(5)
         stiff(13)=0.d0+c3*stri(3)*stri(5)
         stiff(14)=0.d0+c3*stri(4)*stri(5)
         stiff(15)=um+c3*stri(5)*stri(5)-c4/2.d0
         stiff(16)=0.d0+c3*stri(1)*stri(6)
         stiff(17)=0.d0+c3*stri(2)*stri(6)
         stiff(18)=0.d0+c3*stri(3)*stri(6)
         stiff(19)=0.d0+c3*stri(4)*stri(6)
         stiff(20)=0.d0+c3*stri(5)*stri(6)
         stiff(21)=um+c3*stri(6)*stri(6)-c4/2.d0
		 
      do i=1,3
        stiffness(i,i) = am1+c3*stri(i)*stri(i)+c5-c4
        do j=i+1,3
          stiffness(i,j) = al+c3*stri(i)*stri(j)+c5
        enddo
        do j=4,6
          stiffness(i,j) = c3*stri(i)*stri(j)
        enddo
      enddo
      do i=4,6
        stiffness(i,i) = um+c3*stri(i)*stri(i)-c4/2.d0
        do j=i+1,6
          stiffness(i,j) = c3*stri(i)*stri(j)
        enddo
      enddo
      do i=1,6
      do j=i+1,6
        stiffness(j,i)=stiffness(i,j)
      enddo
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

      real(kind=kreal) :: c1,c2,c3,ep0(6),eei(6),e,un,al,am1,um2,ep(6),dg,    &
       ddg,stri(6),p,eeq0,eeq,um,dstri,c4,c5,f,df, eqvs

!
!     state variables
!
      c1 = dsqrt(2.d0/3.d0)
	  c2 = 2.d0/3.d0
      eeq0=extval(1)
      ep0(1:6)=extval(2:7)
!
!     elastic strains
!
      eei(:)=strain(:)-ep0(:)
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
      al=un*um2/(1.d0-2.d0*un)
      am1=al+um2
      um=um2/2.d0

!
!     creep
!
      stri(1)=am1*eei(1)+al*(eei(2)+eei(3))
      stri(2)=am1*eei(2)+al*(eei(1)+eei(3))
      stri(3)=am1*eei(3)+al*(eei(1)+eei(2))
      stri(4)=um2*eei(4)
      stri(5)=um2*eei(5)
      stri(6)=um2*eei(6)
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
      do i=1,6 
         stri(i)=stri(i)/dstri
      enddo
!
      dg=0.d0
	  plstrain = 0.d0
!
!     determination of the consistency parameter
!
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
      eeq=eeq0+c1*dg
      plstrain= plstrain+c1*dg

      do i=1,6
         ep(i)=dg*stri(i)
         eei(i)=eei(i)-ep(i)
         extval(8+i)=ep0(i)+ep(i)
      enddo
!
!     stress values
!
      stress(1)=am1*eei(1)+al*(eei(2)+eei(3))
      stress(2)=am1*eei(2)+al*(eei(1)+eei(3))
      stress(3)=am1*eei(3)+al*(eei(1)+eei(2))
      stress(4)=um2*eei(4)
      stress(5)=um2*eei(5)
      stress(6)=um2*eei(6)
!
!     state variables
!
      extval(8)=eeq
	  
   end subroutine
   
   !> Update viscoplastic state 
   subroutine updateViscoState( gauss )
      use mMechGauss
      type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point
	  
      integer :: i
      do i=1,7
        gauss%fstatus(i) = gauss%fstatus(i+7)
      enddo
      gauss%plstrain =0.d0
   end subroutine
   
end module
