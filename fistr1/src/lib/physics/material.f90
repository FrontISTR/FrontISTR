!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!     Last Update : 2010/08/28                                         !
!        Category : Material                                           !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module summarizes all infomation of material properties
MODULE mMaterial
  USE m_table
  USE Table_DICTS
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: kreal = kind(0.0d0)
  
  
  ! Following algorithm type
      INTEGER, PARAMETER :: INFINITE = 0
      INTEGER, PARAMETER :: TOTALLAG  = 1
      INTEGER, PARAMETER :: UPDATELAG = 2

  ! Following material types. All material type number consists with integer of six digits.
  !   First digit: Indicates physical type
  !                1: mechanical deformation analysis
  !                2: heat conduct analysis
  !                ......
  !   Second digit: 
  !     Mechanical analysis
  !     1: Elastic    
  !     2: Elastoplastic    
  !     3: Hyperelastic
  !     4: Viscoelastic
  !     5: Viscoplastic
  !     Heat conductiovity
  !     ......
  !   Third digit:
  !     For elastic or elastoplastic deformation, elastic 
  !       0: isotropic                      ie. 110000
  !       1: with transever anisotropity        111000
  !     For hyperelastic deformation
  !       0: Neo-Hooke                          130000
  !       1: Mooney-Rivlin                      131000
  !       2: Arruda-Boyce                       132000
  !   Fourth digit
  !     For elastoplastic problem,yield function 
  !       0: isotropic (Mises)
  !       1: Mohr-Coulomb
  !       2: Drucker-Prager
  !   Fifth digit:
  !     For elastoplastic deformation, hardening law
  !       0: Linear hardening           i.e.  120000
  !       1: Multilinear hardening            120010
  !       2: Swift                            
  !       3: Ramberg-Osgood
  !       4: linear kinematic
  !       5: combined (linear kinematic + linear isotropic)
  !   Six digit:
  !     For visco-elastoplastic deformation, visco law
  !       0: Norton                     i.e.  150000
  !       1: Striab                           150001
      INTEGER, PARAMETER :: USERMATERIAL          = 100000
  
      INTEGER, PARAMETER :: ELASTIC               = 110000
      INTEGER, PARAMETER :: MN_ORTHOELASTIC       = 111000
      INTEGER, PARAMETER :: USERELASTIC           = 112000
	  
      INTEGER, PARAMETER :: EPLASTIC              = 120000
	  
      INTEGER, PARAMETER :: NEOHOOKE              = 130000
      INTEGER, PARAMETER :: MOONEYRIVLIN          = 131000
      INTEGER, PARAMETER :: ARRUDABOYCE           = 132000
      INTEGER, PARAMETER :: USERHYPERELASTIC      = 133000

      INTEGER, PARAMETER :: VISCOELASTIC          = 140000
      INTEGER, PARAMETER :: NORTON                = 150000
	  
  ! Following section type
      INTEGER, PARAMETER :: D3            = -1
      INTEGER, PARAMETER :: PlaneStress   = 1
      INTEGER, PARAMETER :: PlaneStrain   = 0
      INTEGER, PARAMETER :: AxisSymetric  = 2
      INTEGER, PARAMETER :: Shell         = 3

  ! Material constants are saved in an array of size 100 and their physical meaning
  ! are conrresponds to their position in the array
      INTEGER, PARAMETER :: M_YOUNGS  = 1
      INTEGER, PARAMETER :: M_POISSON = 2
      INTEGER, PARAMETER :: M_DENSITY = 3
      INTEGER, PARAMETER :: M_THICK   = 4

      ! following plastic constitutive parameter
      INTEGER, PARAMETER :: M_PLCONST1 = 5
      INTEGER, PARAMETER :: M_PLCONST2 = 6
      INTEGER, PARAMETER :: M_PLCONST3 = 7
      INTEGER, PARAMETER :: M_PLCONST4 = 8
      INTEGER, PARAMETER :: M_PLCONST5 = 9
      INTEGER, PARAMETER :: M_KINEHARD = 10

      INTEGER, PARAMETER :: M_EXAPNSION = 20
      
      INTEGER, PARAMETER :: M_ALPHA_OVER_MU = 21

      INTEGER, PARAMETER :: M_BEAM_RADIUS = 22
      INTEGER, PARAMETER :: M_BEAM_ANGLE1 = 23
      INTEGER, PARAMETER :: M_BEAM_ANGLE2 = 24
      INTEGER, PARAMETER :: M_BEAM_ANGLE3 = 25
      INTEGER, PARAMETER :: M_BEAM_ANGLE4 = 26
      INTEGER, PARAMETER :: M_BEAM_ANGLE5 = 27
      INTEGER, PARAMETER :: M_BEAM_ANGLE6 = 28
      
      ! 80~100, hardening parameter
	  
	!<********** Laminated Shell material properties  **********
!		INTEGER, PARAMETER :: SHELL material properties = 101~220 

	!********** Use for i_variables **********
		INTEGER, PARAMETER :: M_SHELL_SOLID = 298
		INTEGER, PARAMETER :: M_TOTAL_LAYER = 299
		INTEGER, PARAMETER :: M_SHELL_MATLTYPE = 300
		
   ! Dictionary constants
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_ISOELASTIC= 'ISOELASTIC'      ! youngs modulus, poisson's ratio
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_ORTHOELASTIC= 'ORTHOELASTIC'  ! ortho elastic modulus
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_YIELD = 'YIELD'               ! plastic strain, yield stress
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_THEMOEXP = 'THEMOEXP'         ! thermo expansion coefficient
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_ORTHOEXP = 'ORTHOEXP'         ! thermo expansion coefficient
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_VISCOELASTIC = 'VISCOELASTIC' ! Prony coeff only curr.
      CHARACTER(len=DICT_KEY_LENGTH) :: MC_NORTON = 'NORTON'             ! NOrton's creep law

  !> Stucture to management all material relates data
  TYPE tMaterial
    integer                    :: nlgeom_flag       !< type of constitutive relation
    integer                    :: mtype             !< material type
    integer                    :: nfstatus          !< number of status variables
    character(len=30)          :: name              !< material name
    real(kind=kreal)           :: variables(300)    !< material properties
    integer                    :: cdsys_ID          !< ID of material coordinate system
    integer                    :: n_table           !< size of table
    REAL(kind=kreal), pointer  :: table(:)=>null()  !< material properties in tables
    type(DICT_STRUCT), pointer :: dict              !< material properties in dictionaried linked list
  END TYPE

  TYPE(tMaterial), ALLOCATABLE :: materials(:)

  CONTAINS

!> Initializer
  SUBROUTINE initMaterial( material )
    TYPE( tMaterial ), INTENT(INOUT) :: material
    material%mtype = -1                  ! not defined yet
    material%nfstatus = 0                ! Default: no status 
    material%nlgeom_flag = INFINITE     ! Default: INFINITE ANALYSIS
    material%variables = 0.d0            ! not defined yet

	call dict_create( material%dict, 'INIT', DICT_NULL )
  END SUBROUTINE

!> Initializer
  SUBROUTINE initializeMatls( nm )
    INTEGER, INTENT(IN) :: nm
    INTEGER :: i
    IF( ALLOCATED(materials) ) DEALLOCATE( materials )
    ALLOCATE( materials( nm ) )
    DO i=1,nm
      CALL initMaterial( materials(i) )
    ENDDO
  END SUBROUTINE

!> Finalizer
  SUBROUTINE finalizeMatls()
    INTEGER :: i
    IF( ALLOCATED( materials ) ) then
      DO i=1,size(materials)
        IF( ASSOCIATED(materials(i)%table) ) DEALLOCATE( materials(i)%table )
      ENDDO
      DEALLOCATE( materials )
    ENDIF
  END SUBROUTINE

!> Set value of variable(m) of material n to v
  SUBROUTINE modifyMatl( n,m,v)
    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(IN)  :: m
    REAL(KIND=kreal), INTENT(IN) :: v
    
    IF( n>SIZE(materials) .OR. m>100 ) RETURN
    materials(n)%variables(m) = v
  END SUBROUTINE

!> Print out the material properties
  SUBROUTINE printMaterial( nfile, material )
    INTEGER, INTENT(IN)           :: nfile
    TYPE( tMaterial ), INTENT(IN) :: material
    INTEGER :: i, nt
    WRITE( nfile, *) "Material type:",material%mtype,material%nlgeom_flag
    DO i=1,100
      if( material%variables(i) /= 0.d0 ) WRITE( nfile, *) i,material%variables(i)
    ENDDO
    IF( ASSOCIATED( material%table ) ) THEN
      nt = SIZE(material%table)
      WRITE( nfile,* ) "--table--"
      DO i=1,nt
        WRITE(nfile,*) i,material%table(i)
      ENDDO
    ENDIF
	call print_TableData( material%dict, nfile )
  END SUBROUTINE

!> Fetch material type
  integer function fetchDigit( npos, cnum )
    integer, intent(in) :: npos
    integer, intent(in) :: cnum
    integer :: i, idum,cdum,dd
    fetchDigit = -1
    cdum = cnum
    if( npos<=0 .or. npos>6) return
    if( cnum<100000 .or. cnum>999999 ) return
    dd = 100000
    do i=1,npos-1
      idum = cdum/dd
      cdum = cdum-idum*dd
      dd = dd/10
    enddo
    fetchDigit = cdum/10**(6-npos)
  end function

!> Modify material type
  subroutine setDigit( npos, ival, mtype )
    integer, intent(in)              :: npos
    integer, intent(in)              :: ival
    integer, INTENT(INOUT)           :: mtype
    integer :: i, idum,cdum, cdum1, dd
    cdum = mtype
    if( npos<=0 .or. npos>6 ) return
    if( ival<0 .or. ival>9 ) return
    dd =100000
    cdum1 = 0
    do i=1,npos-1
      idum = cdum/dd
      cdum1 = cdum1+ idum*dd
      cdum = cdum-idum*dd
      dd=dd/10
    enddo
    cdum1 = cdum1 + ival*dd
    idum = cdum/dd
    cdum = cdum-idum*dd
    dd=dd/10
    do i=npos+1,6
      idum = cdum/dd
      cdum1 = cdum1+ idum*dd
      cdum = cdum-idum*dd
      dd=dd/10
    enddo
    mtype = cdum1
  end subroutine

!> Get elastic type
  integer function getElasticType( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    getElasticType = -1
    itype = fetchDigit( 1, mtype )
    if( itype/=1 ) return  ! not defomration problem
    itype = fetchDigit( 2, mtype )
    if( itype/=1 .and. itype/=2 ) return  ! not defomration problem
    getElasticType = fetchDigit( 3, mtype )
  end function

!> Get type of yield function
  integer function getYieldFunction( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    getYieldFunction = -1
    itype = fetchDigit( 1, mtype )
    if( itype/=1 ) return  ! not defomration problem
    itype = fetchDigit( 2, mtype )
    if( itype/=2 ) return  ! not elstoplastic problem
    getYieldFunction = fetchDigit( 4, mtype )
  end function
  
!> Get type of hardening
  integer function getHardenType( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    getHardenType = -1
    itype = fetchDigit( 1, mtype )
    if( itype/=1 ) return  ! not defomration problem
    itype = fetchDigit( 2, mtype )
    if( itype/=2 ) return  ! not elstoplastic problem
    getHardenType = fetchDigit( 5, mtype )
  end function

!> If it is a kinematic hardening material?
  logical function isKinematicHarden( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    isKinematicHarden = .false.
    itype = fetchDigit( 5, mtype )
    if( itype==4 .or. itype==5 ) isKinematicHarden = .true.
  end function

!> If it is an elastic material?
  logical function isElastic( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    isElastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==1 ) isElastic = .true.
  end function

!> If it is an elastoplastic material?
  logical function isElastoplastic( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    isElastoplastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==2 ) isElastoplastic = .true.
  end function
  
!> If it is an viscoelastic material?
  logical function isViscoelastic( mtype )
    INTEGER, INTENT(IN) :: mtype
    integer :: itype
    isViscoelastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==4 ) isViscoelastic = .true.
  end function

!> Set material type of elastoplastic to elastic
  subroutine ep2e( mtype )
    INTEGER, INTENT(INOUT) :: mtype
    if( .not. isElastoplastic( mtype ) ) return
    call setDigit( 2, 1, mtype )
  end subroutine

END MODULE


  
