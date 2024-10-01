!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module summarizes all information of material properties
module mMaterial
  use hecmw_util
  use m_table
  use Table_DICTS
  implicit none

  ! Following algorithm type
  integer(kind=kint), parameter :: INFINITESIMAL = 0
  integer(kind=kint), parameter :: TOTALLAG  = 1
  integer(kind=kint), parameter :: UPDATELAG = 2

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
  !     6: Incomp newtonian
  !     7: Connector(Spring, Dashpot, Joint etc.)
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
  !     For spring or dashpot, joint etc.
  !       0: spring dof                     ie. 170000
  !       1: spring axial                       171000
  !       2: dashpot dof                        172000
  !       3: dashpot axial                      173000
  !   Fourth digit
  !     For spring_d or dashpot_d, joint_d etc.
  !       k: number of dof param (1st digit)
  !   Fifth digit:
  !     For elastoplastic deformation, hardening law
  !       0: Linear hardening           i.e.  120000
  !       1: Multilinear hardening            120010
  !       2: Swift
  !       3: Ramberg-Osgood
  !       4: linear kinematic
  !       5: combined (linear kinematic + linear isotropic)
  !     For spring_d or dashpot_d, joint_d etc.
  !       k: number of dof param (2nd digit)
  !   Six digit:
  !     For visco-elastoplastic deformation, visco law
  !       0: Norton                     i.e.  150000
  !       1: Striab                           150001
  integer(kind=kint), parameter :: USERMATERIAL          = 100000

  integer(kind=kint), parameter :: ELASTIC               = 110000
  integer(kind=kint), parameter :: MN_ORTHOELASTIC       = 111000
  integer(kind=kint), parameter :: USERELASTIC           = 112000

  integer(kind=kint), parameter :: EPLASTIC              = 120000

  integer(kind=kint), parameter :: NEOHOOKE              = 130000
  integer(kind=kint), parameter :: MOONEYRIVLIN          = 131000
  integer(kind=kint), parameter :: ARRUDABOYCE           = 132000
  integer(kind=kint), parameter :: USERHYPERELASTIC      = 133000
  integer(kind=kint), parameter :: MOONEYRIVLIN_ANISO    = 134000

  integer(kind=kint), parameter :: VISCOELASTIC          = 140000
  integer(kind=kint), parameter :: NORTON                = 150000

  integer(kind=kint), parameter :: INCOMP_NEWTONIAN      = 160000
  integer(kind=kint), parameter :: CONNECTOR             = 170000

  ! Following section type
  integer(kind=kint), parameter :: D3            = -1
  integer(kind=kint), parameter :: PlaneStress   = 1
  integer(kind=kint), parameter :: PlaneStrain   = 0
  integer(kind=kint), parameter :: AxisSymetric  = 2
  integer(kind=kint), parameter :: Shell         = 3

  ! Material constants are saved in an array of size 100 and their physical meaning
  ! correspond to their position in the array
  integer(kind=kint), parameter :: M_YOUNGS  = 1
  integer(kind=kint), parameter :: M_POISSON = 2
  integer(kind=kint), parameter :: M_DENSITY = 3
  integer(kind=kint), parameter :: M_THICK   = 4

  ! following plastic constitutive parameter
  integer(kind=kint), parameter :: M_PLCONST1 = 5
  integer(kind=kint), parameter :: M_PLCONST2 = 6
  integer(kind=kint), parameter :: M_PLCONST3 = 7
  integer(kind=kint), parameter :: M_PLCONST4 = 8
  integer(kind=kint), parameter :: M_PLCONST5 = 9
  integer(kind=kint), parameter :: M_KINEHARD = 10

  integer(kind=kint), parameter :: M_EXAPNSION = 20

  integer(kind=kint), parameter :: M_ALPHA_OVER_MU = 21

  integer(kind=kint), parameter :: M_BEAM_RADIUS = 22
  integer(kind=kint), parameter :: M_BEAM_ANGLE1 = 23
  integer(kind=kint), parameter :: M_BEAM_ANGLE2 = 24
  integer(kind=kint), parameter :: M_BEAM_ANGLE3 = 25
  integer(kind=kint), parameter :: M_BEAM_ANGLE4 = 26
  integer(kind=kint), parameter :: M_BEAM_ANGLE5 = 27
  integer(kind=kint), parameter :: M_BEAM_ANGLE6 = 28

  integer(kind=kint), parameter :: M_VISCOCITY = 29

  ! additional plastic constitutive parameter
  integer(kind=kint), parameter :: M_PLCONST6 = 30
  integer(kind=kint), parameter :: M_PLCONST7 = 31
  integer(kind=kint), parameter :: M_PLCONST8 = 32
  integer(kind=kint), parameter :: M_PLCONST9 = 33
  integer(kind=kint), parameter :: M_PLCONST10 = 34

  integer(kind=kint), parameter :: M_DAMPING_RM = 35
  integer(kind=kint), parameter :: M_DAMPING_RK = 36

  integer(kind=kint), parameter :: M_SPRING_DOF    = 0
  integer(kind=kint), parameter :: M_SPRING_AXIAL  = 1
  integer(kind=kint), parameter :: M_DASHPOT_DOF   = 2
  integer(kind=kint), parameter :: M_DASHPOT_AXIAL = 3

  integer(kind=kint), parameter :: M_SPRING_D_NDOFFSET  = 0
  integer(kind=kint), parameter :: M_SPRING_A_NDOFFSET  = 72
  integer(kind=kint), parameter :: M_DASHPOT_D_NDOFFSET = 73
  integer(kind=kint), parameter :: M_DASHPOT_A_NDOFFSET = 145

  ! Dictionary constants
  character(len=DICT_KEY_LENGTH) :: MC_ISOELASTIC= 'ISOELASTIC'      ! youngs modulus, poisson's ratio
  character(len=DICT_KEY_LENGTH) :: MC_ORTHOELASTIC= 'ORTHOELASTIC'  ! ortho elastic modulus
  character(len=DICT_KEY_LENGTH) :: MC_YIELD = 'YIELD'               ! plastic strain, yield stress
  character(len=DICT_KEY_LENGTH) :: MC_THEMOEXP = 'THEMOEXP'         ! thermo expansion coefficient
  character(len=DICT_KEY_LENGTH) :: MC_ORTHOEXP = 'ORTHOEXP'         ! thermo expansion coefficient
  character(len=DICT_KEY_LENGTH) :: MC_VISCOELASTIC = 'VISCOELASTIC' ! Prony coeff only curr.
  character(len=DICT_KEY_LENGTH) :: MC_NORTON = 'NORTON'             ! NOrton's creep law
  character(len=DICT_KEY_LENGTH) :: MC_INCOMP_NEWTONIAN = 'INCOMP_FLUID' ! viscocity
  character(len=DICT_KEY_LENGTH) :: MC_SPRING= 'SPRING'              ! spring
  character(len=DICT_KEY_LENGTH) :: MC_DASHPOT= 'DASHPOT'            ! dashpot

  type tshellmat
    integer(kind=kint)         :: ortho
    real(kind=kreal)           :: ee
    real(kind=kreal)           :: pp
    real(kind=kreal)           :: ee2
    real(kind=kreal)           :: g12
    real(kind=kreal)           :: g23
    real(kind=kreal)           :: g31
    real(kind=kreal)           :: angle
    real(kind=kreal)           :: rho
    real(kind=kreal)           :: alpha
    real(kind=kreal)           :: alpha_over_mu
    real(kind=kreal)           :: weight
  end type tshellmat

  !> Structure to manage all material related data
  type tMaterial
    integer(kind=kint)         :: nlgeom_flag       !< type of constitutive relation
    integer(kind=kint)         :: mtype             !< material type
    integer(kind=kint)         :: nfstatus          !< number of status variables
    character(len=30)          :: name              !< material name
    real(kind=kreal)           :: variables(200)    !< material properties
    integer(kind=kint)         :: variables_i(200)  !< material properties(integer)
    type(tshellmat), pointer   :: shell_var(:)      !< material properties for shell
    integer(kind=kint)         :: totallyr          !< total layer of element
    integer(kind=kint)         :: cdsys_ID          !< ID of material coordinate system
    integer(kind=kint)         :: n_table           !< size of table
    real(kind=kreal), pointer  :: table(:)=>null()  !< material properties in tables
    type(DICT_STRUCT), pointer :: dict              !< material properties in dictionaried linked list
  end type tMaterial

  type(tMaterial), allocatable :: materials(:)

contains

  !> Initializer
  subroutine initMaterial( material )
    type( tMaterial ), intent(inout) :: material
    material%mtype = -1                  ! not defined yet
    material%nfstatus = 0                ! Default: no status
    material%nlgeom_flag = INFINITESIMAL ! Default: INFINITESIMAL ANALYSIS
    material%variables =  0.d0           ! not defined yet
    material%variables_i =  0            ! not defined yet
    material%totallyr =  0               ! not defined yet

    call dict_create( material%dict, 'INIT', DICT_NULL )
  end subroutine

  !> Finalizer
  subroutine finalizeMaterial( material )
    type( tMaterial ), intent(inout) :: material
    if( associated(material%table) ) deallocate( material%table )
    if( associated(material%dict) ) call dict_destroy( material%dict )
  end subroutine finalizeMaterial

  !> Initializer
  subroutine initializeMatls( nm )
    integer, intent(in) :: nm
    integer :: i
    if( allocated(materials) ) deallocate( materials )
    allocate( materials( nm ) )
    do i=1,nm
      call initMaterial( materials(i) )
    enddo
  end subroutine

  !> Finalizer
  subroutine finalizeMatls()
    integer :: i
    if( allocated( materials ) ) then
      do i=1,size(materials)
        call finalizeMaterial( materials(i) )
      enddo
      deallocate( materials )
    endif
  end subroutine

  !> Set value of variable(m) of material n to v
  subroutine modifyMatl( n,m,v)
    integer, intent(in)  :: n
    integer, intent(in)  :: m
    real(kind=kreal), intent(in) :: v

    if( n>size(materials) .OR. m>100 ) return
    materials(n)%variables(m) = v
  end subroutine

  !> Print out the material properties
  subroutine printMaterial( nfile, material )
    integer, intent(in)           :: nfile
    type( tMaterial ), intent(in) :: material
    integer :: i, nt
    write( nfile, *) "Material type:",material%mtype,material%nlgeom_flag
    do i=1,100
      if( material%variables(i) /= 0.d0 ) write( nfile, *) i,material%variables(i)
    enddo
    if( associated( material%table ) ) then
      nt = size(material%table)
      write( nfile,* ) "--table--"
      do i=1,nt
        write(nfile,*) i,material%table(i)
      enddo
    endif
    call print_TableData( material%dict, nfile )
  end subroutine

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
    integer, intent(inout)           :: mtype
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
    integer, intent(in) :: mtype
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
    integer, intent(in) :: mtype
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
    integer, intent(in) :: mtype
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
    integer, intent(in) :: mtype
    integer :: itype
    isKinematicHarden = .false.
    itype = fetchDigit( 5, mtype )
    if( itype==4 .or. itype==5 ) isKinematicHarden = .true.
  end function

  !> If it is an elastic material?
  logical function isElastic( mtype )
    integer, intent(in) :: mtype
    integer :: itype
    isElastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==1 ) isElastic = .true.
  end function

  !> If it is an elastoplastic material?
  logical function isElastoplastic( mtype )
    integer, intent(in) :: mtype
    integer :: itype
    isElastoplastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==2 ) isElastoplastic = .true.
  end function

  !> If it is a hyperelastic material?
  logical function isHyperelastic( mtype )
    integer, intent(in) :: mtype
    integer :: itype
    isHyperelastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==3 ) isHyperelastic = .true.
  end function

  !> If it is an viscoelastic material?
  logical function isViscoelastic( mtype )
    integer, intent(in) :: mtype
    integer :: itype
    isViscoelastic = .false.
    itype = fetchDigit( 2, mtype )
    if( itype==4 ) isViscoelastic = .true.
  end function

  !> Set material type of elastoplastic to elastic
  subroutine ep2e( mtype )
    integer, intent(inout) :: mtype
    if( .not. isElastoplastic( mtype ) ) return
    call setDigit( 2, 1, mtype )
  end subroutine

  !> Get type of connector
  integer function getConnectorType( mtype )
    integer, intent(in) :: mtype
    integer :: itype
    getConnectorType = -1
    itype = fetchDigit( 1, mtype )
    if( itype/=1 ) return  ! not defomration problem
    itype = fetchDigit( 2, mtype )
    if( itype/=7 ) return  ! not connector
    getConnectorType = fetchDigit( 3, mtype )
  end function

  !> Get number of spring_d parameters 
  integer function getNumOfSpring_dParam( material )
    type( tMaterial ), intent(in) :: material
    getNumOfSpring_dParam = material%variables_i(M_SPRING_D_NDOFFSET+1)
  end function

  !> Get number of spring_a parameters 
  integer function getNumOfSpring_aParam( material )
    type( tMaterial ), intent(in) :: material
    getNumOfSpring_aParam = material%variables_i(M_SPRING_A_NDOFFSET+1)
  end function

  !> Get number of dashpot_d parameters 
  integer function getNumOfDashpot_dParam( material )
    type( tMaterial ), intent(in) :: material
    getNumOfDashpot_dParam = material%variables_i(M_DASHPOT_D_NDOFFSET+1)
  end function

  !> Get number of dashpot_a parameters 
  integer function getNumOfDashpot_aParam( material )
    type( tMaterial ), intent(in) :: material
    getNumOfDashpot_aParam = material%variables_i(M_DASHPOT_A_NDOFFSET+1)
  end function

end module



