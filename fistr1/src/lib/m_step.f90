!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module manages step information
module m_step
  use hecmw
  implicit none

  include 'fstr_ctrl_util_f.inc'

  integer, parameter :: stepStatic = 1
  integer, parameter :: stepVisco  = 2
  integer, parameter :: stepFixedInc = 1
  integer, parameter :: stepAutoInc  = 2

  ! statistics of newton iteration
  integer(kind=kint),parameter :: knstMAXIT  = 1 ! maximum number of newton iteration
  integer(kind=kint),parameter :: knstSUMIT  = 2 ! total number of newton iteration
  integer(kind=kint),parameter :: knstCITER  = 3 ! number of contact iteration
  integer(kind=kint),parameter :: knstDRESN  = 4 ! reason of not to converged

  !> Step control such as active boundary condition, convergent condition etc.
  type step_info
    integer             :: solution           !< solution type; 1: static;  2:visco
    integer             :: inc_type           !< increment type; 1: fixed;  2:auto
    character( len=80 ) :: CONTROL            !< control type, such as arclength etc
    character( len=80 ) :: ConvControl        !< Judgement of convergency, such as nodal force residual
    !< disp increment, energy
    real(kind=kreal)    :: converg            !< value of convergent judgement
    real(kind=kreal)    :: converg_lag        !< value of convergent judgement (Lagrange)
    real(kind=kreal)    :: maxres             !< upper bound of NR residual

    integer :: num_substep                    !< substeps user given
    integer :: max_iter                       !< max number of iteration
    integer :: max_contiter                   !< max number of contact iteration
    integer :: amp_id                         !< id of amplitude definition
    real(kind=kreal) :: initdt                !< time increment
    real(kind=kreal) :: elapsetime            !< elapse time of this step
    real(kind=kreal) :: mindt                 !< lower bound of time increment
    real(kind=kreal) :: maxdt                 !< upper bound of time increment
    real(kind=kreal) :: starttime             !< start time of this step
    integer, pointer :: Boundary(:)=>null()   !< active group of boundary conditions of current step
    integer, pointer :: Load(:)=>null()       !< active group of external load conditions of current step
    integer, pointer :: Contact(:)=>null()    !< active group of contact conditions of current step
    integer :: timepoint_id                   !< id of timepoint
    integer :: AincParam_id                   !< id of auto increment parameter
  end type

  type tParamAutoInc
    character(HECMW_NAME_LEN)  :: name
    real(kind=kreal)     :: ainc_Rs            !< time increment decreasing ratio
    real(kind=kreal)     :: ainc_Rl            !< time increment increasing ratio
    integer( kind=kint ) :: NRbound_s(10)      !< # of NR iteration bound to decrease time increment
    integer( kind=kint ) :: NRbound_l(10)      !< # of NR iteration bound to increase time increment
    integer( kind=kint ) :: NRtimes_s          !< # of times that decreasing condition is satisfied
    integer( kind=kint ) :: NRtimes_l          !< # of times that increasing condition is satisfied
    real(kind=kreal)     :: ainc_Rc            !< time increment decreasing ratio for cutback
    integer( kind=kint ) :: CBbound            !< maximum # of successive cutback
  end type

contains

  !> Initializer
  subroutine init_stepInfo( stepinfo )
    type( step_info ), intent(out) :: stepinfo !< step info
    stepinfo%solution = stepStatic
    stepinfo%inc_type = stepFixedInc
    stepinfo%num_substep = 1
    stepinfo%max_iter = 50
    stepinfo%max_contiter = 10
    stepinfo%amp_id = -1
    stepinfo%initdt = 1.d0
    stepinfo%mindt = 1.d-4
    stepinfo%maxdt = 1.d0
    stepinfo%elapsetime = 1.d0
    stepinfo%starttime = 0.d0
    stepinfo%converg = 1.d-3
    stepinfo%converg_lag = 1.d-4
    stepinfo%maxres = 1.d+10
    stepinfo%timepoint_id = 0
    stepinfo%AincParam_id = 0
  end subroutine

  subroutine setup_stepInfo_starttime( stepinfos )
    type( step_info ), pointer, intent(inout) :: stepinfos(:) !< step info

    integer :: i

    stepinfos(1)%starttime = 0.d0
    do i=1,size(stepinfos)-1
      stepinfos(i+1)%starttime = stepinfos(i)%starttime + stepinfos(i)%elapsetime
    end do
  end subroutine

  !> Is boundary condition in this step active
  logical function isBoundaryActive( bnd, stepinfo )
    integer, intent(in)           :: bnd      !< group number of boundary condition
    type( step_info ), intent(in) :: stepinfo !< current step info
    isBoundaryActive = .false.
    if( .not. associated( stepinfo%Boundary ) ) return
    if( any( stepinfo%Boundary== bnd ) ) isBoundaryActive = .true.
  end function

  !> Is external load in this step active
  logical function isLoadActive( bnd, stepinfo )
    integer, intent(in)           :: bnd      !< group number of boundary condition
    type( step_info ), intent(in) :: stepinfo !< current step info
    isLoadActive = .false.
    if( .not. associated( stepinfo%Load ) ) return
    if( any( stepinfo%Load == bnd ) ) isLoadActive = .true.
  end function

  !> Is contact condition in this step active
  logical function isContactActive( bnd, stepinfo )
    integer, intent(in)           :: bnd      !< group number of boundary condition
    type( step_info ), intent(in) :: stepinfo !< current step info
    isContactActive = .false.
    if( .not. associated( stepinfo%Contact ) ) return
    if( any( stepinfo%Contact== bnd ) ) isContactActive = .true.
  end function

  !> Finalizer
  subroutine free_stepInfo( step )
    type(step_info), intent(inout) :: step  !< step info
    if( associated( step%Boundary ) ) deallocate( step%Boundary )
    if( associated( step%Load ) )     deallocate( step%Load )
    if( associated( step%Contact ) )  deallocate( step%Contact )
  end subroutine

  !> Print out step control
  subroutine fstr_print_steps( nfile, steps )
    integer, intent(in)         :: nfile     !< file number
    type(step_info), intent(in) :: steps(:)  !< step info
    integer :: i, j, nstep, nbc
    nstep = size(steps)

    write( nfile, * ) "-----Information of steps:",nstep

    do i=1,nstep
      write( nfile, * ) "  -----Step:",i
      write(nfile,*) steps(i)%solution, steps(i)%elapsetime, steps(i)%converg, &
        steps(i)%num_substep, steps(i)%max_iter
      if( associated( steps(i)%Boundary ) ) then
        nbc = size( steps(i)%Boundary )
        write(nfile,*) "  Boundary conditions"
        write(nfile,*) ( steps(i)%Boundary(j),j=1,nbc )
      endif
      if( associated( steps(i)%Load ) ) then
        nbc = size( steps(i)%Load )
        write(nfile,*) "  External load conditions"
        write(nfile,*) ( steps(i)%Load(j),j=1,nbc )
      endif
      if( associated( steps(i)%Contact ) ) then
        nbc = size( steps(i)%Contact )
        write(nfile,*) "  Contact conditions"
        write(nfile,*) ( steps(i)%Contact(j),j=1,nbc )
      endif
    enddo
  end subroutine

  !> Initializer
  subroutine init_AincParam( aincparam )
    type( tParamAutoInc ), intent(out) :: aincparam !< auto increment parameter

    aincparam%name      = ''
    aincparam%ainc_Rs   = 0.25d0
    aincparam%ainc_Rl   = 1.25d0
    aincparam%NRbound_s = 0
    aincparam%NRbound_s(knstMAXIT) = 10
    aincparam%NRbound_s(knstSUMIT) = 50
    aincparam%NRbound_s(knstCITER) = 10
    aincparam%NRbound_l = 0
    aincparam%NRbound_l(knstMAXIT) = 1
    aincparam%NRbound_l(knstSUMIT) = 1
    aincparam%NRbound_l(knstCITER) = 1
    aincparam%NRtimes_s = 1
    aincparam%NRtimes_l = 2
    aincparam%ainc_Rc   = 0.25d0
    aincparam%CBbound   = 5
  end subroutine

end module
