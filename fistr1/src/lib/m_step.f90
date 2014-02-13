!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN (AdavanceSoft), K. SATOH (Advancesoft) !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!>   This module manages step infomation
!>
!>  \author     X. YUAN (AdavanceSoft), K. SATOH (Advancesoft)
!>  \date       2009/10/29
!>  \version    0.00
!!
!======================================================================!
module m_step
   use hecmw
   implicit none
   
   include 'fstr_ctrl_util_f.inc'
   
   integer, parameter :: stepStatic = 1
   integer, parameter :: stepVisco  = 2

   !> Step control such as active boundary condition, convergent condition etc.   
   type step_info
      integer             :: solution           !< solution type; 1: static;  2:visco
      character( len=80 ) :: CONTROL            !< control type, such as arclength etc
      character( len=80 ) :: ConvControl        !< Judgement of convergency, such as nodal force residual
                                                !< disp increment, energy
      real(kind=kreal)    :: converg            !< value of convergent judgement
      
      integer :: num_substep                    !< substeps user given
      integer :: max_iter                       !< max number of iteration
      integer :: amp_id                         !< id of amplitude definition
      real(kind=kreal) :: initdt                !< time increment
      real(kind=kreal) :: elapsetime            !< elapse time of this step
      integer, pointer :: Boundary(:)=>null()   !< active group of boundary conditions of current step
      integer, pointer :: Load(:)=>null()       !< active group of external load conditions of current step
      integer, pointer :: Contact(:)=>null()    !< active group of contact conditions of current step
   end type

   contains
   
     !> Initializer
     subroutine init_stepInfo( stepinfo )
        type( step_info ), intent(out) :: stepinfo !< step info
        stepinfo%solution = stepStatic
        stepinfo%num_substep = 1
        stepinfo%max_iter = 50
        stepinfo%amp_id = -1
        stepinfo%initdt = 1.d0
        stepinfo%elapsetime = 1.d0
        stepinfo%converg = 1.d-3
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

end module
