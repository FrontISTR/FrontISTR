!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_RIF_nn
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_contact
  use hecmw_matrix_misc

  private

  public:: hecmw_precond_RIF_nn_setup
  public:: hecmw_precond_RIF_nn_apply
  public:: hecmw_precond_RIF_nn_clear

  integer(4),parameter :: krealp = 8

  integer(kind=kint) :: NPFIU, NPFIL
  integer(kind=kint) :: N
  integer(kind=kint), pointer :: inumFI1L(:) => null()
  integer(kind=kint), pointer :: inumFI1U(:) => null()
  integer(kind=kint), pointer :: FI1L(:) => null()
  integer(kind=kint), pointer :: FI1U(:) => null()

  integer(kind=kint), pointer :: indexL(:) => null()
  integer(kind=kint), pointer :: indexU(:) => null()
  integer(kind=kint), pointer :: itemL(:)  => null()
  integer(kind=kint), pointer :: itemU(:)  => null()
  integer(kind=kint), pointer :: perm(:)   => null()
  integer(kind=kint), pointer :: iperm(:)  => null()
  real(kind=kreal), pointer :: D(:)  => null()
  real(kind=kreal), pointer :: AL(:) => null()
  real(kind=kreal), pointer :: AU(:) => null()

  real(kind=krealp), pointer :: SAINVU(:) => null()
  real(kind=krealp), pointer :: SAINVL(:) => null()
  real(kind=krealp), pointer :: SAINVD(:) => null()
  real(kind=krealp), pointer :: RIFU(:)   => null()
  real(kind=krealp), pointer :: RIFL(:)   => null()
  real(kind=krealp), pointer :: RIFD(:)   => null()

contains

!C***
!C*** hecmw_precond_nn_sainv_setup
!C***
  subroutine hecmw_precond_RIF_nn_setup(hecMAT)
    implicit none
    type(hecmwST_matrix) :: hecMAT

    integer(kind=kint ) :: NP, NPU, NPL
    integer(kind=kint ) :: rcm
    integer(kind=kint ) :: ii, i, j, k, ierror, itr, PRECOND, NDOF, NDOF2

    real(kind=krealp) :: THRESH, FILTER

    N = hecMAT%N
    NDOF  = hecMAT%NDOF
    NDOF2 = NDOF*NDOF
    PRECOND = hecmw_mat_get_precond(hecMAT)

    D => hecMAT%D
    AU=> hecMAT%AU
    AL=> hecMAT%AL
    indexL => hecMAT%indexL
    indexU => hecMAT%indexU
    itemL => hecMAT%itemL
    itemU => hecMAT%itemU

    if (PRECOND.eq.21) call FORM_ILU1_RIF_nn(hecMAT)

    allocate (SAINVD(NDOF2*hecMAT%NP))
    allocate (SAINVL(NDOF2*NPFIU))
    allocate (RIFD(NDOF2*hecMAT%NP))
    allocate (RIFU(NDOF2*NPFIU))
    SAINVD  = 0.0d0
    SAINVL  = 0.0d0
    RIFD  = 0.0d0
    RIFU  = 0.0d0

    FILTER= hecMAT%Rarray(5)

    Write(*,"(a,F15.8)")"### RIF FILTER   :",FILTER

    call hecmw_rif_nn(hecMAT)

    allocate (RIFL(NDOF2*NPFIU))
    RIFL  = 0.0d0

    call hecmw_rif_make_u_nn(hecMAT)

  end subroutine hecmw_precond_RIF_nn_setup

  subroutine hecmw_precond_RIF_nn_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout)  :: ZP(:)
    integer(kind=kint) :: in, i, j, isL, ieL, isU, ieU, k, iold, rcm,NDOF,NDOF2,idof,jdof
    real(kind=kreal) :: SW1, SW2, SW3, X1, X2, X3
    real(kind=kreal) :: SW(NDOF),X(NDOF)
    NDOF2=NDOF*NDOF
      !C-- FORWARD
      do i= 1, N
        do idof = 1, NDOF
          SW(idof) = ZP(NDOF*(i-1)+idof)
        end do 
        isL= inumFI1L(i-1)+1
        ieL= inumFI1L(i)
        do j= isL, ieL
          in= FI1L(j)
          do idof = 1, NDOF
            X(idof) = ZP(NDOF*(in-1)+idof)
          end do 
          do idof = 1, NDOF
            do jdof = 1, NDOF
              SW(idof) = SW(idof) - RIFL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)*X(jdof)
            end do 
          end do 
        enddo
        X = SW
        do idof = 2,NDOF
          do jdof = 1, idof-1
             X(idof) = X(idof) - RIFD(NDOF2*(i-1)+NDOF*(idof-1)+jdof )*X(jdof)
          end do 
        end do
        ZP(NDOF*(i-1)+1:NDOF*(i-1)+NDOF) = X(1:NDOF)
      enddo

      do i=1, N
        do idof=1,NDOF
          ZP(NDOF*(i-1)+idof)= ZP(NDOF*(i-1)+idof)*RIFD(NDOF2*(i-1)+(idof-1)*NDOF+idof)
        end do
      enddo
    !C-- BACKWARD
      do i= N, 1, -1
        do idof = 1, NDOF
          SW(idof) = ZP(NDOF*(i-1)+idof)
        end do 
        isL= inumFI1U(i-1) + 1
        ieL= inumFI1U(i)
        do j= ieL, isL, -1
          in= FI1U(j)
          do idof = 1, NDOF
            X(idof) = ZP(NDOF*(in-1)+idof)
          end do 
          do idof = 1, NDOF
            do jdof = 1, NDOF
              SW(idof) = SW(idof) - RIFU(NDOF2*(j-1)+NDOF*(idof-1)+jdof)*X(jdof)
            end do 
          end do 
        enddo
        X=SW
        do idof = NDOF, 1, -1
          do jdof = NDOF, idof+1, -1
            X(idof) = X(idof) - RIFD(NDOF*NDOF*(i-1)+NDOF*(jdof-1)+idof)*X(jdof)          
          end do 
        end do
        do idof = 1, NDOF
          ZP(NDOF*(i-1)+idof) = X(idof)
        end do 
      enddo
  end subroutine hecmw_precond_RIF_nn_apply


!C***
!C*** hecmw_rif_nn
!C***
    subroutine hecmw_rif_nn(hecMAT)
    implicit none
    type (hecmwST_matrix)     :: hecMAT
    integer(kind=kint) :: i, j, jS, jE, in, itr, iitr, ind(9), PRECOND, NP, NDOF, NDOF2,idof,jdof,k
    real(kind=krealp) :: YV1, YV2, YV3, X1, X2, X3, dd, dd1, dd2, dd3, dtemp(3)
    real(kind=krealp) :: YV(hecMAT%NDOF), X(hecMAT%NDOF)
    real(kind=krealp) :: FILTER, SIGMA_DIAG
    real(kind=krealp), allocatable :: zz(:), vv(:)

    FILTER= hecMAT%Rarray(5)

    NP = hecMAT%NP
    NDOF = hecMAT%NDOF
    NDOF2= NDOF*NDOF
    allocate(vv(NDOF*NP))
    allocate(zz(NDOF*NP))

    do itr=1,NP

    !------------------------------ iitr = 1 ----------------------------------------

    zz(:) = 0.0d0
    vv(:) = 0.0d0

    !{v}=[A]{zi}
    do idof = 1,NDOF
      jdof=1
      zz(NDOF*(itr-1)+idof)= SAINVD(NDOF2*(itr-1)+NDOF*(idof-1)+jdof)
    end do 

    zz(3*itr-2)= 1.0d0

    jS= inumFI1L(itr-1) + 1
    jE= inumFI1L(itr  )
    do j= jS, jE
      in  = FI1L(j)
      do idof = 1, NDOF
        zz(NDOF*(in-1)+idof)=SAINVL(NDOF2*(j-1)+NDOF*(1-1)+idof)
      end do 
!      zz(3*in-2)= SAINVL(9*j-8)
!      zz(3*in-1)= SAINVL(9*j-7)
!      zz(3*in  )= SAINVL(9*J-6)
    enddo

    do i= 1, itr
      do idof = 1,NDOF
        X(idof)=zz(NDOF*(i-1)+idof)
      end do 
      !X(1)= zz(3*i-2)
      !X(2)= zz(3*i-1)
      !X(3)= zz(3*i  )
      !vv(3*i-2) = vv(3*i-2) + D(9*i-8)*X(1) + D(9*i-7)*X(2) + D(9*i-6)*X(3)
      !vv(3*i-1) = vv(3*i-1) + D(9*i-5)*X(1) + D(9*i-4)*X(2) + D(9*i-3)*X(3)
      !vv(3*i  ) = vv(3*i  ) + D(9*i-2)*X(1) + D(9*i-1)*X(2) + D(9*i  )*X(3)
      do idof = 1, NDOF
        do jdof = 1, NDOF
          vv(NDOF*(i-1)+idof) = vv(NDOF*(i-1)+idof) + D(NDOF2*(i-1)+NDOF*(idof-1)+jdof)*X(jdof)
        end do 
      end do 

      jS= indexL(i-1) + 1
      jE= indexL(i  )
      do j=jS,jE
        in = itemL(j)
!        vv(3*in-2)= vv(3*in-2) + AL(9*j-8)*X(1) + AL(9*j-5)*X(2) + AL(9*j-2)*X(3)
!        vv(3*in-1)= vv(3*in-1) + AL(9*j-7)*X(1) + AL(9*j-4)*X(2) + AL(9*j-1)*X(3)
!        vv(3*in  )= vv(3*in  ) + AL(9*j-6)*X(1) + AL(9*j-3)*X(2) + AL(9*j  )*X(3)
        do idof = 1, NDOF
          do jdof = 1, NDOF
            vv(NDOF*(in-1)+idof) = vv(NDOF*(in-1)+idof) + AL(NDOF2*(j-1)+NDOF*(jdof-1)+idof)*X(jdof)
          end do 
        end do 
      enddo

      jS= indexU(i-1) + 1
      jE= indexU(i  )
      do j= jS, jE
        in = itemU(j)
!        vv(3*in-2)= vv(3*in-2) + AU(9*j-8)*X(1) + AU(9*j-5)*X(2) + AU(9*j-2)*X(3)
!        vv(3*in-1)= vv(3*in-1) + AU(9*j-7)*X(1) + AU(9*j-4)*X(2) + AU(9*j-1)*X(3)
!        vv(3*in  )= vv(3*in  ) + AU(9*j-6)*X(1) + AU(9*j-3)*X(2) + AU(9*j  )*X(3)
        do idof = 1, NDOF
          do jdof = 1, NDOF
            vv(NDOF*(in-1)+idof) = vv(NDOF*(in-1)+idof) + AU(NDOF2*(j-1)+NDOF*(jdof-1)+idof)*X(jdof)
          end do 
        end do 
      enddo
    enddo

    !{d}={v^t}{z_j}
    do i= itr,NP
      !SAINVD(9*i-8) = vv(3*i-2)
      !SAINVD(9*i-4) = vv(3*i-2)*SAINVD(9*i-7)   + vv(3*i-1)
      !SAINVD(9*i  ) = vv(3*i-2)*SAINVD(9*i-6)   + vv(3*i-1)*SAINVD(9*i-3)  + vv(3*i)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!ここからまだ未修正！！！






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do idof = 1,NDOF
        SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof) = vv(NDOF*(i-1)+idof)
        !1 5 9 (-8,-4,0) 1 2 3 (2 1 0)
        if (idof==1 .and. (9*i-8)==NDOF2*(i-1)+NDOF*(idof-1)+idof) then 
          write(*,*) idof,"OK"
        else
          write(*,*) idof,"NG"
        end if
        do jdof = 1, idof-1
          !1 5 9         //     3*(2-1-1)+1 (1) 3*(3-1-1)+3(6)
          SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof) = SAINVD(NDOF2*(i-1)+NDOF*(jdof-1)+idof)*vv(NDOF*(i-1)+jdof)
        end do 
      end do
      !write(*,*)NDOF2*(i-1)+NDOF*(2-1)+3,9*i-6,9*i-3
      !1
      !2 1,2 (0*3+2)
      !3 1,6 (1*3+3) 2,3(0*3+3)
      
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      do j= jS, jE
        in  = FI1L(j)
        !X(1)= vv(3*in-2)
        !X(2)= vv(3*in-1)
        !X(3)= vv(3*in  )
        do idof = 1,NDOF
          X(idof)=vv(NDOF*(in-1)+idof)
        end do 
        !SAINVD(9*i-8)= SAINVD(9*i-8) + X(1)*SAINVL(9*j-8) + X(2)*SAINVL(9*j-7) + X(3)*SAINVL(9*j-6)
        !SAINVD(9*i-4)= SAINVD(9*i-4) + X(1)*SAINVL(9*j-5) + X(2)*SAINVL(9*j-4) + X(3)*SAINVL(9*j-3)
        !SAINVD(9*i  )= SAINVD(9*i  ) + X(1)*SAINVL(9*j-2) + X(2)*SAINVL(9*j-1) + X(3)*SAINVL(9*j  )
        do idof = 1, NDOF
          do jdof = 1, NDOF
            SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof) = SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof) + &
                   & SAINVL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)*X(jdof)
          end do 
        end do 
      enddo
    enddo

    !Update D
    !dd = 1.0d0/SAINVD(9*itr-8)
    dd = 1.0d0/SAINVD(NDOF2*(itr-1)+1)
    SAINVD(NDOF2*(itr-1)+NDOF*(1-1)+1)=dd
    do idof = 2, NDOF
      SAINVD(NDOF2*(itr-1)+NDOF*(idof-1)+idof) = SAINVD(NDOF2*(itr-1)+NDOF*(idof-1)+idof)*dd
    end do 
!    SAINVD(9*itr-4) =SAINVD(9*itr-4)*dd
!    SAINVD(9*itr  ) =SAINVD(9*itr  )*dd

    do i =itr+1,NP
      do idof = 1, NDOF
        SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof) = SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof)*dd
      end do 
!      SAINVD(9*i-8) = SAINVD(9*i-8)*dd
!      SAINVD(9*i-4) = SAINVD(9*i-4)*dd
!      SAINVD(9*i  ) = SAINVD(9*i  )*dd
    enddo

    do idof = 1, NDOF
      RIFD(NDOF2*(itr-1)+NDOF*(idof-1)+1) = SAINVD(NDOF2*(itr-1)+NDOF*(idof-1)+idof)
    end do 
!    write(*,*)SAINVD(NDOF2*(itr-1)+NDOF*(idof-1)+idof),dd
    !RIFD(9*itr-8) = dd                       !RIF UPDATE
    !RIFD(9*itr-5) = SAINVD(9*itr-4)          !RIF UPDATE
    !RIFD(9*itr-2) = SAINVD(9*itr  )          !RIF UPDATE

    jS= inumFI1U(itr-1) + 1
    jE= inumFI1U(itr  )
    do j= jS, jE
      do idof = 1, NDOF
        RIFU(NDOF2*(j-1)+NDOF*(1-1)+idof) = SAINVD(NDOF2*(FI1U(j)-1)+NDOF*(idof-1)+idof)
      end do 
!      RIFU(9*j-8) = SAINVD(9*FI1U(j)-8)
!      RIFU(9*j-7) = SAINVD(9*FI1U(j)-4)
!      RIFU(9*j-6) = SAINVD(9*FI1U(j)  )
    enddo

    !Update Z
    do k=2,NDOF
      dd = SAINVD(NDOF2*(itr-1)+NDOF*(k-1)+k)
      if(abs(dd) > FILTER)then
        SAINVD(NDOF2*(itr-1)+(1-1)+k)= SAINVD(NDOF2*(itr-1)+(1-1)+k) - dd*zz(NDOF*(itr-1)+1)
        jS= inumFI1L(itr-1) + 1
        jE= inumFI1L(itr  )
        do j= jS, jE
          in  = FI1L(j)
          do idof = 1, NDOF
            SAINVL(NDOF2*(j-1)+NDOF*(k-1)+idof) = SAINVL(NDOF2*(j-1)+NDOF*(k-1)+idof)-dd*zz(NDOF*(in-1)+idof)
          end do
        enddo
      endif
    end do 
    
    !dd2=SAINVD(9*itr-4)
    !if(abs(dd2) > FILTER)then
    !  SAINVD(9*itr-7)= SAINVD(9*itr-7) - dd2*zz(3*itr-2)
    !  jS= inumFI1L(itr-1) + 1
    !  jE= inumFI1L(itr  )
    !  do j= jS, jE
    !    in  = FI1L(j)
    !    SAINVL(9*j-5) = SAINVL(9*j-5)-dd2*zz(3*in-2)
    !    SAINVL(9*j-4) = SAINVL(9*j-4)-dd2*zz(3*in-1)
    !    SAINVL(9*j-3) = SAINVL(9*j-3)-dd2*zz(3*in  )
    !  enddo
    !endif

   ! dd3=SAINVD(9*itr  )
   ! if(abs(dd3) > FILTER)then
   !   SAINVD(9*itr-6)= SAINVD(9*itr-6) - dd3*zz(3*itr-2)
   !   jS= inumFI1L(itr-1) + 1
   !   jE= inumFI1L(itr  )
   !   do j= jS, jE
   !     in  = FI1L(j)
   !     SAINVL(9*j-2) = SAINVL(9*j-2)-dd3*zz(3*in-2)
   !     SAINVL(9*j-1) = SAINVL(9*j-1)-dd3*zz(3*in-1)
   !     SAINVL(9*j  ) = SAINVL(9*j  )-dd3*zz(3*in  )
   !   enddo
   ! endif
    do i= itr +1,NP
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      do idof = 1, NDOF
        dd = SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof)
        if(abs(dd) > FILTER)then
          do j= jS, jE
            in  = FI1L(j)
            if (in > itr) exit
            do jdof=1,NDOF
            SAINVL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)=SAINVL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)-dd*zz(NDOF*(in-1)+jdof)
            end do 
          enddo
        endif
      end do 
      !dd1=SAINVD(9*i-8)
      !if(abs(dd1) > FILTER)then
      !  do j= jS, jE
      !    in  = FI1L(j)
      !    if (in > itr) exit
      !    SAINVL(9*j-8) = SAINVL(9*j-8)-dd1*zz(3*in-2)
      !    SAINVL(9*j-7) = SAINVL(9*j-7)-dd1*zz(3*in-1)
      !    SAINVL(9*j-6) = SAINVL(9*j-6)-dd1*zz(3*in  )
      !  enddo
      !endif
      !dd2=SAINVD(9*i-4)
      !if(abs(dd2) > FILTER)then
      !  do j= jS, jE
      !    in  = FI1L(j)
      !    if (in > itr) exit
      !    SAINVL(9*j-5) = SAINVL(9*j-5)-dd2*zz(3*in-2)
      !    SAINVL(9*j-4) = SAINVL(9*j-4)-dd2*zz(3*in-1)
      !    SAINVL(9*j-3) = SAINVL(9*j-3)-dd2*zz(3*in  )
      !  enddo
      !endif
      !dd3=SAINVD(9*i  )
      !if(abs(dd3) > FILTER)then
      !  do j= jS, jE
      !    in  = FI1L(j)
      !    if (in > itr) exit
      !    SAINVL(9*j-2) = SAINVL(9*j-2)-dd3*zz(3*in-2)
      !    SAINVL(9*j-1) = SAINVL(9*j-1)-dd3*zz(3*in-1)
      !    SAINVL(9*j  ) = SAINVL(9*j  )-dd3*zz(3*in  )
      !  enddo
      !endif
    enddo

    !------------------------------ iitr = 1 ----------------------------------------

    zz(:) = 0.0d0
    vv(:) = 0.0d0

    !{v}=[A]{zi}

    zz(3*itr-2)= SAINVD(9*itr-7)
    zz(3*itr-1)= SAINVD(9*itr-4)
    zz(3*itr  )= SAINVD(9*itr-1)

    zz(3*itr-1)= 1.0d0

    jS= inumFI1L(itr-1) + 1
    jE= inumFI1L(itr  )
    do j= jS, jE
      in  = FI1L(j)
      zz(3*in-2)= SAINVL(9*j-5)
      zz(3*in-1)= SAINVL(9*j-4)
      zz(3*in  )= SAINVL(9*J-3)
    enddo

    !do i= 1, itr
    !  X(1)= zz(3*i-2)
    !  X(2)= zz(3*i-1)
    !  X(3)= zz(3*i  )
    !  vv(3*i-2) = vv(3*i-2) + D(9*i-8)*X(1) + D(9*i-7)*X(2) + D(9*i-6)*X(3)
    !  vv(3*i-1) = vv(3*i-1) + D(9*i-5)*X(1) + D(9*i-4)*X(2) + D(9*i-3)*X(3)
    !  vv(3*i  ) = vv(3*i  ) + D(9*i-2)*X(1) + D(9*i-1)*X(2) + D(9*i  )*X(3)
    !
    !  jS= indexL(i-1) + 1
    !  jE= indexL(i  )
    !  do j=jS,jE
    !    in = itemL(j)
    !    vv(3*in-2)= vv(3*in-2) + AL(9*j-8)*X(1) + AL(9*j-5)*X(2) + AL(9*j-2)*X(3)
    !    vv(3*in-1)= vv(3*in-1) + AL(9*j-7)*X(1) + AL(9*j-4)*X(2) + AL(9*j-1)*X(3)
    !    vv(3*in  )= vv(3*in  ) + AL(9*j-6)*X(1) + AL(9*j-3)*X(2) + AL(9*j  )*X(3)
    !  enddo
    !
    !  jS= indexU(i-1) + 1
    !  jE= indexU(i  )
    !  do j= jS, jE
    !    in = itemU(j)
    !    vv(3*in-2)= vv(3*in-2) + AU(9*j-8)*X(1) + AU(9*j-5)*X(2) + AU(9*j-2)*X(3)
    !    vv(3*in-1)= vv(3*in-1) + AU(9*j-7)*X(1) + AU(9*j-4)*X(2) + AU(9*j-1)*X(3)
    !    vv(3*in  )= vv(3*in  ) + AU(9*j-6)*X(1) + AU(9*j-3)*X(2) + AU(9*j  )*X(3)
    !  enddo
    !enddo
    do i= 1, itr
      do idof = 1,NDOF
        X(idof)=zz(NDOF*(i-1)+idof)
      end do 
      do idof = 1, NDOF
        do jdof = 1, NDOF
          vv(NDOF*(i-1)+idof) = vv(NDOF*(i-1)+idof) + D(NDOF2*(i-1)+NDOF*(idof-1)+jdof)*X(jdof)
        end do 
      end do 
      jS= indexL(i-1) + 1
      jE= indexL(i  )
      do j=jS,jE
        in = itemL(j)
        do idof = 1, NDOF
          do jdof = 1, NDOF
            vv(NDOF*(in-1)+idof) = vv(NDOF*(in-1)+idof) + AL(NDOF2*(j-1)+NDOF*(jdof-1)+idof)*X(jdof)
          end do 
        end do 
      enddo
      jS= indexU(i-1) + 1
      jE= indexU(i  )
      do j= jS, jE
        in = itemU(j)
        do idof = 1, NDOF
          do jdof = 1, NDOF
            vv(NDOF*(in-1)+idof) = vv(NDOF*(in-1)+idof) + AU(NDOF2*(j-1)+NDOF*(jdof-1)+idof)*X(jdof)
          end do 
        end do 
      enddo
    enddo

    
    
    
    !{d}={v^t}{z_j}
    dtemp(1) = SAINVD(9*itr-8)

    do i=itr,NP
      SAINVD(9*i-8) = vv(3*i-2)
      SAINVD(9*i-4) = vv(3*i-2)*SAINVD(9*i-7)   + vv(3*i-1)
      SAINVD(9*i  ) = vv(3*i-2)*SAINVD(9*i-6)   + vv(3*i-1)*SAINVD(9*i-3)  + vv(3*i)
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      do j= jS, jE
        in  = FI1L(j)
        X(1)= vv(3*in-2)
        X(2)= vv(3*in-1)
        X(3)= vv(3*in  )
        SAINVD(9*i-8)= SAINVD(9*i-8) + X(1)*SAINVL(9*j-8) + X(2)*SAINVL(9*j-7) + X(3)*SAINVL(9*j-6)
        SAINVD(9*i-4)= SAINVD(9*i-4) + X(1)*SAINVL(9*j-5) + X(2)*SAINVL(9*j-4) + X(3)*SAINVL(9*j-3)
        SAINVD(9*i  )= SAINVD(9*i  ) + X(1)*SAINVL(9*j-2) + X(2)*SAINVL(9*j-1) + X(3)*SAINVL(9*j  )
      enddo
    enddo

    !Update D
    dd = 1.0d0/SAINVD(9*itr-4)

    SAINVD(9*itr-8) = dtemp(1)
    SAINVD(9*itr  ) =SAINVD(9*itr  )*dd

    do i =itr+1,NP
      SAINVD(9*i-8) = SAINVD(9*i-8)*dd
      SAINVD(9*i-4) = SAINVD(9*i-4)*dd
      SAINVD(9*i  ) = SAINVD(9*i  )*dd
    enddo

    RIFD(9*itr-4) = dd                       !RIF UPDATE
    RIFD(9*itr-1) = SAINVD(9*itr  )          !RIF UPDATE

    jS= inumFI1U(itr-1) + 1
    jE= inumFI1U(itr  )
    do j= jS, jE
      RIFU(9*j-5) = SAINVD(9*FI1U(j)-8)
      RIFU(9*j-4) = SAINVD(9*FI1U(j)-4)
      RIFU(9*j-3) = SAINVD(9*FI1U(j)  )
    enddo

    !Update Z

    dd3=SAINVD(9*itr  )
    if(abs(dd3) > FILTER)then
      SAINVD(9*itr-6)= SAINVD(9*itr-6) - dd3*zz(3*itr-2)
      SAINVD(9*itr-3)= SAINVD(9*itr-3) - dd3*zz(3*itr-1)

      jS= inumFI1L(itr-1) + 1
      jE= inumFI1L(itr  )
      do j= jS, jE
        in  = FI1L(j)
        SAINVL(9*j-2) = SAINVL(9*j-2)-dd3*zz(3*in-2)
        SAINVL(9*j-1) = SAINVL(9*j-1)-dd3*zz(3*in-1)
        SAINVL(9*j  ) = SAINVL(9*j  )-dd3*zz(3*in  )
      enddo
    endif
    
    
    
    
    do i= itr +1,NP
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      do idof = 1, NDOF
        dd = SAINVD(NDOF2*(i-1)+NDOF*(idof-1)+idof)
        if(abs(dd) > FILTER)then
          do j= jS, jE
            in  = FI1L(j)
            if (in > itr) exit
            do jdof=1,NDOF
            SAINVL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)=SAINVL(NDOF2*(j-1)+NDOF*(idof-1)+jdof)-dd*zz(NDOF*(in-1)+jdof)
            end do 
          enddo
        endif
      end do 
    enddo



    !------------------------------ iitr = 1 ----------------------------------------

    zz(:) = 0.0d0
    vv(:) = 0.0d0

    !{v}=[A]{zi}

    zz(3*itr-2)= SAINVD(9*itr-6)
    zz(3*itr-1)= SAINVD(9*itr-3)
    zz(3*itr  )= SAINVD(9*itr  )

    zz(3*itr  )= 1.0d0

    jS= inumFI1L(itr-1) + 1
    jE= inumFI1L(itr  )
    do j= jS, jE
      in  = FI1L(j)
      zz(3*in-2)= SAINVL(9*j-2)
      zz(3*in-1)= SAINVL(9*j-1)
      zz(3*in  )= SAINVL(9*J  )
    enddo

    do i= 1, itr
      X(1)= zz(3*i-2)
      X(2)= zz(3*i-1)
      X(3)= zz(3*i  )
      vv(3*i-2) = vv(3*i-2) + D(9*i-8)*X(1) + D(9*i-7)*X(2) + D(9*i-6)*X(3)
      vv(3*i-1) = vv(3*i-1) + D(9*i-5)*X(1) + D(9*i-4)*X(2) + D(9*i-3)*X(3)
      vv(3*i  ) = vv(3*i  ) + D(9*i-2)*X(1) + D(9*i-1)*X(2) + D(9*i  )*X(3)

      jS= indexL(i-1) + 1
      jE= indexL(i  )
      do j=jS,jE
        in = itemL(j)
        vv(3*in-2)= vv(3*in-2) + AL(9*j-8)*X(1) + AL(9*j-5)*X(2) + AL(9*j-2)*X(3)
        vv(3*in-1)= vv(3*in-1) + AL(9*j-7)*X(1) + AL(9*j-4)*X(2) + AL(9*j-1)*X(3)
        vv(3*in  )= vv(3*in  ) + AL(9*j-6)*X(1) + AL(9*j-3)*X(2) + AL(9*j  )*X(3)
      enddo

      jS= indexU(i-1) + 1
      jE= indexU(i  )
      do j= jS, jE
        in = itemU(j)
        vv(3*in-2)= vv(3*in-2) + AU(9*j-8)*X(1) + AU(9*j-5)*X(2) + AU(9*j-2)*X(3)
        vv(3*in-1)= vv(3*in-1) + AU(9*j-7)*X(1) + AU(9*j-4)*X(2) + AU(9*j-1)*X(3)
        vv(3*in  )= vv(3*in  ) + AU(9*j-6)*X(1) + AU(9*j-3)*X(2) + AU(9*j  )*X(3)
      enddo
    enddo

    !{d}={v^t}{z_j}

    dtemp(1) = SAINVD(9*itr-8)
    dtemp(2) = SAINVD(9*itr-4)

    do i=itr,NP
      SAINVD(9*i-8) = vv(3*i-2)
      SAINVD(9*i-4) = vv(3*i-2)*SAINVD(9*i-7)   + vv(3*i-1)
      SAINVD(9*i  ) = vv(3*i-2)*SAINVD(9*i-6)   + vv(3*i-1)*SAINVD(9*i-3)  + vv(3*i)
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      do j= jS, jE
        in  = FI1L(j)
        X(1)= vv(3*in-2)
        X(2)= vv(3*in-1)
        X(3)= vv(3*in  )
        SAINVD(9*i-8)= SAINVD(9*i-8) + X(1)*SAINVL(9*j-8) + X(2)*SAINVL(9*j-7) + X(3)*SAINVL(9*j-6)
        SAINVD(9*i-4)= SAINVD(9*i-4) + X(1)*SAINVL(9*j-5) + X(2)*SAINVL(9*j-4) + X(3)*SAINVL(9*j-3)
        SAINVD(9*i  )= SAINVD(9*i  ) + X(1)*SAINVL(9*j-2) + X(2)*SAINVL(9*j-1) + X(3)*SAINVL(9*j  )
      enddo
    enddo

    !Update D
    dd = 1.0d0/SAINVD(9*itr  )

    SAINVD(9*itr-8) = dtemp(1)
    SAINVD(9*itr-4) = dtemp(2)

    do i =itr+1,NP
      SAINVD(9*i-8) = SAINVD(9*i-8)*dd
      SAINVD(9*i-4) = SAINVD(9*i-4)*dd
      SAINVD(9*i  ) = SAINVD(9*i  )*dd
    enddo

    RIFD(9*itr) = dd   !RIF UPDATE

    jS= inumFI1U(itr-1) + 1
    jE= inumFI1U(itr  )
    do j= jS, jE
      RIFU(9*j-2) = SAINVD(9*FI1U(j)-8)
      RIFU(9*j-1) = SAINVD(9*FI1U(j)-4)
      RIFU(9*j  ) = SAINVD(9*FI1U(j)  )
    enddo

    !Update Z

    do i= itr +1,NP
      jS= inumFI1L(i-1) + 1
      jE= inumFI1L(i  )
      dd1=SAINVD(9*i-8)
      if(abs(dd1) > FILTER)then
        do j= jS, jE
          in  = FI1L(j)
          if (in > itr) exit
          SAINVL(9*j-8) = SAINVL(9*j-8)-dd1*zz(3*in-2)
          SAINVL(9*j-7) = SAINVL(9*j-7)-dd1*zz(3*in-1)
          SAINVL(9*j-6) = SAINVL(9*j-6)-dd1*zz(3*in  )
        enddo
      endif
      dd2=SAINVD(9*i-4)
      if(abs(dd2) > FILTER)then
        do j= jS, jE
          in  = FI1L(j)
          if (in > itr) exit
          SAINVL(9*j-5) = SAINVL(9*j-5)-dd2*zz(3*in-2)
          SAINVL(9*j-4) = SAINVL(9*j-4)-dd2*zz(3*in-1)
          SAINVL(9*j-3) = SAINVL(9*j-3)-dd2*zz(3*in  )
        enddo
      endif
      dd3=SAINVD(9*i  )
      if(abs(dd3) > FILTER)then
        do j= jS, jE
          in  = FI1L(j)
          if (in > itr) exit
          SAINVL(9*j-2) = SAINVL(9*j-2)-dd3*zz(3*in-2)
          SAINVL(9*j-1) = SAINVL(9*j-1)-dd3*zz(3*in-1)
          SAINVL(9*j  ) = SAINVL(9*j  )-dd3*zz(3*in  )
        enddo
      endif
    enddo
    enddo
    deallocate(vv)
    deallocate(zz)

  end subroutine hecmw_rif_nn

  subroutine hecmw_rif_make_u_nn(hecMAT)
    implicit none
    type (hecmwST_matrix)     :: hecMAT
    integer(kind=kint) i,j,k,kk,n,m,o,idof,jdof,NDOF,NDOF2
    integer(kind=kint) is,ie,js,je
    NDOF=hecMAT%NDOF
    NDOF2=NDOF*NDOF
    n = 1
    do i= 1, hecMAT%NP
      is=inumFI1L(i-1) + 1
      ie=inumFI1L(i  )
      flag1:do k= is, ie
        m = FI1L(k)
        js=inumFI1U(m-1) + 1
        je=inumFI1U(m  )
        do j= js,je
          o = FI1U(j)
          if (o == i)then
            do idof = 1, NDOF
              do jdof = 1, NDOF
                RIFL(NDOF2*(n-1)+NDOF*(idof-1)+jdof)=RIFU(NDOF2*(j-1)+NDOF*(jdof-1)+idof)
              end do 
            end do 
            n = n + 1
            cycle flag1
          endif
        enddo
      enddo flag1
    enddo
  end subroutine hecmw_rif_make_u_nn

!C***
!C*** FORM_ILU1_nn
!C*** form ILU(1) matrix
  subroutine FORM_ILU0_RIF_nn(hecMAT)
    implicit none
    type(hecmwST_matrix) :: hecMAT

    allocate (inumFI1L(0:hecMAT%NP), inumFI1U(0:hecMAT%NP))
    allocate (FI1L (hecMAT%NPL), FI1U (hecMAT%NPU))

    inumFI1L = 0
    inumFI1U = 0
    FI1L = 0
    FI1U = 0

    inumFI1L(:) = hecMAT%indexL(:)
    inumFI1U(:) = hecMAT%indexU(:)
    FI1L(:) = hecMAT%itemL(:)
    FI1U(:) = hecMAT%itemU(:)

    NPFIU = hecMAT%NPU
    NPFIL = hecMAT%NPL

  end subroutine FORM_ILU0_RIF_nn

!C***
!C*** FORM_ILU1_nn
!C*** form ILU(1) matrix
  subroutine FORM_ILU1_RIF_nn(hecMAT)
    implicit none
    type(hecmwST_matrix) :: hecMAT

    integer(kind=kint),allocatable :: IWsL(:), IWsU(:), IW1(:), IW2(:)
    integer(kind=kint) :: NPLf1,NPUf1
    integer(kind=kint) :: i,jj,jj1,ij0,kk,ik,kk1,kk2,L,iSk,iEk,iSj,iEj
    integer(kind=kint) :: icou,icou0,icouU,icouU1,icouU2,icouU3,icouL,icouL1,icouL2,icouL3
    integer(kind=kint) :: j,k,iSL,iSU
    !C
    !C +--------------+
    !C | find fill-in |
    !C +--------------+
    !C===

    !C
    !C-- count fill-in
    allocate (IW1(hecMAT%NP) , IW2(hecMAT%NP))
    allocate (inumFI1L(0:hecMAT%NP), inumFI1U(0:hecMAT%NP))

    inumFI1L= 0
    inumFI1U= 0

    NPLf1= 0
    NPUf1= 0
    do i= 2, hecMAT%NP
    icou= 0
    IW1= 0
    IW1(i)= 1
    do L= indexL(i-1)+1, indexL(i)
      IW1(itemL(L))= 1
    enddo
    do L= indexU(i-1)+1, indexU(i)
      IW1(itemU(L))= 1
    enddo

      iSk= indexL(i-1) + 1
      iEk= indexL(i)
      do k= iSk, iEk
        kk= itemL(k)
        iSj= indexU(kk-1) + 1
        iEj= indexU(kk  )
        do j= iSj, iEj
          jj= itemU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            inumFI1L(i)= inumFI1L(i)+1
            IW1(jj)= 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            inumFI1U(i)= inumFI1U(i)+1
            IW1(jj)= 1
          endif
        enddo
      enddo
      NPLf1= NPLf1 + inumFI1L(i)
      NPUf1= NPUf1 + inumFI1U(i)
    enddo

    !C
    !C-- specify fill-in
    allocate (IWsL(0:hecMAT%NP), IWsU(0:hecMAT%NP))
    allocate (FI1L (hecMAT%NPL+NPLf1), FI1U (hecMAT%NPU+NPUf1))

  NPFIU = hecMAT%NPU+NPUf1
  NPFIL = hecMAT%NPL+NPLf1

    FI1L= 0
    FI1U= 0

    IWsL= 0
    IWsU= 0
    do i= 1, hecMAT%NP
      IWsL(i)= indexL(i)-indexL(i-1) + inumFI1L(i) + IWsL(i-1)
      IWsU(i)= indexU(i)-indexU(i-1) + inumFI1U(i) + IWsU(i-1)
    enddo

    do i= 2, hecMAT%NP
      icouL= 0
      icouU= 0
      inumFI1L(i)= inumFI1L(i-1) + inumFI1L(i)
      inumFI1U(i)= inumFI1U(i-1) + inumFI1U(i)
      icou= 0
      IW1= 0
      IW1(i)= 1
      do L= indexL(i-1)+1, indexL(i)
        IW1(itemL(L))= 1
      enddo
      do L= indexU(i-1)+1, indexU(i)
        IW1(itemU(L))= 1
      enddo

      iSk= indexL(i-1) + 1
      iEk= indexL(i)
      do k= iSk, iEk
        kk= itemL(k)
        iSj= indexU(kk-1) + 1
        iEj= indexU(kk  )
        do j= iSj, iEj
          jj= itemU(j)
          if (IW1(jj).eq.0 .and. jj.lt.i) then
            icouL           = icouL + 1
            FI1L(icouL+IWsL(i-1)+indexL(i)-indexL(i-1))= jj
            IW1(jj)          = 1
          endif
          if (IW1(jj).eq.0 .and. jj.gt.i) then
            icouU           = icouU + 1
            FI1U(icouU+IWsU(i-1)+indexU(i)-indexU(i-1))= jj
            IW1(jj)          = 1
          endif
        enddo
      enddo
    enddo

    iSL  = 0
    iSU  = 0
    do i= 1, hecMAT%NP
      icouL1= indexL(i) - indexL(i-1)
      icouL2= inumFI1L(i) - inumFI1L(i-1)
      icouL3= icouL1 + icouL2
      icouU1= indexU(i) - indexU(i-1)
      icouU2= inumFI1U(i) - inumFI1U(i-1)
      icouU3= icouU1 + icouU2
      !C
      !C-- LOWER part
      icou0= 0
      do k= indexL(i-1)+1, indexL(i)
        icou0 = icou0 + 1
        IW1(icou0)= itemL(k)
      enddo

      do k= inumFI1L(i-1)+1, inumFI1L(i)
        icou0 = icou0 + 1
        IW1(icou0)= FI1L(icou0+IWsL(i-1))
      enddo

      do k= 1, icouL3
        IW2(k)= k
      enddo
      call RIF_SORT_nn (IW1, IW2, icouL3, hecMAT%NP)

      do k= 1, icouL3
        FI1L (k+isL)= IW1(k)
      enddo
      !C
      !C-- UPPER part
      icou0= 0
      do k= indexU(i-1)+1, indexU(i)
        icou0 = icou0 + 1
        IW1(icou0)= itemU(k)
      enddo

      do k= inumFI1U(i-1)+1, inumFI1U(i)
        icou0 = icou0 + 1
        IW1(icou0)= FI1U(icou0+IWsU(i-1))
      enddo

      do k= 1, icouU3
        IW2(k)= k
      enddo
      call RIF_SORT_nn (IW1, IW2, icouU3, hecMAT%NP)

      do k= 1, icouU3
        FI1U (k+isU)= IW1(k)
      enddo

      iSL= iSL + icouL3
      iSU= iSU + icouU3
    enddo

    !C===
    do i= 1, hecMAT%NP
      inumFI1L(i)= IWsL(i)
      inumFI1U(i)= IWsU(i)
    enddo

    deallocate (IW1, IW2)
    deallocate (IWsL, IWsU)
    !C===
  end subroutine FORM_ILU1_RIF_nn

 !C
  !C***
  !C*** fill_in_S33_SORT
  !C***
  !C
  subroutine RIF_SORT_nn(STEM, INUM, N, NP)
    use hecmw_util
    implicit none
    integer(kind=kint) :: N, NP
    integer(kind=kint), dimension(NP) :: STEM
    integer(kind=kint), dimension(NP) :: INUM
    integer(kind=kint), dimension(:), allocatable :: ISTACK
    integer(kind=kint) :: M,NSTACK,jstack,l,ir,ip,i,j,k,ss,ii,temp,it

    allocate (ISTACK(-NP:+NP))

    M     = 100
    NSTACK= NP

    jstack= 0
    l     = 1
    ir    = N

    ip= 0
1   continue
    ip= ip + 1

    if (ir-l.lt.M) then
      do j= l+1, ir
        ss= STEM(j)
        ii= INUM(j)

        do i= j-1,1,-1
          if (STEM(i).le.ss) goto 2
          STEM(i+1)= STEM(i)
          INUM(i+1)= INUM(i)
        end do
        i= 0

2       continue
        STEM(i+1)= ss
        INUM(i+1)= ii
      end do

      if (jstack.eq.0) then
        deallocate (ISTACK)
        return
      endif

      ir = ISTACK(jstack)
      l = ISTACK(jstack-1)
      jstack= jstack - 2
    else

      k= (l+ir) / 2
      temp = STEM(k)
      STEM(k)  = STEM(l+1)
      STEM(l+1)= temp

      it = INUM(k)
      INUM(k)  = INUM(l+1)
      INUM(l+1)= it

      if (STEM(l+1).gt.STEM(ir)) then
        temp = STEM(l+1)
        STEM(l+1)= STEM(ir)
        STEM(ir )= temp
        it = INUM(l+1)
        INUM(l+1)= INUM(ir)
        INUM(ir )= it
      endif

      if (STEM(l).gt.STEM(ir)) then
        temp = STEM(l)
        STEM(l )= STEM(ir)
        STEM(ir)= temp
        it = INUM(l)
        INUM(l )= INUM(ir)
        INUM(ir)= it
      endif

      if (STEM(l+1).gt.STEM(l)) then
        temp = STEM(l+1)
        STEM(l+1)= STEM(l)
        STEM(l  )= temp
        it = INUM(l+1)
        INUM(l+1)= INUM(l)
        INUM(l  )= it
      endif

      i= l + 1
      j= ir

      ss= STEM(l)
      ii= INUM(l)

3     continue
      i= i + 1
      if (STEM(i).lt.ss) goto 3

4     continue
      j= j - 1
      if (STEM(j).gt.ss) goto 4

      if (j.lt.i)        goto 5

      temp   = STEM(i)
      STEM(i)= STEM(j)
      STEM(j)= temp

      it     = INUM(i)
      INUM(i)= INUM(j)
      INUM(j)= it

      goto 3

5     continue

      STEM(l)= STEM(j)
      STEM(j)= ss
      INUM(l)= INUM(j)
      INUM(j)= ii

      jstack= jstack + 2

      if (jstack.gt.NSTACK) then
        write (*,*) 'NSTACK overflow'
        stop
      endif

      if (ir-i+1.ge.j-1) then
        ISTACK(jstack  )= ir
        ISTACK(jstack-1)= i
        ir= j-1
      else
        ISTACK(jstack  )= j-1
        ISTACK(jstack-1)= l
        l= i
      endif

    endif

    goto 1

  end subroutine RIF_SORT_nn

  subroutine hecmw_precond_RIF_nn_clear()
    implicit none

    if (associated(SAINVD)) deallocate(SAINVD)
    if (associated(SAINVL)) deallocate(SAINVL)
    if (associated(RIFD)) deallocate(RIFD)
    if (associated(RIFU)) deallocate(RIFU)
    if (associated(RIFL)) deallocate(RIFL)
    if (associated(inumFI1L)) deallocate(inumFI1L)
    if (associated(inumFI1U)) deallocate(inumFI1U)
    if (associated(FI1L)) deallocate(FI1L)
    if (associated(FI1U)) deallocate(FI1U)
    nullify(inumFI1L)
    nullify(inumFI1U)
    nullify(FI1L)
    nullify(FI1U)
    nullify(D)
    nullify(AL)
    nullify(AU)
    nullify(indexL)
    nullify(indexU)
    nullify(itemL)
    nullify(itemU)

  end subroutine hecmw_precond_RIF_nn_clear
end module hecmw_precond_RIF_nn
