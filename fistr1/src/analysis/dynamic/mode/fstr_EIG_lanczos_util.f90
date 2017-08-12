!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_fstr_EIG_lanczos_util

  contains

!> Initialize Lanczos iterations
  subroutine lanczos_set_initial_value(hecMESH, hecMAT, fstrEIG, eigvec, p, q, beta)
    use m_fstr
    use hecmw_util
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix    ) :: hecMAT
    type(fstr_eigen        ) :: fstrEIG
    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: i, j
    real(kind=kreal) :: eigvec(:,:), p(:), beta, chk, sigma
    real(kind=kreal), allocatable :: temp(:)
    real(kind=kreal), pointer :: q(:), mass(:), filter(:)

    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    sigma  =  fstrEIG%sigma
    mass   => fstrEIG%mass
    filter => fstrEIG%filter

    allocate(temp(NNDOF))
    temp = 0.0d0

    !> shifting
    !do i = 1,NP
    !  hecMAT%D(9*i-8) = hecMAT%D(9*i-8) + sigma * mass(3*i-2)
    !  hecMAT%D(9*i-4) = hecMAT%D(9*i-4) + sigma * mass(3*i-1)
    !  hecMAT%D(9*i  ) = hecMAT%D(9*i  ) + sigma * mass(3*i  )
    !end do

    call URAND1(NNDOF, temp, hecMESH%my_rank)

    do i=1, NNDOF
      temp(i) = temp(i) * filter(i)
    end do

    !> M-orthogonalization
    do i=1, NNDOF
      eigvec(i,1) = mass(i) * temp(i)
    enddo

    chk = 0.0d0
    do i=1, NNDOF
      chk = chk + temp(i) * eigvec(i,1)
    enddo
    call hecmw_allreduce_R1(hecMESH, chk, hecmw_sum)
    beta = dsqrt(chk)

    if(beta == 0.0d0)then
      call hecmw_finalize()
      stop "Self-orthogonal"
    endif

    chk = 1.0d0/beta
    do i=1, NNDOF
      q(i) = temp(i) * chk
    enddo

    do i=1, NNDOF
      p(i) = mass(i) * q(i)
    enddo
  end subroutine lanczos_set_initial_value

!> Sort eigenvalues
      SUBROUTINE EVSORT(EIG,NEW,NEIG)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION EIG(NEIG),NEW(NEIG)

      DO 10 I=1,NEIG
        NEW(I)=I
   10 CONTINUE

      NM=NEIG-1
      DO 20 I=1,NM
        MINLOC=I
        EMIN=ABS(EIG(NEW(I)))
        IP=I+1
        DO 30 J=IP,NEIG
          IF(ABS(EIG(NEW(J))).LT.EMIN) THEN
            MINLOC=J
            EMIN=ABS(EIG(NEW(J)))
          END IF
   30   CONTINUE
        IBAF=NEW(I)
        NEW(I)=NEW(MINLOC)
        NEW(MINLOC)=IBAF
   20 CONTINUE

      RETURN
      END SUBROUTINE EVSORT



!> Scalar product of two vectors
      SUBROUTINE VECPRO1(Z,X,Y,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)

      DIMENSION X(NN),Y(NN)
        Z=0.0D0
        DO 20 J=1,NN
          Z=Z+X(J)*Y(J)
   20   CONTINUE

      RETURN
      END SUBROUTINE VECPRO1

!> Product of diagonal matrix and vector
      SUBROUTINE MATPRO(Y,A,X,MM,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION Y(MM,NN),A(MM),X(MM,NN)
      DO 10 I=1,NN
        DO 20 J=1,MM
          Y(J,I)=A(J)*X(J,I)
   20   CONTINUE
   10 CONTINUE
      RETURN
      END SUBROUTINE MATPRO




  subroutine URAND1(N, X, SHIFT)
    use hecmw
    implicit none
    REAL(kind=kreal) :: X(N), INVM
    INTEGER(kind=kint), parameter :: MM = 1664501
    INTEGER(kind=kint), parameter :: LAMBDA = 1229
    INTEGER(kind=kint), parameter :: MU = 351750
    INTEGER(kind=kint) :: i, N, IR, SHIFT

    IR = 9999991
    INVM = 1.0D0 / MM
    DO I = 1, SHIFT
      IR = MOD( LAMBDA * IR + MU, MM)
    enddo
    DO I = SHIFT+1, SHIFT+N
      IR = MOD( LAMBDA * IR + MU, MM)
      X(I-SHIFT) = INVM * IR
    enddo
  end subroutine URAND1

end module m_fstr_EIG_lanczos_util

