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
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG
    integer(kind=kint) :: N, NP, NDOF, NDOF2, NNDOF, NPNDOF
    integer(kind=kint) :: i, j
    real(kind=kreal)   :: eigvec(:, :), p(:), beta, chk, sigma
    real(kind=kreal), allocatable :: temp(:)
    real(kind=kreal), pointer     :: q(:), mass(:), filter(:)
    logical :: is_free

    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NDOF2  = NDOF*NDOF
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    sigma  = 0.1d0
    mass   => fstrEIG%mass
    filter => fstrEIG%filter

    allocate(temp(NNDOF))
    temp = 0.0d0

    !> shifting
    if(fstrEIG%is_free)then
      do i = 1, NP
        do j = 1, NDOF
          hecMAT%D(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1) = hecMAT%D(NDOF2*(i-1) + (NDOF+1)*(j-1) + 1) + sigma * mass(NDOF*(i-1) + j)
        enddo
      enddo
    endif

    call URAND1(NNDOF, temp, hecMESH%my_rank)

    do i = 1, NNDOF
      temp(i) = temp(i) * filter(i)
    enddo

    !> M-orthogonalization
    do i = 1, NNDOF
      eigvec(i,1) = mass(i) * temp(i)
    enddo

    chk = 0.0d0
    do i = 1, NNDOF
      chk = chk + temp(i) * eigvec(i,1)
    enddo
    call hecmw_allreduce_R1(hecMESH, chk, hecmw_sum)
    beta = dsqrt(chk)

    if(beta == 0.0d0)then
      call hecmw_finalize()
      stop "Self-orthogonal"
    endif

    chk = 1.0d0/beta
    do i = 1, NNDOF
      q(i) = temp(i) * chk
    enddo

    do i = 1, NNDOF
      p(i) = mass(i) * q(i)
    enddo
  end subroutine lanczos_set_initial_value

  !> Sort eigenvalues
  subroutine EVSORT(EIG, NEW, NEIG)
    use hecmw
    implicit none
    integer(kind=kint) :: i, j, n, ip, minloc, NEIG, IBAF, NEW(NEIG)
    real(kind=kreal) :: EMIN, EIG(NEIG)

    do i = 1, NEIG
      NEW(i) = i
    enddo

    n = NEIG-1
    do i = 1, n
      minloc = i
      EMIN = dabs(EIG(NEW(I)))
      IP = I+1
      do J=IP,NEIG
        if(dabs(EIG(NEW(J))).LT.EMIN)then
          minloc = J
          EMIN = dabs(EIG(NEW(J)))
        endif
      enddo
      IBAF = NEW(I)
      NEW(I) = NEW(minloc)
      NEW(minloc) = IBAF
    enddo
  end subroutine EVSORT

  subroutine URAND1(N, X, SHIFT)
    use hecmw
    implicit none
    real(kind=kreal) :: X(N), INVM
    integer(kind=kint), parameter :: MM = 1664501
    integer(kind=kint), parameter :: LAMBDA = 1229
    integer(kind=kint), parameter :: MU = 351750
    integer(kind=kint) :: i, N, IR, SHIFT

    IR = 0
    INVM = 1.0D0 / MM
    do I = 1, SHIFT
      IR = mod( LAMBDA * IR + MU, MM)
    enddo
    do I = SHIFT+1, SHIFT+N
      IR = mod( LAMBDA * IR + MU, MM)
      X(I-SHIFT) = INVM * IR
    enddo
  end subroutine URAND1

end module m_fstr_EIG_lanczos_util

