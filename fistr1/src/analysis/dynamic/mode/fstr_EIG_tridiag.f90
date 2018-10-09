!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine to find the eigenvalues and eigenvectors
!! of a symmetric tridiagonal matrix by the ql method.
module m_fstr_EIG_tridiag
  use hecmw

  implicit none

  public

  type fstr_eigen_vec
    real(kind=kreal), pointer :: q(:) => null()
  end type fstr_eigen_vec

  type fstr_tri_diag
    real(kind=kreal), allocatable :: alpha(:)
    real(kind=kreal), allocatable :: beta(:)
  end type fstr_tri_diag

contains

  subroutine tridiag(hecMESH, hecMAT, fstrEIG, Q, Tri, iter, is_converge)
    use hecmw
    use m_fstr
    use m_fstr_EIG_lanczos_util
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG
    type(fstr_tri_diag)      :: Tri
    type(fstr_eigen_vec), pointer :: Q(:)

    integer(kind=kint), intent(in) :: iter
    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: i, j, k, in, jn, kn, nget
    integer(kind=kint) :: iter2, ierr, maxiter
    real(kind=kreal)   :: resid, chk, sigma, tolerance, max_eigval
    real(kind=kreal), allocatable :: alpha(:), beta(:), temp(:)
    real(kind=kreal), allocatable :: L(:,:)

    integer(kind=kint), allocatable :: iparm(:)
    real(kind=kreal), pointer       :: eigvec(:,:)
    real(kind=kreal), pointer       :: eigval(:)
    logical :: is_converge

    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF
    eigval => fstrEIG%eigval
    eigvec => fstrEIG%eigvec
    nget      = fstrEIG%nget
    maxiter   = fstrEIG%maxiter
    tolerance = fstrEIG%tolerance

    allocate( iparm(maxiter) )
    allocate( temp(maxiter)  )
    allocate( alpha(iter)    )
    allocate( beta(iter)     )
    allocate( L(iter, iter)  )

    do j=1, iter
      do i = 1,iter
        L(i,j)  = 0.0d0
      enddo
    enddo

    do i=1, iter
      alpha(i) = Tri%alpha(i)
      L(i,i)   = 1.0d0
    enddo

    beta(1) = 0.0d0
    do i=2, iter
      beta(i) = Tri%beta(i)
    enddo

    call QL_decomposition(iter, iter, alpha, beta, L, ierr)

    is_converge = .true.
    chk = 0.0d0
    max_eigval = maxval(alpha)
    do i = 1, min(nget, iter)
      resid = dabs(Tri%beta(iter+1)*L(iter,i))/max_eigval
      chk = max(chk, resid)
      if(tolerance < resid) is_converge = .false.
    enddo
    if(myrank == 0) write(*,"(i8,1pe12.5)")iter, chk

    if(iter < nget) is_converge = .false.
    if(iter == maxiter-1) is_converge = .true.

    if(is_converge)then
      sigma = 0.0d0
      if(fstrEIG%is_free) sigma = 0.1d0
      do i = 1, iter
        if(alpha(i) /= 0.0d0)then
          eigval(i) = 1.0d0/alpha(i) + sigma
        endif
      enddo

      call evsort(eigval, iparm, iter)

      temp = eigval

      eigvec = 0.0d0
      do k=1, iter
        in = iparm(k)
        eigval(k) = temp(in)
        do j=1, iter
          do i=1, NPNDOF
            eigvec(i, k) = eigvec(i, k) + Q(j)%q(i) * L(j, in)
          enddo
        enddo
      enddo

      do j=1, iter
        chk = maxval(eigvec(:,j))
        call hecmw_allreduce_R1(hecMESH, chk, hecmw_max)
        if(chk /= 0.0d0)then
          chk = 1.0d0/chk
          do i = 1, NNDOF
            eigvec(i,j) = eigvec(i,j) * chk
          enddo
        endif
      enddo
    endif

    deallocate(iparm)
    deallocate(temp)
    deallocate(alpha)
    deallocate(beta)
    deallocate(L)
  end subroutine tridiag

  !======================================================================!
  !                       Description                                    !
  !======================================================================!
  !>This subroutine has been adapted from the eispack routine tql2, which is
  !!a translation of the algol procedure tql2, num. math. 11, 293-306(1968)
  !!by bowdler, martin, reinsch and wilkinson.
  !!Handbook for auto. cop., vol.ii-linear algebra, 227-240(1971).
  !!This subroutine finds the eigenvalues and eigenvectors
  !!of a symmetric tridiagonal matrix by the ql method.
  !!the eigenvectors of a full symmetric matrix can also
  !!be found if  tred2  has been used to reduce this
  !!full matrix to tridiagonal form.
  !
  !on input
  !
  !nm must be set to the row dimension of two-dimensional
  !array parameters as declared in the calling program
  !dimension statement.
  !
  !is the order of the matrix.
  !
  !contains the diagonal elements of the input matrix.
  !
  !contains the subdiagonal elements of the input matrix
  !in its last n-1 positions.  e(1) is arbitrary.
  !
  !contains the transformation matrix produced in the
  !reduction by  tred2, if performed.  if the eigenvectors
  !of the tridiagonal matrix are desired, z must contain
  !the identity matrix.
  !
  !on output
  !
  !d
  !contains the eigenvalues in ascending order.  if an
  !error exit is made, the eigenvalues are correct but
  !unordered for indices 1,2,...,ierror-1.
  !-----------------------------------------------------
  !GP
  !du
  !contains the unordered eigenvalues.  if an
  !error exit is made, the eigenvalues are correct and
  !unordered for indices 1,2,...,ierror-1.
  !-----------------------------------------------------
  !e
  !has been destroyed.
  !
  !z
  !contains orthonormal eigenvectors of the symmetric
  !tridiagonal (or full) matrix.  if an error exit is made,
  !z contains the eigenvectors associated with the stored
  !eigenvalues.
  !-----------------------------------------------------
  !zu
  !contains unordered eigenvectors of the symm. tridiag. matrix
  !-----------------------------------------------------
  !ierror is set to
  !  zero       for normal return,
  !  j          if the j-th eigenvalue has not been
  !             determined after 30 iterations.
  !
  !calls a2b2 for  dsqrt(a*a + b*b) .
  !=======================================================================

  !call QL_decomposition(iter, iter, alpha, beta, L, ierr)
  subroutine QL_decomposition(nm, n, d, e, z, ierror)
    use hecmw
    implicit none
    integer(kind=kint) :: i, j, k, l, m, n, ii, l1, l2, nm, mml, ierror
    real(kind=kreal)   :: d(n), e(n), z(nm, n)
    real(kind=kreal)   :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2

    ierror = 0
    if (n .eq. 1) go to 1001

    do i = 2, n
      e(i-1) = e(i)
    enddo

    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0

    do 240 l = 1, n
      j = 0
      h = dabs(d(l)) + dabs(e(l))
      if (tst1 .lt. h) tst1 = h
      !     .......... look for small sub-diagonal element ..........
      bb:do m = l, n
        tst2 = tst1 + dabs(e(m))
        if (tst2 .eq. tst1) exit bb
        !     .......... e(n) is always zero, so there is no exit
        !                through the bottom of the loop ..........
      enddo bb

      if (m .eq. l)  go to 220

      130    if (j .eq. 30) go to 1000
      j = j + 1
      !     .......... form shift ..........
      l1 = l + 1
      l2 = l1 + 1
      g = d(l)
      p = (d(l1) - g) / (2.0d0 * e(l))
      r = a2b2(p,1.0d0)
      d(l) = e(l) / (p + dsign(r,p))
      d(l1) = e(l) * (p + dsign(r,p))
      dl1 = d(l1)
      h = g - d(l)
      if (l2 .gt. n) go to 145

      do i = l2, n
        d(i) = d(i) - h
      enddo

      145    f = f + h
      !     .......... ql transformation ..........
      p = d(m)
      c = 1.0d0
      c2 = c
      el1 = e(l1)
      s = 0.0d0
      s2 = 0.0d0
      c3 = 0.0d0
      mml = m - l
      !     .......... for i=m-1 step -1 until l do -- ..........
      do ii = 1, mml
        c3 = c2
        c2 = c
        s2 = s
        i = m - ii
        g = c * e(i)
        h = c * p
        r = a2b2(p,e(i))
        e(i+1) = s * r
        s = e(i) / r
        c = p / r
        p = c * d(i) - s * g
        d(i+1) = h + s * (c * g + s * d(i))
        !     .......... form vector ..........
        do k = 1, n
          h = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * h
          z(k,i) = c * z(k,i) - s * h
        enddo
      enddo

      p = -s * s2 * c3 * el1 * e(l) / dl1
      e(l) = s * p
      d(l) = c * p
      tst2 = tst1 + dabs(e(l))
      if (tst2 .gt. tst1) go to 130
      220    d(l) = d(l) + f
      240 continue

      !     .......... order eigenvalues and eigenvectors ..........
      !GP: Get unordered eigenvalues and eigenvectors----------------

      do 300 ii = 2, n
        i = ii - 1
        k = i
        p = d(i)

        aa:do j = ii, n
          if (d(j) .ge. p) exit aa
          k = j
          p = d(j)
        enddo aa

        if (k .eq. i) go to 300
        d(k) = d(i)
        d(i) = p

        do j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
        enddo

        300 continue

        go to 1001
        !     .......... set error -- no convergence to an
        !                eigenvalue after 30 iterations ..........
        1000 ierror = l
        1001 return
  end subroutine QL_decomposition

  function a2b2(a,b)
    use hecmw
    implicit none

    real(kind=kreal) :: a2b2
    real(kind=kreal) :: a, b
    real(kind=kreal) :: p, q, r, s, t, u

    p = dmax1(dabs(a), dabs(b))
    if (p /= 0.0d0) then
      r = (dmin1(dabs(a),dabs(b))/p) ** 2
      do
        t = 4.0d0 + r
        if (t == 4.0d0) exit
        s = r/t
        u = 1.0d0 + 2.0d0*s
        p = u * p
        q = s/u
        r = q * q * r
      end do
    end if
    a2b2 = p
    return

  end function a2b2

end module m_fstr_EIG_tridiag
