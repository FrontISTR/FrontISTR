!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Eigen Analysis                                    !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module provides a subroutine to find the eigenvalues and eigenvectors
!! of a symmetric tridiagonal matrix by the ql method.
module m_fstr_EIG_tridiag
contains

      function a2b2(a,b)
      use hecmw
      real(kind=kreal) a,b
      real(kind=kreal) a2b2
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
      real(kind=kreal) p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 a2b2 = p
      return
      end function a2b2

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
      subroutine TRIDIAG(nm,n,d,du,e,z,zu,ierror)
!
      use hecmw
      integer(kind=kint) i,j,k,l,m,n,ii,l1,l2,nm,mml,ierror
      real(kind=kreal) d(n),du(n),e(n),z(nm,n),zu(nm,n)
      real(kind=kreal) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2

      ierror = 0
      if (n .eq. 1) go to 1001
!
      do 100 i = 2, n
        e(i-1) = e(i)
  100 continue
!
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
!
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
  110    continue
!
  120    if (m .eq. l) go to 220
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
!
         do 140 i = l2, n
           d(i) = d(i) - h
  140    continue
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
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
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
!
  200    continue
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
!     .......... order eigenvalues and eigenvectors ..........
!GP: Get unordered eigenvalues and eigenvectors----------------
      do i = 1,n
       du(i) = d(i)
       do j = 1,n
        zu(j,i) = z(j,i)
       end do
      end do
!--------------------------------------------------------------
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
!
  300 continue
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierror = l
 1001 return
      end subroutine TRIDIAG


end module m_fstr_EIG_tridiag

