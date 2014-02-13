!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides aux functions
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/11/16
!>  \version    0.00
!
!======================================================================!
module m_utilities
  implicit none
  
  integer, parameter, private :: kreal = kind(0.0d0)
  real(kind=kreal), parameter, private :: PI=3.14159265358979
  
  contains

  !> Record used memeory
  SUBROUTINE memget(var,dimn,syze)
      INTEGER :: var,dimn,syze,bite
      PARAMETER(bite=1)
      var = var + dimn*syze*bite
  END SUBROUTINE memget

	  
  !> Insert an integer at end of a file name
  subroutine append_int2name( n, fname )
      integer, intent(in)             :: n
      character(len=*), intent(inout) :: fname  
      integer            :: npos, nlen
      character(len=128) :: tmpname, tmp
	  
      npos = scan( fname, '.')
      nlen = len_trim( fname )
      if( nlen>128 ) stop "String too long(>128) in append_int2name"
      if( n>100000 ) stop "Integer too big>100000 in append_int2name"
      tmpname = fname
      if( npos==0 ) then
        write( fname, '(a,i6)') fname(1:nlen),n
      else
        write( tmp, '(i6,a)') n,tmpname(npos:nlen)
        fname = tmpname(1:npos-1) // adjustl(tmp)
      endif
  end subroutine
  
  !> Insert an integer into a integer array
  subroutine insert_int2array( iin, carray )
      integer, intent(in) :: iin  
      integer, pointer :: carray(:)
	  
      integer :: i, oldsize
      integer, pointer :: dumarray(:)
      if( .not. associated(carray) ) then
        allocate( carray(1) )
        carray(1) = iin
      else
        oldsize = size( carray )
        allocate( dumarray(oldsize) )
        do i=1,oldsize
          dumarray(i) = carray(i)
        enddo
        deallocate( carray )
        allocate( carray(oldsize+1) )
        do i=1,oldsize
          carray(i) = dumarray(i)
        enddo
        carray(oldsize+1) = iin
      endif
      if( associated(dumarray) ) deallocate( dumarray )
  end subroutine
   
  !> Given symmetric 3x3 matrix M, compute the eigenvalues
  SUBROUTINE tensor_eigen3( tensor, eigval, eigproj )
    REAL(kind=kreal), INTENT(IN)  :: tensor(6)          !< tensor
    REAL(kind=kreal), INTENT(OUT) :: eigval(3)     !< eigenvalues
    REAL(kind=kreal), INTENT(OUT) :: eigproj(3,3)  !< eigenprojectss
	
    INTEGER  :: i
    REAL(kind=kreal) :: I1,I2,I3,R,sita,Q, X(3,3), XX(3,3), II(3,3)
	
    II(:,:)=0.d0
    II(1,1)=1.d0;  II(2,2)=1.d0;  II(3,3)=1.d0
    X(1,1)=tensor(1); X(2,2)=tensor(2); X(3,3)=tensor(3)
    X(1,2)=tensor(4); X(2,1)=X(1,2)
    X(2,3)=tensor(5); X(3,2)=X(2,3)
    X(3,1)=tensor(6); X(1,3)=X(3,1)

    XX= MATMUL( X,X )
    I1= X(1,1)+X(2,2)+X(3,3)
    I2= 0.5d0*( I1*I1 - (XX(1,1)+XX(2,2)+XX(3,3)) )
    I3= X(1,1)*X(2,2)*X(3,3)+X(2,1)*X(3,2)*X(1,3)+X(3,1)*X(1,2)*X(2,3)    &
     -X(3,1)*X(2,2)*X(1,3)-X(2,1)*X(1,2)*X(3,3)-X(1,1)*X(3,2)*X(2,3)

    R=(-2.d0*I1*I1*I1+9.d0*I1*I2-27.d0*I3)/54.d0
    Q=(I1*I1-3.d0*I2)/9.d0
    sita = acos(R/dsqrt(Q*Q*Q))
 
    eigval(1) = -2.d0*Q*cos(sita/3.d0)+I1/3.d0
    eigval(2) = -2.d0*Q*cos((sita+2.d0*PI)/3.d0)+I1/3.d0
    eigval(3) = -2.d0*Q*cos((sita-2.d0*PI)/3.d0)+I1/3.d0

  END SUBROUTINE
  
  !> Compute eigenvalue and eigenvetor for symmetric 3*3 tensor using
  !> Jacobi iteration adapted from numerical recpies
  SUBROUTINE eigen3 (tensor, eigval, princ)
    real(kind=kreal) :: tensor(6)     !< tensor
    real(kind=kreal) :: eigval(3)     !< vector containing the eigvalches
    real(kind=kreal) :: princ(3, 3)   !< matrix containing the three principal column vectors

    INTEGER, PARAMETER :: msweep = 50
    INTEGER :: i,j, is, ip, iq, ir
    real(kind=kreal) :: fsum, od, theta, t, c, s, tau, g, h, hd, btens(3,3)
	
	btens(1,1)=tensor(1); btens(2,2)=tensor(2); btens(3,3)=tensor(3)
    btens(1,2)=tensor(4); btens(2,1)=btens(1,2)
    btens(2,3)=tensor(5); btens(3,2)=btens(2,3)
    btens(3,1)=tensor(6); btens(1,3)=btens(3,1)
!
!     Initialise princ to the identity
!
      DO i = 1, 3
        DO j = 1, 3
          princ (i, j) = 0.d0
        END DO
        princ (i, i) = 1.d0
        eigval (i) = btens (i, i)
      END DO
!
!     Starts sweeping.
!
      DO is = 1, msweep
        fsum = 0.d0
        DO ip = 1, 2
          DO iq = ip + 1, 3
            fsum = fsum + abs( btens(ip, iq) )
          END DO
        END DO
!
!     If the fsum of off-diagonal terms is zero returns
!
        IF ( fsum < 1.d-10 ) RETURN

!
!     Performs the sweep in three rotations. One per off diagonal term
!
        DO ip = 1, 2
          DO iq = ip + 1, 3
            od = 100.d0 * abs (btens (ip, iq) )
            IF ( (od+abs (eigval (ip) ) /= abs (eigval (ip) )) &
                 .and. (od+abs (eigval (iq) ) /= abs (eigval (iq) ))) then
                hd = eigval (iq) - eigval (ip)
!
!    Evaluates the rotation angle
!
              IF ( abs (hd) + od == abs (hd)  ) then
                t = btens (ip, iq) / hd
              ELSE
                theta = 0.5d0 * hd / btens (ip, iq)
                t = 1.d0 / (abs (theta) + sqrt (1.d0 + theta**2) )
                IF ( theta < 0.d0 ) t = - t
              END IF
!
!     Re-evaluates the diagonal terms
!
              c = 1.d0 / sqrt (1.d0 + t**2)
              s = t * c
              tau = s / (1.d0 + c)
              h = t * btens (ip, iq)
              eigval (ip) = eigval (ip) - h
              eigval (iq) = eigval (iq) + h
!
!     Re-evaluates the remaining off-diagonal terms
!
              ir = 6 - ip - iq
              g = btens (min (ir, ip), max (ir, ip) )
              h = btens (min (ir, iq), max (ir, iq) )
              btens (min (ir, ip), max (ir, ip) ) = g &
                                                  - s * (h + g * tau)
              btens (min (ir, iq), max (ir, iq) ) = h &
                                                  + s * (g - h * tau)
!
!     Rotates the eigenvectors
!
              DO ir = 1, 3
                g = princ (ir, ip)
                h = princ (ir, iq)
                princ (ir, ip) = g - s * (h + g * tau)
                princ (ir, iq) = h + s * (g - h * tau)
              END DO
            END IF
            btens (ip, iq) = 0.d0
          END DO
        END DO
      END DO ! over sweeps
!
!     If convergence is not achieved stops
!
      STOP       ' Jacobi iteration unable to converge'
  END SUBROUTINE eigen3
  
  !> Compute determinant for symmetric 3*3 matrix
  real(kind=kreal) function Determinant( mat )
    real(kind=kreal) :: mat(6)     !< tensor
    real(kind=kreal) :: xj(3,3)
	
	xj(1,1)=mat(1); xj(2,2)=mat(2); xj(3,3)=mat(3)
    xj(1,2)=mat(4); xj(2,1)=xj(1,2)
    xj(2,3)=mat(5); xj(3,2)=xj(2,3)
    xj(3,1)=mat(6); xj(1,3)=xj(3,1)
	
    Determinant=XJ(1,1)*XJ(2,2)*XJ(3,3)               &
           +XJ(2,1)*XJ(3,2)*XJ(1,3)                   &
           +XJ(3,1)*XJ(1,2)*XJ(2,3)                   &
           -XJ(3,1)*XJ(2,2)*XJ(1,3)                   &
           -XJ(2,1)*XJ(1,2)*XJ(3,3)                   &
           -XJ(1,1)*XJ(3,2)*XJ(2,3)
  end function
  
  subroutine fstr_chk_alloc( imsg, sub_name, ierr )
        use hecmw
        character(*) :: sub_name
        integer(kind=kint) :: imsg
        integer(kind=kint) :: ierr

        if( ierr /= 0 ) then
                write(imsg,*) 'Memory overflow at ', sub_name
                write(*,*) 'Memory overflow at ', sub_name
                call hecmw_abort( hecmw_comm_get_comm( ) )
        endif
  end subroutine
		
end module m_utilities
