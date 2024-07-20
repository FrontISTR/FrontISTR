!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides aux functions
module m_utilities
  use hecmw
  implicit none

  real(kind=kreal), parameter, private :: PI=3.14159265358979d0

contains

  !> Record used memory
  subroutine memget(var,dimn,syze)
    integer :: var,dimn,syze,bite
    parameter(bite=1)
    var = var + dimn*syze*bite
  end subroutine memget


  !> Insert an integer at end of a file name
  subroutine append_int2name( n, fname, n1 )
    integer, intent(in)             :: n
    integer, intent(in), optional  ::  n1
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
    if(present(n1).and.n1/=0)then
      write(tmp,'(i8)')n1
      fname = fname(1:len_trim(fname))//'.'//adjustl(tmp)
    endif
  end subroutine

  !> Insert an integer into a integer array
  subroutine insert_int2array( iin, carray )
    integer, intent(in) :: iin
    integer, pointer :: carray(:)

    integer :: i, oldsize
    integer, pointer :: dumarray(:) => null()
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
  subroutine tensor_eigen3( tensor, eigval )
    real(kind=kreal), intent(in)  :: tensor(6)          !< tensor
    real(kind=kreal), intent(out) :: eigval(3)     !< eigenvalues

    real(kind=kreal) :: I1,I2,I3,R,sita,Q, X(3,3), XX(3,3), II(3,3)

    II(:,:)=0.d0
    II(1,1)=1.d0;  II(2,2)=1.d0;  II(3,3)=1.d0
    X(1,1)=tensor(1); X(2,2)=tensor(2); X(3,3)=tensor(3)
    X(1,2)=tensor(4); X(2,1)=X(1,2)
    X(2,3)=tensor(5); X(3,2)=X(2,3)
    X(3,1)=tensor(6); X(1,3)=X(3,1)

    XX= matmul( X,X )
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

  end subroutine

  !> Compute eigenvalue and eigenvetor for symmetric 3*3 tensor using
  !> Jacobi iteration adapted from numerical recipes
  subroutine eigen3 (tensor, eigval, princ)
    real(kind=kreal), intent(in)  :: tensor(6)     !< tensor
    real(kind=kreal), intent(out) :: eigval(3)     !< vector containing the eigvalches
    real(kind=kreal), intent(out) :: princ(3, 3)   !< matrix containing the three principal column vectors

    integer, parameter :: msweep = 50
    integer :: i,j, is, ip, iq, ir
    real(kind=kreal) :: fsum, od, theta, t, c, s, tau, g, h, hd, btens(3,3), factor

    btens(1,1)=tensor(1); btens(2,2)=tensor(2); btens(3,3)=tensor(3)
    btens(1,2)=tensor(4); btens(2,1)=btens(1,2)
    btens(2,3)=tensor(5); btens(3,2)=btens(2,3)
    btens(3,1)=tensor(6); btens(1,3)=btens(3,1)
    !
    !     Initialise princ to the identity
    !
    factor = 0.d0
    do i = 1, 3
      do j = 1, 3
        princ (i, j) = 0.d0
      end do
      princ (i, i) = 1.d0
      eigval (i) = btens (i, i)
      factor = factor + dabs(btens(i,i))
    end do
    !     Scaling and iszero/isnan exception
    if( factor == 0.d0 .or. factor /= factor ) then
      return
    else
      eigval(1:3) = eigval(1:3)/factor
      btens(1:3,1:3) = btens(1:3,1:3)/factor
    end if

    !
    !     Starts sweeping.
    !
    do is = 1, msweep
      fsum = 0.d0
      do ip = 1, 2
        do iq = ip + 1, 3
          fsum = fsum + abs( btens(ip, iq) )
        end do
      end do
      !
      !     If the fsum of off-diagonal terms is zero returns
      !
      if ( fsum < 1.d-10 ) then
        eigval(1:3) = eigval(1:3)*factor
        return
      endif
      !
      !     Performs the sweep in three rotations. One per off diagonal term
      !
      do ip = 1, 2
        do iq = ip + 1, 3
          od = 100.d0 * abs (btens (ip, iq) )
          if ( (od+abs (eigval (ip) ) /= abs (eigval (ip) )) &
              .and. (od+abs (eigval (iq) ) /= abs (eigval (iq) ))) then
            hd = eigval (iq) - eigval (ip)
            !
            !    Evaluates the rotation angle
            !
            if ( abs (hd) + od == abs (hd)  ) then
              t = btens (ip, iq) / hd
            else
              theta = 0.5d0 * hd / btens (ip, iq)
              t = 1.d0 / (abs (theta) + sqrt (1.d0 + theta**2) )
              if ( theta < 0.d0 ) t = - t
            end if
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
            do ir = 1, 3
              g = princ (ir, ip)
              h = princ (ir, iq)
              princ (ir, ip) = g - s * (h + g * tau)
              princ (ir, iq) = h + s * (g - h * tau)
            end do
          end if
          btens (ip, iq) = 0.d0
        end do
      end do
    end do ! over sweeps
    !
    !     If convergence is not achieved stops
    !
    stop       ' Jacobi iteration unable to converge'
  end subroutine eigen3

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
  end function Determinant

  !> Compute determinant for 3*3 matrix
  real(kind=kreal) function Determinant33( XJ )
    real(kind=kreal) :: XJ(3,3)

    Determinant33=XJ(1,1)*XJ(2,2)*XJ(3,3)               &
      +XJ(2,1)*XJ(3,2)*XJ(1,3)                   &
      +XJ(3,1)*XJ(1,2)*XJ(2,3)                   &
      -XJ(3,1)*XJ(2,2)*XJ(1,3)                   &
      -XJ(2,1)*XJ(1,2)*XJ(3,3)                   &
      -XJ(1,1)*XJ(3,2)*XJ(2,3)
  end function Determinant33

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
  end subroutine fstr_chk_alloc

  !> calculate inverse of matrix a
  subroutine calInverse(NN, A)
    integer, intent(in)             :: NN
    real(kind=kreal), intent(inout) :: A(NN,NN)

    integer          :: I, J,K,IW,LR,IP(NN)
    real(kind=kreal) :: W,WMAX,PIVOT,API,EPS,DET
    data EPS/1.0E-35/
    DET=1.d0
    LR = 0.0d0
    do I=1,NN
      IP(I)=I
    enddo
    do K=1,NN
      WMAX=0.d0
      do I=K,NN
        W=dabs(A(I,K))
        if (W.GT.WMAX) then
          WMAX=W
          LR=I
        endif
      enddo
      PIVOT=A(LR,K)
      API=abs(PIVOT)
      if(API.LE.EPS) then
        write(*,'(''PIVOT ERROR AT'',I5)') K
        stop
      end if
      DET=DET*PIVOT
      if (LR.NE.K) then
        DET=-DET
        IW=IP(K)
        IP(K)=IP(LR)
        IP(LR)=IW
        do J=1,NN
          W=A(K,J)
          A(K,J)=A(LR,J)
          A(LR,J)=W
        enddo
      endif
      do I=1,NN
        A(K,I)=A(K,I)/PIVOT
      enddo
      do I=1,NN
        if (I.NE.K) then
          W=A(I,K)
          if (W.NE.0.) then
            do J=1,NN
              if (J.NE.K) A(I,J)=A(I,J)-W*A(K,J)
            enddo
            A(I,K)=-W/PIVOT
          endif
        endif
      enddo
      A(K,K)=1.d0/PIVOT
    enddo

    do I=1,NN
      K=IP(I)
      if (K.NE.I) then
        IW=IP(K)
        IP(K)=IP(I)
        IP(I)=IW
        do J=1,NN
          W=A(J,I)
          A(J,I)=A(J,K)
          A(J,K)=W
        enddo
      endif
    enddo

  end subroutine calInverse

  subroutine cross_product(v1,v2,vn)
    real(kind=kreal),intent(in)  ::  v1(3),v2(3)
    real(kind=kreal),intent(out)  ::  vn(3)

    vn(1) = v1(2)*v2(3) - v1(3)*v2(2)
    vn(2) = v1(3)*v2(1) - v1(1)*v2(3)
    vn(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end subroutine cross_product

  subroutine transformation(jacob, tm)
    real(kind=kreal),intent(in)  ::  jacob(3,3)   !< Jacobian
    real(kind=kreal),intent(out)  ::  tm(6,6)      !< transform matrix

    integer    ::  i,j

    do i=1,3
      do j=1,3
        tm(i,j)= jacob(i,j)*jacob(i,j)
      enddo
      tm(i,4) = jacob(i,1)*jacob(i,2)
      tm(i,5) = jacob(i,2)*jacob(i,3)
      tm(i,6) = jacob(i,3)*jacob(i,1)
    enddo
    tm(4,1) = 2.d0*jacob(1,1)*jacob(2,1)
    tm(5,1) = 2.d0*jacob(2,1)*jacob(3,1)
    tm(6,1) = 2.d0*jacob(3,1)*jacob(1,1)
    tm(4,2) = 2.d0*jacob(1,2)*jacob(2,2)
    tm(5,2) = 2.d0*jacob(2,2)*jacob(3,2)
    tm(6,2) = 2.d0*jacob(3,2)*jacob(1,2)
    tm(4,3) = 2.d0*jacob(1,3)*jacob(2,3)
    tm(5,3) = 2.d0*jacob(2,3)*jacob(3,3)
    tm(6,3) = 2.d0*jacob(3,3)*jacob(1,3)
    tm(4,4) = jacob(1,1)*jacob(2,2) + jacob(1,2)*jacob(2,1)
    tm(5,4) = jacob(2,1)*jacob(3,2) + jacob(2,2)*jacob(3,1)
    tm(6,4) = jacob(3,1)*jacob(1,2) + jacob(3,2)*jacob(1,1)
    tm(4,5) = jacob(1,2)*jacob(2,3) + jacob(1,3)*jacob(2,2)
    tm(5,5) = jacob(2,2)*jacob(3,3) + jacob(2,3)*jacob(3,2)
    tm(6,5) = jacob(3,2)*jacob(1,3) + jacob(3,3)*jacob(1,2)
    tm(4,6) = jacob(1,3)*jacob(2,1) + jacob(1,1)*jacob(2,3)
    tm(5,6) = jacob(2,3)*jacob(3,1) + jacob(2,1)*jacob(3,3)
    tm(6,6) = jacob(3,3)*jacob(1,1) + jacob(3,1)*jacob(1,3)

  end subroutine transformation

  subroutine get_principal (tensor, eigval, princmatrix)

    implicit none
    integer i,j
    real(kind=kreal) :: tensor(1:6)
    real(kind=kreal) :: eigval(3)
    real(kind=kreal) :: princmatrix(3,3)
    real(kind=kreal) :: princnormal(3,3)
    real(kind=kreal) :: tempv(3)
    real(kind=kreal) :: temps

    call eigen3(tensor,eigval,princnormal)

    if (eigval(1)<eigval(2)) then
      temps=eigval(1)
      eigval(1)=eigval(2)
      eigval(2)=temps
      tempv(:)=princnormal(:,1)
      princnormal(:,1)=princnormal(:,2)
      princnormal(:,2)=tempv(:)
    end if
    if (eigval(1)<eigval(3)) then
      temps=eigval(1)
      eigval(1)=eigval(3)
      eigval(3)=temps
      tempv(:)=princnormal(:,1)
      princnormal(:,1)=princnormal(:,3)
      princnormal(:,3)=tempv(:)
    end if
    if (eigval(2)<eigval(3)) then
      temps=eigval(2)
      eigval(2)=eigval(3)
      eigval(3)=temps
      tempv(:)=princnormal(:,2)
      princnormal(:,2)=princnormal(:,3)
      princnormal(:,3)=tempv(:)
    end if

    do j=1,3
      do i=1,3
        princmatrix(i,j) = princnormal(i,j) * eigval(j)
      end do
    end do

  end subroutine get_principal

  subroutine eigen3d (tensor, eigval, princ)
    implicit none

    real(kind=kreal) :: tensor(6)     !< tensor
    real(kind=kreal) :: eigval(3)     !< vector containing the eigvalches
    real(kind=kreal) :: princ(3,3)   !< matrix containing the three principal column vectors

    real(kind=kreal) :: s11, s22, s33, s12, s23, s13, j1, j2, j3
    real(kind=kreal) :: ml,nl
    complex(kind=kreal):: x1,x2,x3
    real(kind=kreal):: rtemp
    integer :: i
    s11 = tensor(1)
    s22 = tensor(2)
    s33 = tensor(3)
    s12 = tensor(4)
    s23 = tensor(5)
    s13 = tensor(6)

    ! invariants of stress tensor
    j1 = s11 + s22 + s33
    j2 = -s11*s22 - s22*s33 - s33*s11 + s12**2 + s23**2 + s13**2
    j3 = s11*s22*s33 + 2*s12*s23*s13 - s11*s23**2 - s22*s13**2 - s33*s12**2
    ! Cardano's method
    ! x^3+ ax^2   + bx  +c =0
    ! s^3 - J1*s^2 -J2s -J3 =0
    call cardano(-j1, -j2, -j3, x1, x2, x3)
    eigval(1)= real(x1)
    eigval(2)= real(x2)
    eigval(3)= real(x3)
    if (eigval(1)<eigval(2)) then
      rtemp=eigval(1)
      eigval(1)=eigval(2)
      eigval(2)=rtemp
    end if
    if (eigval(1)<eigval(3)) then
      rtemp=eigval(1)
      eigval(1)=eigval(3)
      eigval(3)=rtemp
    end if
    if (eigval(2)<eigval(3)) then
      rtemp=eigval(2)
      eigval(2)=eigval(3)
      eigval(3)=rtemp
    end if

    do i=1,3
      if (eigval(i)/(eigval(1)+eigval(2)+eigval(3)) < 1.0d-10 )then
        eigval(i) = 0.0d0
        princ(i,:) = 0.0d0
        exit
      end if
      ml = ( s23*s13 - s12*(s33-eigval(i)) ) / ( -s23**2 + (s22-eigval(i))*(s33-eigval(i)) )
      nl = ( s12**2 - (s22-eigval(i))*(s11-eigval(i)) ) / ( s12*s23 - s13*(s22-eigval(i)) )
      if (abs(ml) >= huge(ml)) then
        ml=0.0d0
      end if
      if (abs(nl) >= huge(nl)) then
        nl=0.0d0
      end if
      princ(i,1) = eigval(i)/sqrt( 1 + ml**2 + nl**2)
      princ(i,2) = ml * princ(i,1)
      princ(i,3) = nl * princ(i,1)
    end do

    write(*,*)
  end subroutine eigen3d

  subroutine cardano(a,b,c,x1,x2,x3)
    real(kind=kreal):: a,b,c
    real(kind=kreal):: p,q,d
    complex(kind=kreal):: w
    complex(kind=kreal):: u,v,y
    complex(kind=kreal):: x1,x2,x3
    w = (-1.0d0 + sqrt(dcmplx(-3.0d0)))/2.0d0
    p = -a**2/9.0d0 + b/3.0d0
    q = 2.0d0/2.7d1*a**3 - a*b/3.0d0 + c
    d = q**2 + 4.0d0*p**3

    u = ((-dcmplx(q) + sqrt(dcmplx(d)))/2.0d0)**(1.0d0/3.0d0)

    if(u.ne.0.0d0) then
      v = -dcmplx(p)/u
      x1 = u + v -dcmplx(a)/3.0d0
      x2 = u*w + v*w**2 -dcmplx(a)/3.0d0
      x3 = u*w**2 + v*w -dcmplx(a)/3.0d0
    else
      y = (-dcmplx(q))**(1.0d0/3.0d0)
      x1 = y -dcmplx(a)/3.0d0
      x2 = y*w -dcmplx(a)/3.0d0
      x3 = y*w**2 -dcmplx(a)/3.0d0
    end if

  end subroutine cardano

  subroutine rotate_3dvector_by_Rodrigues_formula(r,v)
    real(kind=kreal),intent(in)    :: r(3)  !< rotational vector
    real(kind=kreal),intent(inout) :: v(3)  !< vector to be rotated

    real(kind=kreal) :: rotv(3), rv
    real(kind=kreal) :: cosx, sinc(2) !< sinc(1)=sin(x)/x, sinc(2)=(1-cos(x))/x^2
    real(kind=kreal) :: x, x2, x4, x6
    real(kind=kreal), parameter :: c0 = 0.5d0
    real(kind=kreal), parameter :: c2 = -4.166666666666666d-002
    real(kind=kreal), parameter :: c4 = 1.388888888888889d-003
    real(kind=kreal), parameter :: c6 = -2.480158730158730d-005

    x2 = dot_product(r,r)
    x = dsqrt(x2)
    cosx = dcos(x)
    if( x < 1.d-4 ) then
      x4 = x2*x2
      x6 = x4*x2
      sinc(1) = 1.d0-x2/6.d0+x4/120.d0
      sinc(2) = c0+c2*x2+c4*x4+c6*x6
    else
      sinc(1) = dsin(x)/x
      sinc(2) = (1.d0-cosx)/x2
    endif

    ! calc Rot*v
    rv = dot_product(r,v)
    rotv(1:3) = cosx*v(1:3)
    rotv(1:3) = rotv(1:3)+rv*sinc(2)*r(1:3)
    rotv(1) = rotv(1) + (-v(2)*r(3)+v(3)*r(2))*sinc(1)
    rotv(2) = rotv(2) + (-v(3)*r(1)+v(1)*r(3))*sinc(1)
    rotv(3) = rotv(3) + (-v(1)*r(2)+v(2)*r(1))*sinc(1)
    v = rotv

  end subroutine

  !> Compute derivative of a general isotropic tensor function of one tensor
  subroutine deriv_general_iso_tensor_func_3d(dpydpx, dydx, eigv, px, py)
    real(kind=kreal), intent(in) :: dpydpx(3,3)
    real(kind=kreal), intent(out) :: dydx(6,6)
    real(kind=kreal), intent(in) :: eigv(3,3)
    real(kind=kreal), intent(in) :: px(3), py(3)

    real(kind=kreal) :: x(6)
    real(kind=kreal) :: Ep(6,3), dx2dx(6,6)
    integer(kind=kint) :: i, j, m1, m2, m3, ia, ib, ic, ip, jp, k
    real(kind=kreal) :: xmax, dif12, dif23, C1, C2, C3, C4, D1, D2, D3
    real(kind=kreal) :: xa, xc, ya, yc, daa, dac, dca, dcb, dcc
    real(kind=kreal) :: xa_xc, xa_xc2, xa_xc3, ya_yc, xc2
    real(kind=kreal) :: s1, s2, s3, s4, s5, s6
    real(kind=kreal), parameter :: I2(6) = [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    real(kind=kreal), parameter :: Is(6,6) = reshape([&
        1.d0, 0.d0, 0.d0,  0.d0,   0.d0,   0.d0,&
        0.d0, 1.d0, 0.d0,  0.d0,   0.d0,   0.d0,&
        0.d0, 0.d0, 1.d0,  0.d0,   0.d0,   0.d0,&
        0.d0, 0.d0, 0.d0, 0.5d0,   0.d0,   0.d0,&
        0.d0, 0.d0, 0.d0,  0.d0,  0.5d0,   0.d0,&
        0.d0, 0.d0, 0.d0,  0.d0,   0.d0,  0.5d0],&
        [6, 6])
    real(kind=kreal), parameter :: EPS=1.d-6

    do i=1,3
      Ep(1,i) = eigv(1,i)*eigv(1,i)
      Ep(2,i) = eigv(2,i)*eigv(2,i)
      Ep(3,i) = eigv(3,i)*eigv(3,i)
      Ep(4,i) = eigv(1,i)*eigv(2,i)
      Ep(5,i) = eigv(2,i)*eigv(3,i)
      Ep(6,i) = eigv(3,i)*eigv(1,i)
    enddo
    x(:) = 0.d0
    do i=1,3
      do j=1,6
        x(j) = x(j)+px(i)*Ep(j,i)
      enddo
    enddo
    dx2dx = reshape([2.d0*x(1),      0.d0,      0.d0,             x(4),             0.d0,             x(6),&
        &                 0.d0, 2.d0*x(2),      0.d0,             x(4),             x(5),             0.d0,&
        &                 0.d0,      0.d0, 2.d0*x(3),             0.d0,             x(5),             x(6),&
        &                 x(4),      x(4),      0.d0, (x(1)+x(2))/2.d0,        x(6)/2.d0,        x(5)/2.d0,&
        &                 0.d0,      x(5),      x(6),        x(6)/2.d0, (x(2)+x(3))/2.d0,        x(4)/2.d0,&
        &                 x(6),      0.d0,      x(6),        x(5)/2.d0,        x(4)/2.d0, (x(1)+x(3))/2.d0],&
        [6, 6])
    m1 = maxloc(px, 1)
    m3 = minloc(px, 1)
    if( m1 == m3 ) then
      m1 = 1; m2 = 2; m3 = 3
    else
      m2 = 6 - (m1 + m3)
    endif
    xmax = maxval(abs(px), 1)
    if( xmax < EPS ) xmax = 1.d0
    dif12 = abs(px(m1) - px(m2)) / xmax
    dif23 = abs(px(m2) - px(m3)) / xmax
    if( dif12 < EPS .and. dif23 < EPS ) then
      C1 = dpydpx(1,1)-dpydpx(1,2)
      C2 = dpydpx(1,2)
      do j=1,6
        do i=1,6
          dydx(i,j) = C1*Is(i,j)+C2*I2(i)*I2(j)
        enddo
      enddo
    else if( dif12 < EPS .or. dif23 < EPS ) then
      if( dif12 < EPS ) then
        ia = m3; ib = m1; ic = m2
      else
        ia = m1; ib = m2; ic = m3
      endif
      ya = py(ia); yc = py(ic)
      xa = px(ia); xc = px(ic)
      daa = dpydpx(ia,ia)
      dac = dpydpx(ia,ic)
      dca = dpydpx(ic,ia)
      dcb = dpydpx(ic,ib)
      dcc = dpydpx(ic,ic)
      xa_xc = xa-xc
      xa_xc2 = xa_xc*xa_xc
      xa_xc3 = xa_xc2*xa_xc
      ya_yc = ya-yc
      xc2 = xc*xc
      C1 = ya_yc/xa_xc2
      C2 = ya_yc/xa_xc3
      C3 = xc/xa_xc2
      C4 = (xa+xc)/xa_xc
      D1 = dac+dca
      D2 = daa+dcc
      D3 = D1-D2
      s1 = C1+(dcb-dcc)/xa_xc
      s2 = 2.d0*xc*C1+(dcb-dcc)*C4
      s3 = 2.d0*C2+D3/xa_xc2
      s4 = 2.d0*xc*C2+(dac-dcb)/xa_xc+D3*C3
      s5 = 2.d0*xc*C2+(dca-dcb)/xa_xc+D3*C3
      s6 = 2.d0*xc2*C2+(D1*xa-D2*xc)*C3-dcb*C4
      do j=1,6
        do i=1,6
          dydx(i,j) = s1*dx2dx(i,j)-s2*Is(i,j)-s3*x(i)*x(j)+s4*x(i)*I2(j)+s5*I2(i)*x(j)-s6*I2(i)*I2(j)
        enddo
      enddo
    else
      do k=1,3
        select case(k)
        case(1)
          ia = 1; ib = 2; ic = 3
        case(2)
          ia = 2; ib = 3; ic = 1
        case(3)
          ia = 3; ib = 1; ic = 2
        end select
        C1 = py(ia)/((px(ia)-px(ib))*(px(ia)-px(ic)))
        C2 = px(ib)+px(ic)
        C3 = (px(ia)-px(ib))+(px(ia)-px(ic))
        C4 = px(ib)-px(ic)
        do j=1,6
          do i=1,6
            dydx(i,j) = dydx(i,j)+C1*(dx2dx(i,j)-C2*Is(i,j)-C3*Ep(i,ia)*Ep(j,ia)-C4*(Ep(i,ib)*Ep(j,ib)-Ep(i,ic)*Ep(j,ic)))
          enddo
        enddo
      enddo
      do jp=1,3
        do ip=1,3
          C1 = dpydpx(ip,jp)
          do j=1,6
            do i=1,6
              dydx(i,j) = dydx(i,j)+C1*Ep(i,ip)*Ep(j,jp)
            enddo
          enddo
        enddo
      enddo
    endif
  end subroutine deriv_general_iso_tensor_func_3d

end module
