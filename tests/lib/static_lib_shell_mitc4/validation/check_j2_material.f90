program check_j2_material
  use m_ElastoPlastic
  use mMaterial
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  type(tMaterial) :: mat
  type(tTable) :: table
  type(tTable), pointer :: elastic_table
  logical :: missing
  real(kind=kreal) :: deps(6), epsout(6), stress(6), tangent(5,5), zero(6)
  real(kind=kreal) :: state(1), outstate(1), oldstress(6), oldplastic, potential, work
  real(kind=kreal) :: plus(6), minus(6), trial(6), dummy(5,5), fd(5,5), stress0(6)
  real(kind=kreal) :: g, tau, scale, err, h, temperature_table(3,2)
  integer(kind=kint) :: ierr
  integer :: istat, j, stage, component
  integer, parameter :: components(5) = (/1,2,4,5,6/)

  call initMaterial(mat)
  mat%mtype = 120000
  mat%variables(M_YOUNGS) = 210000.0d0
  mat%variables(M_POISSON) = 0.3d0
  mat%variables(M_PLCONST1) = 1000.0d0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_BAD_PROPERTY) stop 2
  call init_table(table, 0, 2, 1, reshape((/210000.0d0,0.3d0/),(/2,1/)))
  call dict_add_key(mat%dict, MC_ISOELASTIC, table)
  call fetch_Table(MC_ISOELASTIC, mat%dict, elastic_table, missing)
  oldstress = 0.0d0
  zero = 0.0d0
  oldplastic = 0.0d0
  state = 0.0d0
  g = 210000.0d0/(2.0d0*1.3d0)

  ! Prescribed uniaxial stress in the elastic range: eps22=eps33=-nu*eps11.
  deps = (/1.0d-4,-3.0d-5,0.0d0,0.0d0,0.0d0,0.0d0/)
  call response(deps, stress, tangent)
  call check('elastic thickness strain', epsout(3), -3.0d-5, 1.0d-12)
  call check('uniaxial stress', stress(1), 21.0d0, 1.0d-10)
  call check('lateral stress', stress(2), 0.0d0, 1.0d-10)

  do component = 4, 6
    scale = 1.0d0
    if(component > 4) scale = sqrt(5.0d0/6.0d0)
    deps = 0.0d0
    deps(component) = 1.0d-4
    call response(deps, stress, tangent)
    call check('elastic shear stiffness', stress(component)/deps(component), scale**2*g, 1.0d-8)
    deps(component) = 0.03d0
    call response(deps, stress, tangent)
    tau = scale*1000.0d0/sqrt(3.0d0)
    call check('scaled shear yield', stress(component), tau, 1.0d-8)
    call check('shear work', work, tau*deps(component)-tau*tau/(2.0d0*scale**2*g), 1.0d-8)
  enddo

  h = 1.0d-8
  do stage = 1, 3
    if(stage == 1) deps = (/1.0d-4, -2.0d-5, 0.0d0, 3.0d-5, 2.0d-5, -1.0d-5/)
    if(stage == 2) deps = (/0.012d0, 0.001d0, 0.0d0, 0.002d0, 0.001d0, -0.001d0/)
    if(stage == 3) then
      oldstress = stress0
      oldplastic = outstate(1)
      state(1) = oldplastic
      deps = -0.02d0*deps
    endif
    call response(deps, stress, tangent)
    if(stage == 2) then
      call check('plastic incompressibility', sum(epsout(1:3)), &
        (1.0d0-0.6d0)/210000.0d0*sum(stress(1:3)), 1.0d-11)
    endif
    do j = 1, 5
      trial = deps
      trial(components(j)) = trial(components(j))+h
      call response(trial, plus, dummy)
      trial(components(j)) = trial(components(j))-2.0d0*h
      call response(trial, minus, dummy)
      fd(:,j) = (plus(components)-minus(components))/(2.0d0*h)
    enddo
    err = sqrt(sum((fd-tangent)**2)/sum(tangent**2))
    print *, 'FD relative error, stage', stage, err
    if(err > 1.0d-6) stop 1
    call response(deps, stress0, dummy)
  enddo
  call check('elastic unloading preserves plastic strain', outstate(1), oldplastic, 1.0d-12)
  call check_stress_units()

  elastic_table%tbval(1,1) = -1.0d0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_BAD_PROPERTY) stop 2
  elastic_table%tbval(1,1) = 210000.0d0
  elastic_table%tbval(2,1) = 0.5d0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_BAD_PROPERTY) stop 2
  elastic_table%tbval(2,1) = 0.3d0
  mat%variables(M_PLCONST1) = 0.0d0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_BAD_PROPERTY) stop 2
  mat%variables(M_PLCONST1) = 1000.0d0
  mat%variables(M_PLCONST2) = 10.0d0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_UNSUPPORTED) stop 2
  mat%variables(M_PLCONST2) = 0.0d0
  mat%mtype = 121000
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_UNSUPPORTED) stop 2
  mat%mtype = 120000
  allocate(mat%shell_var(2))
  mat%shell_var%ortho = 0
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= 0) stop 2
  mat%shell_var(2)%ortho = 1
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_UNSUPPORTED) stop 2
  deallocate(mat%shell_var)
  call Update_PlaneStressJ2(mat, zero, zero, 0.0d0, state(1:0), stress, epsout, tangent, &
    istat, outstate, potential, work, ierr)
  if(ierr /= J2_BAD_STATE) stop 2
  deps = zero
  deps(1) = ieee_value(1.0d0, ieee_quiet_nan)
  call Update_PlaneStressJ2(mat, deps, zero, 0.0d0, state, stress, epsout, tangent, &
    istat, outstate, potential, work, ierr)
  if(ierr /= J2_NONFINITE) stop 2
  temperature_table(:,1) = (/210000.0d0,0.3d0,0.0d0/)
  temperature_table(:,2) = (/200000.0d0,0.3d0,100.0d0/)
  call init_table(table, 1, 3, 2, temperature_table)
  call dict_add_key(mat%dict, MC_ISOELASTIC, table)
  call checkPlaneStressJ2Material(mat, ierr)
  if(ierr /= J2_TEMPERATURE) stop 2
  print *, 'PASS thickness strain, transverse shear, work, tangents, stress units, input errors'
  call finalizeMaterial(mat)
contains
  subroutine check_stress_units()
    real(kind=kreal), parameter :: factors(5) = (/1.0d0, 1.0d-12, 1.0d-6, 1.0d6, 1.0d12/)
    real(kind=kreal) :: ref_stress(6,4), ref_strain(6,4), ref_tangent(5,5,4), ref_state(4), ref_work(4)
    real(kind=kreal) :: factor, error
    integer :: unit, step

    do unit = 1, size(factors)
      factor = factors(unit)
      mat%variables(M_YOUNGS) = factor*210000.0d0
      mat%variables(M_PLCONST1) = factor*1000.0d0
      elastic_table%tbval(1,1) = mat%variables(M_YOUNGS)
      oldstress = 0.0d0
      oldplastic = 0.0d0
      state = 0.0d0
      error = 0.0d0
      do step = 1, 4
        select case(step)
        case(1)
          deps = 0.0d0
        case(2)
          deps = (/1.0d-4, -2.0d-5, 0.0d0, 3.0d-5, 2.0d-5, -1.0d-5/)
        case(3)
          deps = (/0.012d0, 0.001d0, 0.0d0, 0.002d0, 0.001d0, -0.001d0/)
        case(4)
          deps = -0.02d0*deps
        end select
        call response(deps, stress, tangent)
        if(unit == 1) then
          ref_stress(:,step) = stress
          ref_strain(:,step) = epsout
          ref_tangent(:,:,step) = tangent
          ref_state(step) = outstate(1)
          ref_work(step) = work
        else
          error = max(error, maxval(abs(stress/factor-ref_stress(:,step)))/ &
            max(1.0d0, maxval(abs(ref_stress(:,step)))))
          call check('unit-invariant stress', error, 0.0d0, 1.0d-9)
          call check('unit-invariant strain', maxval(abs(epsout-ref_strain(:,step))), 0.0d0, 1.0d-11)
          call check('unit-invariant tangent', maxval(abs(tangent/factor-ref_tangent(:,:,step)))/ &
            maxval(abs(ref_tangent(:,:,step))), 0.0d0, 1.0d-9)
          call check('unit-invariant plastic strain', outstate(1), ref_state(step), 1.0d-11)
          call check('unit-invariant work', work/factor, ref_work(step), 1.0d-9)
        endif
        oldstress = stress
        oldplastic = outstate(1)
        state = outstate
      enddo
      print *, 'Unit-scaled response error, factor', factor, error
    enddo
    mat%variables(M_YOUNGS) = 210000.0d0
    mat%variables(M_PLCONST1) = 1000.0d0
    elastic_table%tbval(1,1) = mat%variables(M_YOUNGS)
  end subroutine check_stress_units

  subroutine response(increment, s, d)
    real(kind=kreal), intent(in) :: increment(6)
    real(kind=kreal), intent(out) :: s(6), d(5,5)
    call Update_PlaneStressJ2(mat, increment, oldstress, oldplastic, state, &
      s, epsout, d, istat, outstate, potential, work, ierr)
    if(ierr /= 0) then
      print *, trim(PlaneStressJ2ErrorMessage(ierr))
      stop 3
    endif
    call check('plane stress', s(3), 0.0d0, 1.0d-8)
  end subroutine response

  subroutine check(label, actual, expected, tolerance)
    character(len=*), intent(in) :: label
    real(kind=kreal), intent(in) :: actual, expected, tolerance
    if(abs(actual-expected) > tolerance) then
      print *, 'FAIL ', label, actual, expected
      stop 4
    endif
  end subroutine check
end program check_j2_material
