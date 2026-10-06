program check_j2_potential
  use m_static_LIB_shell, only: UPDATE_Shell_MITC, STF_Shell_MITC
  use mMechGauss
  use m_ElastoPlastic, only: updateEPState
  implicit none
  type(tMaterial), target :: material
  type(tTable) :: table
  type(tElement) :: element, saved
  real(kind=kreal) :: coords(3,4), base(6,4), trial(6,4), perturb(6,4)
  real(kind=kreal) :: force(24), dummy(24), fd(24), potential, plus, minus, repeated, error, h
  real(kind=kreal) :: stiffness(24,24), fd_stiffness(24,24), plus_force(24), minus_force(24)
  integer :: ig, stage, j, node, component, ih

  call initMaterial(material)
  material%mtype = 120000
  material%variables(M_YOUNGS) = 210000.0d0
  material%variables(M_POISSON) = 0.3d0
  material%variables(M_PLCONST1) = 1000.0d0
  material%totallyr = 2
  allocate(material%shell_var(2))
  material%shell_var%ortho = 0
  material%shell_var%ee = 210000.0d0
  material%shell_var%pp = 0.3d0
  material%shell_var%alpha_over_mu = 0.1d0
  material%shell_var(1)%weight = 0.25d0
  material%shell_var(2)%weight = 0.75d0
  call init_table(table, 0, 2, 1, reshape((/210000.0d0,0.3d0/),(/2,1/)))
  call dict_add_key(material%dict, MC_ISOELASTIC, table)
  element%etype = 741
  allocate(element%gausses(4))
  do ig = 1, 4
    element%gausses(ig)%pMaterial => material
    call fstr_init_gauss(element%gausses(ig))
  enddo
  call fstr_init_shell_layer_gausses(element, 4, 2, 2)
  saved%etype = element%etype
  allocate(saved%gausses(4))
  do ig = 1, 4
    saved%gausses(ig)%pMaterial => material
    call fstr_init_gauss(saved%gausses(ig))
  enddo
  call fstr_init_shell_layer_gausses(saved, 4, 2, 2)
  print *, 'History bytes per point (descriptor and status payload, excluding allocator overhead):', &
    storage_size(element%shell_layer_gausses(1))/8 &
    +size(element%shell_layer_gausses(1)%fstatus)*storage_size(0.0_kreal)/8 &
    +size(element%shell_layer_gausses(1)%istatus)*storage_size(0)/8
  ! Nonuniform surface Jacobian and unequal layer thicknesses exercise both weights.
  coords(:,1) = (/0.0d0,0.0d0,0.0d0/)
  coords(:,2) = (/2.0d0,0.0d0,0.0d0/)
  coords(:,3) = (/1.0d0,1.0d0,0.0d0/)
  coords(:,4) = (/0.0d0,1.0d0,0.0d0/)
  base = 0.0d0

  do stage = 1, 3
    if(stage == 1) then
      do node = 1, 4
        trial(:,node) = (/1.0d-4*coords(1,node), -2.0d-5*coords(2,node), &
          2.0d-4*coords(1,node)*coords(2,node), 3.0d-4*coords(2,node), &
          -2.0d-4*coords(1,node), 1.0d-3/)
      enddo
    else if(stage == 2) then
      trial = 100.0d0*base
    else
      trial = 0.98d0*base
    endif
    call response(trial, potential, force)
    if(stage == 1) then
      if(abs(potential-0.5d0*dot_product(force,reshape(trial,(/24/)))) > 1.0d-11) stop 1
    endif
    call response(trial, repeated, dummy)
    if(abs(repeated-potential) > 1.0d-12 .or. maxval(abs(dummy-force)) > 1.0d-12) stop 2
    call STF_Shell_MITC(741, 4, 6, coords, element%gausses, stiffness, 0.2d0, 0, &
      nddisp=trial, element=element)
    error = maxval(abs(stiffness-transpose(stiffness)))/maxval(abs(stiffness))
    if(error > 1.0d-12) stop 'FAIL tangent symmetry'
    do ih = 1, 2
      h = 10.0d0**(-6-ih)
      do j = 1, 24
        node = (j-1)/6+1
        component = mod(j-1,6)+1
        perturb = trial
        perturb(component,node) = perturb(component,node)+h
        call response(perturb, plus, plus_force)
        perturb(component,node) = trial(component,node)-h
        call response(perturb, minus, minus_force)
        fd(j) = (plus-minus)/(2.0d0*h)
        fd_stiffness(:,j) = (plus_force-minus_force)/(2.0d0*h)
      enddo
      error = maxval(abs(fd-force)/max(1.0d0,abs(force)))
      print *, 'Potential gradient error, stage, h', stage, h, error
      if(error > 1.0d-6) stop 3
      error = sqrt(sum((fd_stiffness-stiffness)**2)/sum(stiffness**2))
      print *, 'Element tangent FD relative error, stage, h', stage, h, error
      if(error > 1.0d-6) stop 'FAIL element tangent'
    enddo
    call response(trial, potential, force)
    if(stage == 2 .and. maxval(element%shell_layer_gausses%plpotential) >= 0.0d0) stop 4
    do ig = 1, size(element%shell_layer_gausses)
      call updateEPState(element%shell_layer_gausses(ig))
      element%shell_layer_gausses(ig)%strain_bak = element%shell_layer_gausses(ig)%strain
      element%shell_layer_gausses(ig)%stress_bak = element%shell_layer_gausses(ig)%stress
      element%shell_layer_gausses(ig)%strain_energy_bak = element%shell_layer_gausses(ig)%strain_energy
    enddo
    base = trial
    ! A copied history must not alias the live status arrays during rollback.
    call fstr_copy_shell_layer_gausses(element, saved)
    if(associated(saved%shell_layer_gausses(1)%fstatus, element%shell_layer_gausses(1)%fstatus)) &
      stop 'FAIL history alias'
    call response(1.5d0*trial, repeated, dummy)
    call fstr_copy_shell_layer_gausses(saved, element)
    call response(trial, repeated, dummy)
    if(abs(repeated-potential) > 1.0d-10 .or. maxval(abs(dummy-force)) > 1.0d-10) &
      stop 'FAIL restored history response'
  enddo
  call fstr_finalize_shell_layer_gausses(saved)
  if(associated(saved%shell_layer_gausses)) stop 'FAIL history deallocation'
  do ig = 1, size(saved%gausses)
    call fstr_finalize_gauss(saved%gausses(ig))
  enddo
  deallocate(saved%gausses)
  do ig = 1, size(element%shell_layer_gausses)
    call fstr_finalize_gauss(element%shell_layer_gausses(ig))
  enddo
  do ig = 1, size(element%gausses)
    call fstr_finalize_gauss(element%gausses(ig))
  enddo
  deallocate(element%gausses, element%shell_layer_gausses, material%shell_var)
  call finalizeMaterial(material)
  print *, 'PASS potential, element tangent, repeated trials, plastic loading/unloading, history copy/restore'
contains
  subroutine response(displacement, energy, qforce)
    real(kind=kreal), intent(in) :: displacement(6,4)
    real(kind=kreal), intent(out) :: energy, qforce(24)
    call UPDATE_Shell_MITC(741, 4, 6, coords, base, displacement-base, element%gausses, &
      qforce, 0.2d0, 0, element=element)
    energy = sum(element%gausses%strain_energy)
  end subroutine response
end program check_j2_potential
