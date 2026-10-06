program check_j2_output
  use m_fstr_NodalStress, only: NodalStress_ShellJ2
  use mMechGauss, only: fstr_element_average_plstrain
  use m_fstr
  use elementInfo, only: getQuadPoint, getNodalNaturalCoord
  implicit none
  type(tElement) :: element
  type(tMaterial), target :: material
  type(fstr_solid) :: solid
  type(fstr_solid_physic_val), pointer :: layer
  integer :: ig, ilayer, ithick, ishell, side
  real(kind=kreal) :: strain(4,6), stress(4,6), estrain(6), estress(6), coord(2), nodes(4,2)
  real(kind=kreal) :: zeta, x, centers(2), fractions(2)

  element%etype = 741
  element%shell_nlayer = 2
  element%shell_nthick = 2
  material%totallyr = 2
  allocate(material%shell_var(2))
  fractions = (/0.25d0, 0.75d0/)
  centers = (/-0.75d0, 0.25d0/)
  material%shell_var(1)%weight = fractions(1)
  material%shell_var(2)%weight = fractions(2)
  allocate(element%gausses(4), element%shell_layer_gausses(16))
  element%gausses(1)%pMaterial => material
  allocate(solid%SHELL)
  allocate(solid%SHELL%LAYER(2))
  do ilayer = 1, 2
    allocate(solid%SHELL%LAYER(ilayer)%PLUS, solid%SHELL%LAYER(ilayer)%MINUS)
    call allocate_layer(solid%SHELL%LAYER(ilayer)%PLUS)
    call allocate_layer(solid%SHELL%LAYER(ilayer)%MINUS)
  enddo
  do ig = 1, 4
    call getQuadPoint(741, ig, coord)
    x = coord(1)+2.0d0*coord(2)
    do ilayer = 1, 2
      do ithick = 1, 2
        zeta = centers(ilayer)+fractions(ilayer)*(2*ithick-3)/sqrt(3.0d0)
        ishell = ((ig-1)*2+ilayer-1)*2+ithick
        element%shell_layer_gausses(ishell)%strain_out = 2.0d0+x+3.0d0*zeta
        element%shell_layer_gausses(ishell)%stress_out = 4.0d0+x+5.0d0*zeta*zeta
        element%shell_layer_gausses(ishell)%plstrain = 0.02d0+0.01d0*(x+zeta*zeta)
      enddo
    enddo
  enddo
  call NodalStress_ShellJ2(element, 4, solid, 1, (/1,2,3,4/), strain, stress, estrain, estress)
  ! Retain the existing history-point mean, not a physical-volume average.
  call check('history strain mean', maxval(abs(estrain-1.25d0)), 0.0d0)
  call check('history stress mean', maxval(abs(estress-73.0d0/12.0d0)), 0.0d0)
  call check('history plastic strain mean', fstr_element_average_plstrain(element), 0.02d0+0.01d0*5.0d0/12.0d0)
  call getNodalNaturalCoord(741, nodes)
  call check('nodal strain recovery', maxval(abs(strain(:,1)-(1.25d0+nodes(:,1)+2.0d0*nodes(:,2)))), 0.0d0)
  do ilayer = 1, 2
    do side = 1, 2
      layer => solid%SHELL%LAYER(ilayer)%MINUS
      if(side == 2) layer => solid%SHELL%LAYER(ilayer)%PLUS
      zeta = centers(ilayer)+fractions(ilayer)*(2*side-3)/sqrt(3.0d0)
      call check('layer strain mean', layer%ESTRAIN(1), 2.0d0+3.0d0*zeta)
      call check('layer stress mean', layer%ESTRESS(1), 4.0d0+5.0d0*zeta*zeta)
      call check('layer plastic strain mean', layer%EPLSTRAIN(1), 0.02d0+0.01d0*zeta*zeta)
    enddo
  enddo
  print *, 'PASS stored-history output and nodal recovery with existing averaging convention'
contains
  subroutine allocate_layer(phys)
    type(fstr_solid_physic_val), intent(inout) :: phys
    allocate(phys%ESTRAIN(6), phys%ESTRESS(6), phys%EPLSTRAIN(1), phys%STRAIN(24), phys%STRESS(24))
    phys%STRAIN = 0.0d0
    phys%STRESS = 0.0d0
  end subroutine allocate_layer
  subroutine check(label, actual, expected)
    character(len=*), intent(in) :: label
    real(kind=kreal), intent(in) :: actual, expected
    if(abs(actual-expected) > 1.0d-11) then
      print *, 'FAIL ', label, actual, expected
      stop 1
    endif
  end subroutine check
end program check_j2_output
