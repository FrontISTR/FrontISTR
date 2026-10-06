program check_shell_output
  use m_static_LIB_shell, only: ElementStress_Shell_MITC
  use mMechGauss
  use elementInfo, only: NumOfQuadPoints, getNodalNaturalCoord, NumOfShellThicknessQuadPoints
  implicit none
  type(tMaterial), target :: material
  type(tTable) :: table
  type(tElement) :: element
  integer :: itype, etype, nn, ng, ig, node, layer, thick, ierr
  integer, parameter :: types(3) = (/731, 741, 743/), nodes(3) = (/3, 4, 9/)
  real(kind=kreal), allocatable :: coords(:,:), disp(:,:), natural(:,:), strain(:,:), stress(:,:)
  real(kind=kreal) :: one_strain(1,6), one_stress(1,6), zeta, weight, x, y, error, max_error

  call initMaterial(material)
  material%mtype = ELASTIC
  material%nlgeom_flag = INFINITESIMAL
  material%variables(M_YOUNGS) = 210000.0d0
  material%variables(M_POISSON) = 0.3d0
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
  max_error = 0.0d0

  do itype = 1, size(types)
    etype = types(itype)
    nn = nodes(itype)
    ng = NumOfQuadPoints(etype)
    element%etype = etype
    allocate(element%gausses(ng), coords(3,nn), disp(6,nn), natural(nn,2), strain(ng,6), stress(ng,6))
    do ig = 1, ng
      element%gausses(ig)%pMaterial => material
      call fstr_init_gauss(element%gausses(ig))
    enddo
    call getNodalNaturalCoord(etype, natural)
    ! Skewed geometry, rotated out of the xy plane, with nonuniform displacement.
    do node = 1, nn
      x = 2.0d0*natural(node,1)+0.3d0*natural(node,2)
      y = natural(node,2)
      coords(:,node) = (/0.8d0*x, y, -0.6d0*x/)
      disp(:,node) = 1.0d-3*(/x*y, x, y, x, y, x-y/)
    enddo
    if(fstr_shell_output_point_count(element) /= ng*2*NumOfShellThicknessQuadPoints(etype)) stop 1
    do layer = 1, 2
      do thick = 1, NumOfShellThicknessQuadPoints(etype)
        call fstr_shell_thickness_quadrature(etype, thick, zeta, weight, ierr)
        if(ierr /= 0) stop 2
        call ElementStress_Shell_MITC(etype, nn, 6, coords, element%gausses, disp, &
          strain, stress, 0.2d0, zeta, layer, surface_gauss_points=.true.)
        do ig = 1, ng
          call ElementStress_Shell_MITC(etype, nn, 6, coords, element%gausses, disp, &
            one_strain, one_stress, 0.2d0, zeta, layer, surface_gauss_index=ig)
          error = max(maxval(abs(one_strain(1,:)-strain(ig,:))), maxval(abs(one_stress(1,:)-stress(ig,:))))
          max_error = max(max_error, error)
          if(error > 1.0d-11) stop 3
        enddo
      enddo
    enddo
    do ig = 1, ng
      call fstr_finalize_gauss(element%gausses(ig))
    enddo
    deallocate(element%gausses, coords, disp, natural, strain, stress)
  enddo
  deallocate(material%shell_var)
  call finalizeMaterial(material)
  print *, 'PASS MITC3/4/9 single-point and all-point output, two layers; max error', max_error
end program check_shell_output
