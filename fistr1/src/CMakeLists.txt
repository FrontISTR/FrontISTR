cmake_minimum_required(VERSION 2.8.11)

if(NOT WITH_PARACON)
  set(Fistr_MAIN_SRCS main/fistr_main.f90)
else()
  set(Fistr_MAIN_SRCS main/fistr_main_contact.f90)
endif()

set(Fistr_SRCS
lib/fstr_lczparm.f90
lib/m_step.f90
lib/m_Makrose.f90
lib/m_MakrosePartMesh.f90
lib/m_out.f90
lib/m_common_struct.f90
lib/m_fstr.f90
lib/m_fstr_freqdata.f90
lib/GaussM.f90
lib/static_LIB_1d.f90
lib/static_LIB_2d.f90
lib/static_LIB_3d.f90
lib/static_LIB_C3D8.f90
lib/static_LIB_3dIC.f90
lib/static_LIB_beam.f90
lib/static_LIB_shell.f90
lib/heat_LIB_CAPACITY.f90
lib/heat_LIB_DFLUX.f90
lib/heat_LIB_FILM.f90
lib/heat_LIB_NEUTRAL.f90
lib/heat_LIB_RADIATE.f90
lib/heat_LIB_THERMAL.f90
lib/eigen_LIB_2d1mass.f90
lib/eigen_LIB_2d2mass.f90
lib/eigen_LIB_3d1mass.f90
lib/eigen_LIB_3d2mass.f90
lib/precheck_LIB_2d.f90
lib/precheck_LIB_3d.f90
lib/precheck_LIB_shell.f90
lib/solve_LINEQ.f90
lib/static_LIB.f90
lib/eigen_LIB.f90
lib/heat_LIB.f90
lib/m_fstr_para_contact.F90
lib/utilities/utilities.f90
lib/utilities/ttable.f90
lib/user/umat.f90
lib/user/uyield.f90
lib/user/uelastic.f90
lib/user/uload.f90
lib/element/line2n.f90
lib/element/line3n.f90
lib/element/tri3n.f90
lib/element/tri6n.f90
lib/element/quad4n.f90
lib/element/quad8n.f90
lib/element/quad9n.f90
lib/element/hex8n.f90
lib/element/hex20n.f90
lib/element/tet4n.f90
lib/element/tet10n.f90
lib/element/prism6n.f90
lib/element/prism15n.f90
lib/element/quadrature.f90
lib/element/element.f90
lib/physics/material.f90
lib/physics/mechgauss.f90
lib/physics/ElasticLinear.f90
lib/physics/Hyperelastic.f90
lib/physics/Elastoplastic.f90
lib/physics/Viscoelastic.f90
lib/physics/creep.f90
lib/physics/calMatMatrix.f90
lib/contact/surf_ele.f90
lib/contact/contact_lib.f90
lib/contact/fstr_contact_def.F90
common/hecmw2fstr_connect_conv.c
common/fstr_ctrl_util.c
common/fstr_sort_index.c
common/hecmw2fstr_mesh_conv.f90
common/fstr_setup_util.f90
common/fstr_ctrl_modifier.f90
common/fstr_contact.f90
common/fstr_ctrl_common.f90
common/fstr_ctrl_material.f90
common/fstr_ctrl_static.f90
common/fstr_ctrl_heat.f90
common/fstr_ctrl_eigen.f90
common/fstr_ctrl_dynamic.f90
common/fstr_get_prop.f90
common/fstr_setup.f90
common/fstr_debug_dump.f90
common/fstr_precheck.f90
common/fstr_rcap_io.F90
analysis/static/fstr_Spring.f90
analysis/static/fstr_Restart.f90
analysis/static/fstr_NodalStress.f90
analysis/static/fstr_StiffMatrix.f90
analysis/static/fstr_Update.f90
analysis/static/readtemp.f90
analysis/static/static_echo.f90
analysis/static/static_make_result.f90
analysis/static/fstr_ass_load.f90
analysis/static/static_output.f90
analysis/static/static_mat_ass_main.f90
analysis/static/fstr_mat_con_contact.f90
analysis/static/fstr_Residual.f90
analysis/static/fstr_AddContactStiff.f90
analysis/static/fstr_AddBC.f90
analysis/static/static_mat_ass.f90
analysis/static/fstr_solve_LINEAR.f90
analysis/static/fstr_mat_resid_contact.f90
analysis/static/sparse_matrix_contact.f90
analysis/static/solve_LINEQ_MUMPS_contact.f90
analysis/static/set_arrays_DirectSolver.f90
analysis/static/solve_LINEQ_mkl.F90
analysis/static/solve_LINEQ_direct_serial_lag.f90
analysis/static/solve_LINEQ_iter_contact.f90
analysis/static/solve_LINEQ_contact.f90
analysis/static/fstr_solve_NonLinear.f90
analysis/static/fstr_solve_NLGEOM.f90
analysis/dynamic/mode/fstr_EIG_lczeigenm.f90
analysis/dynamic/mode/fstr_EIG_getamat.f90
analysis/dynamic/mode/fstr_EIG_lanczos.f90
analysis/dynamic/mode/fstr_EIG_matmult.f90
analysis/dynamic/mode/fstr_EIG_mgs1.f90
analysis/dynamic/mode/fstr_EIG_setMASS.f90
analysis/dynamic/mode/fstr_EIG_tridiag.f90
analysis/dynamic/mode/fstr_EIG_output.f90
analysis/dynamic/mode/fstr_solve_eigen.f90
analysis/dynamic/freq/fstr_frequency_analysis.f90
analysis/dynamic/transit/table_dyn.f90
analysis/dynamic/transit/dynamic_mat_ass_load.f90
analysis/dynamic/transit/dynamic_mat_ass_bc.f90
analysis/dynamic/transit/dynamic_mat_ass_bc_vl.f90
analysis/dynamic/transit/dynamic_mat_ass_bc_ac.f90
analysis/dynamic/transit/dynamic_mat_ass_couple.f90
analysis/dynamic/transit/dynamic_var_init.f90
analysis/dynamic/transit/dynamic_make_result.f90
analysis/dynamic/transit/dynamic_output.f90
analysis/dynamic/transit/fstr_dynamic_explicit.f90
analysis/dynamic/transit/fstr_dynamic_implicit.f90
analysis/dynamic/transit/fstr_dynamic_nlexplicit.f90
analysis/dynamic/transit/fstr_dynamic_nlimplicit.f90
analysis/dynamic/transit/fstr_solve_dynamic.f90
analysis/heat/heat_init.f90
analysis/heat/heat_get_amplitude.f90
analysis/heat/heat_echo.f90
analysis/heat/heat_mat_ass_bc_CFLUX.f90
analysis/heat/heat_mat_ass_bc_DFLUX.f90
analysis/heat/heat_mat_ass_bc_FILM.f90
analysis/heat/heat_mat_ass_bc_FIXT.f90
analysis/heat/heat_mat_ass_bc_RADIATE.f90
analysis/heat/heat_mat_ass_boundary.f90
analysis/heat/heat_mat_ass_capacity.f90
analysis/heat/heat_mat_ass_conductivity.f90
analysis/heat/heat_solve_SS.f90
analysis/heat/heat_solve_TRAN.f90
analysis/heat/fstr_solve_heat.f90
)
set(Fistr_WITH_PARACON_SRCS
lib/m_Metis403API_f.F90
lib/m_Metis510API_f.F90
)
set(Fistr_WITHOUT_PARACON_SRCS
lib/m_Metis403API_f_dumm.f90
lib/m_Metis510API_f_dumm.f90
)

if(WITH_PARACON)
  set(Fistr_SRCS ${Fistr_SRCS} ${Fistr_WITH_PARACON_SRCS})
else()
  set(Fistr_SRCS ${Fistr_SRCS} ${Fistr_WITHOUT_PARACON_SRCS})
endif()

add_library(fistr ${Fistr_SRCS})
add_executable(fistr1 ${Fistr_MAIN_SRCS})

if(WITH_ML)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()

if(WITH_MUMPS)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${MUMPS_LIBRARIES} ${SCALAPACK_LIBRARIES})
endif()

if(WITH_LAPACK)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if(WITH_METIS)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${METIS_LIBRARIES})
endif()

if(WITH_REVOCAP)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${REVOCAP_LIBRARIES})
endif()

if(WITH_REFINER)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${REFINER_LIBRARIES})
endif()

if(WITH_MPI)
  set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_Fortran_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()

target_link_libraries(fistr1 fistr hecmw ${EXTRA_LIBS})

install(TARGETS fistr1 DESTINATION bin CONFIGURATIONS ${CMAKE_BUILD_TYPE})
