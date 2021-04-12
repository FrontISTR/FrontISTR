###############################################################################
# Copyright (c) 2021 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

add_subdirectory(contact)
add_subdirectory(element)
add_subdirectory(physics)
add_subdirectory(user)
add_subdirectory(utilities)

target_sources(fistr
	PRIVATE
    GaussM.f90
    dynamic_mass.f90
    eigen_LIB.f90
    heat_LIB.f90
    heat_LIB_CAPACITY.f90
    heat_LIB_CONDUCTIVITY.f90
    heat_LIB_DFLUX.f90
    heat_LIB_FILM.f90
    heat_LIB_NEUTRAL.f90
    heat_LIB_RADIATE.f90
    m_common_struct.f90
    m_fstr.f90
    m_timepoint.f90
    m_out.f90
    m_step.f90
    precheck_LIB_2d.f90
    precheck_LIB_3d.f90
    precheck_LIB_shell.f90
    solve_LINEQ.f90
    static_LIB.f90
    static_LIB_1d.f90
    static_LIB_2d.f90
    static_LIB_3d.f90
    static_LIB_3dIC.f90
    static_LIB_3d_vp.f90
    static_LIB_C3D8.f90
    static_LIB_Fbar.f90
    static_LIB_beam.f90
    static_LIB_shell.f90
)
