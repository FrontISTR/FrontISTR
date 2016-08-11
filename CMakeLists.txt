cmake_minimum_required(VERSION 2.8.11)

project(FrontISTR Fortran C CXX)

###############################################################################
# Windows (MinGW) specific setting
###############################################################################
if(CMAKE_GENERATOR MATCHES "(MinGW|MSYS)")
	add_definitions(-D_WINDOWS)
	set(CMAKE_Fortran_FLAGS "-fno-range-check")
	set(_WINDOWS ON)
endif()

###############################################################################
# Options
###############################################################################
# -DWITH_TOOLS       : compile tools
# -DWITH_MPI         : for parallel environment with MPI
# -DWITH_OPENMP      : for multi-(core|processor) environment
# -DWITH_REFINER     : compile with REVOCAP_Refiner
# -DWITH_REVOCAP     : compile with REVOCAP_Coupler
# -DWITH_METIS       : compile with METIS graph partitioning package
#   -DMETIS_VER_4    : uses METIS-4
# -DWITH_MUMPS       : compile with MUMPS direct solver
# -DWITH_LAPACK      : compile with LAPACK
# -DWITH_ML          : compile with Trilinos multigrid preconditioning package
# -DBUILD_DOC        : generate API documents

# add module path for FindPackage.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# -DWITH_TOOLS
option(WITH_TOOLS "compile tools" OFF)

# -DWITH_MPI
option(WITH_MPI "for parallel environment with MPI" OFF)
if(NOT WITH_MPI)
 	set(HECMW_SERIAL ON)
	add_definitions(-DHECMW_SERIAL)
else()
	find_package(MPI REQUIRED)
	if(MPI_Fortran_FOUND AND MPI_CXX_FOUND AND MPI_C_FOUND)
		include_directories(${MPI_C_INCLUDE_PATH})
		set(CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
		set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
		set(CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortra_LINK_FLAGS})
	endif()
endif()

# -DWITH_OPENMP
option(WITH_OPENMP "for multi-(core|processor) environment" OFF)
if(WITH_OPENMP)
	find_package(OpenMP)
	if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
	else()
    set(WITH_OPENMP OFF)
	endif()
endif()

# -DWITH_REFINER
option(WITH_REFINER "compile with REVOCAP_Refiner" OFF)
if(WITH_REFINER)
	find_package(Refiner)
	if(REFINER_FOUND)
		include_directories(${REFINER_INCLUDE_PATH})
    add_definitions(-DHECMW_WITH_REFINER)
		set(HECMW_WITH_REFINER ON)
	else()
		set(WITH_REFINER OFF)
		set(HECMW_WITH_REFINER OFF)
	endif()
endif()

# -DWITH_REVOCAP
option(WITH_REVOCAP "compile with REVOCAP_Coupler" OFF)
if(WITH_REVOCAP)
	find_package(Revocap)
	if(REVOCAP_FOUND)
		include_directories(${REVOCAP_INCLUDE_PATH})
	else()
		set(WITH_REVOCAP OFF)
	endif()
endif()

# -DWITH_METIS, -DMETIS_VER_4
option(WITH_METIS "compile with METIS graph partitioning package" OFF)
if(WITH_METIS)
	option(METIS_VER_4 "uses METIS-4" OFF)
	find_package(Metis)
	if(METIS_FOUND)
		include_directories(${METIS_INCLUDE_PATH})
		add_definitions(-DHECMW_WITH_METIS)
		if(NOT METIS_VER_4)
			add_definitions(-DHECMW_METIS_VER=5)
		else()
			add_definitions(-DHECMW_METIS_VER=4)
		endif()
	else()
		set(WITH_METIS OFF)
	endif()
endif()

# -DWITH_MUMPS
option(WITH_MUMPS "compile with MUMPS direct solver" OFF)
if(WITH_MUMPS)
	find_package(Mumps)
	find_package(Scalapack)
	if(MUMPS_FOUND AND SCALAPACK_FOUND AND WITH_METIS AND WITH_MPI)
		include_directories(${MUMPS_INCLUDE_PATH})
		add_definitions(-DWITH_MUMPS)
	else()
		message(WARNING "Please turn on WITH_MPI and install scalapack")
		set(WITH_MUMPS OFF)
	endif()
endif()

# -DWITH_LAPACK
option(WITH_LAPACK "compile with LAPACK" OFF)
if(WITH_LAPACK)
	find_package(LAPACK)
	if(LAPACK_FOUND AND BLAS_FOUND)
		add_definitions(-DWITH_LAPACK)
	else()
    set(WITH_LAPACK OFF)
    set(BLAS_LIBRARIES "BLAS_LIBRARIES-NOTFOUND" CACHE FILEPATH "BLAS Library")
    set(LAPACK_LIBRARIES "LAPACK_LIBRARIES-NOTFOUND" CACHE FILEPATH "LAPACK Library")
    message(WARNING "Please specify -DBLAS_LIBRARIES and -DLAPACK_LIBRARIES manually.")
	endif()
endif()

# -DWITH_ML
option(WITH_ML "compile with Trilinos multigrid preconditioning package" OFF)
if(WITH_ML)
	find_package(Trilinos)
	if(Trilinos_FOUND AND WITH_MPI)
		include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
		set(HECMW_WITH_ML ON)
		add_definitions(-DHECMW_WITH_ML)
	else()
		message(WARNING "Please turn on WITH_MPI")
		set(WITH_ML OFF)
	endif()
endif()

# -DBUILD_DOC
option(BUILD_DOC "Generate API documents." OFF)
if(BUILD_DOC)
	file(COPY
		${CMAKE_SOURCE_DIR}/doc/fistr.css
		DESTINATION
		${CMAKE_BINARY_DIR}/doc/
	)
	file(COPY
		${CMAKE_SOURCE_DIR}/doc/fistr_logo.png
		DESTINATION
		${CMAKE_BINARY_DIR}/doc/
	)
  set(doxyfile_in ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
	set(doxyfile ${CMAKE_BINARY_DIR}/doc/Doxyfile)
	configure_file(${doxyfile_in} ${doxyfile} @ONLY)

	find_package(Doxygen)
	if(DOXYGEN_FOUND)
		message(STATUS "Generate API documents. Please type \"make doc\"")
		set(doxy_main_page "${CMAKE_SOURCE_DIR}/INSTALL.md")
		add_custom_target(doc
			COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			COMMENT "Generate API documents."
			VERBATIM
		)
	else()
		message(WARNING "doxygen not found.")
	endif()
endif()

add_subdirectory(hecmw1)
