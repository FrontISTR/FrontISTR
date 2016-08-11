cmake_minimum_required(VERSION 2.8.11)

project(Hecmw C Fortran)

include_directories(
	${Hecmw_SOURCE_DIR}/src/common
	${Hecmw_SOURCE_DIR}/src/hecmw
	${Hecmw_SOURCE_DIR}/src/visualizer
	${Hecmw_BINARY_DIR}
)

add_subdirectory(src)
add_subdirectory(tools)
