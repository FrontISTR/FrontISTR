###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

#
# Create executable targets first
#
add_executable(hec2rcap)
add_executable(hecmw_part1)
add_executable(rmerge)
add_executable(rconv)
add_executable(hecmw_vis1)

#
# Configure common properties for all targets
#
foreach(target hec2rcap hecmw_part1 rmerge rconv hecmw_vis1)
  set_target_properties(${target} PROPERTIES LINKER_LANGUAGE "CXX")
  target_link_libraries(${target} PRIVATE hecmw)
  target_include_directories(${target} PUBLIC ${FrontISTR_INCLUDE_DIRS})
  target_compile_definitions(${target} PUBLIC ${FrontISTR_DEFINITIONS})
endforeach()

#
# Add subdirectories that will contribute sources to the respective targets
# Use second parameter for hec2rcap to avoid directory-target name conflict
#
add_subdirectory(hec2rcap hec2rcap_dir)
add_subdirectory(partitioner)
add_subdirectory(result_file_merger)
add_subdirectory(result_type_converter)
add_subdirectory(visualizer)

#
# Set installing files
#
install(TARGETS hec2rcap hecmw_part1 rmerge rconv hecmw_vis1
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)
