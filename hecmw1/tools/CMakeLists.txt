###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

#
# Build executable hec2rcap
#
set(hec2rcap_SOURCES)
include(hec2rcap/CMakeLists.txt)
add_executable(hec2rcap ${hec2rcap_SORUCES})
set_target_properties(hec2rcap PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(hec2rcap hecmw ${hecmw_LIBRARIES})
target_include_directories(hec2rcap PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(hec2rcap PUBLIC ${hecmw_DEFINITIONS})

#
# Build executable hecmw_part1
#
set(hecmw_part1_SOURCES)
include(partitioner/CMakeLists.txt)
add_executable(hecmw_part1 ${hecmw_part1_SOURCES})
set_target_properties(hecmw_part1 PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(hecmw_part1 hecmw ${hecmw_LIBRARIES})
target_include_directories(hecmw_part1 PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(hecmw_part1 PUBLIC ${hecmw_DEFINITIONS})

#
# Build executable rmerge
#
set(rmerge_SOURCES)
include(result_file_merger/CMakeLists.txt)
add_executable(rmerge ${rmerge_SOURCES})
set_target_properties(rmerge PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(rmerge hecmw ${hecmw_LIBRARIES})
target_include_directories(rmerge PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(rmerge PUBLIC ${hecmw_DEFINITIONS})

#
# Build executable rconv
#
set(rconv_SOURCES)
include(result_type_converter/CMakeLists.txt)
add_executable(rconv ${rconv_SOURCES})
set_target_properties(rconv PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(rconv hecmw ${hecmw_LIBRARIES})
target_include_directories(rconv PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(rconv PUBLIC ${hecmw_DEFINITIONS})

#
# Build executable hecmw_vis1
#
set(hecmw_vis1_SOURCES)
include(visualizer/CMakeLists.txt)
add_executable(hecmw_vis1 ${hecmw_vis1_SOURCES})
set_target_properties(hecmw_vis1 PROPERTIES LINKER_LANGUAGE "CXX")
target_link_libraries(hecmw_vis1 hecmw ${hecmw_LIBRARIES})
target_include_directories(hecmw_vis1 PUBLIC ${hecmw_INCLUDE_DIRS})
target_compile_definitions(hecmw_vis1 PUBLIC ${hecmw_DEFINITIONS})

#
# Set installing files
#
install(TARGETS hec2rcap hecmw_part1 rmerge rconv hecmw_vis1
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)
