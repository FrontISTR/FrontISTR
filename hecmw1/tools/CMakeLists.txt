###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

# hec2rcap
set(hec2rcap_SOURCES)
include(hec2rcap/CMakeLists.txt)
add_executable(hec2rcap ${hec2rcap_SORUCES})
target_link_libraries(hec2rcap hecmw ${EXTRA_LIBS})

# hecmw_part1
set(hecmw_part1_SOURCES)
include(partitioner/CMakeLists.txt)
add_executable(hecmw_part1 ${hecmw_part1_SOURCES})
target_link_libraries(hecmw_part1 hecmw ${EXTRA_LIBS})

# rmerge
set(rmerge_SOURCES)
include(result_file_merger/CMakeLists.txt)
add_executable(rmerge ${rmerge_SOURCES})
target_link_libraries(rmerge hecmw ${EXTRA_LIBS})

# rconv
set(rconv_SOURCES)
include(result_type_converter/CMakeLists.txt)
add_executable(rconv ${rconv_SOURCES})
target_link_libraries(rconv hecmw ${EXTRA_LIBS})

# hecmw_vis1
set(hecmw_vis1_SOURCES)
include(visualizer/CMakeLists.txt)
add_executable(hecmw_vis1 ${hecmw_vis1_SOURCES})
target_link_libraries(hecmw_vis1 hecmw ${EXTRA_LIBS})

if(WITH_REVOCAP OR WITH_REFINER)
  set_target_properties(hec2rcap PROPERTIES LINKER_LANGUAGE "CXX")
  set_target_properties(hecmw_part1 PROPERTIES LINKER_LANGUAGE "CXX")
  set_target_properties(rmerge PROPERTIES LINKER_LANGUAGE "CXX")
  set_target_properties(rconv PROPERTIES LINKER_LANGUAGE "CXX")
  set_target_properties(hecmw_vis1 PROPERTIES LINKER_LANGUAGE "CXX")
endif()

install(TARGETS hec2rcap hecmw_part1 rmerge rconv hecmw_vis1
  DESTINATION bin
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)
