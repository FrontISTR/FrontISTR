# Build excutable hec2rcap
add_subdirectory(hec2rcap)

# Build executable hecmw_part1
add_subdirectory(partitioner)

# Build executable rmerge
add_subdirectory(result_file_merger)

# Build executable rconv
add_subdirectory(result_type_converter)

# Build executable hecmw_vis1
add_subdirectory(visualizer)

# Set intalling files
#install(TARGETS hec2rcap hecmw_part1 rmerge rconv hecmw_vis1
#  DESTINATION bin
#  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
#)
