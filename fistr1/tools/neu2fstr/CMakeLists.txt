add_executable(neu2fstr)

target_sources(neu2fstr
  PRIVATE
  neu2fstr.cpp
)
target_include_directories(neu2fstr
  PRIVATE
  converter
  NFD
  HECD
  )
add_subdirectory(HECD)
add_subdirectory(NFD)
add_subdirectory(converter)
