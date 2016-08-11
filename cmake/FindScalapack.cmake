#
# SCALAPACK_FOUND
# SCALAPACK_LIBRARIES
#
find_library(SCALAPACK_LIBRARIES
  NAMES scalapack
  PATHS $ENV{SCALAPACK_ROOT}
    ${CMAKE_SOURCE_DIR}/../scalapack-2.0.2
    $ENV{HOME}/local/lib
    $ENV{HOME}/.local/lib
    /usr/local/lib
    /usr/lib
    NO_DEFAULT_PATH
)
if(SCALAPACK_LIBRARIES)
  set(SCALAPACK_FOUND ON)
else()
  set(SCALAPACK_LIBRARIES "SCALAPACK_LIBRARIES-NOTFOUND" CACHE FILEPATH "libscalapack")
endif()
