###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

include("${CMAKE_CURRENT_LIST_DIR}/communication/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/matrix/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/init/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/sparse_matrix/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_direct_mumps/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_direct_clustermkl/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_11/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_22/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_33/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_44/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_66/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_direct/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_direct_parallel/CMakeLists.txt")
include("${CMAKE_CURRENT_LIST_DIR}/solver_direct_lag/CMakeLists.txt")
