###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(fistr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/material.f90
  ${CMAKE_CURRENT_LIST_DIR}/mechgauss.f90
  ${CMAKE_CURRENT_LIST_DIR}/ElasticLinear.f90
  ${CMAKE_CURRENT_LIST_DIR}/Hyperelastic.f90
  ${CMAKE_CURRENT_LIST_DIR}/Elastoplastic.f90
  ${CMAKE_CURRENT_LIST_DIR}/Viscoelastic.f90
  ${CMAKE_CURRENT_LIST_DIR}/creep.f90
  ${CMAKE_CURRENT_LIST_DIR}/calMatMatrix.f90
)
