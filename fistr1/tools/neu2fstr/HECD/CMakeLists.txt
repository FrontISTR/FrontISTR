###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

target_sources(neu2fstr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/hecd_util.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECData.cpp
  ${CMAKE_CURRENT_LIST_DIR}/HECDB_util.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Amplitude.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_EGroup.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Element.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Header.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Material.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_NGroup.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Node.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Section.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_SGroup.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Visual.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CHECDB_Zero.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRData.cpp
  ${CMAKE_CURRENT_LIST_DIR}/FSTRDB_util.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Boundary.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_CFlux.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_CLoad.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_DFlux.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_DLoad.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Echo.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Eigen.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Film.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Fixtemp.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Heat.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Radiate.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Reftemp.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_SFilm.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_SFlux.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Solution.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Solver.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_SRadiate.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Static.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Step.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Temperature.cpp
  ${CMAKE_CURRENT_LIST_DIR}/CFSTRDB_Write.cpp
)

target_include_directories(neu2fstr PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
)

