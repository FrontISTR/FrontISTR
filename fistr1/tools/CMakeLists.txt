###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

include_directories(
  neu2fstr/
  neu2fstr/NFD/
  neu2fstr/HECD/
  neu2fstr/converter/
)

set(NEU2FSTR_SRCS
neu2fstr/neu2fstr.cpp
neu2fstr/NFD/CNFData.cpp
neu2fstr/NFD/CNFDB_100.cpp
neu2fstr/NFD/CNFDB_402.cpp
neu2fstr/NFD/CNFDB_403.cpp
neu2fstr/NFD/CNFDB_404.cpp
neu2fstr/NFD/CNFDB_405.cpp
neu2fstr/NFD/CNFDB_408.cpp
neu2fstr/NFD/CNFDB_506.cpp
neu2fstr/NFD/CNFDB_507.cpp
neu2fstr/NFD/CNFDB_601.cpp
neu2fstr/NFD/CNFMessage.cpp
neu2fstr/HECD/hecd_util.cpp
neu2fstr/HECD/CHECData.cpp
neu2fstr/HECD/HECDB_util.cpp
neu2fstr/HECD/CHECDB_Amplitude.cpp
neu2fstr/HECD/CHECDB_EGroup.cpp
neu2fstr/HECD/CHECDB_Element.cpp
neu2fstr/HECD/CHECDB_Header.cpp
neu2fstr/HECD/CHECDB_Material.cpp
neu2fstr/HECD/CHECDB_NGroup.cpp
neu2fstr/HECD/CHECDB_Node.cpp
neu2fstr/HECD/CHECDB_Section.cpp
neu2fstr/HECD/CHECDB_SGroup.cpp
neu2fstr/HECD/CHECDB_Visual.cpp
neu2fstr/HECD/CHECDB_Zero.cpp
neu2fstr/HECD/CFSTRData.cpp
neu2fstr/HECD/FSTRDB_util.cpp
neu2fstr/HECD/CFSTRDB_Boundary.cpp
neu2fstr/HECD/CFSTRDB_CFlux.cpp
neu2fstr/HECD/CFSTRDB_CLoad.cpp
neu2fstr/HECD/CFSTRDB_DFlux.cpp
neu2fstr/HECD/CFSTRDB_DLoad.cpp
neu2fstr/HECD/CFSTRDB_Echo.cpp
neu2fstr/HECD/CFSTRDB_Eigen.cpp
neu2fstr/HECD/CFSTRDB_Film.cpp
neu2fstr/HECD/CFSTRDB_Fixtemp.cpp
neu2fstr/HECD/CFSTRDB_Heat.cpp
neu2fstr/HECD/CFSTRDB_Radiate.cpp
neu2fstr/HECD/CFSTRDB_Reftemp.cpp
neu2fstr/HECD/CFSTRDB_SFilm.cpp
neu2fstr/HECD/CFSTRDB_SFlux.cpp
neu2fstr/HECD/CFSTRDB_Solution.cpp
neu2fstr/HECD/CFSTRDB_Solver.cpp
neu2fstr/HECD/CFSTRDB_SRadiate.cpp
neu2fstr/HECD/CFSTRDB_Static.cpp
neu2fstr/HECD/CFSTRDB_Step.cpp
neu2fstr/HECD/CFSTRDB_Temperature.cpp
neu2fstr/HECD/CFSTRDB_Write.cpp
neu2fstr/converter/cconv_mat.cpp
neu2fstr/converter/CConvMessage.cpp
neu2fstr/converter/conv_neu2hec.cpp
neu2fstr/converter/conv_neu2fstr_static.cpp
neu2fstr/converter/conv_neu2fstr_dynamic.cpp
neu2fstr/converter/conv_neu2fstr_heat.cpp
)

# neu2fstr
add_executable(neu2fstr ${NEU2FSTR_SRCS})
