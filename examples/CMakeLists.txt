###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################

cmake_minimum_required(VERSION 2.8.11)

project(Fistr_Test)

find_program(BASH_PROGRAM bash)

add_test(
 NAME Static_exA_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exA
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exB_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exB
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exC_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exC
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exD_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exD
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exE_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exE
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exF_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exF
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exG_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exG
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exI_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/exI
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_spring_boundary_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src static/spring_boundary
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Eigen_exJ_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src eigen/exJ
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Eigen_exK_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src eigen/exK
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exM_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exM
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exN_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exN
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exO_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exO
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exP_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exP
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exQ_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exQ
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exR_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exR
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exS_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exS
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exT_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exT
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exU_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exU
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exV_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src heat/exV
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Dynamic_exW_Test
 COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src dynamic/exW
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

#add_test(
# NAME Dynamic_exX_Test
# COMMAND ruby ./test_FrontISTR.rb ${Fistr_BINARY_DIR}/src dynamic/exX
# WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#)
