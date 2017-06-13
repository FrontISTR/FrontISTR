###############################################################################
# Copyright (c) 2016 FrontISTR Forum, The University of Tokyo
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

project(Fistr_Test)

find_program(BASH_PROGRAM bash)

add_test(
 NAME Static_exA_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exA
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exB_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exB
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exC_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exC
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exD_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exD
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exE_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exE
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exF_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exF
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exG_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exG
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_exI_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/exI
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_restart_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/restart
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_spring_boundary_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/spring_boundary
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_thermal_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/thermal
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Static_torque_rot_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 static/torque_rot
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Eigen_exJ_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 eigen/exJ
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Eigen_exK_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 eigen/exK
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exM_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exM
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exN_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exN
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exO_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exO
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exP_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exP
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exQ_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exQ
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exR_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exR
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exS_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exS
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exT_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exT
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exU_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exU
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Heat_exV_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 heat/exV
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Dynamic_exW_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 dynamic/exW
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(
 NAME Dynamic_exX_Test
 COMMAND ruby ./test_FrontISTR.rb ${CMAKE_BINARY_DIR}/fistr1 dynamic/exX
 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
