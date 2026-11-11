###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

file(READ "${INPUT}" mesh)
foreach(mode combined mesh_only)
  set(dir "${WORK_DIR}/${mode}")
  file(MAKE_DIRECTORY "${dir}")
  # DLOAD shares the data reader; its last record is followed by blanks and EOF.
  file(WRITE "${dir}/model.inp"
    "${mesh}*DLOAD\n\n1, P1, 1.0\n \n2, P1, 1.0\n\t\n")
  file(WRITE "${dir}/model.cnt" "!END\n")
  if(mode STREQUAL combined)
    set(control "./model.inp")
  else()
    set(control "model.cnt")
  endif()
  file(WRITE "${dir}/hecmw_ctrl.dat"
    "!MESH, NAME=part_in, TYPE=INP\nmodel.inp\n!MESH, NAME=part_out, TYPE=HECMW-DIST\nmodel.inp\n!CONTROL, NAME=analysis_input\n${control}\n")
  file(WRITE "${dir}/hecmw_part_ctrl.dat"
    "!PARTITION, TYPE=NODE-BASED, METHOD=KMETIS, DOMAIN=2, DEPTH=1\n")
  execute_process(COMMAND "${PROGRAM}" WORKING_DIRECTORY "${dir}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
  file(WRITE "${dir}/partition.log" "${output}${errors}")
  if(NOT "${status}" STREQUAL "0" OR
      NOT "${output}${errors}" MATCHES "Domain decomposition done")
    message(FATAL_ERROR "${mode}: partition failed (${status}); see ${dir}/partition.log")
  endif()
  foreach(pattern "HECMW-IO-INP-W0099" "Auto-generated cards" "NGRP=")
    string(FIND "${output}${errors}" "${pattern}" position)
    if(mode STREQUAL combined AND NOT position EQUAL -1)
      message(FATAL_ERROR "Unexpected ${pattern} in combined input")
    elseif(mode STREQUAL mesh_only AND position EQUAL -1)
      message(FATAL_ERROR "Missing ${pattern} in mesh-only input")
    endif()
  endforeach()
endforeach()
