###############################################################################
# Copyright (c) 2026 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################

foreach(case data_line solver_twice parameters unknown_card static_frequency frequency_static
    frequency_matrixstorage frequency_no_data frequency_zero step_warnings nonlinear_solver step_amplitude time_delay
    nlgeom_later load_overlap load_cancelled static_missing op_new_later dload_label dload_element dload_time_delay
    dload_amplitude_steps dload_overlap dload_centrif_absolute dload_plane_thickness dload_centrif_plane_steps
    dsload_label surface_face surface_node_name)
  add_test(NAME test_serial_unit/input/${case}
    COMMAND ${CMAKE_COMMAND}
      -DPROGRAM=$<TARGET_FILE:fistr1>
      -DINPUT=${FrontISTR_SOURCE_DIR}/tests/analysis/inp/standalone_static/standalone_static.inp
      -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${case}
      -DCASE=${case}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/check_inp_input.cmake)
  set_tests_properties(test_serial_unit/input/${case} PROPERTIES
    LABELS unit/input,serial)
endforeach()
