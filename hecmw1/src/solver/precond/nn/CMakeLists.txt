
target_sources(hecmw
  PRIVATE
  hecmw_ML_helper_nn_c.c
  hecmw_ML_wrapper.c
  hecmw_precond_BILU_nn.f90
  hecmw_precond_DIAG_nn.f90
  hecmw_precond_SSOR_nn.f90
  hecmw_precond_SAINV_nn.f90
  hecmw_precond_RIF_nn.f90
  hecmw_ML_helper.f90
  hecmw_ML_helper_nn_f.f90
  hecmw_precond_ML_nn.f90
  hecmw_precond_nn.f90
)

