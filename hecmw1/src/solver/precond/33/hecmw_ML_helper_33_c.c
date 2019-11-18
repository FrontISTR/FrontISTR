/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifdef HECMW_WITH_ML

#include "hecmw_ML_helper_33.h"

int hecmw_ML_getrow_33(ML_Operator *mat_in, int N_requested_rows,
                       int requested_rows[], int allocated_space,
                       int cols[], double values[], int row_lengths[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_getrow_33_(id, &N_requested_rows, requested_rows, &allocated_space,
                      cols, values, row_lengths, &ierr);
  return ierr;
}

int hecmw_ML_matvec_33(ML_Operator *mat_in, int in_length, double p[],
                       int out_length, double ap[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_matvec_33_(id, &in_length, p, &out_length, ap, &ierr);
  return ierr;
}

int hecmw_ML_comm_33(double x[], void *A_data) {
  int *id, ierr;
  id = (int *)A_data;
  hecmw_ml_comm_33_(id, x, &ierr);
  return ierr;
}

int hecmw_ML_smoother_diag_apply_33(ML_Smoother *data, int x_length, double x[],
                                    int rhs_length, double rhs[]) {
  int *id, ierr;
  id = (int *)ML_Get_MySmootherData(data);
  hecmw_ml_smoother_diag_apply_33_(id, &x_length, x, &rhs_length, rhs, &ierr);
  return ierr;
}

int hecmw_ML_smoother_ssor_apply_33(ML_Smoother *data, int x_length, double x[],
                                    int rhs_length, double rhs[]) {
  int *id, ierr;
  id = (int *)ML_Get_MySmootherData(data);
  hecmw_ml_smoother_ssor_apply_33_(id, &x_length, x, &rhs_length, rhs, &ierr);
  return ierr;
}

#endif /* HECMW_WITH_ML */
