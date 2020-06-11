/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#ifndef HECMW_ML_HELPER_NN_INCLUDED
#define HECMW_ML_HELPER_NN_INCLUDED

#ifdef HECMW_WITH_ML

#include "ml_include.h"

/*
 * prototype of helper functions in hecmw_ML_helper_f_nn.f90
 */
extern void hecmw_ml_getrow_nn_(int *id, int *n_requested_rows,
                                int *requested_rows, int *allocated_space,
                                int *cols, double *values, int *row_lengths,
                                int *ierr);
extern void hecmw_ml_matvec_nn_(int *id, int *in_length, double *p,
                                int *out_length, double *ap, int *ierr);
extern void hecmw_ml_comm_nn_(int *id, double *x, int *ierr);
extern void hecmw_ml_smoother_diag_setup_nn_(int *id, int *ierr);
extern void hecmw_ml_smoother_diag_apply_nn_(int *id, int *x_length, double x[],
                                             int *rhs_length, double rhs[], int *ierr);
extern void hecmw_ml_smoother_diag_clear_nn_(int *id, int *ierr);
extern void hecmw_ml_smoother_ssor_setup_nn_(int *id, int *ierr);
extern void hecmw_ml_smoother_ssor_apply_nn_(int *id, int *x_length, double x[],
                                             int *rhs_length, double rhs[], int *ierr);
extern void hecmw_ml_smoother_ssor_clear_nn_(int *id, int *ierr);

/*
 * prototype of helper functions in hecmw_ML_helper_nn.c
 */

extern int hecmw_ML_getrow_nn(ML_Operator *mat_in, int N_requested_rows,
                              int requested_rows[], int allocated_space,
                              int cols[], double values[], int row_lengths[]);
extern int hecmw_ML_matvec_nn(ML_Operator *mat_in, int in_length, double p[],
                              int out_length, double ap[]);
extern int hecmw_ML_comm_nn(double x[], void *A_data);
extern int hecmw_ML_smoother_diag_apply_nn(ML_Smoother *data, int x_length, double x[],
                                           int rhs_length, double rhs[]);
extern int hecmw_ML_smoother_ssor_apply_nn(ML_Smoother *data, int x_length, double x[],
                                           int rhs_length, double rhs[]);

#endif /* HECMW_WITH_ML */

#endif /* HECMW_ML_HELPER_NN_INCLUDED */
