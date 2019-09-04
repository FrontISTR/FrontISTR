/*****************************************************************************
 * Copyright (c) 2019 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_util.h"

#ifdef HECMW_WITH_ML

#include "ml_include.h"

extern void hecmw_ml_getrow_(int *id, int *n_requested_rows,
                             int *requested_rows, int *allocated_space,
                             int *cols, double *values, int *row_lengths,
                             int *ierr);
extern void hecmw_ml_matvec_(int *id, int *in_length, double *p,
                             int *out_length, double *ap, int *ierr);
extern void hecmw_ml_comm_(int *id, double *x, int *ierr);
extern void hecmw_ml_get_nlocal_(int *id, int *nlocal, int *nlocal_allcolumns,
                                 int *ierr);
extern void hecmw_ml_get_coord_(int *id, double x[], double y[], double z[],
                                int *ierr);

static int hecmw_ML_getrow(ML_Operator *mat_in, int N_requested_rows,
                           int requested_rows[], int allocated_space,
                           int cols[], double values[], int row_lengths[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_getrow_(id, &N_requested_rows, requested_rows, &allocated_space,
                   cols, values, row_lengths, &ierr);
  return ierr;
}

static int hecmw_ML_matvec(ML_Operator *mat_in, int in_length, double p[],
                           int out_length, double ap[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_matvec_(id, &in_length, p, &out_length, ap, &ierr);
  return ierr;
}

static int hecmw_ML_comm(double x[], void *A_data) {
  int *id, ierr;
  id = (int *)A_data;
  hecmw_ml_comm_(id, x, &ierr);
  return ierr;
}

struct ml_info {
  ML *ml_object;
  ML_Aggregate *agg_object;
};

#define MAX_MI 8

static struct ml_info mlinfo[MAX_MI];

void hecmw_ML_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  int loglevel;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  ML *ml_object;
  ML_Aggregate *agg_object;

  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  hecmw_ml_get_loglevel_(id, &loglevel);
  ML_Set_PrintLevel(loglevel);

  /* ML object */
  N_grids = 4;
  ML_Create(&ml_object, N_grids);
  hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
  if (*ierr != HECMW_SUCCESS) return;
  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, id);
  ML_Set_Amatrix_Getrow(ml_object, 0, hecmw_ML_getrow, hecmw_ML_comm,
                        nlocal_allcolumns);
  ML_Set_Amatrix_Matvec(ml_object, 0, hecmw_ML_matvec);

  /* if (!(*sym)) ML_Set_Symmetrize(ml_object, ML_YES); */

  /* Aggregate */
  ML_Aggregate_Create(&agg_object);

  /* Null Space (Rigid Body Mode) */
  {
    int num_PDE_eqns = *Ndof;
    int null_dim     = 6;
    double *null_vect;
    int leng  = nlocal;
    null_vect = (double *)HECMW_malloc(sizeof(double) * null_dim * leng);
    if (!null_vect) {
      HECMW_set_error(errno, "");
      abort();
    }
    hecmw_ml_get_rbm_(id, null_vect, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    ML_Aggregate_Set_NullSpace(agg_object, num_PDE_eqns, null_dim, null_vect,
                               leng);
    HECMW_free(null_vect);
  }
  /* ML_Aggregate_Set_MaxCoarseSize(agg_object,1); */ /* default is 32 */

  /* options */
  /* CoarsenScheme : default is UncoupledMIS */
  /*  others are: Coupled, Uncoupled, MIS, UncoupledCoupled, */
  /*              DD, METIS, ParMETIS, Zoltan, User */
  /* ML_Aggregate_Set_CoarsenScheme_Coupled(agg_object); */
  /* ML_Aggregate_Set_CoarsenScheme_MIS(agg_object); */
  /* ML_Aggregate_Set_CoarsenScheme_METIS(agg_object); */
  /* ML_Aggregate_Set_Threshold(agg_object, threshold); */
  /* ML_Aggregate_Set_DampingFactor(agg_object, dampingfactor); */

  if (*sym)
    ML_Set_SpectralNormScheme_Calc(ml_object); /* default is PowerMethod */
  /* ML_Set_SpectralNorm_Iterations(); */      /* ddfault is ten */

  /* repartitioning */
  ML_Repartition_Activate(ml_object);
  /* ML_Repartition_Set_LargestMinMaxRatio(ml_object, 1.3); */ /* default: 1.3
                                                                  */
  /* ML_Repartition_Set_MinPerProc(ml_object, 512); */ /* default: 512 */
  /* ML_Repartition_Set_PutOnSingleProc(ml_object, i); */
  /* ML_Repartition_Set_StartLevel(ml_object, 1); */ /* default: 1 */
  /* ML_Repartition_Set_Partitioner(ml_object, ML_USEPARMETIS); */ /* default:
                                                                      ML_USEZOLTAN
                                                                      */
  /* ML_Aggregate_Set_Dimensions(agg_object, 3); */

  /* Generate MultiGrid */
  /* N_levels = ML_Gen_MGHierarchy_UsingAggregation( */
  /* 	ml_object, 0, ML_INCREASING, agg_object); */
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(
      ml_object, 0, ML_INCREASING, agg_object);
  /* fprintf(stderr, "DEBUG: N_levels = %d\n", N_levels); */

  /* Smoother */
  /* ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, */
  /* 		       ML_PRESMOOTHER, 1, ML_DEFAULT); */
  /* ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, */
  /* 		       ML_BOTH, 1, ML_DEFAULT); */
  /* ML_Gen_Smoother_SymGaussSeidel(ml_object, ML_ALL_LEVELS, */
  /* 		       ML_BOTH, 1, ML_DEFAULT); */
  /* ML_Gen_Smoother_SymGaussSeidel(ml_object, ML_ALL_LEVELS, */
  /* 		       ML_BOTH, 1, 1.0); */
  /* ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, ML_ALL_LEVELS, */
  /* 				    ML_BOTH, 1, ML_DEFAULT, Ndof); */
  /* ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, ML_ALL_LEVELS, */
  /* 				    ML_BOTH, 1, 1.0, Ndof); */
  {
    int level;
    int coarsest_level = N_levels - 1;
    for (level = 0; level < coarsest_level; level++) {
      ML_Gen_Smoother_Cheby(ml_object, level, ML_BOTH, 20.0, 2);
    }
    /* ML_Gen_Smoother_SymGaussSeidel(ml_object, coarsest_level, */
    /* 			       ML_BOTH, 3, ML_DEFAULT); */
    ML_Gen_Smoother_Cheby(ml_object, coarsest_level, ML_BOTH, 20.0, 2);
    /* ML_Gen_Smoother_Amesos(ml_object, coarsest_level, */
    /* 		       ML_AMESOS_KLU, 1, 0.0); */
    /* ML_Gen_Smoother_Amesos(ml_object, coarsest_level, */
    /* 		       ML_AMESOS_MUMPS, 1, 0.0); */
  }

  /* Solver */
  ML_Gen_Solver(ml_object, ML_MGV, 0, N_levels - 1);

  /* Save objects */
  mlinfo[*id - 1].ml_object  = ml_object;
  mlinfo[*id - 1].agg_object = agg_object;
}

void hecmw_ML_wrapper_apply(int *id, double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  double *sol;
  int i;
  ML *ml_object;
  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }
  ml_object = mlinfo[*id - 1].ml_object;
  hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
  if (*ierr != HECMW_SUCCESS) return;
  sol = (double *)HECMW_malloc(sizeof(double) * nlocal_allcolumns);
  if (!sol) {
    HECMW_set_error(errno, "");
    abort();
  }
  /* MultiGrid V-cycle */
  ML_Solve_MGV(ml_object, rhs, sol);
  for (i = 0; i < nlocal; i++) {
    rhs[i] = sol[i];
  }
  HECMW_free(sol);
}

void hecmw_ML_wrapper_clear(int *id, int *ierr) {
  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }
  ML_Aggregate_Destroy(&(mlinfo[*id - 1].agg_object));
  ML_Destroy(&(mlinfo[*id - 1].ml_object));
}

#else /* WITH_ML */

void hecmw_ML_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}
void hecmw_ML_wrapper_apply(int *id, double rhs[], int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}
void hecmw_ML_wrapper_clear(int *id, int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}

#endif /* WITH_ML */
/* Fortran interface */

void hecmw_ml_wrapper_setup_(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_ML_wrapper_setup(id, sym, ndof, ierr);
}
void hecmw_ml_wrapper_setup__(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_ML_wrapper_setup(id, sym, ndof, ierr);
}
void HECMW_ML_WRAPPER_SETUP(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_ML_wrapper_setup(id, sym, ndof, ierr);
}

void hecmw_ml_wrapper_apply_(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply(id, rhs, ierr);
}
void hecmw_ml_wrapper_apply__(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply(id, rhs, ierr);
}
void HECMW_ML_WRAPPER_APPLY(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply(id, rhs, ierr);
}

void hecmw_ml_wrapper_clear_(int *id, int *ierr) {
  hecmw_ML_wrapper_clear(id, ierr);
}
void hecmw_ml_wrapper_clear__(int *id, int *ierr) {
  hecmw_ML_wrapper_clear(id, ierr);
}
void HECMW_ML_WRAPPER_CLEAR(int *id, int *ierr) {
  hecmw_ML_wrapper_clear(id, ierr);
}
