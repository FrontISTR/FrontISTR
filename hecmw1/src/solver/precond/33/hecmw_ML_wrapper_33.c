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

/*
 * prototype of helper functions in hecmw_ML_helper_33.f90
 */
extern void hecmw_ml_getrow_33_(int *id, int *n_requested_rows,
                                int *requested_rows, int *allocated_space,
                                int *cols, double *values, int *row_lengths,
                                int *ierr);
extern void hecmw_ml_matvec_33_(int *id, int *in_length, double *p,
                                int *out_length, double *ap, int *ierr);
extern void hecmw_ml_comm_33_(int *id, double *x, int *ierr);
extern void hecmw_ml_get_nlocal_33_(int *id, int *nlocal,
                                    int *nlocal_allcolumns, int *ierr);
extern void hecmw_ml_get_coord_33_(int *id, double x[], double y[], double z[],
                                   int *ierr);
extern void hecmw_ml_get_rbm_33_(int *id, double rbm[], int *ierr);
extern void hecmw_ml_get_loglevel_33_(int *id, int *level);
extern void hecmw_ml_smoother_setup_33_(int *id, int *ierr);
extern void hecmw_ml_smoother_apply_33_(int *id, int *x_length, double x[],
                                        int *rhs_length, double rhs[], int *ierr);
extern void hecmw_ml_smoother_clear_33_(int *id, int *ierr);
extern void hecmw_ml_get_opt1_(int *id, int *opt1, int *ierr);
extern void hecmw_ml_get_opt2_(int *id, int *opt2, int *ierr);
extern void hecmw_ml_get_opt3_(int *id, int *opt3, int *ierr);

/*
 * static functions
 */

static int hecmw_ML_getrow_33(ML_Operator *mat_in, int N_requested_rows,
                              int requested_rows[], int allocated_space,
                              int cols[], double values[], int row_lengths[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_getrow_33_(id, &N_requested_rows, requested_rows, &allocated_space,
                      cols, values, row_lengths, &ierr);
  return ierr;
}

static int hecmw_ML_matvec_33(ML_Operator *mat_in, int in_length, double p[],
                              int out_length, double ap[]) {
  int *id, ierr;
  id = (int *)ML_Get_MyGetrowData(mat_in);
  hecmw_ml_matvec_33_(id, &in_length, p, &out_length, ap, &ierr);
  return ierr;
}

static int hecmw_ML_comm_33(double x[], void *A_data) {
  int *id, ierr;
  id = (int *)A_data;
  hecmw_ml_comm_33_(id, x, &ierr);
  return ierr;
}

static int hecmw_ML_smoother_apply_33(ML_Smoother *data, int x_length, double x[],
                                      int rhs_length, double rhs[]) {
  int *id, ierr;
  id = (int *)ML_Get_MySmootherData(data);
  hecmw_ml_smoother_apply_33_(id, &x_length, x, &rhs_length, rhs, &ierr);
  return ierr;
}


/*
 * static variable
 */

struct ml_info {
  ML *ml_object;
  ML_Aggregate *agg_object;
};

#define MAX_MI 8

static struct ml_info MLInfo[MAX_MI];


/*
 * Settings
 */

/* Whether coasest level is solved with direct solver
 * If not, a few sweeps of smoother is applied.
 * (Note: Trilinos must be build with Amesos package enabled)
 */
static int FlgDirectSolveCoarsest = 1;

/* Direct solver for the coarsest level
 *  available types: KLU, MUMPS
 * (KLU is a serial direct solver that comes with Trilinos/Amesos)
 */
enum direct_solver {KLU, MUMPS};
static enum direct_solver DirectSolver = KLU;

/* Smoother type
 *  available types: Jacobi, GaussSeidel, BlockGaussSeidel, SymGaussSeidel,
 *                   SymBlockGaussSeidel, Cheby, Amesos, etc.
 * However, the following three types are currently available from this interface
 */
enum smoother_type {SymGaussSeidel, SymBlockGaussSeidel, Cheby};
static enum smoother_type SmootherType = Cheby;

/* Whether HEC-MW smoother is used at finest level when SmootherType is SymBlockGaussSeidel
 */
static int FlgUseHECMWSmoother = 1;

/* Solver cycle
 *  available types: ML_MGV (V-cycle), ML_MGW (W-cycle), ML_MGFULLV (Full V-Cycle), etc.
 */
static int MGType = ML_MGW;


/*
 * public functions
 */

void hecmw_ML_wrapper_setup_33(int *id, int *sym, int *ierr) {
  int loglevel;
  int N_grids, N_levels;
  int nlocal, nlocal_allcolumns;
  ML *ml_object;
  ML_Aggregate *agg_object;
  int Ndof = 3;

  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  hecmw_ml_get_loglevel_33_(id, &loglevel);
  ML_Set_PrintLevel(loglevel);

  /* Get options */
  {
    int opt1, opt2, opt3;

    hecmw_ml_get_opt1_(id, &opt1, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    switch (opt1) {
    case 1:
      FlgDirectSolveCoarsest = 0;
      if (loglevel > 0) fprintf(stderr, "INFO: ML coarse solver is smoother\n");
      break;
    case 2:
      FlgDirectSolveCoarsest = 1;
      DirectSolver = KLU;
      if (loglevel > 0) fprintf(stderr, "INFO: ML coarse solver is KLU\n");
      break;
    case 3:
      FlgDirectSolveCoarsest = 1;
      DirectSolver = MUMPS;
      if (loglevel > 0) fprintf(stderr, "INFO: ML coarse solver is MUMPS\n");
      break;
    default:
      fprintf(stderr, "WARNING: invalid solver_opt1=%d (ignored)\n", opt1);
    }

    hecmw_ml_get_opt2_(id, &opt2, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    switch (opt2) {
    case 1:
      SmootherType = Cheby;
      if (loglevel > 0) fprintf(stderr, "INFO: ML smoother is Cheby\n");
      break;
    case 2:
      SmootherType = SymBlockGaussSeidel;
      FlgUseHECMWSmoother = 1;
      if (loglevel > 0) fprintf(stderr, "INFO: ML smoother is SymBlockGaussSeidel\n");
      break;
    case 3:
      SmootherType = SymGaussSeidel;
      if (loglevel > 0) fprintf(stderr, "INFO: ML smoother is SymGaussSeidel\n");
      break;
    default:
      fprintf(stderr, "WARNING: invalid solver_opt2=%d (ignored)\n", opt2);
    }

    hecmw_ml_get_opt3_(id, &opt3, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    switch (opt3) {
    case 1:
      MGType = ML_MGV;
      if (loglevel > 0) fprintf(stderr, "INFO: ML multigrid type is V-cycle\n");
      break;
    case 2:
      MGType = ML_MGW;
      if (loglevel > 0) fprintf(stderr, "INFO: ML multigrid type is W-cycle\n");
      break;
    case 3:
      MGType = ML_MGFULLV;
      if (loglevel > 0) fprintf(stderr, "INFO: ML multigrid type is Full-V-cycle\n");
      break;
    default:
      fprintf(stderr, "WARNING: invalid solver_opt3=%d (ignored)\n", opt3);
    }
  }

  /* ML object */
  N_grids = 4;
  ML_Create(&ml_object, N_grids);
  hecmw_ml_get_nlocal_33_(id, &nlocal, &nlocal_allcolumns, ierr);
  if (*ierr != HECMW_SUCCESS) return;
  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, id);
  ML_Set_Amatrix_Getrow(ml_object, 0, hecmw_ML_getrow_33, hecmw_ML_comm_33,
                        nlocal_allcolumns);
  ML_Set_Amatrix_Matvec(ml_object, 0, hecmw_ML_matvec_33);

  /* if (!(*sym)) ML_Set_Symmetrize(ml_object, ML_YES); */

  /* Aggregate */
  ML_Aggregate_Create(&agg_object);

  /* Null Space (Rigid Body Mode) */
  {
    int num_PDE_eqns = Ndof;
    int null_dim     = 6;
    double *null_vect;
    int leng  = nlocal;
    null_vect = (double *)HECMW_malloc(sizeof(double) * null_dim * leng);
    if (!null_vect) {
      HECMW_set_error(errno, "");
      abort();
    }
    hecmw_ml_get_rbm_33_(id, null_vect, ierr);
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
  /*     ml_object, 0, ML_INCREASING, agg_object); */
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(
      ml_object, 0, ML_INCREASING, agg_object);
  /* fprintf(stderr, "DEBUG: N_levels = %d\n", N_levels); */

  /* Smoother */
  /*
   * Set pre- and post-smoother for each level
   *  level      : num in (0, N_levels-1) or ML_ALL_LEVELS
   *  pre-or-post: ML_PRESMOOTHER, ML_POSTSMOOTHER or ML_BOTH
   *  omega      : damping factor for Jacobi, GaussSeidel, etc. (ML_DEFAULT=1.0)
   */
  {
    int level;
    int coarsest_level = N_levels - 1;
    /*
     * levels other than the coarsest level
     */
    if (SmootherType == SymGaussSeidel) {
      for (level = 0; level < coarsest_level; level++) {
        ML_Gen_Smoother_SymGaussSeidel(ml_object, level,
                                       ML_BOTH, 1, ML_DEFAULT);
      }
    } else if (SmootherType == SymBlockGaussSeidel) {
      level = 0;
      if (FlgUseHECMWSmoother) {
        /* use HEC-MW's Block-SSOR preconditioner at the finest level */
        hecmw_ml_smoother_setup_33_(id, ierr);
        if (*ierr != HECMW_SUCCESS) return;
        ML_Set_Smoother(ml_object, 0, ML_BOTH, id, hecmw_ML_smoother_apply_33, "HEC-MW");
        level++;
      }
      /* use ML's smoother at other levels */
      for (; level < coarsest_level; level++) {
        ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, level,
                                            ML_BOTH, 1, ML_DEFAULT, Ndof);
      }
    } else /* if (SmootherType == Cheby) */ {
      for (level = 0; level < coarsest_level; level++) {
        ML_Gen_Smoother_Cheby(ml_object, level, ML_BOTH, 20.0, 2);
      }
    }
    /*
     * coarsest level
     */
    if (FlgDirectSolveCoarsest) {
      if (DirectSolver == MUMPS) {
        ML_Gen_Smoother_Amesos(ml_object, coarsest_level,
                               ML_AMESOS_MUMPS, 1, 0.0);
      } else /* if (DirectSolver == KLU) */ {
        ML_Gen_Smoother_Amesos(ml_object, coarsest_level,
                               ML_AMESOS_KLU, 1, 0.0);
      }
    } else {
      if (SmootherType == SymGaussSeidel) {
        ML_Gen_Smoother_SymGaussSeidel(ml_object, coarsest_level,
                                       ML_BOTH, 3, ML_DEFAULT);
      } else if (SmootherType == SymBlockGaussSeidel) {
        ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, coarsest_level,
                                            ML_BOTH, 3, ML_DEFAULT, Ndof);
      } else /* if (SmootherType == Cheby) */ {
        ML_Gen_Smoother_Cheby(ml_object, coarsest_level, ML_BOTH, 20.0, 2);
      }
    }
  }

  /* Solver */
  ML_Gen_Solver(ml_object, MGType, 0, N_levels - 1);

  /* Save objects */
  MLInfo[*id - 1].ml_object  = ml_object;
  MLInfo[*id - 1].agg_object = agg_object;
}

void hecmw_ML_wrapper_apply_33(int *id, double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  double *sol;
  int i;
  ML *ml_object;
  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }
  ml_object = MLInfo[*id - 1].ml_object;
  hecmw_ml_get_nlocal_33_(id, &nlocal, &nlocal_allcolumns, ierr);
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

void hecmw_ML_wrapper_clear_33(int *id, int *ierr) {
  if (*id <= 0 && MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }
  ML_Aggregate_Destroy(&(MLInfo[*id - 1].agg_object));
  ML_Destroy(&(MLInfo[*id - 1].ml_object));

  if (SmootherType == SymBlockGaussSeidel) {
    if (FlgUseHECMWSmoother) {
      hecmw_ml_smoother_clear_33_(id, ierr);
    }
  }
}

#else /* WITH_ML */

void hecmw_ML_wrapper_setup_33(int *id, int *sym, int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}
void hecmw_ML_wrapper_apply_33(int *id, double rhs[], int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}
void hecmw_ML_wrapper_clear_33(int *id, int *ierr) {
  fprintf(stderr, "ERROR: ML not enabled\n");
  *ierr = HECMW_ERROR;
}

#endif /* WITH_ML */
/* Fortran interface */

void hecmw_ml_wrapper_setup_33_(int *id, int *sym, int *ierr) {
  hecmw_ML_wrapper_setup_33(id, sym, ierr);
}
void hecmw_ml_wrapper_setup_33__(int *id, int *sym, int *ierr) {
  hecmw_ML_wrapper_setup_33(id, sym, ierr);
}
void HECMW_ML_WRAPPER_SETUP_33(int *id, int *sym, int *ierr) {
  hecmw_ML_wrapper_setup_33(id, sym, ierr);
}

void hecmw_ml_wrapper_apply_33_(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply_33(id, rhs, ierr);
}
void hecmw_ml_wrapper_apply_33__(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply_33(id, rhs, ierr);
}
void HECMW_ML_WRAPPER_APPLY_33(int *id, double rhs[], int *ierr) {
  hecmw_ML_wrapper_apply_33(id, rhs, ierr);
}

void hecmw_ml_wrapper_clear_33_(int *id, int *ierr) {
  hecmw_ML_wrapper_clear_33(id, ierr);
}
void hecmw_ml_wrapper_clear_33__(int *id, int *ierr) {
  hecmw_ML_wrapper_clear_33(id, ierr);
}
void HECMW_ML_WRAPPER_CLEAR_33(int *id, int *ierr) {
  hecmw_ML_wrapper_clear_33(id, ierr);
}
