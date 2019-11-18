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
#include "hecmw_ML_helper.h"
#include "hecmw_ML_helper_33.h"
#include "hecmw_ML_helper_nn.h"

/*
 * static variable
 */

struct ml_info {
  ML *ml_object;
  ML_Aggregate *agg_object;
  int ndof;
};

#define MAX_MI 8

static struct ml_info MLInfo[MAX_MI];


/*
 * Settings
 */

/* Whether coasest level is solved with direct solver
 * If not, a few sweeps of smoother is applied.
 * (Note: Trilinos must be built with Amesos package enabled)
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
enum smoother_type {Jacobi, SymBlockGaussSeidel, Cheby};
static enum smoother_type SmootherType = Cheby;

/* Whether HEC-MW smoother is used at finest level when SmootherType is SymBlockGaussSeidel
 */
static int FlgUseHECMWSmoother = 1;

/* Solver cycle
 *  available types: ML_MGV (V-cycle), ML_MGW (W-cycle), ML_MGFULLV (Full V-Cycle), etc.
 */
static int MGType = ML_MGW;

/* Max num of levels
 */
static int MaxLevels = 10;

/* Coarsening scheme
 *  available types: Coupled, Uncoupled, MIS, UncoupledCoupled, UncoupledMIS,
 *                   DD, METIS, ParMETIS, Zoltan, User
 * However, the following three types are currently available from this interface
 */
enum coarsen_scheme {UncoupledMIS, METIS, ParMETIS, Zoltan, DD};
static enum coarsen_scheme CoarsenScheme = UncoupledMIS;

/* Num of smoother sweeps
 */
static int NumSweeps = 2;


/*
 * public functions
 */

void hecmw_ML_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank;
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

  HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);

  /* Get options */
  {
    int opt1, opt2, opt3, opt4, opt5, opt6;

    hecmw_ml_get_opt1_(id, &opt1, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    switch (opt1) {
    case 0: /* default */
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default coarse solver\n");
    case 1:
      FlgDirectSolveCoarsest = 0;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarse solver is smoother\n");
      break;
    case 2:
      FlgDirectSolveCoarsest = 1;
      DirectSolver = KLU;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarse solver is KLU\n");
      break;
    case 3:
      FlgDirectSolveCoarsest = 1;
      DirectSolver = MUMPS;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarse solver is MUMPS\n");
      break;
    default:
      if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt1=%d (ignored)\n", opt1);
    }

    hecmw_ml_get_opt2_(id, &opt2, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    switch (opt2) {
    case 0: /* default */
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default smoother\n");
    case 1:
      SmootherType = Cheby;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML smoother is Cheby\n");
      break;
    case 2:
      SmootherType = SymBlockGaussSeidel;
      FlgUseHECMWSmoother = 1;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML smoother is SymBlockGaussSeidel\n");
      break;
    case 3:
      SmootherType = Jacobi;
      FlgUseHECMWSmoother = 1;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML smoother is Jacobi\n");
      break;
    default:
      if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt2=%d (ignored)\n", opt2);
    }

    hecmw_ml_get_opt3_(id, &opt3, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    switch (opt3) {
    case 0: /* default */
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default multigrid type\n");
    case 1:
      MGType = ML_MGV;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML multigrid type is V-cycle\n");
      break;
    case 2:
      MGType = ML_MGW;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML multigrid type is W-cycle\n");
      break;
    case 3:
      MGType = ML_MGFULLV;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML multigrid type is Full-V-cycle\n");
      break;
    default:
      if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt3=%d (ignored)\n", opt3);
    }

    hecmw_ml_get_opt4_(id, &opt4, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    if (opt4 > 0) {
      MaxLevels = opt4;
    } else {
      if (opt4 < 0) {
        if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt4=%d (ignored)\n", opt4);
      }
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default MaxLevels\n");
      MaxLevels = 10; /* default */
    }
    if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML num of max levels is %d\n", MaxLevels);

    hecmw_ml_get_opt5_(id, &opt5, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    switch (opt5) {
    case 0: /* default */
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default coarsening scheme\n");
    case 1:
      CoarsenScheme = UncoupledMIS;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarsening scheme is UncoupledMIS\n");
      break;
    case 2:
      CoarsenScheme = METIS;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarsening scheme is METIS\n");
      break;
    case 3:
      CoarsenScheme = ParMETIS;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarsening scheme is ParMETIS\n");
      break;
    case 4:
      CoarsenScheme = Zoltan;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarsening scheme is Zoltan\n");
      break;
    case 5:
      CoarsenScheme = DD;
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML coarsening scheme is DD\n");
      break;
    default:
      if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt5=%d (ignored)\n", opt5);
    }

    hecmw_ml_get_opt6_(id, &opt6, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    if (opt6 > 0) {
      NumSweeps = opt6;
    } else {
      if (opt6 < 0) {
        if (myrank == 0) fprintf(stderr, "WARNING: invalid solver_opt6=%d (ignored)\n", opt6);
      }
      if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML using default num of sweeps\n");
      NumSweeps = 2; /* default */
      hecmw_ml_set_opt6_(id, &NumSweeps, ierr);
      if (*ierr != HECMW_SUCCESS) return;
    }
    if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML num of smoother sweeps is %d\n", NumSweeps);

  }

  /* ML object */
  N_grids = MaxLevels;
  ML_Create(&ml_object, N_grids);
  hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
  if (*ierr != HECMW_SUCCESS) return;
  ML_Init_Amatrix(ml_object, 0, nlocal, nlocal, id);
  if (*Ndof == 3) {
    ML_Set_Amatrix_Getrow(ml_object, 0, hecmw_ML_getrow_33, hecmw_ML_comm_33,
                          nlocal_allcolumns);
    ML_Set_Amatrix_Matvec(ml_object, 0, hecmw_ML_matvec_33);
  } else {
    ML_Set_Amatrix_Getrow(ml_object, 0, hecmw_ML_getrow_nn, hecmw_ML_comm_nn,
                          nlocal_allcolumns);
    ML_Set_Amatrix_Matvec(ml_object, 0, hecmw_ML_matvec_nn);
  }

  /* if (!(*sym)) ML_Set_Symmetrize(ml_object, ML_YES); */

  /* Aggregate */
  ML_Aggregate_Create(&agg_object);

  /* Null Space (Rigid Body Mode) */
  {
    int num_PDE_eqns = *Ndof;
    int null_dim;
    double *null_vect;
    int leng  = nlocal;
    if (*Ndof == 1) {
      null_dim = 1;
    } else if (*Ndof == 2) {
      null_dim = 3;
    } else {
      null_dim = 6;
    }
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
  /* ML_Aggregate_Set_MaxCoarseSize(agg_object, 128); */ /* default: 128 */

  /* options */
  /* CoarsenScheme */
  {
    if (CoarsenScheme == UncoupledMIS) {
      ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);
    } else if (CoarsenScheme == METIS) {
      ML_Aggregate_Set_CoarsenScheme_METIS(agg_object);
    } else if (CoarsenScheme == ParMETIS) {
      ML_Aggregate_Set_CoarsenScheme_ParMETIS(agg_object);
    } else if (CoarsenScheme == Zoltan) {
      ML_Aggregate_Set_CoarsenScheme_Zoltan(agg_object);
    } else if (CoarsenScheme == DD) {
      ML_Aggregate_Set_CoarsenScheme_DD(agg_object);
    }
    /*
    if (MaxLevels == 2) {
      ML_Aggregate_Set_LocalNumber(ml_object, agg_object, ML_ALL_LEVELS, 1);
    } else if (MaxLevels == 3) {
      ML_Aggregate_Set_NodesPerAggr(ml_object, agg_object, ML_ALL_LEVELS, 512);
      ML_Aggregate_Set_ReqLocalCoarseSize(ml_object->ML_num_levels, agg_object, ML_ALL_LEVELS, 128);
    }
    */
  }
  /* ML_Aggregate_Set_Threshold(agg_object, threshold); */
  /* ML_Aggregate_Set_DampingFactor(agg_object, dampingfactor); */

  /* eigen-analysis */
  /* ML_Set_SpectralNormScheme_PowerMethod(ml_object); */ /* power-method (default) */
  if (*sym) {
    ML_Set_SpectralNormScheme_Calc(ml_object); /* cg */
  }
  /* ML_Set_SpectralNorm_Iterations(ml_object, 10); */ /* default: 10 */

  /* repartitioning */
  /* ML_Repartition_Activate(ml_object); */
  /* ML_Repartition_Set_LargestMinMaxRatio(ml_object, 1.3); */ /* default: 1.3 */
  /* ML_Repartition_Set_MinPerProc(ml_object, 512); */ /* default: 512 */
  /* ML_Repartition_Set_PutOnSingleProc(ml_object, i); */
  /* ML_Repartition_Set_StartLevel(ml_object, 1); */ /* default: 1 */
  /* ML_Repartition_Set_Partitioner(ml_object, ML_USEPARMETIS); */ /* default: ML_USEZOLTAN */

  ML_Aggregate_Set_Dimensions(agg_object, *Ndof);

  /* Generate MultiGrid */
  /* N_levels = ML_Gen_MGHierarchy_UsingAggregation( */
  /*     ml_object, 0, ML_INCREASING, agg_object); */
  N_levels = ML_Gen_MultiLevelHierarchy_UsingAggregation(
      ml_object, 0, ML_INCREASING, agg_object);
  if (loglevel > 0 && myrank == 0) fprintf(stderr, "INFO: ML generated num of levels is %d\n", N_levels);

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
    if (SmootherType == Jacobi) {
      level = 0;
      if (FlgUseHECMWSmoother) {
        /* use HEC-MW's Block-Diag preconditioner at the finest level */
        if (*Ndof == 3) {
          hecmw_ml_smoother_diag_setup_33_(id, ierr);
          if (*ierr != HECMW_SUCCESS) return;
          ML_Set_Smoother(ml_object, 0, ML_BOTH, id, hecmw_ML_smoother_diag_apply_33, "HEC-MW");
        } else {
          hecmw_ml_smoother_diag_setup_nn_(id, ierr);
          if (*ierr != HECMW_SUCCESS) return;
          ML_Set_Smoother(ml_object, 0, ML_BOTH, id, hecmw_ML_smoother_diag_apply_nn, "HEC-MW");
        }
        level++;
      }
      /* use ML's smoother at other levels */
      for (; level < coarsest_level; level++) {
        ML_Gen_Smoother_Jacobi(ml_object, level,
                               ML_BOTH, NumSweeps, ML_DEFAULT);
      }
    } else if (SmootherType == SymBlockGaussSeidel) {
      level = 0;
      if (FlgUseHECMWSmoother) {
        /* use HEC-MW's Block-SSOR preconditioner at the finest level */
        if (*Ndof == 3) {
          hecmw_ml_smoother_ssor_setup_33_(id, ierr);
          if (*ierr != HECMW_SUCCESS) return;
          ML_Set_Smoother(ml_object, 0, ML_BOTH, id, hecmw_ML_smoother_ssor_apply_33, "HEC-MW");
        } else {
          hecmw_ml_smoother_ssor_setup_nn_(id, ierr);
          if (*ierr != HECMW_SUCCESS) return;
          ML_Set_Smoother(ml_object, 0, ML_BOTH, id, hecmw_ML_smoother_ssor_apply_nn, "HEC-MW");
        }
        level++;
      }
      /* use ML's smoother at other levels */
      for (; level < coarsest_level; level++) {
        ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, level,
                                            ML_BOTH, NumSweeps, ML_DEFAULT, *Ndof);
      }
    } else /* if (SmootherType == Cheby) */ {
      for (level = 0; level < coarsest_level; level++) {
        ML_Gen_Smoother_Cheby(ml_object, level, ML_BOTH, 20.0, NumSweeps);
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
      if (SmootherType == Jacobi) {
        ML_Gen_Smoother_Jacobi(ml_object, coarsest_level,
                               ML_BOTH, 3, ML_DEFAULT);
      } else if (SmootherType == SymBlockGaussSeidel) {
        ML_Gen_Smoother_SymBlockGaussSeidel(ml_object, coarsest_level,
                                            ML_BOTH, 3, ML_DEFAULT, *Ndof);
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
  MLInfo[*id - 1].ndof = *Ndof;
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
  ml_object = MLInfo[*id - 1].ml_object;
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
  ML_Aggregate_Destroy(&(MLInfo[*id - 1].agg_object));
  ML_Destroy(&(MLInfo[*id - 1].ml_object));

  if (FlgUseHECMWSmoother) {
    if (SmootherType == Jacobi) {
      if (MLInfo[*id - 1].ndof == 3) {
        hecmw_ml_smoother_diag_clear_33_(id, ierr);
      } else {
        hecmw_ml_smoother_diag_clear_nn_(id, ierr);
      }
    } else if (SmootherType == SymBlockGaussSeidel) {
      if (MLInfo[*id - 1].ndof == 3) {
        hecmw_ml_smoother_ssor_clear_33_(id, ierr);
      } else {
        hecmw_ml_smoother_ssor_clear_nn_(id, ierr);
      }
    }
  }
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
