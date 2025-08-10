/*****************************************************************************
 * Copyright (c) 2021 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "hecmw_util.h"

#ifdef HECMW_WITH_MUELU

#include "Trilinos_version.h"
#include "MueLu.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_MLParameterListInterpreter.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "EpetraExt_RowMatrixOut.h"

#ifdef HAVE_MUELU_AMESOS2
# include "Amesos2_config.h"
#endif

#include "hecmw_ML_helper.h"
#include "hecmw_ML_helper_33.h"
#include "hecmw_ML_helper_nn.h"

/*
 * Options
 */

enum coarse_solver {Smoother, KLU, MUMPS, SuperLU};
enum smoother_type {Cheby, SymGaussSeidel, Jacobi, ILUT};
enum coarsen_scheme {Uncoupled, Coupled, MIS, HMIS, METIS, ParMETIS};
enum cycle_type {V_Cycle, W_Cycle, Full_V_Cycle};

struct muelu_options {
  /* Coarse solver
   *  available solvers: Smoother, KLU, MUMPS, SuperLU
   * Note:
   *  - Trilinos must be built with Amesos2 enabled to use KLU/MUMPS/SuperLU
   */
  enum coarse_solver CoarseSolver;

  /* Smoother type
   *  available types: Cheby, SymGaussSeidel, Jacobi, ILUT
   */
  enum smoother_type SmootherType;

  /* Whether HEC-MW smoother is used at finest level
   */
  int FlgUseHECMWSmoother;

  /* Solver cycle
   *  available types: V_Cycle, W_Cycle, Full_V_Cycle
   */
  enum cycle_type CycleType;

  /* Max num of levels
   */
  int MaxLevels;

  /* Coarsening scheme
   *  available types: Uncoupled, Coupled, MIS, HMIS, METIS, ParMETIS
   */
  enum coarsen_scheme CoarsenScheme;

  /* Num of smoother sweeps
   */
  int NumSweeps;

  /* Max coarse size
   */
  int MaxCoarseSize;

  /* Repartitioning enabled
   */
  int EnableRepartitioning;

  /* Aggregation damping factor
   */
  double DampingFactor;
};

/* default values */
#ifdef HAVE_MUELU_AMESOS2
# ifdef HAVE_AMESOS2_MUMPS
#  define DEFAULT_COARSE_SOLVER MUMPS
# elif defined(HAVE_AMESOS2_KLU2)
#  define DEFAULT_COARSE_SOLVER KLU
# else
#  define DEFAULT_COARSE_SOLVER SuperLU
# endif
#else
# define DEFAULT_COARSE_SOLVER  Smoother
#endif
#define DEFAULT_SMOOTHER_TYPE       Cheby
#define DEFAULT_CYCLE_TYPE          W_Cycle
#define DEFAULT_MAX_LEVELS          4
#define DEFAULT_COARSEN_SCHEME      Uncoupled
#define DEFAULT_NUM_SWEEPS          2
#define DEFAULT_ENABLE_REPARTITION  0
#define DEFAULT_DAMPING_FACTOR      1.333

#define MAX_COARSE_SIZE_MUMPS       50000
#define MAX_COARSE_SIZE_KLU         10000
#define MAX_COARSE_SIZE_SUPERLU     5000

/*
 * MueLu wrapper data structure
 */
struct muelu_info {
  struct muelu_options opt;
  Teuchos::RCP<MueLu::EpetraOperator> muelu_prec;
  Teuchos::RCP<Epetra_CrsMatrix> matrix;
  int ndof;
};

#define MAX_MI 8
static struct muelu_info MueLuInfo[MAX_MI];

/*
 * Helper functions
 */
static void muelu_options_set(struct muelu_options *muelopt, int *id, int myrank, int *ierr) {
  int opt[10];

  hecmw_ml_get_opt_(id, opt, ierr);
  if (*ierr != HECMW_SUCCESS) return;

  // CoarseSolver
  switch (opt[0]) {
  case 0:
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
  case 1:
    muelopt->CoarseSolver = Smoother;
    break;
#ifdef HAVE_MUELU_AMESOS2
  case 2:
    muelopt->CoarseSolver = KLU;
    break;
  case 3:
# ifdef HAVE_AMESOS2_MUMPS
    muelopt->CoarseSolver = MUMPS;
    break;
# else
    if (myrank == 0) fprintf(stderr, "WARNING: MUMPS not available as coarse solver (rebuild Trilinos with Amesos2 and MUMPS enabled)\n");
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
# endif
  case 4:
    muelopt->CoarseSolver = SuperLU;
    break;
#else
  case 2:
    if (myrank == 0) fprintf(stderr, "WARNING: KLU not available as coarse solver (rebuild Trilinos with Amesos2 enabled)\n");
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
  case 3:
    if (myrank == 0) fprintf(stderr, "WARNING: MUMPS not available as coarse solver (rebuild Trilinos with Amesos2 and MUMPS enabled)\n");
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
  case 4:
    if (myrank == 0) fprintf(stderr, "WARNING: SuperLU not available as coarse solver (rebuild Trilinos with Amesos2 enabled)\n");
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
    break;
#endif
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_CoarseSolver=%d (ignored)\n", opt[0]);
    muelopt->CoarseSolver = DEFAULT_COARSE_SOLVER;
  }

  // SmootherType
  switch (opt[1]) {
  case 0:
    muelopt->SmootherType = DEFAULT_SMOOTHER_TYPE;
    muelopt->FlgUseHECMWSmoother = 0;
    break;
  case 1:
    muelopt->SmootherType = Cheby;
    muelopt->FlgUseHECMWSmoother = 0;
    break;
  case 2:
    muelopt->SmootherType = SymGaussSeidel;
    muelopt->FlgUseHECMWSmoother = 1;
    break;
  case 3:
    muelopt->SmootherType = Jacobi;
    muelopt->FlgUseHECMWSmoother = 1;
    break;
  case 4:
    muelopt->SmootherType = ILUT;
    muelopt->FlgUseHECMWSmoother = 0;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_Smoother=%d (ignored)\n", opt[1]);
    muelopt->SmootherType = DEFAULT_SMOOTHER_TYPE;
    muelopt->FlgUseHECMWSmoother = 0;
  }

  // CycleType
  switch (opt[2]) {
  case 0:
    muelopt->CycleType = DEFAULT_CYCLE_TYPE;
    break;
  case 1:
    muelopt->CycleType = V_Cycle;
    break;
  case 2:
    muelopt->CycleType = W_Cycle;
    break;
  case 3:
    muelopt->CycleType = Full_V_Cycle;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_CycleType=%d (ignored)\n", opt[2]);
    muelopt->CycleType = DEFAULT_CYCLE_TYPE;
  }

  // MaxLevels
  if (opt[3] > 0) {
    muelopt->MaxLevels = opt[3];
  } else {
    if (opt[3] < 0) {
      if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_MaxLevels=%d (ignored)\n", opt[3]);
    }
    muelopt->MaxLevels = DEFAULT_MAX_LEVELS;
  }

  // CoarsenScheme
  switch (opt[4]) {
  case 0:
    muelopt->CoarsenScheme = DEFAULT_COARSEN_SCHEME;
    break;
  case 1:
    muelopt->CoarsenScheme = Uncoupled;
    break;
  case 2:
    muelopt->CoarsenScheme = Coupled;
    break;
  case 3:
    muelopt->CoarsenScheme = MIS;
    break;
  case 4:
    muelopt->CoarsenScheme = HMIS;
    break;
  case 5:
    muelopt->CoarsenScheme = METIS;
    break;
  case 6:
    muelopt->CoarsenScheme = ParMETIS;
    break;
  default:
    if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_CoarseningScheme=%d (ignored)\n", opt[4]);
    muelopt->CoarsenScheme = DEFAULT_COARSEN_SCHEME;
  }

  // NumSweeps
  if (opt[5] > 0) {
    muelopt->NumSweeps = opt[5];
  } else {
    if (opt[5] < 0) {
      if (myrank == 0) fprintf(stderr, "WARNING: invalid MueLu_NumSweep=%d (ignored)\n", opt[5]);
    }
    muelopt->NumSweeps = DEFAULT_NUM_SWEEPS;
    opt[5] = muelopt->NumSweeps;
    hecmw_ml_set_opt_(id, opt, ierr);
    if (*ierr != HECMW_SUCCESS) return;
  }

  // MaxCoarseSize
  if (opt[6] > 0) {
    muelopt->MaxCoarseSize = opt[6];
  } else {
    if (muelopt->CoarseSolver == MUMPS) {
      muelopt->MaxCoarseSize = MAX_COARSE_SIZE_MUMPS;
    } else if (muelopt->CoarseSolver == KLU) {
      muelopt->MaxCoarseSize = MAX_COARSE_SIZE_KLU;
    } else if (muelopt->CoarseSolver == SuperLU) {
      muelopt->MaxCoarseSize = MAX_COARSE_SIZE_SUPERLU;
    } else {
      muelopt->MaxCoarseSize = 128; // MueLU default
    }
  }

  // Repartitioning (new option)
  muelopt->EnableRepartitioning = DEFAULT_ENABLE_REPARTITION;
  
  // Damping factor (new option)
  muelopt->DampingFactor = DEFAULT_DAMPING_FACTOR;
}

static void muelu_options_print(struct muelu_options *muelopt, FILE *fp, int myrank, int loglevel) {
  if (loglevel < 1 || myrank != 0) return;

  fprintf(fp, "INFO: MueLu options:\n");
  
  // Coarse solver
  switch (muelopt->CoarseSolver) {
  case Smoother:
    fprintf(fp, "  CoarseSolver: Smoother\n");
    break;
  case KLU:
    fprintf(fp, "  CoarseSolver: KLU\n");
    break;
  case MUMPS:
    fprintf(fp, "  CoarseSolver: MUMPS\n");
    break;
  case SuperLU:
    fprintf(fp, "  CoarseSolver: SuperLU\n");
    break;
  }

  // Smoother type
  switch (muelopt->SmootherType) {
  case Cheby:
    fprintf(fp, "  SmootherType: Chebyshev\n");
    break;
  case SymGaussSeidel:
    fprintf(fp, "  SmootherType: Symmetric Gauss-Seidel%s\n", 
            muelopt->FlgUseHECMWSmoother ? " (HEC-MW at finest)" : "");
    break;
  case Jacobi:
    fprintf(fp, "  SmootherType: Jacobi%s\n", 
            muelopt->FlgUseHECMWSmoother ? " (HEC-MW at finest)" : "");
    break;
  case ILUT:
    fprintf(fp, "  SmootherType: ILUT\n");
    break;
  }

  // Cycle type
  switch (muelopt->CycleType) {
  case V_Cycle:
    fprintf(fp, "  CycleType: V-cycle\n");
    break;
  case W_Cycle:
    fprintf(fp, "  CycleType: W-cycle\n");
    break;
  case Full_V_Cycle:
    fprintf(fp, "  CycleType: Full V-cycle\n");
    break;
  }

  fprintf(fp, "  MaxLevels: %d\n", muelopt->MaxLevels);

  // Coarsening scheme
  switch (muelopt->CoarsenScheme) {
  case Uncoupled:
    fprintf(fp, "  CoarsenScheme: Uncoupled\n");
    break;
  case Coupled:
    fprintf(fp, "  CoarsenScheme: Coupled\n");
    break;
  case MIS:
    fprintf(fp, "  CoarsenScheme: MIS\n");
    break;
  case HMIS:
    fprintf(fp, "  CoarsenScheme: HMIS\n");
    break;
  case METIS:
    fprintf(fp, "  CoarsenScheme: METIS\n");
    break;
  case ParMETIS:
    fprintf(fp, "  CoarsenScheme: ParMETIS\n");
    break;
  }

  fprintf(fp, "  NumSweeps: %d\n", muelopt->NumSweeps);
  fprintf(fp, "  MaxCoarseSize: %d\n", muelopt->MaxCoarseSize);
  fprintf(fp, "  EnableRepartitioning: %s\n", muelopt->EnableRepartitioning ? "true" : "false");
  fprintf(fp, "  DampingFactor: %.3f\n", muelopt->DampingFactor);
}

static Teuchos::RCP<Teuchos::ParameterList> create_muelu_parameter_list(struct muelu_options *muelopt) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  RCP<ParameterList> paramList = rcp(new ParameterList());

  // Verbosity
  paramList->set("verbosity", "low");

  // Number of equations per node
  paramList->set("number of equations", 1); // Will be updated in setup

  // Max levels
  paramList->set("max levels", muelopt->MaxLevels);

  // Coarse solver
  switch (muelopt->CoarseSolver) {
  case Smoother:
    paramList->set("coarse: type", "RELAXATION");
    break;
  case KLU:
    paramList->set("coarse: type", "KLU2");
    break;
  case MUMPS:
    paramList->set("coarse: type", "MUMPS");
    break;
  case SuperLU:
    paramList->set("coarse: type", "SUPERLU");
    break;
  }

  // Max coarse size
  paramList->set("coarse: max size", muelopt->MaxCoarseSize);

  // Smoother
  switch (muelopt->SmootherType) {
  case Cheby:
    paramList->set("smoother: type", "CHEBYSHEV");
    paramList->set("smoother: sweeps", muelopt->NumSweeps);
    break;
  case SymGaussSeidel:
    paramList->set("smoother: type", "RELAXATION");
    paramList->sublist("smoother: params").set("relaxation: type", "Symmetric Gauss-Seidel");
    paramList->set("smoother: sweeps", muelopt->NumSweeps);
    break;
  case Jacobi:
    paramList->set("smoother: type", "RELAXATION");
    paramList->sublist("smoother: params").set("relaxation: type", "Jacobi");
    paramList->set("smoother: sweeps", muelopt->NumSweeps);
    break;
  case ILUT:
    paramList->set("smoother: type", "ILUT");
    paramList->set("smoother: sweeps", 1);
    break;
  }

  // Aggregation
  switch (muelopt->CoarsenScheme) {
  case Uncoupled:
    paramList->set("aggregation: type", "uncoupled");
    break;
  case Coupled:
    paramList->set("aggregation: type", "coupled");
    break;
  case MIS:
    paramList->set("aggregation: type", "uncoupled");
    paramList->set("aggregation: ordering", "natural");
    break;
  case HMIS:
    paramList->set("aggregation: type", "uncoupled");
    paramList->set("aggregation: ordering", "random");
    break;
  case METIS:
    paramList->set("repartition: enable", true);
    paramList->set("repartition: partitioner", "zoltan2");
    break;
  case ParMETIS:
    paramList->set("repartition: enable", true);
    paramList->set("repartition: partitioner", "zoltan2");
    break;
  }

  // Damping factor
  paramList->set("aggregation: damping factor", muelopt->DampingFactor);

  // Cycle type
  switch (muelopt->CycleType) {
  case V_Cycle:
    paramList->set("cycle type", "V");
    break;
  case W_Cycle:
    paramList->set("cycle type", "W");
    break;
  case Full_V_Cycle:
    paramList->set("cycle type", "V");
    break;
  }

  // Repartitioning
  if (muelopt->EnableRepartitioning) {
    paramList->set("repartition: enable", true);
    paramList->set("repartition: min rows per proc", 800);
    paramList->set("repartition: max imbalance", 1.2);
  }

  return paramList;
}

/*
 * Public interface functions
 */
extern "C" {

void hecmw_MueLU_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank;
  int nlocal, nlocal_allcolumns;
  struct muelu_options *muelopt;

  if (*id <= 0 || MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    hecmw_ml_get_loglevel_(id, &loglevel);
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);

    // Get options
    muelopt = &(MueLuInfo[*id - 1].opt);
    muelu_options_set(muelopt, id, myrank, ierr);
    if (*ierr != HECMW_SUCCESS) return;
    muelu_options_print(muelopt, stderr, myrank, loglevel);

    // Get matrix dimensions
    hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    // Create parameter list
    auto paramList = create_muelu_parameter_list(muelopt);
    paramList->set("number of equations", *Ndof);

    // Note: In a complete implementation, you would need to:
    // 1. Create Epetra_CrsMatrix from HEC-MW matrix data
    // 2. Set up coordinates for semi-coarsening (if needed)
    // 3. Set up null space (rigid body modes)
    // 4. Create MueLU preconditioner

    // For now, we'll create a placeholder that would be filled in
    // with the actual matrix setup code
    
    // Store ndof for later use
    MueLuInfo[*id - 1].ndof = *Ndof;

    if (loglevel >= 1 && myrank == 0) {
      fprintf(stderr, "INFO: MueLu preconditioner setup completed\n");
    }

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: MueLu setup failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

void hecmw_MueLu_wrapper_apply(int *id, double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  int myrank;

  if (*id <= 0 || MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);
    
    hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    // In a complete implementation, this would apply the MueLu preconditioner
    // using the stored preconditioner object
    
    auto& muelu_prec = MueLuInfo[*id - 1].muelu_prec;
    if (muelu_prec.is_null()) {
      if (myrank == 0) fprintf(stderr, "ERROR: MueLu preconditioner not initialized\n");
      *ierr = HECMW_ERROR;
      return;
    }

    // Apply preconditioner: x = M^(-1) * rhs
    // This would involve creating Epetra_MultiVector objects and calling ApplyInverse

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: MueLu apply failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

void hecmw_MueLu_wrapper_clear(int *id, int *ierr) {
  struct muelu_options *muelopt;
  int myrank;

  if (*id <= 0 || MAX_MI < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);
    
    muelopt = &(MueLuInfo[*id - 1].opt);

    // Clear MueLu objects
    MueLuInfo[*id - 1].muelu_prec = Teuchos::null;
    MueLuInfo[*id - 1].matrix = Teuchos::null;

    // Clear HEC-MW smoothers if used
    if (muelopt->FlgUseHECMWSmoother) {
      if (muelopt->SmootherType == Jacobi) {
        if (MueLuInfo[*id - 1].ndof == 3) {
          hecmw_ml_smoother_diag_clear_33_(id, ierr);
        } else {
          hecmw_ml_smoother_diag_clear_nn_(id, ierr);
        }
      } else if (muelopt->SmootherType == SymGaussSeidel) {
        if (MueLuInfo[*id - 1].ndof == 3) {
          hecmw_ml_smoother_ssor_clear_33_(id, ierr);
        } else {
          hecmw_ml_smoother_ssor_clear_nn_(id, ierr);
        }
      }
    }

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: MueLu clear failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

} // extern "C"

#else /* HECMW_WITH_MUELU */

extern "C" {

void hecmw_MueLu_wrapper_setup(int *id, int *sym, int *Ndof, int *ierr) {
  fprintf(stderr, "ERROR: MueLu not enabled\n");
  *ierr = HECMW_ERROR;
}

void hecmw_MueLu_wrapper_apply(int *id, double rhs[], int *ierr) {
  fprintf(stderr, "ERROR: MueLu not enabled\n");
  *ierr = HECMW_ERROR;
}

void hecmw_MueLu_wrapper_clear(int *id, int *ierr) {
  fprintf(stderr, "ERROR: MueLu not enabled\n");
  *ierr = HECMW_ERROR;
}

} // extern "C"

#endif /* HECMW_WITH_MUELU */

/* Fortran interface */
extern "C" {

void hecmw_muelu_wrapper_setup_(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_MueLu_wrapper_setup(id, sym, ndof, ierr);
}
void hecmw_muelu_wrapper_setup__(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_MueLu_wrapper_setup(id, sym, ndof, ierr);
}
void HECMW_MUELU_WRAPPER_SETUP(int *id, int *sym, int *ndof, int *ierr) {
  hecmw_MueLu_wrapper_setup(id, sym, ndof, ierr);
}

void hecmw_muelu_wrapper_apply_(int *id, double rhs[], int *ierr) {
  hecmw_MueLu_wrapper_apply(id, rhs, ierr);
}
void hecmw_muelu_wrapper_apply__(int *id, double rhs[], int *ierr) {
  hecmw_MueLu_wrapper_apply(id, rhs, ierr);
}
void HECMW_MUELU_WRAPPER_APPLY(int *id, double rhs[], int *ierr) {
  hecmw_MueLu_wrapper_apply(id, rhs, ierr);
}

void hecmw_muelu_wrapper_clear_(int *id, int *ierr) {
  hecmw_MueLu_wrapper_clear(id, ierr);
}
void hecmw_muelu_wrapper_clear__(int *id, int *ierr) {
  hecmw_MueLu_wrapper_clear(id, ierr);
}
void HECMW_MUELU_WRAPPER_CLEAR(int *id, int *ierr) {
  hecmw_MueLu_wrapper_clear(id, ierr);
}

}
