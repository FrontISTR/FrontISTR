/*****************************************************************************
 * Copyright (c) 2021 FrontISTR Commons
 * This software is released under the MIT License, see LICENSE.txt
 *****************************************************************************/

/*
 * Complete MueLu wrapper implementation with matrix setup
 * This file provides a more complete implementation showing how to
 * integrate MueLu with HEC-MW matrix data structures.
 */

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
#include "Epetra_MpiComm.h"
#include "EpetraExt_RowMatrixOut.h"

#include "hecmw_ML_helper.h"
#include "hecmw_ML_helper_33.h"
#include "hecmw_ML_helper_nn.h"

using Teuchos::RCP;
using Teuchos::rcp;

/*
 * Matrix assembly helper class
 */
class HECMWToEpetraMatrix {
private:
  int id_;
  int ndof_;
  int nlocal_;
  int nlocal_allcolumns_;
  RCP<Epetra_Comm> comm_;
  RCP<Epetra_Map> map_;
  RCP<Epetra_CrsMatrix> matrix_;

public:
  HECMWToEpetraMatrix(int id, int ndof) : id_(id), ndof_(ndof) {
    // Get MPI communicator
    MPI_Comm mpi_comm = HECMW_comm_get_comm();
    comm_ = rcp(new Epetra_MpiComm(mpi_comm));
    
    // Get matrix dimensions
    int ierr;
    hecmw_ml_get_nlocal_(&id_, &nlocal_, &nlocal_allcolumns_, &ierr);
    if (ierr != HECMW_SUCCESS) {
      throw std::runtime_error("Failed to get matrix dimensions");
    }

    // Create map
    map_ = rcp(new Epetra_Map(nlocal_, 0, *comm_));
  }

  void assembleMatrix() {
    // Estimate number of non-zeros per row
    int max_nnz_per_row = 27 * ndof_; // Conservative estimate for 3D problems
    
    // Create matrix
    matrix_ = rcp(new Epetra_CrsMatrix(Copy, *map_, max_nnz_per_row));

    // Assembly loop - get rows from HEC-MW and insert into Epetra matrix
    int max_entries = 1000; // Adjust based on your problem
    std::vector<int> indices(max_entries);
    std::vector<double> values(max_entries);
    std::vector<int> row_lengths(1);

    for (int local_row = 0; local_row < nlocal_; local_row++) {
      int global_row = map_->GID(local_row);
      int allocated_space = max_entries;
      int ierr;

      // Get row data from HEC-MW
      if (ndof_ == 3) {
        hecmw_ml_getrow_33_(&id_, &local_row, &local_row, &allocated_space,
                           indices.data(), values.data(), row_lengths.data(), &ierr);
      } else {
        hecmw_ml_getrow_nn_(&id_, &local_row, &local_row, &allocated_space,
                           indices.data(), values.data(), row_lengths.data(), &ierr);
      }

      if (ierr != HECMW_SUCCESS) {
        throw std::runtime_error("Failed to get matrix row");
      }

      // Insert row into Epetra matrix
      int num_entries = row_lengths[0];
      matrix_->InsertGlobalValues(global_row, num_entries, 
                                 values.data(), indices.data());
    }

    // Complete assembly
    matrix_->FillComplete();
  }

  RCP<Epetra_CrsMatrix> getMatrix() { return matrix_; }
  RCP<Epetra_Map> getMap() { return map_; }
  RCP<Epetra_Comm> getComm() { return comm_; }
};

/*
 * Null space (rigid body modes) setup
 */
class RigidBodyModes {
private:
  int id_;
  int ndof_;
  int nlocal_;
  RCP<Epetra_MultiVector> null_space_;

public:
  RigidBodyModes(int id, int ndof, int nlocal, RCP<Epetra_Map> map) 
    : id_(id), ndof_(ndof), nlocal_(nlocal) {
    
    int null_dim;
    if (ndof == 1) {
      null_dim = 1;
    } else if (ndof == 2) {
      null_dim = 3;
    } else {
      null_dim = 6;
    }

    // Create MultiVector for null space
    null_space_ = rcp(new Epetra_MultiVector(*map, null_dim, true));

    // Get rigid body modes from HEC-MW
    std::vector<double> rbm_data(null_dim * nlocal);
    int ierr;
    hecmw_ml_get_rbm_(&id_, rbm_data.data(), &ierr);
    if (ierr != HECMW_SUCCESS) {
      throw std::runtime_error("Failed to get rigid body modes");
    }

    // Copy data to MultiVector
    for (int vec = 0; vec < null_dim; vec++) {
      for (int node = 0; node < nlocal; node++) {
        (*null_space_)[vec][node] = rbm_data[vec * nlocal + node];
      }
    }
  }

  RCP<Epetra_MultiVector> getNullSpace() { return null_space_; }
};

/*
 * Complete MueLu setup with matrix assembly
 */
extern "C" {

static struct {
  RCP<HECMWToEpetraMatrix> matrix_helper;
  RCP<MueLu::EpetraOperator> muelu_prec;
  RCP<RigidBodyModes> rbm_helper;
  int ndof;
  bool is_setup;
} MueLuData[8]; // MAX_MI = 8

void hecmw_MueLu_wrapper_setup_complete(int *id, int *sym, int *Ndof, int *ierr) {
  int loglevel, myrank;
  
  if (*id <= 0 || 8 < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    hecmw_ml_get_loglevel_(id, &loglevel);
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);

    int idx = *id - 1;
    
    // Create matrix helper and assemble matrix
    MueLuData[idx].matrix_helper = rcp(new HECMWToEpetraMatrix(*id, *Ndof));
    MueLuData[idx].matrix_helper->assembleMatrix();
    
    auto matrix = MueLuData[idx].matrix_helper->getMatrix();
    auto map = MueLuData[idx].matrix_helper->getMap();

    // Setup rigid body modes
    int nlocal = map->NumMyElements();
    MueLuData[idx].rbm_helper = rcp(new RigidBodyModes(*id, *Ndof, nlocal, map));
    
    // Create parameter list
    RCP<Teuchos::ParameterList> paramList = rcp(new Teuchos::ParameterList());
    
    // Basic MueLu parameters
    paramList->set("verbosity", "low");
    paramList->set("number of equations", *Ndof);
    paramList->set("max levels", 4);
    
    // Use MUMPS as coarse solver if available
#ifdef HAVE_MUELU_AMESOS2
# ifdef HAVE_AMESOS2_MUMPS
    paramList->set("coarse: type", "MUMPS");
# else
    paramList->set("coarse: type", "KLU2");
# endif
#else
    paramList->set("coarse: type", "RELAXATION");
#endif

    // Smoother settings
    paramList->set("smoother: type", "CHEBYSHEV");
    paramList->set("smoother: sweeps", 2);
    
    // Aggregation settings
    paramList->set("aggregation: type", "uncoupled");
    paramList->set("aggregation: damping factor", 1.333);
    
    // Create preconditioner
    MueLuData[idx].muelu_prec = MueLu::CreateEpetraPreconditioner(
      matrix, *paramList, MueLuData[idx].rbm_helper->getNullSpace());
    
    MueLuData[idx].ndof = *Ndof;
    MueLuData[idx].is_setup = true;

    if (loglevel >= 1 && myrank == 0) {
      fprintf(stderr, "INFO: Complete MueLu preconditioner setup finished\n");
    }

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: Complete MueLu setup failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

void hecmw_MueLu_wrapper_apply_complete(int *id, double rhs[], int *ierr) {
  int nlocal, nlocal_allcolumns;
  int myrank;

  if (*id <= 0 || 8 < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);
    
    int idx = *id - 1;
    if (!MueLuData[idx].is_setup) {
      if (myrank == 0) fprintf(stderr, "ERROR: MueLu not setup\n");
      *ierr = HECMW_ERROR;
      return;
    }

    hecmw_ml_get_nlocal_(id, &nlocal, &nlocal_allcolumns, ierr);
    if (*ierr != HECMW_SUCCESS) return;

    auto map = MueLuData[idx].matrix_helper->getMap();
    
    // Create Epetra_MultiVector for RHS and solution
    Epetra_MultiVector rhs_vec(*map, 1, false);
    Epetra_MultiVector sol_vec(*map, 1, true); // Initialize to zero

    // Copy RHS data
    for (int i = 0; i < nlocal; i++) {
      rhs_vec[0][i] = rhs[i];
    }

    // Apply preconditioner: sol = M^(-1) * rhs
    int result = MueLuData[idx].muelu_prec->ApplyInverse(rhs_vec, sol_vec);
    if (result != 0) {
      if (myrank == 0) fprintf(stderr, "ERROR: MueLu ApplyInverse failed\n");
      *ierr = HECMW_ERROR;
      return;
    }

    // Copy solution back
    for (int i = 0; i < nlocal; i++) {
      rhs[i] = sol_vec[0][i];
    }

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: Complete MueLu apply failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

void hecmw_MueLu_wrapper_clear_complete(int *id, int *ierr) {
  int myrank;

  if (*id <= 0 || 8 < *id) {
    *ierr = HECMW_ERROR;
    return;
  }

  try {
    HECMW_Comm_rank(HECMW_comm_get_comm(), &myrank);
    
    int idx = *id - 1;
    
    // Clear all objects
    MueLuData[idx].muelu_prec = Teuchos::null;
    MueLuData[idx].matrix_helper = RCP<HECMWToEpetraMatrix>();
    MueLuData[idx].rbm_helper = RCP<RigidBodyModes>();
    MueLuData[idx].is_setup = false;

  } catch (std::exception& e) {
    if (myrank == 0) fprintf(stderr, "ERROR: Complete MueLu clear failed: %s\n", e.what());
    *ierr = HECMW_ERROR;
  }
}

} // extern "C"

#endif /* HECMW_WITH_MUELU */
