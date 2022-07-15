#include <stdio.h>
#include <stdlib.h>
#include "lobpcg.h"
#include "interpreter.h"
#include "multivector.h"
#include "multi_vector.h"
#include "pcg_multi.h"

typedef struct {
  mv_InterfaceInterpreter ii;
} data_struct;

BlopexInt dsygv_ (
  BlopexInt *itype, char *jobz, char *uplo, BlopexInt *
  n, double *a, BlopexInt *lda, double *b, BlopexInt *ldb,
  double *w, double *work, BlopexInt *lwork, BlopexInt *info);

BlopexInt dpotrf_ (
  char *uplo, BlopexInt *n, double *a, BlopexInt *
  lda, BlopexInt *info);

void (*callback_matmult_opA)(void *, void*, void*);
void (*callback_matmult_opB)(void *, void*, void*);
void (*callback_matmult_opT)(void *, void*, void*);

void OperatorAMultiVector(void *data, void *x, void *y) {
  serial_Multi_Vector* a;
  serial_Multi_Vector* b;
  a = (serial_Multi_Vector*)x;
  b = (serial_Multi_Vector*)y;
  callback_matmult_opA(data, a->data, b->data);
}

void OperatorBMultiVector(void * data, void * x, void * y) {
  serial_Multi_Vector* a;
  serial_Multi_Vector* b;
  a = (serial_Multi_Vector*)x;
  b = (serial_Multi_Vector*)y;
  callback_matmult_opB(data, a->data, b->data);
}

void OperatorTMultiVector(void *data, void *x, void *y) {
  serial_Multi_Vector* a;
  serial_Multi_Vector* b;
  a = (serial_Multi_Vector*)x;
  b = (serial_Multi_Vector*)y;
  callback_matmult_opT(data, a->data, b->data);
}

void blopex_lobpcg_solve_c_(
    int     *n_eigs,       /* number of eigenvalues             */
    int     *maxit,
    double  *tol,
    int     *mat_n,
    void    *matmult_opA, /* Fortran routine for operator A    */
    void    *matmult_opB, /* Fortran routine for operator B    */
    void    *matmult_opT, /* Fortran routine for operator T    */
    int     *is_B,
    int     *is_T,
    int     *log_level,
    double* eigval,
    double* eigvec
  ) {
  int                        ierr;         /* for PETSc return code        */
  mv_MultiVectorPtr          eigenvectors; /* the eigenvectors             */
  double*                    resid;        /* the residuals                */
  int                        iterations;   /* number of iterations         */
  lobpcg_Tolerance           lobpcg_tol;   /* residual tolerance           */
  mv_InterfaceInterpreter    ii;           /* Interface Interpreter        */
  lobpcg_BLASLAPACKFunctions blap_fn;      /* BLAS functions               */
  data_struct                aux_data;     /* auxillary data               */
  serial_Multi_Vector* x;

  iterations = 0;

  callback_matmult_opA = matmult_opA;
  callback_matmult_opB = matmult_opB;
  callback_matmult_opT = matmult_opT;

  void *ptr_opB = NULL;
  void *ptr_opT = NULL;
  if(*is_B == 1){ptr_opB = OperatorBMultiVector;};
  if(*is_T == 1){ptr_opT = OperatorTMultiVector;};

  blap_fn.dpotrf = dpotrf_;
  blap_fn.dsygv = dsygv_;

  lobpcg_tol.absolute = *tol;
  lobpcg_tol.relative = 1.0e-12;

  resid = (double *)malloc(sizeof(double)*(*n_eigs));

  x = serial_Multi_VectorCreate(*mat_n, *n_eigs);
  serial_Multi_VectorInitialize(x);
  serial_Multi_VectorSetRandomValues(x, 1);

  SerialSetupInterpreter(&ii);
  eigenvectors = mv_MultiVectorWrap(&ii, x, 0);

  ierr = lobpcg_solve_double(
    eigenvectors,     /*input eigen vectors */
    &aux_data,        /*input dummy */
    OperatorAMultiVector, /*input operation A */
    &aux_data,        /*input dummy */
    ptr_opB,          /*input operation B */
    &aux_data,        /*input dummy */
    ptr_opT,          /*input operation precond */
    NULL,
    blap_fn,          /*input-lapack functions */
    lobpcg_tol,       /*input-tolerances */
    *maxit,           /*input-max iterations */
    *log_level,       /*input-verbosity level */
    &iterations,      /*output-actual iterations */
    eigval,           /*output-eigenvalues */
    NULL,             /*output-eigenvalues history */
    0,                /*output-history global height */
    resid,            /*output-residual norms */
    NULL,             /*output-residual norms history */
    0                 /*output-history global height  */
  );

  /* eigen vector output */
  for (int i = 0; i < *n_eigs; ++i){
    for (int j = 0; j < *mat_n; ++j){
      eigvec[i*(*n_eigs) + j] = x->data[i*(*n_eigs) + j];
    }
  }

  return;
}
