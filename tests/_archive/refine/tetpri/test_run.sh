#!/bin/sh
FSTR=../../../../fistr1/bin/fistr1
NPROCS=1

exec 2> test_run.log

error_exit() {
  echo $1
  exit 1
}

cp hecmw_ctrl.dat.0 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr without refinement failed."

cp hecmw_ctrl.dat.1 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with single-step refinement failed"

cp hecmw_ctrl.dat.2 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with double-step refinement failed"

cp hecmw_ctrl.dat.3 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with triple-step refinement failed"
