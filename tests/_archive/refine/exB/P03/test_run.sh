#!/bin/sh
PART=../../../../../bin/hecmw_part
FSTR=../../../../../bin/fstr
NPROCS=3

exec 2> test_run.log

error_exit() {
  echo $1
  exit 1
}

#cp hecmw_ctrl.dat.part hecmw_ctrl.dat
#$PART || error_exit "Error: partitioning failed."

cp hecmw_ctrl.dat.0 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr without refinement failed."

cp hecmw_ctrl.dat.1 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with single-step refinement failed"

cp hecmw_ctrl.dat.2 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with double-step refinement failed"

cp hecmw_ctrl.dat.3 hecmw_ctrl.dat
mpirun -np $NPROCS $FSTR 1>&2 || error_exit "Error: fstr with triple-step refinement failed"
