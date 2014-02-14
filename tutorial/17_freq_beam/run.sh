#!/bin/sh
cp hecmw_ctrl_eigen.dat hecmw_ctrl.dat
../../fistr1/bin/fistr1
mv 0.log eigen_0.log
cp hecmw_ctrl_freq.dat hecmw_ctrl.dat
../../fistr1/bin/fistr1
