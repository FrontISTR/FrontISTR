#!/bin/sh
clean_dir() {
  cwd=$PWD
  cd $1
  rm -f hecmw_ctrl.dat FSTR.* *.log *.bmp *.inp *.ini
  cd $cwd
}

for i in . P02 P03
do
  clean_dir $i
done
