#!/bin/bash
#export PATH=$HOME/src/autofistr/INTEL_IMPI/FrontISTR/bin:$PATH
export PATH=$HOME/work/FrontISTR/fistr1/bin:$HOME/work/FrontISTR/hecmw1/bin:$PATH
for np in 1 2
do
	echo Partitioning to ${np} domains
	cd P${np}
	mpirun -np 1 hecmw_part1 -v > hecmw_part1.out 2>&1
	cd ..
	echo Running FrontISTR
	for i in P${np} P${np}R1 P${np}R2 P${np}R3
	do
		echo $i
		cd $i
		mpirun -np ${np} fistr1 > fistr1.out 2>&1
		cd ..
	done
done
