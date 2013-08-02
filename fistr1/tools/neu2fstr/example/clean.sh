#!/bin/sh

DIRS="\
	A \
	B \
	C \
	D \
	heat \
	."


for i in ${DIRS}
do
	echo "clean up ${i}"
	rm -f $i/hecmw_ctrl.dat
	rm -f $i/*.msh
	rm -f $i/*.cnt
	rm -f $i/*.log
	rm -f $i/*.msg
	rm -f $i/*.sta
	rm -f $i/*.bmp
	rm -f $i/*.res.*
	rm -f $i/*.dbg.*
	rm -f $i/*.neu
	rm -f $i/*.out
	rm -f $i/*.ini
done

