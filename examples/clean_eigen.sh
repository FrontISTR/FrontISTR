#!/bin/sh

DIRS="\
	eigen\
	eigen/exJ\
	eigen/exK"

for i in ${DIRS}
do
	echo "clean up ${i}"
	rm -f $i/*.dat
	rm -f $i/*.log
	rm -f $i/*.msg
	rm -f $i/*.sta
	rm -f $i/*.bmp
	rm -f $i/*.res.*
	rm -f $i/*.dbg.*
	rm -f $i/*.neu
	rm -f $i/*.out
	rm -f $i/*.ini
	rm -f $i/*.inp
done

