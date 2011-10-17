#!/bin/sh

test=test_exR
prg=../test_heat_sub.sh
test_log="test.log"

run_2d () {
	${prg} ${test_log} R231 R230
	${prg} ${test_log} R232 R230
	${prg} ${test_log} R241 R240
	${prg} ${test_log} R242 R240
}


run_3d () {
	${prg} ${test_log} R341 R340
	${prg} ${test_log} R342 R340
	${prg} ${test_log} R351 R350
	${prg} ${test_log} R352 R350
	${prg} ${test_log} R361 R360
	${prg} ${test_log} R362 R360
}

run_shell() {
	${prg} ${test_log} R731 R700
	${prg} ${test_log} R741 R700
}


help () {
	echo "FrontSTR executing test"
	echo "[usage] ${test} (options)"
	echo " -h      : help (this message)"
	echo "  2d     : 2 dimentional model"
	echo "  3d     : 3 dimentional model"
	echo "  shell  : shell model"
	echo "  all or no options : all model" 
}


############################# MAIN ################################

echo "Max/Min Temperature" > ${test_log}

if [ $# -lt 1 -o "${1}" = "all" ]; then
	run_2d
	run_3d
	run_shell
	exit
fi

for i in $*
do
	if   [ "${i}" = "-h"    ]; then
		help
		exit
	elif [ "${i}" = "-l"    ]; then
		list_up
		exit
	elif [ "${i}" = "2d"    ]; then run_2d
	elif [ "${i}" = "3d"    ]; then run_3d
	elif [ "${i}" = "shell" ]; then run_shell
	else
		echo "## Error in ${0}: unknown parameter ${i}"
		echo "   show help with -h"
		exit
	fi
done
