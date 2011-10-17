#!/bin/sh

test=test_exO
prg=../test_heat_sub.sh
test_log="test.log"

run_2d () {
	${prg} ${test_log} O231 O200
	${prg} ${test_log} O232 O200
	${prg} ${test_log} O241 O200
	${prg} ${test_log} O242 O200
}


run_3d () {
	${prg} ${test_log} O341 O300
	${prg} ${test_log} O342 O300
	${prg} ${test_log} O351 O300
	${prg} ${test_log} O352 O300
	${prg} ${test_log} O361 O300
	${prg} ${test_log} O362 O300
}

run_shell() {
	${prg} ${test_log} O731 O700
	${prg} ${test_log} O741 O700
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
