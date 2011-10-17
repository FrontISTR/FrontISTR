#!/bin/sh

test=test_exQ
prg=../test_heat_sub.sh
test_log="test.log"

run_2d () {
	${prg} ${test_log} Q231 Q230
	${prg} ${test_log} Q232 Q230
	${prg} ${test_log} Q241 Q240
	${prg} ${test_log} Q242 Q240
}


run_3d () {
	${prg} ${test_log} Q341 Q340
	${prg} ${test_log} Q342 Q340
	${prg} ${test_log} Q351 Q350
	${prg} ${test_log} Q352 Q350
	${prg} ${test_log} Q361 Q360
	${prg} ${test_log} Q362 Q360
}

run_shell() {
	${prg} ${test_log} Q731 Q700
	${prg} ${test_log} Q741 Q700
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
