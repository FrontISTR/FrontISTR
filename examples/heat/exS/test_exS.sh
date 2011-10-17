#!/bin/sh

test=test_exS
prg=../test_heat_sub.sh
test_log="test.log"

run_2d () {
	${prg} ${test_log} S231 S
	${prg} ${test_log} S232 S
	${prg} ${test_log} S241 S
	${prg} ${test_log} S242 S
}


run_3d () {
	${prg} ${test_log} S341 S
	${prg} ${test_log} S342 S
	${prg} ${test_log} S351 S
	${prg} ${test_log} S352 S
	${prg} ${test_log} S361 S
	${prg} ${test_log} S362 S
}

run_shell() {
	${prg} ${test_log} S731 S
	${prg} ${test_log} S741 S
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
