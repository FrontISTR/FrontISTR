#!/bin/sh

test=test_exN
prg=../test_heat_sub.sh
test_log="test.log"

run_2d () {
	${prg} ${test_log} N231 N
	${prg} ${test_log} N232 N
	${prg} ${test_log} N241 N
	${prg} ${test_log} N242 N
}


run_3d () {
	${prg} ${test_log} N341 N
	${prg} ${test_log} N342 N
	${prg} ${test_log} N351 N
	${prg} ${test_log} N352 N
	${prg} ${test_log} N361 N
	${prg} ${test_log} N362 N
}

run_shell() {
	${prg} ${test_log} N731 N
	${prg} ${test_log} N741 N
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
