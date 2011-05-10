#!/bin/sh
# test_neu_heat_sub.sh ver.1.0
# 2006.04.20 by N.Imai
# ----------------------------
# FrontSTR Test for NEU Examples

fstr=../../../../bin/fstr
launcher=../../../launcher/fstr_srun.sh

test_log="test.log"

print_u () {
	echo "---------------------------------------------------"
	cat $1 | while read s1 s2
	do
		if [ "${s1}" == "Maximum" ]; then
			echo "Maximum ${s2}"
		elif [ "${s1}" == "Minimum" ]; then
			echo "Minimum ${s2}"
		fi
	done
}


run () {
	echo ${1}
	echo "===================================================" >> ${test_log}
	echo ${1} >> ${test_log}
	mesh=${1}.msh
	cnt=${1}.cnt
	log=${1}.log
	res=${1}.res
	vis=${1}.vis
	rm -f 0.log
	${launcher} -f ${fstr} -l ${log} ${mesh} ${cnt} ${res} ${vis}
	if [ -e "0.log" ]; then
			print_u 0.log >> ${test_log}
	fi
}


run_2d () {
	for i in ${model_2d}
	do
		run ${i}
	done
}

run_3d () {
	for i in ${model_3d}
	do
		run ${i}
	done
}

run_shell () {
	for i in ${model_shell}
	do
		run ${i}
	done
}

list_up () {
	echo "2d model   : ${model_2d}"
	echo "3d model   : ${model_3d}"
	echo "shell model: ${model_shell}"
}


help () {
	echo "FrontSTR executing test"
	echo "[usage] test.sh (options)"
	echo " -h      : help (this message)"
	echo " -l      : list up models"
	echo "  2d     : 2 dimentional model"
	echo "  3d     : 3 dimentional model"
	echo "  shell  : shell model"
	echo "  all or no options : all model" 
}

############################# MAIN ################################

echo "Max/Min Temperature" > ${test_log}

if [ $# -lt 1 -o "${1}" == "all" ]; then
	run_2d
	run_3d
	run_shell
	exit
fi

for i in $*
do
	if   [ "${i}" == "-h"    ]; then
		help
		exit
	elif [ "${i}" == "-l"    ]; then
		list_up
		exit
	elif [ "${i}" == "2d"    ]; then run_2d
	elif [ "${i}" == "3d"    ]; then run_3d
	elif [ "${i}" == "shell" ]; then run_shell
	else
		echo "## Error in ${0}: unknown parameter ${i}"
		echo "   show help with -h"
		exit
	fi
done
