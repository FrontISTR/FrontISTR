#!/bin/sh
# test_sub.sh ver.1.0
# 2006.04.05 by N.Imai
# ----------------------------
# FrontSTR Test for Examples

PATH=../../../fistr1/bin:$PATH
fstr=fistr1
launcher=../../../fistr1/tools/launcher/fstr_srun.sh

test_log="test.log"


print_u () {
	fg=1
	cat $1 | while read s
	do
		if [ "${s}" = "##### Global Summary :Max/Min####" ]; then
			if [ fg ]; then
				echo "---------------------------------------------------"
			fi
			fg=0
			# U1,U2,U3
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			# E11-S13
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			# SMS
			read s
			echo ${s}
			#
		fi
	done
}


print_shell_u () {
	fg=1
	cat $1 | while read s
	do
		if [ "${s}" = "##### Global Summary :Max/Min####" ]; then
			if [ fg ]; then
				echo "---------------------------------------------------"
			fi
			fg=0
			# U1,U2,U3,R1,R2,R3
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			# E11(+)-S13(+)
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			read s
			# SMS(+)
			read s
			echo ${s}
			read s
			read s
			read s
			read s
			read s
			# SMS(-)
			read s
			echo ${s}
			#
		fi
	done
}


run () {
	echo ${1}
	echo "===================================================" >> ${test_log}
	echo ${1} >> ${test_log}
	mesh=${1}.msh
	log=${1}.log
	res=${1}.res
	vis=${1}.vis
	fg_shell=${2}
	rm 0.log
	${launcher} -f ${fstr} -l ${log} ${mesh} ${cnt} ${res} ${vis}
	if [ -e "0.log" ]; then
		if [ "${fg_shell}" = "true" ]; then
			print_shell_u 0.log >> ${test_log}
		else
			print_u 0.log >> ${test_log}
		fi
	else
		echo "### ERROR ###" >> ${test_log}
	fi
}

run_2d () {
	cnt=${cnt_2d}
	for i in ${model_2d}
	do
		run ${i} false
	done
}

run_3d () {
	cnt=${cnt_3d}
	for i in ${model_3d}
	do
		run ${i} false
	done
}

run_shell () {
	cnt=${cnt_shell}
	for i in ${model_shell}
	do
		run ${i} true
	done
}

run_shell33 () {
	cnt=${cnt_shell33}
	for i in ${model_shell33}
	do
		run ${i} true
	done
}

run_beam33 () {
	cnt=${cnt_beam33}
	for i in ${model_beam33}
	do
		run ${i} true
	done
}

list_up () {
	echo "2d model   : ${model_2d}"
	echo "2d ctrl    : ${cnt_2d}"
	echo "3d model   : ${model_3d}"
	echo "3d ctrl    : ${cnt_3d}"
	echo "shell model: ${model_shell}"
	echo "shell ctrl : ${cnt_shell}"
	echo "shell33 model: ${model_shell33}"
	echo "shell33 ctrl : ${cnt_shell33}"
	echo "beam33 model: ${model_beam33}"
	echo "beam33 ctrl : ${cnt_beam33}"
}


help () {
	echo "FrontSTR executing test"
	echo "[usage] test.sh (options)"
	echo " -h      : help (this message)"
	echo " -l      : list up models"
	echo "  2d     : 2 dimentional model"
	echo "  3d     : 3 dimentional model"
	echo "  shell  : shell model"
	echo "  shell33: shell33 model"
	echo "  beam33: beam33 model"
	echo "  all or no options : all model" 
}


############################# MAIN ################################

echo "Max/Min Displacement" > ${test_log}

if [ $# -lt 1 -o "${1}" = "all" ]; then
	run_2d
	run_3d
	run_shell
	run_shell33
	run_beam33
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
	elif [ "${i}" = "shell33" ]; then run_shell33
	elif [ "${i}" = "beam33" ]; then run_beam33
	else
		echo "## Error in ${0}: unknown parameter ${i}"
		echo "   show help with -h"
		exit
	fi
done


