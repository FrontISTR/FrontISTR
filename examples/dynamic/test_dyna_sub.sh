#!/bin/sh
# test_sub.sh ver.1.0
# 2010.03.26 Xi YUAN
# ----------------------------
# FrontSTR Test for Examples

PATH=../../../bin:$PATH
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
	cnt=${1}.cnt
	log=${1}.log
	res=${1}.res
	vis=${1}.vis
	displog=${1}.disp.txt
	fg_shell=${2}
	rm 0.log
	rm dyna_disp_p1.txt
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
	if [ -e "dyna_disp_p1.txt" ]; then
		cp dyna_disp_p1.txt ${displog}
	fi
	
}

run_342 () {
	mesh=${model_342}
	for i in ${cnt_342}
	do
		run ${i} false
	done
}

run_361 () {
	mesh=${model_361}
	for i in ${cnt_361}
	do
		run ${i} false
	done
}

run_741 () {
	mesh=${model_741}
	for i in ${cnt_741}
	do
		run ${i} false
	done
}

run_731 () {
	mesh=${model_731}
	for i in ${cnt_731}
	do
		run ${i} false
	done
}


list_up () {
	echo "2d model   : ${model_2d}"
	echo "2d ctrl    : ${cnt_2d}"
	echo "3d model   : ${model_3d}"
	echo "3d ctrl    : ${cnt_3d}"
	echo "shell model: ${model_shell}"
	echo "shell ctrl : ${cnt_shell}"
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

echo "Max/Min Displacement" > ${test_log}

if [ $# -lt 1 -o "${1}" = "all" ]; then
	run_342
	run_361
	run_741
	run_731
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
	elif [ "${i}" = "342"    ]; then run_342
	elif [ "${i}" = "361"    ]; then run_361
	elif [ "${i}" = "shell" ]; then run_shell
	else
		echo "## Error in ${0}: unknown parameter ${i}"
		echo "   show help with -h"
		exit
	fi
done


