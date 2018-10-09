#!/bin/sh
# test_sub.sh ver.1.0
# 2006.04.05 by N.Imai
# --------------------------------
# FrontSTR Test for Eigen Examples

PATH=../../../fistr1/bin:$PATH
fstr=fistr1
launcher=../../../fistr1/tools/launcher/fstr_srun.sh

test_log="test.log"


print_u () {
	fg=1
	OLDIFS=$IFS
	IFS=
	cat $1 | while read s
	do
		if [ "${s}" = "*RESULT OF EIGEN VALUE ANALYSIS*" ]; then
			if [ fg ]; then
				echo "---------------------------------------------------"
			fi
			fg=0
			# skip lines
			read s
			read s
			read s
			echo ${s}
			read s
			read s
			# eigen values and etc.
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
			read s
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
			read s
			echo ${s}
		fi
	done
	IFS=$OLDIFS
}


run () {
	echo ${1}
	echo "===================================================" >> ${test_log}
	echo ${1} >> ${test_log}
	mesh=${1}.msh
	log=${1}.log
	res=${1}.res
	vis=${1}.vis
	rm 0.log
	${launcher} -f ${fstr} -l ${log} ${mesh} ${cnt} ${res} ${vis}
	if [ -e "0.log" ]; then
		print_u 0.log >> ${test_log}
	else
		echo "### ERROR ###" >> ${test_log}
	fi
}

run_2d () {
	if [ "${model_2d}" = "" ]; then
		echo "NO 2D MODEL"
		return
	fi
	cnt=${cnt_2d}
	for i in ${model_2d}
	do
		run ${i}
	done
}

run_3d () {
	if [ "${model_3d}" = "" ]; then
		echo "NO 3D MODEL"
		return
	fi
	cnt=${cnt_3d}
	for i in ${model_3d}
	do
		run ${i}
	done
}

run_shell () {
	if [ "${model_shell}" = "" ]; then
		echo "NO SHELL MODEL"
		return
	fi
	cnt=${cnt_shell}
	for i in ${model_shell}
	do
		run ${i}
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


