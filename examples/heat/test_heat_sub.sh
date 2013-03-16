#!/bin/sh
# test_sub2.sh ver.1.0
# 2006.04.05 by N.Imai
# ----------------------------
# FrontSTR Test for Examples

PATH=../../../bin:$PATH
fstr=fistr1
launcher=../../../fistr1/tools/launcher/fstr_srun.sh

test_log=${1}

print_u () {
	echo "---------------------------------------------------"
	cat $1 | while read s1 s2
	do
		if [ "${s1}" = "Maximum" ]; then
			echo "Maximum ${s2}"
		elif [ "${s1}" = "Minimum" ]; then
			echo "Minimum ${s2}"
		fi
	done
}


# run [mesh] [cnt]
run () {
	echo ${1}
	echo "===================================================" >> ${test_log}
	echo ${1} >> ${test_log}
	mesh=${1}.msh
	log=${1}.log
	res=${1}.res
	vis=${1}.vis
	cnt=${2}.cnt
	rm -f 0.log
	fg_shell=${2}
	${launcher} -f ${fstr} -l ${log} ${mesh} ${cnt} ${res} ${vis}
	if [ -e "0.log" ]; then
			print_u 0.log >> ${test_log}
	fi
}


############################# MAIN ################################


run $2 $3


