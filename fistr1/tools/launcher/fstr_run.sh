#!/bin/sh
############################################################
# FrontSTR Launcher Ver.1.0
# 2006.04.12 by N.Imai
# ----------------------------------------------------------
# This B shell script automatically generates hecmw_ctrl.dat
# with command parameters, and launch FrontSTR
# ----------------------------------------------------------
# Copyright (C) 2006 RSS21 Project (Univ. of Tokyo)
############################################################

fstr=fistr1
this=fstr_srun.sh
fg_log=0
fg_p=0

help () {
	echo "FrontSTR launcher (serial)"
	echo "[usage] ./${this} (options) [mesh] [cnt] ([res] [vis])"
	echo "   -h        : help (this message)"
	echo "   -l  [log] : rename 0.log to [log]"
	echo "   -f  [fstr]: FrontSTR command [fstr] (default:${fstr})"
	echo "   -p        : parallel execution"
	echo "   [mesh]    : mesh file"
	echo "   [cnt]     : fstr control file"
	echo "   [res]     : result file (default:[mesh].res)"
	echo "   [vis]     : visual file (default:[mesh].vis)"
}

#
# PARSE COMMAND LINE PARAMETERS
#

if [ "$1" = "" ] ; then
	help
	exit
fi

while [ "$1" = "-h" -o "$1" = "--help" -o "$1" = "-l"  -o "$1" = "-f" -o "$1" = "-p" ]
do
	if [ "$1" = "" -o "$1" = "-h" -o "$1" = "--help" ] ;then
		help
		exit
	elif [ "$1" = "-l" ]; then
		fg_log=1
		log=$2
		if [ "${log}" = "" ]; then
			echo "log file name required"
			exit
		fi
		shift
	elif [ "$1" = "-f" ]; then
		fstr=$2
		shift
	fi
	elif [ "$1" = "-p" ]; then
		fg_p=1
	shift
done

if [ $# -lt 2 ]; then
	echo "[mesh] and [cnt] must be required (see help:-h)"
	exit
fi

#
# SETUP FILE NAMES
#

hec_ctrl="hecmw_ctrl.dat"
mesh=$1
cnt=$2
res=$3
vis=$4

if [ "${res}" = "" ]; then
	res=${mesh}.res
fi
if [ "${vis}" = "" ]; then
	vis=${mesh}.vis
fi

#
# GENERATE HEC-MW CONTROL FILE
#

echo "##"                                      >  $hec_ctrl
echo "## HEC-MW control file for FrontSTR"     >> $hec_ctrl
echo "## Auto created by "${this}              >> $hec_ctrl
echo "## "`date`                               >> $hec_ctrl
echo "##"                                      >> $hec_ctrl

if [ "${fg_p}" = "0" ]; then
	echo "!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE"   >> $hec_ctrl
else
	echo "!MESH, NAME=fstrMSH,TYPE=HECMW-DIST"     >> $hec_ctrl
fi
echo ${mesh}                                   >> $hec_ctrl
echo "!CONTROL,NAME=fstrCNT"                   >> $hec_ctrl
echo ${cnt}                                    >> $hec_ctrl
echo "!RESULT,NAME=fstrRES,IO=OUT"             >> $hec_ctrl
echo ${res}                                    >> $hec_ctrl
echo "!RESULT,NAME=vis_out,IO=OUT"             >> $hec_ctrl
echo ${vis}                                    >> $hec_ctrl

#
# LAUNCH FrontSTR & COPY LOG FILE
#

echo "---------------- "${hec_ctrl}" ------------------"
cat $hec_ctrl
echo "--------------------------------------------------"
echo
echo launch ${fstr}

${fstr}

if [ $fg_log -eq 1 ]; then
	echo cp 0.log ${log}
	cp 0.log ${log}
fi
