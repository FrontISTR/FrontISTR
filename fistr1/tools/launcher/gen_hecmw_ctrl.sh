#!/bin/sh
############################################################
# hecmw_ctrl.dat generator
# 2006.07.27 by N.Imai
# ----------------------------------------------------------
# This B shell script automatically generates hecmw_ctrl.dat
# with command parameters
# ----------------------------------------------------------
# Copyright (C) 2006 RSS21 Project (Univ. of Tokyo)
############################################################

this=gen_hecmw_ctrl.sh
fg_p=0
fg_part=0

help () {
	echo "hecmw_ctrl.dat generator"
	echo "[usage] ./${this} (options) [mesh] [cnt] ([res] [vis])"
	echo "   -h        : help (this message)"
	echo "   -p        : for parallel execution"
	echo " -part [dist]: for partitioner execution"
	echo "   [mesh]    : mesh file"
	echo "   [cnt]     : fstr control file"
	echo "   [res]     : result file (default:[mesh].res)"
	echo "   [vis]     : visual file (default:[mesh].vis)"
}

#
# PARSE COMMAND LINE PARAMETERS
#

if [ "$1" == "" ] ; then
	help
	exit
fi

while [ "$1" == "-h" -o "$1" == "--help" -o "$1" == "-p"  -o "$1" == "-part" ]
do
	if [ "$1" == "" -o "$1" == "-h" -o "$1" == "--help" ] ;then
		help
		exit
	elif [ "$1" == "-p" ]; then
		fg_p=1
	elif [ "$1" == "-part" ]; then
		fg_part=1
		dist=$2
		if [ "${dist}" == "" ]; then
			echo "distribution mesh name required"
			exit
		fi
		shift
	fi
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

if [ "${res}" == "" ]; then
	res=${mesh}.res
fi
if [ "${vis}" == "" ]; then
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
if [ "${fg_p}" == "0" ]; then
	echo "!MESH, NAME=fstrMSH,TYPE=HECMW-ENTIRE"   >> $hec_ctrl
	echo ${mesh}                                   >> $hec_ctrl
else
	echo "!MESH, NAME=fstrMSH,TYPE=HECMW-DIST"   >> $hec_ctrl
	if [ "${fg_part}" == "1" ]; then
		echo ${dist}                         >> $hec_ctrl
	else
		echo ${mesh}                         >> $hec_ctrl
	fi
fi
echo "!CONTROL,NAME=fstrCNT"                   >> $hec_ctrl
echo ${cnt}                                    >> $hec_ctrl
echo "!RESULT,NAME=fstrRES,IO=OUT"             >> $hec_ctrl
echo ${res}                                    >> $hec_ctrl
echo "!RESULT,NAME=vis_out,IO=OUT"             >> $hec_ctrl
echo ${vis}                                    >> $hec_ctrl

if [ "${fg_part}" == "1" ]; then
	echo "!MESH, NAME=part_in, TYPE=HECMW-ENTIRE" >> $hec_ctrl
	echo ${mesh}                                  >> $hec_ctrl
	echo "!MESH, NAME=part_out,TYPE=HECMW-DIST"   >> $hec_ctrl
	echo ${dist}                                  >> $hec_ctrl
fi

cat $hec_ctrl
