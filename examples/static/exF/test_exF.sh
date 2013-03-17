#!/bin/sh

export model_2d="\
	F231\
	F232\
	F241\
	F242"

export model_3d="\
	F341\
	F342\
	F351\
	F352\
	F361\
	F362"

export model_shell=

export cnt_2d=F200.cnt
export cnt_3d=F300.cnt
export cnt_shell=

../test_static_sub.sh $*

