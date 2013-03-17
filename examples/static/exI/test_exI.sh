#!/bin/sh

export model_2d=""

export model_3d="\
	A341\
	A342\
	A351\
	A352\
	A361\
	A362"

export model_shell=

export cnt_2d=I200.cnt
export cnt_3d=I300.cnt
export cnt_shell=

../test_static_sub.sh $*

