#!/bin/sh

export model_2d="\
	A231\
    A232\
	A241\
	A242"

export model_3d="\
	A341\
	A342\
	A351\
	A352\
	A361\
	A362"

export model_shell="\
	A731\
	A741"

export cnt_2d=A200.cnt
export cnt_3d=A300.cnt
export cnt_shell=A700.cnt

../test_static_sub.sh $*

