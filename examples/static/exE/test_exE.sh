#!/bin/sh

export model_2d="\
	E231\
	E232\
	E241\
	E242"

export model_3d="\
	E341\
	E342\
	E351\
	E352\
	E361\
	E362"

export model_shell="\
	E731\
	E741"

export cnt_2d=E200.cnt
export cnt_3d=E300.cnt
export cnt_shell=E700.cnt

../test_static_sub.sh $*

