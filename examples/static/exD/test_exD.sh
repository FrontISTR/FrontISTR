#!/bin/sh

export model_2d="\
	D231\
    D232\
	D241\
	D242"

export model_3d="\
	D341\
	D342\
	D351\
	D352\
	D361\
	D362"

export model_shell="\
	D731\
	D741"

export cnt_2d=D200.cnt
export cnt_3d=D300.cnt
export cnt_shell=D700.cnt

../test_static_sub.sh $*

