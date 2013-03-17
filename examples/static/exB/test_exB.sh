#!/bin/sh

export model_2d="\
	B231\
    B232\
	B241\
	B242"

export model_3d="\
	B341\
	B342\
	B351\
	B352\
	B361\
	B362"

export model_shell="\
	B731\
	B741"

../test_static_sub2.sh $*

