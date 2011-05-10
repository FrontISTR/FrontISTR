#!/bin/sh

export model_2d=""

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

../test_neu_static_sub.sh $*



