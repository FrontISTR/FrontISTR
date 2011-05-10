#!/bin/sh

export model_2d=""

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

../test_neu_static_sub.sh $*



