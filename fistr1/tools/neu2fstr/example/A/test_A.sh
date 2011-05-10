#!/bin/sh

export model_2d=""

export model_3d="\
	A341\
	A342\
	A351\
	A352\
	A361\
	A362"

export model_shell="\
	A731"

../test_neu_static_sub.sh $*



