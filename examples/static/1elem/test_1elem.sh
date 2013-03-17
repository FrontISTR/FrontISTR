#!/bin/sh

export model_2d=""

export model_3d="\
	neohooke\
	rivlin\
	arruda\
	mises\
	mohr\
	drucker\
	swift\
	ramberg\
	viscoe\
	creep\
	relax"

export model_shell=""

../test_static_sub2.sh $*

