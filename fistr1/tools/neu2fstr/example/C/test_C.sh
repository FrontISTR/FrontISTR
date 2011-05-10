#!/bin/sh

export model_2d=""

export model_3d="\
	C341\
	C342\
	C351\
	C352\
	C361\
	C362"

export model_shell="\
	C731\
	C741"

../test_neu_static_sub.sh $*



