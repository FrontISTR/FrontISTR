#!/bin/sh

export model_2d="\
    C231\
    C232\
    C241\
    C242"

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

export cnt_2d=C200.cnt
export cnt_3d=C300.cnt
export cnt_shell=C700.cnt

../test_static_sub.sh $*

