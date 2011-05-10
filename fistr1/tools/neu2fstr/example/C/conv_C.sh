#!/bin/sh

export model="\
	C341\
	C342\
	C351\
	C352\
	C361\
	C362"

export model_shell="\
	C731\
	C741"

export analysis="s"
../conv_sub.sh ${1}


