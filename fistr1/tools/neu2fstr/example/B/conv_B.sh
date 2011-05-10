#!/bin/sh

export model="\
	B341\
	B342\
	B351\
	B352\
	B361\
	B362"

export model_shell="\
	B731\
	B741"

export analysis="s"
../conv_sub.sh ${1}


