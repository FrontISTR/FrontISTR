#!/bin/sh

export model="\
	A341\
	A342\
	A351\
	A352\
	A361\
	A362"

export model_shell="\
	A731"

export analysis="s"
../conv_sub.sh ${1}
