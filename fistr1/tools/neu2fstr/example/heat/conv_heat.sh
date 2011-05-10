#!/bin/sh

export model="\
	MB361_CFlux\
	MC361\
	MD361\
	ME361\
	MF361"

export model_shell=""

export analysis="h"
../conv_sub.sh ${1}


