#!/bin/sh

export model_2d=""

export model_3d="\
	MB361_CFlux\
	MC361\
	MD361\
	ME361\
	MF361"

export model_shell=""

../test_neu_heat_sub.sh $*



