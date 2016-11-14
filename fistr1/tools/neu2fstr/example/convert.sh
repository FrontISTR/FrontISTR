#!/bin/sh
#-------------------------------------------------------------------------------
# Copyright (c) 2016 The University of Tokyo
# This software is released under the MIT License, see LICENSE.txt
#-------------------------------------------------------------------------------

test_log="test_neu.log"

path="\
	A\
	B\
	C\
	D\
	heat"

for i in ${path}
do
	cd ${i}
	./conv_${i}.sh $*
	cd ..
done

