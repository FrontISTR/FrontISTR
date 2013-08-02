#!/bin/sh
# convert.sh ver.1.0
# 2006.04.20 by N.Imai
# ----------------------------
# FrontSTR Test for NEU Examples

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

