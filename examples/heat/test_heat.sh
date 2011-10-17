#!/bin/sh
# test_heat.sh ver.1.0
# 2006.04.07 by N.Imai
# ----------------------------
# FrontSTR Test for Examples

test_log="test_heat.log"

path="\
	exM\
	exN\
	exO\
	exP\
	exQ\
	exR\
	exS\
	exT"

rm -f ${test_log}

echo "FrontSTR test result for heat analysis" >> ${test_log}
date >> ${test_log}
echo "test option : $*"  >> ${test_log}

for i in ${path}
do
	echo "Test of ${i}"
	echo "***************************************************" >> ${test_log}
	echo "${i}" >> ${test_log}
	cd ${i}
	./test_${i}.sh $*
	cat test.log >> ../${test_log}
	cd ..
done

echo "end of test"

