#!/bin/sh
# test_neu.sh ver.1.0
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

########################## MAIN ############################


rm -f ${test_log}

echo "FrontSTR test result for analysis" >> ${test_log}
date >> ${test_log}
echo "test option : $*"  >> ${test_log}

for i in ${path}
do
	echo "Test of ${i}"
	echo "***************************************************" >> ${test_log}
	echo "${i}" >> ${test_log}
	cd ${i}
	test_${i}.sh $*
	cat test.log >> ../${test_log}
	cd ..
done

echo "end of test"

