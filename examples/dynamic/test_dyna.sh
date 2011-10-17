#!/bin/sh
# test_eigen.sh ver.1.0
# 2010.03.28 Xi YUAN
# ----------------------------
# FrontISTR Test for Examples

test_log="test_dyna.log"

path="\
	exW\
	exX"

rm -f ${test_log}

echo "FrontSTR test result for dynamic analysis" >> ${test_log}
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

