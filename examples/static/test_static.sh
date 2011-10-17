#!/bin/sh
# test_static.sh ver.1.0
# 2006.04.07 by N.Imai
# ----------------------------
# FrontSTR Test for Examples

test_log="test_static.log"

path="\
  exA\
  exB\
  exC\
  exD\
  exE\
  exF\
  exG\
  exX\
  1elem"

rm -f ${test_log}

echo "FrontSTR test result for static analysis" >> ${test_log}
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

