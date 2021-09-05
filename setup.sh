date >> setup.log
echo "  setup option : " $@ >> setup.log

cp ./Makefile.conf ./hecmw1
cd ./hecmw1
echo "Executing hecmw1/setup_hecmw.sh"
./setup_hecmw.sh $@
cd ..

cp ./Makefile.conf ./fistr1
cd ./fistr1
echo "Executing fistr1/setup_fistr.sh"
./setup_fistr.sh $@
cd ..

echo "Executing /setup_fistr.sh"
./setup_fistr.sh $@

echo "" > ./fistr1/src/main/FrontISTRConfig.h
VERSION_MAJOR=$(cat VERSION | sed -r "s/^[^0-9]*(([0-9]+)\\.?([0-9]+)?\\.?([0-9]+)?).*$/\\1/"|cut -d. -f1|sed "s/^$/0/")
VERSION_MINOR=$(cat VERSION | sed -r "s/^[^0-9]*(([0-9]+)\\.?([0-9]+)?\\.?([0-9]+)?).*$/\\1/"|cut -d. -f2|sed "s/^$/0/")
VERSION_PATCH=$(cat VERSION | sed -r "s/^[^0-9]*(([0-9]+)\\.?([0-9]+)?\\.?([0-9]+)?).*$/\\1/"|cut -d. -f3|sed "s/^$/0/")
echo "#define VERSION_MAJOR $VERSION_MAJOR" >> ./fistr1/src/main/FrontISTRConfig.h
echo "#define VERSION_MINOR $VERSION_MINOR" >> ./fistr1/src/main/FrontISTRConfig.h
echo "#define VERSION_PATCH $VERSION_PATCH" >> ./fistr1/src/main/FrontISTRConfig.h
echo '#define BUILD_DATE "'`date '+%Y-%m-%dT%H:%M:%S%z'`'"' >> ./fistr1/src/main/FrontISTRConfig.h

if type "git" > /dev/null 2>&1
then
    hash=`git rev-parse HEAD 2> /dev/null`
else
    hash=0
fi
echo "#define GIT_HASH" \"$hash\" >> ./fistr1/src/main/FrontISTRConfig.h
echo "#define OPENMP_UNKNOWN" >> ./fistr1/src/main/FrontISTRConfig.h

for i in $*
do
	if [ "\"$i\"" = "\"-p\"" -o "\"$i\"" = "\"-parallel\"" -o "\"$i\"" = "\"--parallel\"" ]; then
    echo "#define WITH_MPI" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-metis\"" -o "\"$i\"" = "\"--with-metis\"" ]; then
    echo "#define WITH_METIS" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-parmetis\"" -o "\"$i\"" = "\"--with-parmetis\"" ]; then
    echo "#define WITH_PARMETIS" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-tools\"" -o "\"$i\"" = "\"--with-tools\"" ]; then
    echo "#define WITH_TOOLS" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-revocap\"" -o "\"$i\"" = "\"--with-revocap\"" ]; then
    echo "#define WITH_REVOCAP" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-refiner\"" -o "\"$i\"" = "\"--with-refiner\"" ]; then
    echo "#define WITH_REFINER" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-mkl\"" -o "\"$i\"" = "\"--with-mkl\"" ]; then
    echo "#define WITH_MKL" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-mumps\"" -o "\"$i\"" = "\"--with-mumps\"" ]; then
    echo "#define WITH_MUMPS" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-ml\"" -o "\"$i\"" = "\"--with-ml\"" ]; then
    echo "#define WITH_ML" >> ./fistr1/src/main/FrontISTRConfig.h
	elif [ "\"$i\"" = "\"-with-lapack\"" -o "\"$i\"" = "\"--with-lapack\"" ]; then
    echo "#define WITH_LAPACK" >> ./fistr1/src/main/FrontISTRConfig.h
  fi
done


