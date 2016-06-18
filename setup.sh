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
