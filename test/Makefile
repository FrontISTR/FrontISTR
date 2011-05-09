#==========================================================
#
# Software Name :HEC middleware Ver. 3.0beta
#
#   Makelfile for MW3
#
#                     Written by T.Takeda,    2010/06/01
#                                K.Goto,      2010/01/12
#                                K.Matsubara, 2010/06/01
#
#   Contact address : IIS, The University of Tokyo CISS
#
#==========================================================


include ../Makefile.in

OBJ_MAIN=MW_Str.o

default: $(EXE)

$(EXE): $(OBJ_MAIN) $(LIB)
	$(CPC) $(LOPT) $(LLIB) -o $(EXE) $(OBJ_MAIN)

testall: test1 test2 test3 test4

test1:
	echo '=begin' > mw3.cnt
	echo 'this is test data1.' >> mw3.cnt
	echo '=end' >> mw3.cnt
	echo 'MeshFileName' >> mw3.cnt
	echo 'beam' >> mw3.cnt
	echo 'End' >> mw3.cnt
	$(RUN) -np 1 $(EXE)

test2:
	echo '=begin' > mw3.cnt
	echo 'this is test data1.' >> mw3.cnt
	echo '=end' >> mw3.cnt
	echo 'MeshFileName' >> mw3.cnt
	echo 'beamMPC' >> mw3.cnt
	echo 'End' >> mw3.cnt
	$(RUN) -np 1 $(EXE)

test3:
	echo '=begin' > mw3.cnt
	echo 'this is test data1.' >> mw3.cnt
	echo '=end' >> mw3.cnt
	echo 'MeshFileName' >> mw3.cnt
	echo 'beamPara2' >> mw3.cnt
	echo 'End' >> mw3.cnt
	$(RUN) -np 2 $(EXE)

test4:
	echo '=begin' > mw3.cnt
	echo 'this is test data1.' >> mw3.cnt
	echo '=end' >> mw3.cnt
	echo 'MeshFileName' >> mw3.cnt
	echo 'beamPara4' >> mw3.cnt
	echo 'End' >> mw3.cnt
	$(RUN) -np 4 $(EXE)

run1:
	$(RUN) -np 1 $(EXE)

run2:
	$(RUN) -np 2 $(EXE)

run4:
	$(RUN) -np 4 $(EXE)

run8:
	$(RUN) -np 8 $(EXE)

clean:
	rm -rf *~ *.o *.out *UCD*.inp hec*log $(EXE)
