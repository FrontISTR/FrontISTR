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


include ./Makefile.in

default:
	(cd src ; make )
	(cd test ; make ) 

check:
	(cd test ; make testall )

clean:
	rm -f core *~ $(MW3_HOME)/lib/libmw3.a
	(cd src ; make clean )
	(cd test ; make clean )

