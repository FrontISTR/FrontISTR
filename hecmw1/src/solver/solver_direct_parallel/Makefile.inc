#======================================================================!
#                                                                      !
#   Software Name : HEC-MW Library for PC-cluster                      !
#         Version : 1.00                                               !
#                                                                      !
#     Last Update : 2007/11/21                                         !
#        Category : Makefile
#                                                                      !
#            Written by Takeshi Kitayama  (Univ. of Tokyo)             !
#                                                                      !
#     Contact address :  IIS,The University of Tokyo RSS21 project     !
#                                                                      !
#     "Structural Analysis System for General-purpose Coupling         !
#      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
#                                                                      !
#======================================================================!

# define compiler settings
# included by all Makefile

#FC=pgf90

# options for intel compilar
FFLAGS= -quiet -static -O3
#FFLAGS= -quiet -static -check all -traceback -g
FC=mpif90 $(FFLAGS) # compilar is defined via $MPI

.SUFFIXES: .f .f90
.f.o:
	$(FC) -c $<
.f90.o:
	$(FC) -c $<

