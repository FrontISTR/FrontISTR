MPIDIR                  = @mpidir@
MPIBINDIR               = @mpibindir@
MPILIBDIR               = @mpilibdir@
MPIINCDIR               = @mpiincdir@
MPILIBS                 = @mpilibs@

HECMWDIR                = ../../../hecmw/src_dist
HECMWBINDIR             = $(HECMWDIR)/src
HECMWLIBDIR             = $(HECMWDIR)/lib
HECMWINCDIR             = $(HECMWDIR)/src
HECMWLIBS               = @hecmwlibs@

BOOST_HOME     = /C/boost_1_44

FSTRBINDIR              = ../../bin
FSTRLIBDIR              = ../../lib
FSTRINCDIR              = ../../include
FSTRLIBS                = @fstrlibs@
FSTRLIB_TARGETFILE      = @fstrlib_targetfile@
FSTRLIB_F90TARGETFILE   = @fstrlib_f90targetfile@

REVOCAPDIR              = /C/MyPrograms/REVOCAP_Refiner-1.0.10
REVOCAPINCDIR           = $(REVOCAPDIR)/Refiner
REVOCAPLIBDIR           = $(REVOCAPDIR)/lib/i86_64
REVOCAPLIBS             = -lRcapRefiner

PREFIX                  = @prefix@
BINDIR                  = @bindir@

DEBUG_OPTFLAGS          = @debug_optflags@
NODEBUG_OPTFLAGS        = @nodebug_optflags@

CC                      = @cc@
CPP                     = @cpp@
CFLAGS                  = @cflags@
BASE_CFLAGS             = @base_cflags@
MPI_CFLAGS              = @mpi_cflags@
HECMW_CFLAGS            = @hecmw_cflags@ -I$(BOOST_HOME)
FSTR_CFLAGS             = @fstr_cflags@
ALL_CFLAGS              = $(BASE_CFLAGS) $(FSTR_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                 = @ldflags@
MPI_LDFLAGS             = @mpi_ldflags@
HECMW_LDFLAGS           = @hecmw_ldflags@
FSTR_LDFLAGS            = @fstr_ldflags@
ALL_LDFLAGS             = $(FSTR_LDFLAGS) $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS                = @optflags@

RCAP_F90FLAGS           = @rcap_f90flags@
RCAP_F90LDFLAGS         = @rcap_f90ldflags@

F90                     = @f90@
F90FLAGS                = @f90flags@
BASE_F90FLAGS           = @base_f90flags@
MPI_F90FLAGS            = @mpi_f90flags@
HECMW_F90FLAGS          = @hecmw_f90flags@
FSTR_F90FLAGS           = @fstr_f90flags@
DIRECTSOLVER_F90FLAGS   = @directsolver_f90flags@
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(F90FLAGS) \
                          $(MPI_F90FLAGS) $(HECMW_F90FLAGS) $(FSTR_F90FLAGS) $(RCAP_F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = $(FSTR_F90LDFLAGS) $(HECMW_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS             = @f90optflags@

AR                      = @ar@
CP                      = @cp@
RM                      = @rm@
RANLIB                  = @ranlib@
MKDIR                   = @mkdir@

#------------------------------------------------------------------------------#
F90MODULEPOSTFIX        = @f90modulepostfix@
COBJFILEPOSTFIX         = @cobjfilepostfix@
F90OBJFILEPOSTFIX       = @f90objfilepostfix@

TARGET                  = @fstrexec_targetfile@

#------------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@cppobjfilepostfix@ .@f90objfilepostfix@ .f90 .f .cpp

.c.@cobjfilepostfix@:
	$(CC) -c $(ALL_CFLAGS) $(OPTFLAGS) $< -o $@

.cpp.@cppobjfilepostfix@:
	$(CPP) -c $(ALL_CFLAGS) $(OPTFLAGS) $< -o $@

.f90.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

.f.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

#------------------------------------------------------------------------------#
OBJS = 

OBJSF = \
	istr.@f90objfilepostfix@ \
	../lib/*.@f90objfilepostfix@ \
	../common/*.@f90objfilepostfix@ \
	../analysis/static/*.@f90objfilepostfix@ \
	../analysis/dynamic/mode/*.@f90objfilepostfix@ \
	../analysis/dynamic/transit/*.@f90objfilepostfix@ \
	../analysis/heat/*.@f90objfilepostfix@

HEADERS = 

#------------------------------------------------------------------------------#
all: @build_target@

build-default: $(OBJSF) $(OBJS)
	$(F90) $(ALL_F90FLAGS) $(F90OPTFLAGS) \
	-o $(TARGET) \
	$(OBJSF) $(OBJS) $(ALL_F90LDFLAGS) -lstdc++
	$(CP) $(TARGET) $(FSTRBINDIR)/.

build-with-rcap: $(OBJSF)
	$(F90) $(ALL_F90FLAGS) $(F90OPTFLAGS) \
	-o $(TARGET) \
	$(OBJSF) $(OBJS) $(ALL_F90LDFLAGS) $(RCAP_F90LDFLAGS) -lstdc++
	$(CP) $(TARGET) $(FSTRBINDIR)/.

$(OBJS): $(HEADERS)

$(OBJSF):


install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)

maintainer-clean: distclean

