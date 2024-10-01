MPIDIR            = @mpidir@
MPIBINDIR         = @mpibindir@
MPILIBDIR         = @mpilibdir@
MPIINCDIR         = @mpiincdir@
MPILIBS           = @mpilibs@

PREFIX            = @prefix@
BINDIR            = @bindir@
LIBDIR            = @libdir@
INCLUDEDIR        = @includedir@

HECMWBINDIR       = ../../../bin
HECMWLIBDIR       = ../../../lib
HECMWINCDIR       = ../../../include
HECMWLIBS         = @hecmwlibs@

REFINERDIR        = @refinerdir@
REFINERINCDIR     = @refinerincdir@
REFINERLIBDIR     = @refinerlibdir@
REFINERLIBS       = @refinerlibs@

CC                = @cc@
CFLAGS            = @cflags@
BASE_CFLAGS       = @base_cflags@
MPI_CFLAGS        = @mpi_cflags@
HECMW_CFLAGS      = @hecmw_cflags@
ALL_CFLAGS        = $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS           = @ldflags@
MPI_LDFLAGS       = @mpi_ldflags@
HECMW_LDFLAGS     = @hecmw_ldflags@
ALL_LDFLAGS       = $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS          = @optflags@

F90               = @f90@
F90FLAGS          = @f90flags@
BASE_F90FLAGS     = @base_f90flags@
MPI_F90FLAGS      = @mpi_f90flags@
HECMW_F90FLAGS    = @hecmw_f90flags@
ALL_F90FLAGS      = $(BASE_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS)
F90LDFLAGS        = @f90ldflags@
MPI_F90LDFLAGS    = @mpi_f90ldflags@
HECMW_F90LDFLAGS  = @hecmw_f90ldflags@
ALL_F90LDFLAGS    = $(HECMW_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS       = @f90optflags@

AR                = @ar@
CP                = @cp@
RM                = @rm@
RANLIB            = @ranlib@
MKDIR             = @mkdir@
MAKE = @make@

#-----------------------------------------------------------------------------#
TARGET            = @hecmwlib_targetfile@
F90TARGET         = @hecmwlib_f90targetfile@

F90MODULEPOSTFIX  = @f90modulepostfix@
COBJFILEPOSTFIX   = @cobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ @f90objfilepostfix@ .c .f90 .f

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

.f90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

.f.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS =

OBJSF = \
	hecmw_precond.@f90objfilepostfix@
  
HEADERS =

#-----------------------------------------------------------------------------#
all: @build_target@

@build_target@:
	@cd ssor   && $(MAKE) @build_target@
	@cd diag   && $(MAKE) @build_target@
	@cd bilu   && $(MAKE) @build_target@
	@cd sainv  && $(MAKE) @build_target@
	@cd rif    && $(MAKE) @build_target@
	@cd ml     && $(MAKE) @build_target@
	$(MAKE) BUILD_CURR

build-without-f: $(TARGET)

BUILD_CURR: $(TARGET) $(F90TARGET) install-module

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSF):

install-header:
	@cd ml     && $(MAKE) install-header

install-module: $(OBJSF)
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:
	@cd ssor   && $(MAKE) install
	@cd diag   && $(MAKE) install
	@cd bilu   && $(MAKE) install
	@cd sainv  && $(MAKE) install
	@cd rif    && $(MAKE) install
	@cd ml     && $(MAKE) install

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)
	@cd ssor   && $(MAKE) clean
	@cd diag   && $(MAKE) clean
	@cd bilu   && $(MAKE) clean
	@cd sainv  && $(MAKE) clean
	@cd rif    && $(MAKE) clean
	@cd ml     && $(MAKE) clean
  
distclean: clean
	@cd ssor   && $(MAKE) distclean
	@cd diag   && $(MAKE) distclean
	@cd bilu   && $(MAKE) distclean
	@cd sainv  && $(MAKE) distclean
	@cd rif    && $(MAKE) distclean
	@cd ml     && $(MAKE) distclean

maintainer-clean:
	@cd ssor   && $(MAKE) maintainer-clean
	@cd diag   && $(MAKE) maintainer-clean
	@cd bilu   && $(MAKE) maintainer-clean
	@cd sainv  && $(MAKE) maintainer-clean
	@cd rif    && $(MAKE) maintainer-clean
	@cd ml     && $(MAKE) maintainer-clean

include .dependf90

depend .dependf90:
	makedepf90 *.[fF]90 > .dependf90
	@cd ssor   && $(MAKE) depend
	@cd diag   && $(MAKE) depend
	@cd bilu   && $(MAKE) depend
	@cd sainv  && $(MAKE) depend
	@cd rif    && $(MAKE) depend
	@cd ml     && $(MAKE) depend
