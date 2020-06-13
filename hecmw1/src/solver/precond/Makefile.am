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

build-default: @build_target@_subdir $(TARGET) $(F90TARGET) install-module

build-serial: $(TARGET) $(F90TARGET) install-module

build-without-f: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSF):

install-header:
	@cd nn     && $(MAKE) install-header
	@cd 33     && $(MAKE) install-header

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

@build_target@_subdir:
	@cd nn     && $(MAKE) @build_target@
	@cd 11     && $(MAKE) @build_target@
	@cd 22     && $(MAKE) @build_target@
	@cd 33     && $(MAKE) @build_target@
	@cd 44     && $(MAKE) @build_target@
	@cd 66     && $(MAKE) @build_target@
  
install:
	@cd nn     && $(MAKE) install
	@cd 11     && $(MAKE) install
	@cd 22     && $(MAKE) install
	@cd 33     && $(MAKE) install
	@cd 44     && $(MAKE) install
	@cd 66     && $(MAKE) install

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)
	@cd nn     && $(MAKE) clean
	@cd 11     && $(MAKE) clean
	@cd 22     && $(MAKE) clean
	@cd 33     && $(MAKE) clean
	@cd 44     && $(MAKE) clean
	@cd 66     && $(MAKE) clean
  
distclean: clean
	@cd nn     && $(MAKE) distclean
	@cd 11     && $(MAKE) distclean
	@cd 22     && $(MAKE) distclean
	@cd 33     && $(MAKE) distclean
	@cd 44     && $(MAKE) distclean
	@cd 66     && $(MAKE) distclean

maintainer-clean:
	@cd nn     && $(MAKE) maintainer-clean
	@cd 11     && $(MAKE) maintainer-clean
	@cd 22     && $(MAKE) maintainer-clean
	@cd 33     && $(MAKE) maintainer-clean
	@cd 44     && $(MAKE) maintainer-clean
	@cd 66     && $(MAKE) maintainer-clean
