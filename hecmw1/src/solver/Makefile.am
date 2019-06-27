MPIDIR            = @mpidir@
MPIBINDIR         = @mpibindir@
MPILIBDIR         = @mpilibdir@
MPIINCDIR         = @mpiincdir@
MPILIBS           = @mpilibs@

PREFIX            = @prefix@
BINDIR            = @bindir@
LIBDIR            = @libdir@
INCLUDEDIR        = @includedir@

HECMWBINDIR       = ../../bin
HECMWLIBDIR       = ../../lib
HECMWINCDIR       = ../../include
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
	hecmw_solver.@f90objfilepostfix@
  
HEADERS = .

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


install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

@build_target@_subdir:
	@cd communication && $(MAKE) @build_target@
	@cd matrix        && $(MAKE) @build_target@
	@cd init          && $(MAKE) @build_target@
	@cd sparse_matrix && $(MAKE) @build_target@
	@cd las           && $(MAKE) @build_target@
	@cd mumps         && $(MAKE) @build_target@
	@cd mkl           && $(MAKE) @build_target@
	@cd clustermkl    && $(MAKE) @build_target@
	@cd precond       && $(MAKE) @build_target@
	@cd iterative     && $(MAKE) @build_target@
	@cd solver_direct && $(MAKE) @build_target@
	@cd solver_direct_parallel && $(MAKE) @build_target@
	@cd solver_direct_lag && $(MAKE) @build_target@
	@cd mpc           && $(MAKE) @build_target@

install:
	@cd communication && $(MAKE) install
	@cd matrix        && $(MAKE) install
	@cd init          && $(MAKE) install
	@cd sparse_matrix && $(MAKE) install
	@cd las           && $(MAKE) install
	@cd mumps         && $(MAKE) install
	@cd mkl           && $(MAKE) install
	@cd clustermkl    && $(MAKE) install
	@cd precond       && $(MAKE) install
	@cd iterative     && $(MAKE) install
	@cd solver_direct && $(MAKE) install
	@cd solver_direct_parallel && $(MAKE) install
	@cd solver_direct_lag && $(MAKE) install
	@cd mpc           && $(MAKE) install

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)
	@cd communication && $(MAKE) clean
	@cd matrix        && $(MAKE) clean
	@cd init          && $(MAKE) clean
	@cd sparse_matrix && $(MAKE) clean
	@cd las           && $(MAKE) clean
	@cd mumps         && $(MAKE) clean
	@cd mkl           && $(MAKE) clean
	@cd clustermkl    && $(MAKE) clean
	@cd precond       && $(MAKE) clean
	@cd iterative     && $(MAKE) clean
	@cd solver_direct && $(MAKE) clean
	@cd solver_direct_parallel && $(MAKE) clean
	@cd solver_direct_lag && $(MAKE) clean
	@cd mpc           && $(MAKE) clean

distclean: clean
	@cd communication && $(MAKE) distclean
	@cd matrix        && $(MAKE) distclean
	@cd init          && $(MAKE) distclean
	@cd sparse_matrix && $(MAKE) distclean
	@cd las           && $(MAKE) distclean
	@cd mumps         && $(MAKE) distclean
	@cd mkl           && $(MAKE) distclean
	@cd clustermkl    && $(MAKE) distclean
	@cd precond       && $(MAKE) distclean
	@cd iterative     && $(MAKE) distclean
	@cd solver_direct && $(MAKE) distclean
	@cd solver_direct_parallel && $(MAKE) distclean
	@cd solver_direct_lag && $(MAKE) distclean
	@cd mpc           && $(MAKE) distclean

maintainer-clean:
	@cd communication && $(MAKE) maintainer-clean
	@cd matrix        && $(MAKE) maintainer-clean
	@cd init          && $(MAKE) maintainer-clean
	@cd sparse_matrix && $(MAKE) maintainer_clean
	@cd las           && $(MAKE) maintainer-clean
	@cd mumps         && $(MAKE) maintainer_clean
	@cd mkl           && $(MAKE) maintainer_clean
	@cd clustermkl    && $(MAKE) maintainer_clean
	@cd precond       && $(MAKE) maintainer-clean
	@cd iterative     && $(MAKE) maintainer-clean
	@cd solver_direct && $(MAKE) maintainer-clean
	@cd solver_direct_parallel && $(MAKE) maintainer-clean
	@cd solver_direct_lag && $(MAKE) maintainer-clean
	@cd mpc           && $(MAKE) maintainer-clean
