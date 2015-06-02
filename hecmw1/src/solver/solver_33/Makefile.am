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

MLDIR             = @mldir@
MLLIBDIR          = @mllibdir@
MLINCDIR          = @mlincdir@
MLLIBS            = @mllibs@
ML_CFLAGS         = @ml_cflags@
ML_LDFLAGS        = @ml_ldflags@
ML_F90FLAGS       = @ml_f90flags@
ML_F90LDFLAGS     = @ml_f90ldflags@

CC                = @cc@
CFLAGS            = @cflags@
BASE_CFLAGS       = @base_cflags@
MPI_CFLAGS        = @mpi_cflags@
HECMW_CFLAGS      = @hecmw_cflags@
ALL_CFLAGS        = $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(ML_CFLAGS) $(CFLAGS)
LDFLAGS           = @ldflags@
MPI_LDFLAGS       = @mpi_ldflags@
HECMW_LDFLAGS     = @hecmw_ldflags@
ALL_LDFLAGS       = $(HECMW_LDFLAGS) $(ML_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS          = @optflags@

F90               = @f90@
F90FLAGS          = @f90flags@
BASE_F90FLAGS     = @base_f90flags@
MPI_F90FLAGS      = @mpi_f90flags@
HECMW_F90FLAGS    = @hecmw_f90flags@
ALL_F90FLAGS      = $(BASE_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(ML_F90FLAGS) $(F90FLAGS)
F90LDFLAGS        = @f90ldflags@
MPI_F90LDFLAGS    = @mpi_f90ldflags@
HECMW_F90LDFLAGS  = @hecmw_f90ldflags@
ALL_F90LDFLAGS    = $(HECMW_F90LDFLAGS) $(ML_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS       = @f90optflags@

AR                = @ar@
CP                = @cp@
RM                = @rm@
RANLIB            = @ranlib@
MKDIR             = @mkdir@

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
OBJS = \
	hecmw_ML_wrapper.@cobjfilepostfix@ \

OBJSF = \
	hecmw_jadm.@f90objfilepostfix@ \
	hecmw_tuning_fx.@f90objfilepostfix@ \
	hecmw_solver_las_33.@f90objfilepostfix@ \
	hecmw_solver_scaling_33.@f90objfilepostfix@ \
	hecmw_precond_BILU_33.@f90objfilepostfix@ \
	hecmw_precond_DIAG_33.@f90objfilepostfix@ \
	hecmw_precond_SSOR_33.@f90objfilepostfix@ \
	hecmw_precond_SAINV_33.@f90objfilepostfix@ \
	hecmw_precond_RIF_33.@f90objfilepostfix@ \
	hecmw_mat_id.@f90objfilepostfix@ \
	hecmw_ML_helper.@f90objfilepostfix@ \
	hecmw_precond_ML_33.@f90objfilepostfix@ \
	hecmw_precond_33.@f90objfilepostfix@ \
	hecmw_solver_CG_33.@f90objfilepostfix@ \
	hecmw_solver_BiCGSTAB_33.@f90objfilepostfix@ \
	hecmw_solver_GMRES_33.@f90objfilepostfix@ \
	hecmw_solver_GPBiCG_33.@f90objfilepostfix@ \
	hecmw_solver_33.@f90objfilepostfix@

OBJSF_LAPACK = hecmw_estimate_condition.@f90objfilepostfix@

OBJSF_LAPACK_DUMM = hecmw_estimate_condition_dumm.@f90objfilepostfix@

HEADERS =

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-module

build-serial: $(TARGET) $(F90TARGET) install-module

build-with-lapack: $(TARGET) F90TARGET_LAPACK install-module

build-serial-with-lapack: $(TARGET) F90TARGET_LAPACK install-module

build-without-f: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF_LAPACK_DUMM) $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF_LAPACK_DUMM) $(OBJSF)

F90TARGET_LAPACK: $(OBJSF_LAPACK) $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$(F90TARGET) $(OBJSF_LAPACK) $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSF):

install-header: $(HEADERS)
	$(CP) $(HEADERS) $(HECMWINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
