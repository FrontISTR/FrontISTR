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
HECMWBINDIR2      = ../../../hecmw1/bin
HECMWLIBDIR2      = ../../../hecmw1/lib
HECMWINCDIR2      = ../../../hecmw1/include
HECMWLIBS2        = @hecmwlibs2@

REFINERDIR        = @refinerdir@
REFINERINCDIR     = @refinerincdir@
REFINERLIBDIR     = @refinerlibdir@
REFINERLIBS       = @refinerlibs@

CC                = @cc@
CFLAGS            = @cflags@
BASE_CFLAGS       = @base_cflags@
MPI_CFLAGS        = @mpi_cflags@
HECMW_CFLAGS      = @hecmw_cflags@
HECMW_CFLAGS2     = @hecmw_cflags2@
ALL_CFLAGS        = $(BASE_CFLAGS) $(HECMW_CFLAGS) $(HECMW_CFLAGS2) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS           = @ldflags@
MPI_LDFLAGS       = @mpi_ldflags@
HECMW_LDFLAGS     = @hecmw_ldflags@
HECMW_LDFLAGS2    = @hecmw_ldflags2@
ALL_LDFLAGS       = $(HECMW_LDFLAGS) $(HECMW_LDFLAGS2) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS          = @optflags@

CPP               = @cpp@
CPPFLAGS          = @cppflags@
BASE_CPPFLAGS     = @base_cppflags@
HECMW_CPPFLAGS    = @hecmw_cppflags@
ALL_CPPFLAGS      = $(BASE_CPPFLAGS) $(HECMW_CPPFLAGS) $(MPI_CFLAGS) $(CPPFLAGS)
CPPLDFLAGS        = @cppldflags@
HECMW_CPPLDFLAGS  = @hecmw_cppldflags@
ALL_CPPLDFLAGS    = $(HECMW_CPPLDFLAGS) $(MPI_LDFLAGS) $(CPPLDFLAGS)
CPPOPTFLAGS       = @cppoptflags@

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

#-----------------------------------------------------------------------------#
TARGET            = @hecmwlib_targetfile@
F90TARGET         = @hecmwlib_f90targetfile@

F90MODULEPOSTFIX  = @f90modulepostfix@
COBJFILEPOSTFIX   = @cobjfilepostfix@
CPPOBJFILEPOSTFIX = @cppobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@cppobjfilepostfix@ .@f90objfilepostfix@ .c .cpp .f90

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

.cpp.@cppobjfilepostfix@:
	$(CPP) -c $< $(ALL_CPPFLAGS) $(CPPOPTFLAGS)

.f90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS = \
	MW3_visualizer.@cobjfilepostfix@

OBJSP =

OBJSF =

HEADERS = \
	MW3_visualizer.h

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-header

build-serial: $(TARGET) $(F90TARGET) install-header

build-without-f: $(TARGET) install-header

$(TARGET): $(OBJS) $(OBJSP)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSP)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSP): $(HEADERS)

$(OBJSF):

install-header: $(HEADERS)
	$(CP) $(HEADERS) $(HECMWINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(CPPOBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
