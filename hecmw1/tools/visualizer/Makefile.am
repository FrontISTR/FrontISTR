MPIDIR                 = @mpidir@
MPIBINDIR              = @mpibindir@
MPILIBDIR              = @mpilibdir@
MPIINCDIR              = @mpiincdir@
MPILIBS                = @mpilibs@

PREFIX                 = @prefix@
BINDIR                 = @bindir@
LIBDIR                 = @libdir@
INCLUDEDIR             = @includedir@

HECMWBINDIR            = ../../bin
HECMWLIBDIR            = ../../lib
HECMWINCDIR            = ../../include
HECMWLIBS              = @hecmwlibs@

REFINERDIR             = @refinerdir@
REFINERINCDIR          = @refinerincdir@
REFINERLIBDIR          = @refinerlibdir@
REFINERLIBS            = @refinerlibs@

CC                     = @cc@
CFLAGS                 = @cflags@
BASE_CFLAGS            = @base_cflags@
MPI_CFLAGS             = @mpi_cflags@
HECMW_CFLAGS           = @hecmw_cflags@
VISUALIZER_CFLAGS      = @visualizer_cflags@
ALL_CFLAGS             = $(VISUALIZER_CFLAGS) $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                = @ldflags@
MPI_LDFLAGS            = @mpi_ldflags@
HECMW_LDFLAGS          = @hecmw_ldflags@
VISUALIZER_LDFLAGS     = @visualizer_ldflags@
ALL_LDFLAGS            = $(VISUALIZER_LDFLAGS) $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS               = @optflags@
CLINKER                = @clinker@

F90                    = @f90@
F90FLAGS               = @f90flags@
BASE_F90FLAGS          = @base_f90flags@
MPI_F90FLAGS           = @mpi_f90flags@
HECMW_F90FLAGS         = @hecmw_f90flags@
ALL_F90FLAGS           = $(BASE_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS)
F90LDFLAGS             = @f90ldflags@
MPI_F90LDFLAGS         = @mpi_f90ldflags@
HECMW_F90LDFLAGS       = @hecmw_f90ldflags@
ALL_F90LDFLAGS         = $(HECMW_F90LDLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS            = @f90optflags@

AR                     = @ar@
CP                     = @cp@
RM                     = @rm@
RANLIB                 = @ranlib@
MKDIR                  = @mkdir@

#-----------------------------------------------------------------------------#
TARGET                 = @visualizer_targetfile@

COBJFILEPOSTFIX        = @cobjfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .c

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS = \
	visualizer.@cobjfilepostfix@

HEADERS =

#-----------------------------------------------------------------------------#
all: @build_target@

@build_target@: $(TARGET)

$(TARGET): $(OBJS) $(HECMWLIBDIR)/*.*
	$(CLINKER) $(ALL_CFLAGS) $(OPTFLAGS) \
	-o $@ $(OBJS) $(ALL_LDFLAGS)
	$(CP) $(TARGET) $(HECMWBINDIR)/.

$(OBJS): $(HEADERS)

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)
