MPIDIR            = @mpidir@
MPIBINDIR         = @mpibindir@
MPILIBDIR         = @mpilibdir@
MPIINCDIR         = @mpiincdir@
MPILIBS           = @mpilibs@

PREFIX            = @prefix@
BINDIR            = @bindir@
LIBDIR            = @libdir@
INCLUDEDIR        = @includedir@

NEUINCDIR	      = ../include
NEULIBDIR	      = ../lib
NEULIBS           = -lnfd -lhecd -lconv
NEU_CPPFLAGS      = -I$(NEUINCDIR)
NEU_CPPLDFLAGS    = -L$(NEULIBDIR) $(NEULIBS)

CC                = @cc@
CFLAGS            = @cflags@
BASE_CFLAGS       = @base_cflags@
MPI_CFLAGS        = @mpi_cflags@
ALL_CFLAGS        = $(BASE_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS           = @ldflags@
MPI_LDFLAGS       = @mpi_ldflags@
ALL_LDFLAGS       = $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS          = @optflags@

CPP               = @cpp@
CPPFLAGS          = @cppflags@
BASE_CPPFLAGS     = @base_cppflags@
MPI_CPPFLAGS      = @mpi_f90flags@
ALL_CPPFLAGS      = $(BASE_CPPFLAGS) $(MPI_CPPFLAGS) $(CPPFLAGS)
CPPLDFLAGS        = @cppldflags@
MPI_CPPLDFLAGS    = @mpi_f90ldflags@
ALL_CPPLDFLAGS    = $(MPI_CPPLDFLAGS) $(CPPLDFLAGS)
CPPOPTFLAGS       = @cppoptflags@

F90               = @f90@
F90FLAGS          = @f90flags@
BASE_F90FLAGS     = @base_f90flags@
MPI_F90FLAGS      = @mpi_f90flags@
ALL_F90FLAGS      = $(BASE_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS)
F90LDFLAGS        = @f90ldflags@
MPI_F90LDFLAGS    = @mpi_f90ldflags@
ALL_F90LDFLAGS    = $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS       = @f90optflags@

AR                = @ar@
CP                = @cp@
RM                = @rm@
RANLIB            = @ranlib@
MKDIR             = @mkdir@

#-----------------------------------------------------------------------------#
TARGET            = libconv.a

F90MODULEPOSTFIX  = @f90modulepostfix@
COBJFILEPOSTFIX   = @cobjfilepostfix@
CPPOBJFILEPOSTFIX = @cobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@cppobjfilepostfix@ .@f90objfilepostfix@ .c .cpp .f90

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

.cpp.@cppobjfilepostfix@:
	$(CPP) -c $< $(ALL_CPPFLAGS) $(CPPOPTFLAGS) $(NEU_CPPFLAGS)

.f90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS =

OBJSP = \
	cconv_mat.@cppobjfilepostfix@ \
	CConvMessage.@cppobjfilepostfix@ \
	conv_neu2hec.@cppobjfilepostfix@ \
	conv_neu2fstr_static.@cppobjfilepostfix@ \
	conv_neu2fstr_dynamic.@cppobjfilepostfix@ \
	conv_neu2fstr_heat.@cppobjfilepostfix@

OBJSF =

HEADERS = \
	cconv_mat.h \
	CConvMessage.h \
	conv_neu2hec.h \
	conv_neu2fstr_static.h \
	conv_neu2fstr_dynamic.h \
	conv_neu2fstr_heat.h \
	conv_util.h

#-----------------------------------------------------------------------------#
all: @build_target@

@build_target@: $(TARGET) install-header

$(TARGET): $(OBJS) $(OBJSP)
	$(AR) $(NEULIBDIR)/$@ $(OBJS)
	$(AR) $(NEULIBDIR)/$@ $(OBJSP)

$(F90TARGET): $(OBJSF)
	$(AR) $(NEULIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSP): $(HEADERS)

$(OBJSF):

install-header: $(HEADERS)
	$(CP) $(HEADERS) $(NEUINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(MEUINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(CPPOBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
