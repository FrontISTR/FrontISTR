## Caution) Opt-flags are removed

MPIDIR                  = @mpidir@
MPIBINDIR               = @mpibindir@
MPILIBDIR               = @mpilibdir@
MPIINCDIR               = @mpiincdir@
MPILIBS                 = @mpilibs@

HECMWDIR                = ../../hecmw
HECMWBINDIR             = $(HECMWDIR)/bin
HECMWLIBDIR             = $(HECMWDIR)/lib
HECMWINCDIR             = $(HECMWDIR)/include
HECMWLIBS               = @hecmwlibs@

FSTRBINDIR              = ../../bin
FSTRLIBDIR              = ../../lib
FSTRINCDIR              = ../../include
FSTRLIBS                = @fstrlibs@
FSTRLIB_TARGETFILE      = @fstrlib_targetfile@
FSTRLIB_F90TARGETFILE   = @fstrlib_f90targetfile@

PREFIX                  = @prefix@
BINDIR                  = @bindir@

DEBUG_OPTFLAGS          = @debug_optflags@
NODEBUG_OPTFLAGS        = @nodebug_optflags@

CC                      = @cc@
CFLAGS                  = @cflags@
BASE_CFLAGS             = @base_cflags@
MPI_CFLAGS              = @mpi_cflags@
HECMW_CFLAGS            = @hecmw_cflags@
FSTR_CFLAGS             = @fstr_cflags@
ALL_CFLAGS              = $(BASE_CFLAGS) $(FSTR_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                 = @ldflags@
MPI_LDFLAGS             = @mpi_ldflags@
HECMW_LDFLAGS           = @hecmw_ldflags@
FSTR_LDFLAGS            = @fstr_ldflags@
ALL_LDFLAGS             = 
OPTFLAGS                = @optflags@

F90                     = @f90@
F90FLAGS                = @f90flags@
BASE_F90FLAGS           = @base_f90flags@
MPI_F90FLAGS            = @mpi_f90flags@
HECMW_F90FLAGS          = @hecmw_f90flags@
FSTR_F90FLAGS           = @fstr_f90flags@
DIRECTSOLVER_F90FLAGS   = @directsolver_f90flags@
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(FSTR_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = 
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

TARGET                  = @fstrlib_targetfile@
F90TARGET               = @fstrlib_f90targetfile@

#------------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@f90objfilepostfix@ .f90 .f

.c.@cobjfilepostfix@:
	$(CC) -c $(ALL_CFLAGS) $(OPTFLAGS) $< -o $@

.f90.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

.f.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

#------------------------------------------------------------------------------#
OBJS = 

OBJSF = \
        fstr_lczparm.@f90objfilepostfix@ \
        m_step.@f90objfilepostfix@ \
        m_out.@f90objfilepostfix@ \
        m_fstr.@f90objfilepostfix@ \
        GaussM.@f90objfilepostfix@ \
        static_LIB_2d.@f90objfilepostfix@ \
        static_LIB_3d.@f90objfilepostfix@ \
        static_LIB_C3D8.@f90objfilepostfix@ \
        static_LIB_3dIC.@f90objfilepostfix@ \
        static_LIB_shell.@f90objfilepostfix@ \
        heat_LIB_CAPACITY.@f90objfilepostfix@ \
        heat_LIB_DFLUX.@f90objfilepostfix@ \
        heat_LIB_FILM.@f90objfilepostfix@ \
        heat_LIB_RADIATE.@f90objfilepostfix@ \
        heat_LIB_THERMAL.@f90objfilepostfix@ \
        eigen_LIB_2d1mass.@f90objfilepostfix@ \
        eigen_LIB_2d2mass.@f90objfilepostfix@ \
        eigen_LIB_3d1mass.@f90objfilepostfix@ \
        eigen_LIB_3d2mass.@f90objfilepostfix@ \
        precheck_LIB_2d.@f90objfilepostfix@ \
        precheck_LIB_3d.@f90objfilepostfix@ \
        precheck_LIB_shell.@f90objfilepostfix@ \
        solve_LINEQ.@f90objfilepostfix@ \
        static_LIB.@f90objfilepostfix@ \
        eigen_LIB.@f90objfilepostfix@ \
        heat_LIB.@f90objfilepostfix@

HEADERS = 

#------------------------------------------------------------------------------#
all: @build_target@

@build_curr@: $(TARGET) $(F90TARGET) install-header install-module

@build_target@: 
	@cd utilities && $(MAKE) @build_target@
	@cd user      && $(MAKE) @build_target@
	@cd element   && $(MAKE) @build_target@
	@cd physics   && $(MAKE) @build_target@
	@cd contact   && $(MAKE) @build_target@
	$(MAKE) @build_curr@

#$(TARGET): $(OBJS)
#	$(AR) $(FSTRLIBDIR)/$@ $(OBJS)
$(TARGET):

$(F90TARGET): $(OBJSF)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJSF)

#install-header:
#	$(CP) $(HEADERS) $(FSTRINCDIR)/.
install-header:

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(FSTRINCDIR)/.

#$(OBJS): $(HEADERS)
$(OBJS):

$(OBJSF):


install:

clean:
	@cd utilities  && $(MAKE) clean
	@cd user       && $(MAKE) clean
	@cd element && $(MAKE) clean
	@cd physics && $(MAKE) clean
	@cd contact && $(MAKE) clean
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean

maintainer-clean: distclean

