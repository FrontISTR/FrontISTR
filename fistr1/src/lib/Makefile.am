MPIDIR                  = @mpidir@
MPIBINDIR               = @mpibindir@
MPILIBDIR               = @mpilibdir@
MPIINCDIR               = @mpiincdir@
MPILIBS                 = @mpilibs@

HECMWDIR                = ../../../hecmw1
HECMWBINDIR             = $(HECMWDIR)/bin
HECMWLIBDIR             = $(HECMWDIR)/lib
HECMWINCDIR             = $(HECMWDIR)/include
HECMWLIBS               = @hecmwlibs@

FSTRBINDIR              = ../../bin
FSTRLIBDIR              = ../../lib
FSTRINCDIR              = ../../include
FSTRLIBS                = @fstrlibs@

METISDIR                = @metisdir@
METISLIBDIR             = @metislibdir@
METISINCDIR             = @metisincdir@
METISLIBS               = @metislibs@
METIS_CFLAGS            = @metis_cflags@
METIS_LDFLAGS           = @metis_ldflags@
METIS_F90FLAGS          = @metis_f90flags@
METIS_F90LDFLAGS        = @metis_f90ldflags@
HECMW_METIS_VER         = -DHECMW_METIS_VER=@hecmw_metis_ver@

MUMPSDIR                = @mumpsdir@
MUMPSLIBDIR             = @mumpslibdir@
MUMPSINCDIR             = @mumpsincdir@
MUMPSLIBS               = @mumpslibs@
MUMPS_CFLAGS            = @mumps_cflags@
MUMPS_LDFLAGS           = @mumps_ldflags@
MUMPS_F90FLAGS          = @mumps_f90flags@
MUMPS_F90LDFLAGS        = @mumps_f90ldflags@

REVOCAPDIR              = @revocapdir@
REVOCAPINCDIR           = @revocapincdir@
REVOCAPLIBDIR           = @revocaplibdir@
REVOCAPLIBS             = @revocaplibs@

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
ALL_LDFLAGS             = $(FSTR_LDFLAGS) $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS                = @optflags@

REVOCAP_F90FLAGS        = @revocap_f90flags@
REVOCAP_F90LDFLAGS      = @revocap_f90ldflags@

F90                     = @f90@
F90FLAGS                = @f90flags@
BASE_F90FLAGS           = @base_f90flags@
MPI_F90FLAGS            = @mpi_f90flags@
HECMW_F90FLAGS          = @hecmw_f90flags@
FSTR_F90FLAGS           = @fstr_f90flags@
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(FSTR_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS) $(REVOCAP_F90FLAGS) $(MUMPS_F90FLAGS) $(METIS_F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = $(FSTR_F90LDFLAGS) $(HECMW_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS) $(MUMPS_F90LDFLAGS)
F90FPP                  = @f90fpp@
F90OPTFLAGS             = @f90optflags@

AR                      = @ar@
CP                      = @cp@
RM                      = @rm@
RANLIB                  = @ranlib@
MKDIR                   = @mkdir@

#------------------------------------------------------------------------------#
TARGET                  = @fstrlib_targetfile@
F90TARGET               = @fstrlib_f90targetfile@

F90MODULEPOSTFIX        = @f90modulepostfix@
COBJFILEPOSTFIX         = @cobjfilepostfix@
F90OBJFILEPOSTFIX       = @f90objfilepostfix@

#------------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@f90objfilepostfix@ .c .f90 .F90 .f

.c.@cobjfilepostfix@:
	$(CC) -c $(ALL_CFLAGS) $(OPTFLAGS) $< -o $@

.f90.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

.F90.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(HECMW_METIS_VER) $(F90FPP) $(F90OPTFLAGS) $< -o $@

.f.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

#------------------------------------------------------------------------------#
OBJS =

OBJSF = \
        m_step.@f90objfilepostfix@ \
        m_timepoint.@f90objfilepostfix@ \
        m_out.@f90objfilepostfix@ \
        m_common_struct.@f90objfilepostfix@ \
        m_fstr.@f90objfilepostfix@ \
        GaussM.@f90objfilepostfix@ \
        static_LIB_1d.@f90objfilepostfix@ \
        static_LIB_2d.@f90objfilepostfix@ \
        static_LIB_3d.@f90objfilepostfix@ \
        static_LIB_3d_vp.@f90objfilepostfix@ \
        static_LIB_C3D8.@f90objfilepostfix@ \
        static_LIB_Fbar.@f90objfilepostfix@ \
        static_LIB_3dIC.@f90objfilepostfix@ \
        static_LIB_beam.@f90objfilepostfix@ \
        static_LIB_shell.@f90objfilepostfix@ \
        dynamic_mass.@f90objfilepostfix@ \
        heat_LIB_CAPACITY.@f90objfilepostfix@ \
        heat_LIB_DFLUX.@f90objfilepostfix@ \
        heat_LIB_FILM.@f90objfilepostfix@ \
        heat_LIB_NEUTRAL.@f90objfilepostfix@ \
        heat_LIB_RADIATE.@f90objfilepostfix@ \
        heat_LIB_THERMAL.@f90objfilepostfix@ \
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

@build_target@:
	@cd utilities && $(MAKE) @build_target@
	@cd user      && $(MAKE) @build_target@
	@cd element   && $(MAKE) @build_target@
	@cd physics   && $(MAKE) @build_target@
	@cd contact   && $(MAKE) @build_target@
	$(MAKE) BUILD_CURR

build-with-rcap: @build_target@

BUILD_CURR: $(TARGET) $(F90TARGET) install-module

$(TARGET): $(OBJS)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJSF)
	F90="$(F90)" F90FLAGS="$(F90FLAGS)" F90TARGET="$(FSTRLIBDIR)/$(F90TARGET)" AR="$(AR)" RM="$(RM)" sh test_flush.sh

install-header:
	$(CP) $(HEADERS) $(FSTRINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(FSTRINCDIR)/.

$(OBJS): $(HEADERS)

$(OBJSF):

install:

clean:
	@cd utilities  && $(MAKE) clean
	@cd user       && $(MAKE) clean
	@cd element    && $(MAKE) clean
	@cd physics    && $(MAKE) clean
	@cd contact    && $(MAKE) clean
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)
	${RM} test_flush.log

distclean: clean
