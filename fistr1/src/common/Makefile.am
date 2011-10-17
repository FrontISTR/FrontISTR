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
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(FSTR_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(F90FLAGS) $(REVOCAP_F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = $(FSTR_F90LDFLAGS) $(HECMW_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS) $(REVOCAP_F90LDFLAGS)
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
.SUFFIXES: .@cobjfilepostfix@ .@f90objfilepostfix@ .c .f90 .f

.c.@cobjfilepostfix@:
	$(CC) -c $(ALL_CFLAGS) $(OPTFLAGS) $< -o $@

.f90.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

.f.@f90objfilepostfix@:
	$(F90) -c $(ALL_F90FLAGS) $(F90OPTFLAGS) $< -o $@

#------------------------------------------------------------------------------#
OBJS = \
	hecmw2fstr_connect_conv.@cobjfilepostfix@ \
	fstr_ctrl_util.@cobjfilepostfix@ \
	fstr_sort_index.@cobjfilepostfix@

OBJSF = \
	hecmw2fstr_mesh_conv.@f90objfilepostfix@ \
	fstr_setup_util.@f90objfilepostfix@ \
	fstr_ctrl_modifier.@f90objfilepostfix@ \
	fstr_contact.@f90objfilepostfix@ \
	fstr_ctrl_common.@f90objfilepostfix@ \
	fstr_ctrl_material.@f90objfilepostfix@ \
	fstr_ctrl_static.@f90objfilepostfix@ \
	fstr_ctrl_heat.@f90objfilepostfix@ \
	fstr_ctrl_eigen.@f90objfilepostfix@ \
	fstr_ctrl_dynamic.@f90objfilepostfix@ \
	fstr_get_prop.@f90objfilepostfix@ \
	fstr_setup.@f90objfilepostfix@ \
	fstr_debug_dump.@f90objfilepostfix@ \
	fstr_precheck.@f90objfilepostfix@

OBJSF_RCAP = \
	fstr_rcap_io.@f90objfilepostfix@

OBJSF_RCAP_DUMM = \
	fstr_rcap_io_dummy.@f90objfilepostfix@

INCLUDESF = \
	fstr_ctrl_util_f.inc

HEADERS = \
	hecmw2fstr_connect_conv.h \
	fstr_ctrl_util.h \
	fstr_sort_index.h

#------------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-module

build-with-rcap: $(TARGET) F90TARGET_RCAP install-module

$(TARGET): $(OBJS)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF) $(OBJSF_RCAP_DUMM)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJSF) $(OBJSF_RCAP_DUMM)

F90TARGET_RCAP: $(OBJSF) $(OBJSF_RCAP)
	$(AR) $(FSTRLIBDIR)/$(F90TARGET) $(OBJSF) $(OBJSF_RCAP)

install-header:
	$(CP) $(HEADERS) $(FSTRINCDIR)/.
	$(CP) $(INCLUDESF) $(FSTRINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(FSTRINCDIR)/.

$(OBJS): $(HEADERS)

$(OBJSF): $(INCLUDESF)

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
