MPIDIR                  = @mpidir@
MPIBINDIR               = @mpibindir@
MPILIBDIR               = @mpilibdir@
MPIINCDIR               = @mpiincdir@
MPILIBS                 = @mpilibs@

HECMWDIR                = ../../../../../hecmw2
HECMWBINDIR             = $(HECMWDIR)/bin
HECMWLIBDIR             = $(HECMWDIR)/lib
HECMWINCDIR             = $(HECMWDIR)/include
HECMWLIBS               = @hecmwlibs@
HECMWDIR2               = ../../../../../hecmw1
HECMWBINDIR2            = $(HECMWDIR2)/bin
HECMWLIBDIR2            = $(HECMWDIR2)/lib
HECMWINCDIR2            = $(HECMWDIR2)/include
HECMWLIBS2              = @hecmwlibs2@

FSTRBINDIR              = ../../../../bin
FSTRLIBDIR              = ../../../../lib
FSTRINCDIR              = ../../../../include
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
HECMW_CFLAGS2           = @hecmw_cflags2@
FSTR_CFLAGS             = @fstr_cflags@
ALL_CFLAGS              = $(BASE_CFLAGS) $(FSTR_CFLAGS) $(HECMW_CFLAGS) $(HECMW_CFLAGS2) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                 = @ldflags@
MPI_LDFLAGS             = @mpi_ldflags@
HECMW_LDFLAGS           = @hecmw_ldflags@
HECMW_LDFLAGS2          = @hecmw_ldflags2@
FSTR_LDFLAGS            = @fstr_ldflags@
ALL_LDFLAGS             = $(FSTR_LDFLAGS) $(HECMW_LDFLAGS) $(HECMW_LDFLAGS2) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS                = @optflags@

REVOCAP_F90FLAGS        = @revocap_f90flags@
REVOCAP_F90LDFLAGS      = @revocap_f90ldflags@

F90                     = @f90@
F90FLAGS                = @f90flags@
BASE_F90FLAGS           = @base_f90flags@
MPI_F90FLAGS            = @mpi_f90flags@
HECMW_F90FLAGS          = @hecmw_f90flags@
HECMW_F90FLAGS2         = @hecmw_f90flags2@
FSTR_F90FLAGS           = @fstr_f90flags@
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(FSTR_F90FLAGS) $(HECMW_F90FLAGS) $(HECMW_F90FLAGS2) $(MPI_F90FLAGS) $(F90FLAGS) $(REVOCAP_F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
HECMW_F90LDFLAGS2       = @hecmw_f90ldflags2@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = $(FSTR_F90LDFLAGS) $(HECMW_F90LDFLAGS) $(HECMW_F90LDFLAGS2) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
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
OBJS =

OBJSF = \
	table_dyn.@f90objfilepostfix@ \
	dynamic_mat_ass_load.@f90objfilepostfix@ \
	dynamic_mat_ass_bc.@f90objfilepostfix@ \
	dynamic_mat_ass_bc_vl.@f90objfilepostfix@ \
	dynamic_mat_ass_bc_ac.@f90objfilepostfix@ \
	dynamic_mat_ass_couple.@f90objfilepostfix@ \
	dynamic_post.@f90objfilepostfix@ \
	dynamic_output.@f90objfilepostfix@ \
	fstr_dynamic_explicit.@f90objfilepostfix@ \
	fstr_dynamic_implicit.@f90objfilepostfix@ \
	fstr_solve_dynamic.@f90objfilepostfix@

HEADERS =

#------------------------------------------------------------------------------#
all: @build_target@

@build_target@: $(TARGET) $(F90TARGET) install-module

$(TARGET): $(OBJS)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(FSTRLIBDIR)/$@ $(OBJSF)

install-header:
	$(CP) $(HEADERS) $(FSTRINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(FSTRINCDIR)/.

$(OBJS): $(HEADERS)

$(OBJSF):

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
