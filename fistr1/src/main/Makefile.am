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

PARMETISDIR             = @parmetisdir@
PARMETISLIBDIR          = @parmetislibdir@
PARMETISINCDIR          = @parmetisincdir@
PARMETISLIBS            = @parmetislibs@
PARMETIS_CFLAGS         = @parmetis_cflags@
PARMETIS_LDFLAGS        = @parmetis_ldflags@
PARMETIS_F90FLAGS       = @parmetis_f90flags@
PARMETIS_F90LDFLAGS     = @parmetis_f90ldflags@

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

REFINERDIR              = @refinerdir@
REFINERINCDIR           = @refinerincdir@
REFINERLIBDIR           = @refinerlibdir@
REFINERLIBS             = @refinerlibs@

CC                      = @cc@
CFLAGS                  = @cflags@
BASE_CFLAGS             = @base_cflags@
MPI_CFLAGS              = @mpi_cflags@
HECMW_CFLAGS            = @hecmw_cflags@
FSTR_CFLAGS             = @fstr_cflags@
ALL_CFLAGS              = $(BASE_CFLAGS) $(FSTR_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(METIS_CFLAGS) $(PARMETIS_CFLAGS) $(MUMPS_CFLAGS) $(CFLAGS)
LDFLAGS                 = @ldflags@
MPI_LDFLAGS             = @mpi_ldflags@
HECMW_LDFLAGS           = @hecmw_ldflags@
FSTR_LDFLAGS            = @fstr_ldflags@
ALL_LDFLAGS             = $(FSTR_LDFLAGS) $(HECMW_LDFLAGS) $(MUMPS_LDFLAGS) $(METIS_LDFLAGS) $(PARMETIS_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS                = @optflags@

REVOCAP_F90FLAGS        = @revocap_f90flags@
REVOCAP_F90LDFLAGS      = @revocap_f90ldflags@

F90                     = @f90@
F90FLAGS                = @f90flags@
BASE_F90FLAGS           = @base_f90flags@
MPI_F90FLAGS            = @mpi_f90flags@
HECMW_F90FLAGS          = @hecmw_f90flags@
FSTR_F90FLAGS           = @fstr_f90flags@
ALL_F90FLAGS            = $(BASE_F90FLAGS) $(FSTR_F90FLAGS) $(HECMW_F90FLAGS) $(MPI_F90FLAGS) $(METIS_F90FLAGS) $(PARMETIS_F90FLAGS) $(MUMPS_F90FLAGS) $(F90FLAGS) $(REVOCAP_F90FLAGS)
F90LDFLAGS              = @f90ldflags@
MPI_F90LDFLAGS          = @mpi_f90ldflags@
HECMW_F90LDFLAGS        = @hecmw_f90ldflags@
FSTR_F90LDFLAGS         = @fstr_f90ldflags@
ALL_F90LDFLAGS          = $(FSTR_F90LDFLAGS) $(HECMW_F90LDFLAGS) $(MUMPS_F90LDFLAGS) $(METIS_F90LDFLAGS) $(PARMETIS_F90LDFLAGS) $(REVOCAP_F90LDFLAGS) $(MPI_F90LDFLAGS) $(F90LDFLAGS)
F90OPTFLAGS             = @f90optflags@
F90LINKER               = @f90linker@

AR                      = @ar@
CP                      = @cp@
RM                      = @rm@
RANLIB                  = @ranlib@
MKDIR                   = @mkdir@

#------------------------------------------------------------------------------#
TARGET                  = @fstrexec_targetfile@

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
	fistr_main.@f90objfilepostfix@

OBJSFP = \
	fistr_main_contact.@f90objfilepostfix@


HEADERS =

#------------------------------------------------------------------------------#
all: @build_target_paracon@

build-default: $(OBJSF) $(HECMWLIBDIR) $(FSTRLIBDIR)
	$(F90LINKER) $(ALL_F90FLAGS) $(F90OPTFLAGS) \
	-o $(TARGET) \
	$(OBJSF) $(ALL_F90LDFLAGS)
	$(CP) $(TARGET) $(FSTRBINDIR)/.

build-with-paracon: $(OBJSFP) $(HECMWLIBDIR) $(FSTRLIBDIR)
	$(F90LINKER) $(ALL_F90FLAGS) $(F90OPTFLAGS) \
	-o $(TARGET) \
	$(OBJSFP) $(ALL_F90LDFLAGS)
	$(CP) $(TARGET) $(FSTRBINDIR)/.

$(OBJS): $(HEADERS)

$(OBJSF):

$(OBJSFP):

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)
