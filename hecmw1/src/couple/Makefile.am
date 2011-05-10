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
	hecmw_couple_comm.@cobjfilepostfix@ \
	hecmw_couple_copy_c2f.@cobjfilepostfix@ \
	hecmw_couple_copy_f2c.@cobjfilepostfix@ \
	hecmw_couple_control.@cobjfilepostfix@ \
	hecmw_couple_info.@cobjfilepostfix@ \
	hecmw_couple_info_if.@cobjfilepostfix@ \
	hecmw_couple_get_mesh.@cobjfilepostfix@ \
	hecmw_couple_boundary_info.@cobjfilepostfix@ \
	hecmw_couple_bounding_box.@cobjfilepostfix@ \
	hecmw_couple_background_cell.@cobjfilepostfix@ \
	hecmw_couple_mapped_point.@cobjfilepostfix@ \
	hecmw_couple_judge.@cobjfilepostfix@ \
	hecmw_couple_inter_iftable.@cobjfilepostfix@ \
	hecmw_couple_intra_iftable.@cobjfilepostfix@ \
	hecmw_couple_weight.@cobjfilepostfix@ \
	hecmw_couple_s2n_average.@cobjfilepostfix@ \
	hecmw_couple_s2n_with_area.@cobjfilepostfix@ \
	hecmw_couple_s2n_dist_node.@cobjfilepostfix@ \
	hecmw_couple_s2n_dist_surf.@cobjfilepostfix@ \
	hecmw_couple_n2s_average.@cobjfilepostfix@ \
	hecmw_couple_n2s_with_area.@cobjfilepostfix@ \
	hecmw_couple_interpolate_info.@cobjfilepostfix@ \
	hecmw_couple_init.@cobjfilepostfix@ \
	hecmw_couple_startup.@cobjfilepostfix@ \
	hecmw_couple_finalize.@cobjfilepostfix@ \
	hecmw_couple.@cobjfilepostfix@ \
	hecmw_couple_init_if.@cobjfilepostfix@ \
	hecmw_couple_startup_if.@cobjfilepostfix@ \
	hecmw_couple_if.@cobjfilepostfix@ \
	hecmw_couple_finalize_if.@cobjfilepostfix@

OBJSF = \
    hecmw_couple_define_f.@f90objfilepostfix@ \
    hecmw_couple_struct_f.@f90objfilepostfix@ \
	hecmw_couple_copy_c2f_f.@f90objfilepostfix@ \
	hecmw_couple_copy_f2c_f.@f90objfilepostfix@ \
	hecmw_couple_info_f.@f90objfilepostfix@ \
	hecmw_couple_get_mesh_f.@f90objfilepostfix@ \
	hecmw_couple_init_f.@f90objfilepostfix@ \
	hecmw_couple_startup_f.@f90objfilepostfix@ \
	hecmw_couple_cleanup_f.@f90objfilepostfix@ \
	hecmw_couple_f.@f90objfilepostfix@ \
	hecmw_couple_finalize_f.@f90objfilepostfix@

HEADERS = \
	hecmw_couple_define.h \
	hecmw_couple_struct.h \
	hecmw_couple_comm.h \
	hecmw_couple_control.h \
	hecmw_couple_info.h \
	hecmw_couple_boundary_info.h \
	hecmw_couple_bounding_box.h \
	hecmw_couple_background_cell.h \
	hecmw_couple_mapped_point.h \
	hecmw_couple_judge.h \
	hecmw_couple_inter_iftable.h \
	hecmw_couple_intra_iftable.h \
	hecmw_couple_weight.h \
	hecmw_couple_s2n_average.h \
	hecmw_couple_s2n_with_area.h \
	hecmw_couple_s2n_dist_node.h \
	hecmw_couple_s2n_dist_surf.h \
	hecmw_couple_n2s_average.h \
	hecmw_couple_n2s_with_area.h \
	hecmw_couple_interpolate_info.h \
	hecmw_couple_init.h \
	hecmw_couple_startup.h \
	hecmw_couple.h \
	hecmw_couple_finalize.h \
	hecmw_couple_copy_c2f.h \
	hecmw_couple_copy_f2c.h

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-module

build-serial: $(TARGET) $(F90TARGET) install-module

build-without-f: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

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
