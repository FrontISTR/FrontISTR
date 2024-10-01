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
F90FPP            = @f90fpp@
F90OPTFLAGS       = @f90optflags@

AR                = @ar@
CP                = @cp@
RM                = @rm@
RANLIB            = @ranlib@
MKDIR             = @mkdir@

#-----------------------------------------------------------------------------#
LEX               = @lex@

TARGET            = @hecmwlib_targetfile@
F90TARGET         = @hecmwlib_f90targetfile@

F90MODULEPOSTFIX  = @f90modulepostfix@
COBJFILEPOSTFIX   = @cobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@f90objfilepostfix@ .c .f90 .F90 .f

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

.f90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

.F90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90FPP) $(F90OPTFLAGS)

.f.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

#-----------------------------------------------------------------------------#
LEXSRC = \
	hecmw_ctrllex.c \
	hecmw_gflex.c \
	hecmw_heclex.c \
	hecmw_ablex.c

MSGSRC = \
	hecmw_msgno.h \
	hecmw_msgno_f.f90 \
	hecmw_msg_table.c

OBJS = \
	hecmw_hash.@cobjfilepostfix@ \
	hecmw_init.@cobjfilepostfix@ \
	hecmw_finalize.@cobjfilepostfix@ \
	hecmw_comm.@cobjfilepostfix@ \
	hecmw_conn_conv.@cobjfilepostfix@ \
	hecmw_error.@cobjfilepostfix@ \
	hecmw_util.@cobjfilepostfix@ \
	hecmw_malloc.@cobjfilepostfix@ \
	hecmw_log.@cobjfilepostfix@ \
	hecmw_msg.@cobjfilepostfix@ \
	hecmw_msg_table.@cobjfilepostfix@ \
	hecmw_lib_fc.@cobjfilepostfix@ \
	hecmw_control.@cobjfilepostfix@ \
	hecmw_ctrllex.@cobjfilepostfix@ \
	hecmw_bit_array.@cobjfilepostfix@ \
	hecmw_map_int.@cobjfilepostfix@ \
	hecmw_set_int.@cobjfilepostfix@ \
	hecmw_varray_int.@cobjfilepostfix@ \
	hecmw_etype.@cobjfilepostfix@ \
	hecmw_debug_write_dist.@cobjfilepostfix@ \
	hecmw_reorder.@cobjfilepostfix@ \
	hecmw_gflex.@cobjfilepostfix@ \
	hecmw_heclex.@cobjfilepostfix@ \
	hecmw_ablex.@cobjfilepostfix@ \
	hecmw_io_hec.@cobjfilepostfix@ \
	hecmw_io_abaqus.@cobjfilepostfix@ \
	hecmw_io_nastran_dummy.@cobjfilepostfix@ \
	hecmw_io_geofem.@cobjfilepostfix@ \
	hecmw_io_dist.@cobjfilepostfix@ \
	hecmw_io_mesh.@cobjfilepostfix@ \
	hecmw_io_get_mesh.@cobjfilepostfix@ \
	hecmw_io_get_mesh_if.@cobjfilepostfix@ \
	hecmw_io_put_mesh.@cobjfilepostfix@ \
	hecmw_dist.@cobjfilepostfix@ \
	hecmw_dist_alloc.@cobjfilepostfix@ \
	hecmw_dist_print.@cobjfilepostfix@ \
	hecmw_dist_free.@cobjfilepostfix@ \
	hecmw_dist_copy_f2c.@cobjfilepostfix@ \
	hecmw_dist_copy_c2f.@cobjfilepostfix@ \
	hecmw_dist_refine.@cobjfilepostfix@ \
	hecmw_ucd_print.@cobjfilepostfix@ \
	hecmw_geometric.@cobjfilepostfix@ \
	hecmw_path.@cobjfilepostfix@ \
	hecmw_put_mesh_if.@cobjfilepostfix@ \
	hecmw_restart.@cobjfilepostfix@ \
	hecmw_bin_io.@cobjfilepostfix@ \
	hecmw_result.@cobjfilepostfix@ \
	hecmw_result_copy_c2f.@cobjfilepostfix@ \
	hecmw_result_copy_f2c.@cobjfilepostfix@ \
	hecmw_result_io.@cobjfilepostfix@ \
	hecmw_result_io_bin.@cobjfilepostfix@ \
	hecmw_result_io_txt.@cobjfilepostfix@ \
	hecmw_visual_if.@cobjfilepostfix@ \
	hecmw_system.@cobjfilepostfix@ \
	hecmw_time.@cobjfilepostfix@

OBJSF = \
	hecmw_util_f.@f90objfilepostfix@ \
	hecmw_etype_f.@f90objfilepostfix@ \
	hecmw_logging.@f90objfilepostfix@ \
	hecmw_msgno_f.@f90objfilepostfix@ \
	hecmw_msg_f.@f90objfilepostfix@ \
	hecmw_control_f.@f90objfilepostfix@ \
	hecmw_result_f.@f90objfilepostfix@ \
	hecmw_restart_f.@f90objfilepostfix@ \
	hecmw_dist_copy_f2c_f.@f90objfilepostfix@ \
	hecmw_dist_copy_c2f_f.@f90objfilepostfix@ \
	hecmw_dist_free_f.@f90objfilepostfix@ \
	hecmw_dist_print_f.@f90objfilepostfix@ \
	hecmw_visualizer_f.@f90objfilepostfix@ \
	hecmw_varray_int_f.@f90objfilepostfix@ \
	hecmw_io.@f90objfilepostfix@

HEADERS = \
	hecmw_hash.h\
	hecmw_struct.h \
	hecmw_util.h \
	hecmw_comm.h \
	hecmw_conn_conv.h \
	hecmw_error.h \
	hecmw_malloc.h \
	hecmw_log.h \
	hecmw_msg.h \
	hecmw_msgno.h \
	hecmw_config.h \
	hecmw_init.h \
	hecmw_finalize.h \
	hecmw_lib_fc.h \
	hecmw_control.h \
	hecmw_ctrllex.h \
	hecmw_bit_array.h \
	hecmw_map_int.h \
	hecmw_set_int.h \
	hecmw_varray_int.h \
	hecmw_common.h \
	hecmw_common_define.h \
	hecmw_etype.h \
	hecmw_debug_write_dist.h \
	hecmw_reorder.h \
	hecmw_dist.h \
	hecmw_dist_alloc.h \
	hecmw_dist_print.h \
	hecmw_dist_copy_f2c.h \
	hecmw_dist_copy_c2f.h \
	hecmw_dist_free.h \
	hecmw_dist_refine.h \
	hecmw_geometric.h \
	hecmw_ablex.h \
	hecmw_gflex.h \
	hecmw_heclex.h \
	hecmw_io.h \
	hecmw_io_dist.h \
	hecmw_io_geofem.h \
	hecmw_io_nastran.h \
	hecmw_io_get_mesh.h \
	hecmw_io_hec.h \
	hecmw_io_abaqus.h \
	hecmw_io_mesh.h \
	hecmw_io_put_mesh.h \
	hecmw_io_struct.h \
	hecmw_path.h \
	hecmw_restart.h \
	hecmw_result.h \
	hecmw_ucd_print.h \
	hecmw_result_copy_c2f.h \
	hecmw_result_copy_f2c.h \
	hecmw_result_io.h \
	hecmw_result_io_bin.h \
	hecmw_result_io_txt.h \
	hecmw_system.h \
	hecmw_time.h

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-module

build-without-f: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSF):

install-header: $(HEADERS)
	$(CP) $(HEADERS) $(HECMWINCDIR)/.

install-module: $(OBJSF)
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean

include .dependf90

depend .dependf90:
	makedepf90 *.[fF]90 > .dependf90
