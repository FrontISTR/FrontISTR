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

METISDIR               = @metisdir@
METISLIBDIR            = @metislibdir@
METISINCDIR            = @metisincdir@
METISLIBS              = @metislibs@
METIS_CFLAGS           = @metis_cflags@
METIS_LDFLAGS          = @metis_ldflags@

REFINERDIR             = @refinerdir@
REFINERINCDIR          = @refinerincdir@
REFINERLIBDIR          = @refinerlibdir@
REFINERLIBS            = @refinerlibs@

CC                     = @cc@
CFLAGS                 = @cflags@
BASE_CFLAGS            = @base_cflags@
MPI_CFLAGS             = @mpi_cflags@
HECMW_CFLAGS           = @hecmw_cflags@
METIS_CFLAGS           = @metis_cflags@
PARTITIONER_CFLAGS     = @partitioner_cflags@
ALL_CFLAGS             = $(PARTITIONER_CFLAGS) $(METIS_CFLAGS) $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                = @ldflags@
MPI_LDFLAGS            = @mpi_ldflags@
HECMW_LDFLAGS          = @hecmw_ldflags@
METIS_LDFLAGS          = @metis_ldflags@
PARTITIONER_LDFLAGS    = @partitioner_ldflags@
ALL_LDFLAGS            = $(PARTITIONER_LDFLAGS) $(METIS_LDFLAGS) $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS               = @optflags@
PARTITIONER_OPTFLAGS   = @partitioner_optflags@
ALL_OPTFLAGS           = $(OPTFLAGS) $(PARTITIONER_OPTFLAGS)
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
TARGET                 = @partitioner_targetfile@

COBJFILEPOSTFIX        = @cobjfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .c

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(ALL_OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS = \
	hecmw_partlex.@cobjfilepostfix@ \
	hecmw_part_log.@cobjfilepostfix@ \
	hecmw_mesh_hash_sort.@cobjfilepostfix@ \
	hecmw_mesh_edge_info.@cobjfilepostfix@ \
	hecmw_part_get_control.@cobjfilepostfix@ \
	hecmw_init_for_partition.@cobjfilepostfix@ \
	hecmw_graph.@cobjfilepostfix@ \
	hecmw_partition.@cobjfilepostfix@ \
	hecmw_partitioner.@cobjfilepostfix@

HEADERS = \
	hecmw_part_define.h \
	hecmw_partlex.h \
	hecmw_part_log.h \
	hecmw_mesh_hash_sort.h \
	hecmw_mesh_edge_info.h \
	hecmw_part_get_control.h \
	hecmw_init_for_partition.h \
	hecmw_graph.h \
	hecmw_partition.h

LEXSRC = \
	hecmw_partlex.c

#-----------------------------------------------------------------------------#
all: @build_target@

@build_target@: $(TARGET)

$(TARGET): $(OBJS) $(HECMWLIBDIR)/*.*
	$(CLINKER) $(ALL_CFLAGS) $(ALL_OPTFLAGS) \
	-o $@ $(OBJS) $(ALL_LDFLAGS)
	$(CP) $(TARGET) $(HECMWBINDIR)/.

$(OBJS): $(HEADERS)

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)
