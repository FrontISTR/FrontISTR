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
HECMWBINDIR2           = ../../../hecmw1/bin
HECMWLIBDIR2           = ../../../hecmw1/lib
HECMWINCDIR2           = ../../../hecmw1/include
HECMWLIBS2             = -lfhecmw -lhecmw

REFINERDIR             = @refinerdir@
REFINERINCDIR          = @refinerincdir@
REFINERLIBDIR          = @refinerlibdir@
REFINERLIBS            = @refinerlibs@

CC                     = @cc@
CFLAGS                 = @cflags@
BASE_CFLAGS            = @base_cflags@
MPI_CFLAGS             = @mpi_cflags@
HECMW_CFLAGS           = @hecmw_cflags@
ALL_CFLAGS             = $(BASE_CFLAGS) $(HECMW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)
LDFLAGS                = @ldflags@
MPI_LDFLAGS            = @mpi_ldflags@
HECMW_LDFLAGS          = @hecmw_ldflags@
ALL_LDFLAGS            = $(HECMW_LDFLAGS) $(MPI_LDFLAGS) $(LDFLAGS)
OPTFLAGS               = @optflags@

CPP                    = @cpp@
CPPFLAGS               = @cppflags@
BASE_CPPFLAGS          = @base_cppflags@
HECMW_CPPFLAGS         = @hecmw_cppflags@
ALL_CPPFLAGS           = $(BASE_CPPFLAGS) $(HECMW_CPPFLAGS) $(MPI_CFLAGS) $(CPPFLAGS)
CPPLDFLAGS             = @cppldflags@
HECMW_CPPLDFLAGS       = @hecmw_cppldflags@
ALL_CPPLDFLAGS         = $(HECMW_CPPLDFLAGS) $(MPI_LDFLAGS) $(CPPLDFLAGS)
CPPOPTFLAGS            = @cppoptflags@

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
TARGET                 = @converter_targetfile@

F90MODULEPOSTFIX       = @f90modulepostfix@
COBJFILEPOSTFIX        = @cobjfilepostfix@
CPPOBJFILEPOSTFIX      = @cppobjfilepostfix@
F90OBJFILEPOSTFIX      = @f90objfilepostfix@

#-----------------------------------------------------------------------------#
.SUFFIXES:
.SUFFIXES: .@cobjfilepostfix@ .@cppobjfilepostfix@ .@f90objfilepostfix@ .c .cpp .f90

.c.@cobjfilepostfix@:
	$(CC) -c $< $(ALL_CFLAGS) $(OPTFLAGS)

.cpp.@cppobjfilepostfix@:
	$(CPP) -c $< $(ALL_CPPFLAGS) $(CPPOPTFLAGS)

.f90.@f90objfilepostfix@:
	$(F90) -c $< $(ALL_F90FLAGS) $(F90OPTFLAGS)

#-----------------------------------------------------------------------------#
OBJS =

OBJSP = \
	AssyModel.@cppobjfilepostfix@ \
	Beam.@cppobjfilepostfix@ \
	Beam2.@cppobjfilepostfix@ \
	ConPair.@cppobjfilepostfix@ \
	ConvMain.@cppobjfilepostfix@ \
	Egroup.@cppobjfilepostfix@ \
	Element.@cppobjfilepostfix@ \
	FileBase.@cppobjfilepostfix@ \
	FileReader.@cppobjfilepostfix@ \
	FileWriter.@cppobjfilepostfix@ \
	Group.@cppobjfilepostfix@ \
	GroupPair.@cppobjfilepostfix@ \
	Hexa.@cppobjfilepostfix@ \
	Hexa2.@cppobjfilepostfix@ \
	IFaceQuad.@cppobjfilepostfix@ \
	IFaceQuad2.@cppobjfilepostfix@ \
	Lgroup.@cppobjfilepostfix@ \
	Mesh.@cppobjfilepostfix@ \
	Message.@cppobjfilepostfix@ \
	MPCPair.@cppobjfilepostfix@ \
	Ngroup.@cppobjfilepostfix@ \
	Node.@cppobjfilepostfix@ \
	Prism.@cppobjfilepostfix@ \
	Prism2.@cppobjfilepostfix@ \
	Quad.@cppobjfilepostfix@ \
	Quad2.@cppobjfilepostfix@ \
	QuadShell.@cppobjfilepostfix@ \
	QuadShell2.@cppobjfilepostfix@ \
	Sgroup.@cppobjfilepostfix@ \
	Tetra.@cppobjfilepostfix@ \
	Tetra2.@cppobjfilepostfix@ \
	Triangle.@cppobjfilepostfix@ \
	Triangle2.@cppobjfilepostfix@ \
	TriShell.@cppobjfilepostfix@ \
	TriShell2.@cppobjfilepostfix@ \
	main.@cppobjfilepostfix@

OBJSF =

HEADERS = \
	AssyModel.h \
	Beam.h \
	Beam2.h \
	ConPair.h \
	ConvHeader.h \
	ConvMain.h \
	Egroup.h \
	Element.h \
	FileBase.h \
	FileReader.h \
	FileWriter.h \
	FrontISTR_Type.h \
	Group.h \
	GroupPair.h \
	Hexa.h \
	Hexa2.h \
	IFaceQuad.h \
	IFaceQuad2.h \
	Lgroup.h \
	Mesh.h \
	Message.h \
	MPCPair.h \
	Ngroup.h \
	Node.h \
	Prism.h \
	Prism2.h \
	Quad.h \
	Quad2.h \
	QuadShell.h \
	QuadShell2.h \
	Sgroup.h \
	Tetra.h \
	Tetra2.h \
	Triangle.h \
	Triangle2.h \
	TriShell.h \
	TriShell2.h

#-----------------------------------------------------------------------------#
all: @build_target@

@build_target@: $(TARGET)

$(TARGET): $(OBJSP) $(HECMWLIBDIR)/*.* $(HECMWLIBDIR2)/*.*
	$(CPP) $(ALL_CPPFLAGS) $(CPPOPTFLAGS) \
	-o $@ $(OBJSP) $(ALL_CPPLDFLAGS)
	$(CP) $(TARGET) $(HECMWBINDIR)/.

$(OBJS): $(HEADERS)

$(OBJSP): $(HEADERS)

$(OBJSF):

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(CPPOBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
	$(RM) $(TARGET)
