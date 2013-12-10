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

CPP               = @cpp@
CPPFLAGS          = @cppflags@
BASE_CPPFLAGS     = @base_cppflags@
HECMW_CPPFLAGS    = @hecmw_cppflags@
ALL_CPPFLAGS      = $(BASE_CPPFLAGS) $(HECMW_CPPFLAGS) $(MPI_CFLAGS) $(CPPFLAGS)
CPPLDFLAGS        = @cppldflags@
HECMW_CPPLDFLAGS  = @hecmw_cppldflags@
ALL_CPPLDFLAGS    = $(HECMW_CPPLDFLAGS) $(MPI_LDFLAGS) $(CPPLDFLAGS)
CPPOPTFLAGS       = @cppoptflags@

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
CPPOBJFILEPOSTFIX = @cppobjfilepostfix@
F90OBJFILEPOSTFIX = @f90objfilepostfix@

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
	AggregateElement.@cppobjfilepostfix@ \
	AggregateNode.@cppobjfilepostfix@ \
	API_Fortran.@cppobjfilepostfix@ \
	AssyMatrix.@cppobjfilepostfix@ \
	AssyModel.@cppobjfilepostfix@ \
	AssyVector.@cppobjfilepostfix@ \
	Beam.@cppobjfilepostfix@ \
	Beam2.@cppobjfilepostfix@ \
	BeamBase.@cppobjfilepostfix@ \
	BndVertex.@cppobjfilepostfix@ \
	BNodeMeshGrp.@cppobjfilepostfix@ \
	BoundaryEdge.@cppobjfilepostfix@ \
	BoundaryEdgeMesh.@cppobjfilepostfix@ \
	BoundaryFace.@cppobjfilepostfix@ \
	BoundaryFaceMesh.@cppobjfilepostfix@ \
	BoundaryHexa.@cppobjfilepostfix@ \
	BoundaryMesh.@cppobjfilepostfix@ \
	BoundaryNode.@cppobjfilepostfix@ \
	BoundaryNodeMesh.@cppobjfilepostfix@ \
	BoundaryParts.@cppobjfilepostfix@ \
	BoundaryPrism.@cppobjfilepostfix@ \
	BoundarySBNode.@cppobjfilepostfix@ \
	BoundaryTetra.@cppobjfilepostfix@ \
	BoundaryVolume.@cppobjfilepostfix@ \
	BoundaryVolumeMesh.@cppobjfilepostfix@ \
	BoundingBox.@cppobjfilepostfix@ \
	Calc.@cppobjfilepostfix@ \
	CommBeam.@cppobjfilepostfix@ \
	CommElement.@cppobjfilepostfix@ \
	CommFace.@cppobjfilepostfix@ \
	CommHexa.@cppobjfilepostfix@ \
	CommMesh2.@cppobjfilepostfix@ \
	CommNode.@cppobjfilepostfix@ \
	CommPrism.@cppobjfilepostfix@ \
	CommQuad.@cppobjfilepostfix@ \
	CommTetra.@cppobjfilepostfix@ \
	CommTriangle.@cppobjfilepostfix@ \
	CommunicationMesh.@cppobjfilepostfix@ \
	ContactMatrix.@cppobjfilepostfix@ \
	ContactMesh.@cppobjfilepostfix@ \
	ContactNode.@cppobjfilepostfix@ \
	DiscreteVolume.@cppobjfilepostfix@ \
	Edge2ISTR.@cppobjfilepostfix@ \
	EdgeFaceTree.@cppobjfilepostfix@ \
	EdgeTree.@cppobjfilepostfix@ \
	Element.@cppobjfilepostfix@ \
	ElementGroup.@cppobjfilepostfix@ \
	FaceProgFaceTree.@cppobjfilepostfix@ \
	FaceTree.@cppobjfilepostfix@ \
	FileIO.@cppobjfilepostfix@ \
	FileIOManage.@cppobjfilepostfix@ \
	FileReader.@cppobjfilepostfix@ \
	FileReaderAlgebra.@cppobjfilepostfix@ \
	FileReaderAssyModel.@cppobjfilepostfix@ \
	FileReaderBinCheck.@cppobjfilepostfix@ \
	FileReaderBoundaryEdge.@cppobjfilepostfix@ \
	FileReaderBoundaryEdgeMesh.@cppobjfilepostfix@ \
	FileReaderBoundaryFace.@cppobjfilepostfix@ \
	FileReaderBoundaryFaceMesh.@cppobjfilepostfix@ \
	FileReaderBoundaryNode.@cppobjfilepostfix@ \
	FileReaderBoundaryNodeMesh.@cppobjfilepostfix@ \
	FileReaderBoundaryVolume.@cppobjfilepostfix@ \
	FileReaderBoundaryVolumeMesh.@cppobjfilepostfix@ \
	FileReaderChunk.@cppobjfilepostfix@ \
	FileReaderCnt.@cppobjfilepostfix@ \
	FileReaderCommElement.@cppobjfilepostfix@ \
	FileReaderCommFace.@cppobjfilepostfix@ \
	FileReaderCommMesh.@cppobjfilepostfix@ \
	FileReaderCommMesh2.@cppobjfilepostfix@ \
	FileReaderCommNode.@cppobjfilepostfix@ \
	FileReaderCommNode_CM2.@cppobjfilepostfix@ \
	FileReaderContactMesh.@cppobjfilepostfix@ \
	FileReaderElement.@cppobjfilepostfix@ \
	FileReaderElementGroup.@cppobjfilepostfix@ \
	FileReaderElementGroupEntity.@cppobjfilepostfix@ \
	FileReaderMaterial.@cppobjfilepostfix@ \
	FileReaderNode.@cppobjfilepostfix@ \
	FileReaderRefine.@cppobjfilepostfix@ \
	FileReaderRes.@cppobjfilepostfix@ \
	FileReaderVersion.@cppobjfilepostfix@ \
	FileWriter.@cppobjfilepostfix@ \
	FileWriterAVS.@cppobjfilepostfix@ \
	FileWriterBoundaryEdge.@cppobjfilepostfix@ \
	FileWriterBoundaryFace.@cppobjfilepostfix@ \
	FileWriterBoundaryNode.@cppobjfilepostfix@ \
	FileWriterBoundaryVolume.@cppobjfilepostfix@ \
	FileWriterChunk.@cppobjfilepostfix@ \
	FileWriterCommMesh2.@cppobjfilepostfix@ \
	FileWriterContactMesh.@cppobjfilepostfix@ \
	FileWriterElement.@cppobjfilepostfix@ \
	FileWriterElementGroup.@cppobjfilepostfix@ \
	FileWriterNode.@cppobjfilepostfix@ \
	FileWriterRes.@cppobjfilepostfix@ \
	FileWriterUNS.@cppobjfilepostfix@ \
	FileWriterVTK.@cppobjfilepostfix@ \
	Film.@cppobjfilepostfix@ \
	GMGModel.@cppobjfilepostfix@ \
	HEC_MPI.@cppobjfilepostfix@ \
	HEC_MW3.@cppobjfilepostfix@ \
	Hexa.@cppobjfilepostfix@ \
	Hexa2.@cppobjfilepostfix@ \
	IndexBucket.@cppobjfilepostfix@ \
	IndexBucketMesh.@cppobjfilepostfix@ \
	ISTR2Edge.@cppobjfilepostfix@ \
	Jacobian.@cppobjfilepostfix@ \
	Logger.@cppobjfilepostfix@ \
	MasterFace.@cppobjfilepostfix@ \
	Material.@cppobjfilepostfix@ \
	Matrix.@cppobjfilepostfix@ \
	MatrixBCRS.@cppobjfilepostfix@ \
	Mesh.@cppobjfilepostfix@ \
	MeshFactory.@cppobjfilepostfix@ \
	MPCMatrix.@cppobjfilepostfix@ \
	Node.@cppobjfilepostfix@ \
	NodeConnectNodeTree.@cppobjfilepostfix@ \
	OctreeKnot.@cppobjfilepostfix@ \
	Poland.@cppobjfilepostfix@ \
	Prism.@cppobjfilepostfix@ \
	Prism2.@cppobjfilepostfix@ \
	ProgElementTree.@cppobjfilepostfix@ \
	Quad.@cppobjfilepostfix@ \
	Quad2.@cppobjfilepostfix@ \
	ScalarNode.@cppobjfilepostfix@ \
	ScalarVectorNode.@cppobjfilepostfix@ \
	ShapeFunctionBase.@cppobjfilepostfix@ \
	ShapeFunctionCatalog.@cppobjfilepostfix@ \
	ShapeHexa.@cppobjfilepostfix@ \
	ShapeHexaNic.@cppobjfilepostfix@ \
	ShapeLine.@cppobjfilepostfix@ \
	ShapePrism.@cppobjfilepostfix@ \
	ShapeQuad.@cppobjfilepostfix@ \
	ShapeTetra.@cppobjfilepostfix@ \
	ShapeTriangle.@cppobjfilepostfix@ \
	ShellBase.@cppobjfilepostfix@ \
	SkinFace.@cppobjfilepostfix@ \
	SolidBase.@cppobjfilepostfix@ \
	Solver.@cppobjfilepostfix@ \
	SolverBiCGSTAB.@cppobjfilepostfix@ \
	SolverCG.@cppobjfilepostfix@ \
	SolverGMRES.@cppobjfilepostfix@ \
	SolverGPBiCG.@cppobjfilepostfix@ \
	SolverMG.@cppobjfilepostfix@ \
	Tetra.@cppobjfilepostfix@ \
	Tetra2.@cppobjfilepostfix@ \
	Triangle.@cppobjfilepostfix@ \
	Triangle2.@cppobjfilepostfix@ \
	Vector.@cppobjfilepostfix@ \
	VectorNode.@cppobjfilepostfix@ \
	Vertex.@cppobjfilepostfix@

OBJSF =

HEADERS = \
	AggregateElement.h \
	AggregateNode.h \
	API_Fortran.h \
	API_Fortran.hxx \
	AssyMatrix.h \
	AssyModel.h \
	AssyVector.h \
	Beam.h \
	Beam2.h \
	BeamBase.h \
	BndVertex.h \
	BNodeMeshGrp.h \
	BoundaryEdge.h \
	BoundaryEdgeMesh.h \
	BoundaryFace.h \
	BoundaryFaceMesh.h \
	BoundaryGroup.h \
	BoundaryHexa.h \
	BoundaryMesh.h \
	BoundaryNode.h \
	BoundaryNodeMesh.h \
	BoundaryParts.h \
	BoundaryPrism.h \
	BoundarySBNode.h \
	BoundaryTetra.h \
	BoundaryType.h \
	BoundaryVolume.h \
	BoundaryVolumeMesh.h \
	BoundingBox.h \
	Calc.h \
	CodeType.h \
	CommBeam.h \
	CommElement.h \
	CommFace.h \
	CommHexa.h \
	CommMesh2.h \
	CommNode.h \
	CommonFile.h \
	CommonStd.h \
	CommPrism.h \
	CommQuad.h \
	CommTetra.h \
	CommTriangle.h \
	CommunicationMesh.h \
	ContactMatrix.h \
	ContactMesh.h \
	ContactNode.h \
	DeleteObject.h \
	DiscreteVolume.h \
	Edge2ISTR.h \
	EdgeFaceTree.h \
	EdgeTree.h \
	Element.h \
	Element4Factory.h \
	ElementGroup.h \
	ElementProperty.h \
	ElementType.h \
	EqnTerm.h \
	Equation.h \
	FaceProgFaceTree.h \
	FaceTree.h \
	FileBlockName.h \
	FileBlockNameMW2.h \
	FileIO.h \
	FileIOManage.h \
	FileReader.h \
	FileReaderAlgebra.h \
	FileReaderAssyModel.h \
	FileReaderBinCheck.h \
	FileReaderBoundaryEdge.h \
	FileReaderBoundaryEdgeMesh.h \
	FileReaderBoundaryFace.h \
	FileReaderBoundaryFaceMesh.h \
	FileReaderBoundaryNode.h \
	FileReaderBoundaryNodeMesh.h \
	FileReaderBoundaryVolume.h \
	FileReaderBoundaryVolumeMesh.h \
	FileReaderChunk.h \
	FileReaderCnt.h \
	FileReaderCommElement.h \
	FileReaderCommFace.h \
	FileReaderCommMesh.h \
	FileReaderCommMesh2.h \
	FileReaderCommNode.h \
	FileReaderCommNode_CM2.h \
	FileReaderContactMesh.h \
	FileReaderElement.h \
	FileReaderElementGroup.h \
	FileReaderElementGroupEntity.h \
	FileReaderMaterial.h \
	FileReaderNode.h \
	FileReaderRefine.h \
	FileReaderRes.h \
	FileReaderVersion.h \
	FileWriter.h \
	FileWriterAVS.h \
	FileWriterBoundaryEdge.h \
	FileWriterBoundaryFace.h \
	FileWriterBoundaryNode.h \
	FileWriterBoundaryVolume.h \
	FileWriterChunk.h \
	FileWriterCommMesh2.h \
	FileWriterContactMesh.h \
	FileWriterElement.h \
	FileWriterElementGroup.h \
	FileWriterNode.h \
	FileWriterRes.h \
	FileWriterUNS.h \
	FileWriterVTK.h \
	Film.h \
	GMGModel.h \
	HEC_MPI.h \
	HEC_MW3.h \
	HEC_MW3.hxx \
	HEC_MW3_For.h \
	Hexa.h \
	Hexa2.h \
	IndexBucket.h \
	IndexBucketMesh.h \
	IndexValuePair.h \
	Info.h \
	ISTR2Edge.h \
	Jacobian.h \
	Logger.h \
	LoggerType.h \
	MasterFace.h \
	Material.h \
	MatrialPropType.h \
	Matrix.h \
	MatrixBCRS.h \
	Mesh.h \
	MeshFactory.h \
	MPCMatrix.h \
	MPCValueType.h \
	Node.h \
	NodeConnectNodeTree.h \
	NodeType.h \
	OctreeKnot.h \
	Poland.h \
	Prism.h \
	Prism2.h \
	ProgElementTree.h \
	Quad.h \
	Quad2.h \
	QuickSort.h \
	ScalarNode.h \
	ScalarVectorNode.h \
	ShapeFunctionBase.h \
	ShapeFunctionCatalog.h \
	ShapeHexa.h \
	ShapeHexaNic.h \
	ShapeLine.h \
	ShapePrism.h \
	ShapeQuad.h \
	ShapeTetra.h \
	ShapeTriangle.h \
	ShellBase.h \
	SkinFace.h \
	SolidBase.h \
	SolutionType.h \
	Solver.h \
	SolverBiCGSTAB.h \
	SolverCG.h \
	SolverGMRES.h \
	SolverGPBiCG.h \
	SolverMG.h \
	Tetra.h \
	Tetra2.h \
	Triangle.h \
	Triangle2.h \
	TypeDef.h \
	Vector.h \
	VectorNode.h \
	Vertex.h

#-----------------------------------------------------------------------------#
all: @build_target@

build-default: $(TARGET) $(F90TARGET) install-header

build-serial: $(TARGET) $(F90TARGET) install-header

build-without-f: $(TARGET) install-header

$(TARGET): $(OBJS) $(OBJSP)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJS)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSP)

$(F90TARGET): $(OBJSF)
	$(AR) $(HECMWLIBDIR)/$@ $(OBJSF)

$(OBJS): $(HEADERS)

$(OBJSP): $(HEADERS)

$(OBJSF):

install-header: $(HEADERS)
	$(CP) $(HEADERS) $(HECMWINCDIR)/.

install-module:
	$(CP) *.$(F90MODULEPOSTFIX) $(HECMWINCDIR)/.

install:

clean:
	$(RM) *.$(COBJFILEPOSTFIX) *.$(CPPOBJFILEPOSTFIX) *.$(F90OBJFILEPOSTFIX) *.$(F90MODULEPOSTFIX)

distclean: clean
