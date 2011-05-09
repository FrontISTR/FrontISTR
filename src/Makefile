
include ../Makefile.in
  
OBJS= \
AggregateElement.o AggregateNode.o AssyMatrix.o AssyModel.o \
AssyVector.o BeamBase.o Beam.o BoundaryFace.o BoundaryNode.o \
BoundaryVolume.o BoundingBox.o CommBeam.o CommElement.o \
CommFace.o CommHexa.o CommMesh2.o CommNode.o CommPrism.o \
CommPyramid.o CommQuad.o CommTetra.o CommTriangle.o CommunicationMesh.o \
ContactMatrix.o ContactMesh.o ContactNode.o Edge2ISTR.o EdgeFaceTree.o \
EdgeTree.o Element.o FaceTree.o FileIO.o FileReaderAssyModel.o \
FileReaderBoundaryFace.o FileReaderBoundaryNode.o FileReaderBoundaryVolume.o \
FileReaderChunk.o FileReaderCnt.o FileReaderCommElement.o \
FileReaderCommFace.o FileReaderCommMesh2.o FileReaderCommMesh.o \
FileReaderCommNode_CM2.o FileReaderCommNode.o FileReaderContactMesh.o \
FileReader.o FileReaderElement.o FileReaderMaterial.o FileReaderNode.o \
FileReaderRefine.o FileWriterBoundaryNode.o FileWriterChunk.o \
FileWriterCommMesh2.o FileWriterContactMesh.o FileWriter.o \
FileWriterElement.o FileWriterNode.o Film.o GeneralBoundary.o \
GMGModel.o HEC_MPI.o HEC_MW3.o Hexa.o IndexBucket.o IndexBucketMesh.o \
ISTR2Edge.o Jacobian.o Logger.o MasterFace.o Material.o MatrixBCRS.o \
Matrix.o Mesh.o MeshFactory.o MPCMatrix.o NodeConnectNodeTree.o Node.o \
OctreeKnot.o Prism.o ProgElementTree.o Pyramid.o Quad.o ScalarNode.o \
ScalarVectorNode.o ShapeFunctionBase.o ShapeFunctionCatalog.o \
ShapeHexa.o ShapeHexaNic.o ShapeLine.o ShapePrism.o ShapeQuad.o \
ShapeTetra.o ShapeTriangle.o ShellBase.o SkinFace.o SolidBase.o \
SolverBiCGSTAB.o SolverCG.o Solver.o SolverGMRES.o SolverGPBiCG.o \
SolverMG.o Tetra.o Triangle.o Vector.o VectorNode.o Vertex.o

default: $(OBJS)
	$(AR) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

clean:
	rm -rf *.o *~ core
