/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/HEC_MW3.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "HEC_MPI.h"
#include "FileIO.h"
#include "Logger.h"
#include "SolutionType.h"
#include "MeshFactory.h"
#include "GMGModel.h"
#include "AssyMatrix.h"
#include "ShapeHexa.h"
#include "ShapeHexaNic.h"
#include "ShapeTetra.h"
#include "ShapePrism.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeLine.h"
#include "ShapeFunctionCatalog.h"
#include "ISTR2Edge.h"
#include "Edge2ISTR.h"
#include "Jacobian.h"
#include <cstdarg>
#include <iomanip>
#include <cstring>
typedef pair<uiint,uiint> integPair;
namespace pmw{
#ifndef PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#define PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
class CMW
{
public:
    static CMW*  Instance(){
        static CMW moMW;
        return &moMW;
    }
private:
     CMW(void);
public:
     virtual ~CMW(void);
protected:
    FileIO::CFileIO   *mpFileIO;
    Utility::CLogger  *mpLogger;
    uiint mnSolutionType;
    bool mb_file;
    string msInputFileName;  
    string msOutputFileName; 
    string msResFileName;    
    string msRltFileName;    
    string msPartOutFileName;
    CGMGModel     *mpGMGModel;
    CMeshFactory  *mpFactory;
    CHecMPI       *mpMPI;
    CAssyModel  *mpAssy;
    CAssyMatrix *mpAssyMatrix;
    CAssyVector *mpRHSAssyVector;
    CAssyVector *mpSolAssyVector;
    CMesh       *mpMesh;
    CElement    *mpElement;
    CCommMesh2  *mpComMesh;
    CShapeHexa     *mpShapeHexa;
    CShapeHexaNic  *mpShapeHexaNic;
    CShapeTetra    *mpShapeTetra;
    CShapePrism    *mpShapePrism;
    CShapeQuad     *mpShapeQuad;
    CShapeTriangle *mpShapeTriangle;
    CShapeLine     *mpShapeLine;
    CShapeFunctionCatalog *mpShapeCatalog;
    CISTR2Edge *mpISTR2Edge;
    CEdge2ISTR *mpEdge2ISTR;
    vvvdouble mvdNdx;
    vuint mvAlgebraEquation;
    vuint mv_ItemL;
    vuint mv_ItemU;
    uiint BaseName_BCast(int& nLength, string& sName, int nType);
public:
    uiint Initialize( int argc, char** argv);
    uiint Initialize_fstr( int argc, char** argv, string& ctrlname);
    uiint Finalize();  
    void  Banner_Display();
    uiint FileRead(string& basename, bool bBinary);
    uiint FileRead_fstr(bool bBinary);             
    uiint FileDebugWrite();          
    uiint FileWriteRes(const uiint& nStep, bool bBinary);
    uiint SetRestart(const uiint& nStep, bool bBinary);
    void PrintRlt_Start(const uiint& nStep, bool bBinary);
    void PrintRlt_P(const uiint& width, const char* format, ... );
    void PrintRlt_R(const uiint& width, const char* format, ... );
    void PrintRlt_End();
    void PrintMicroAVS_Basis();
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void PrintMicroAVS_FEM();  
    string getFstr_FileName_Mesh();
    string getFstr_FileName_Control();
    string getFstr_FileName_Result();
    string getFstr_FileName_Restart();
    string getFstr_FileName_PartIn();
    string getFstr_FileName_PartOut();
    string getFstr_FileName_VisMesh();
    string getFstr_FileName_VisIn();
    string getFstr_FileName_VisOut();
    void GeneLinearAlgebra(const uiint& nNumOfAlgebra, uiint* vDOF);
    void SelectAlgebra(const uiint& iequ);
    uiint Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    uiint Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double *NodalMatrix);
    void Matrix_Clear(const uiint& iMesh);
    void Vector_Clear(const uiint& iMesh);
    uiint Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& sol_value);
    uiint Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& rhs_value); 
    uiint Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                        
    uiint Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                        
    void  Sample_Set_BC(uiint iMesh);
    uiint Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition);
    void GetSolution_Vector(double* buf, const uiint& imesh);
    void GetSolution_AssyVector(double* buf);
    void GetRHS_Vector(double* buf, const uiint& imesh);
    void GetRHS_AssyVector(double* buf);
    double& GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    double& GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    uiint& GetSolutionVector_DOF();
    uiint& GetRHSVector_DOF();
    void multVector(double* vX, double* vB);
    uiint  Refine(const uiint& nNumOfRefine);
    uiint GetNumOfAssembleModel();
    void SelectAssembleModel(const uiint& mgLevel);
    uiint GetNumOfMeshPart();
    void SelectMeshPart_ID(const uiint& mesh_id);
    void SelectMeshPart_IX(const uiint& index);
    void SelectElement_ID(const uiint& elem_id);
    void SelectElement_IX(const uiint& index);
    uiint GetElementType();
    uiint GetNumOfElementVert();
    void GetElementVertNodeID(iint* vNodeID);
    uiint GetNumOfElementEdge();
    void GetElementEdgeNodeID(iint* vNodeID);
    void GetNodeCoord(const uiint& node_id, double& x, double& y, double& z);
    uiint GetNumOfDOF(const uiint& node_id);
    uiint GetNumOfScalar(const uiint& node_id);
    uiint GetNumOfVector(const uiint& node_id);
    uiint& GetNodeType(const uiint& node_id);
    uiint  getNodeSize(){ return mpMesh->getNodeSize();}
    uiint  getElementSize(){ return mpMesh->getElementSize();}
    uiint  getNodeSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    uiint  getElementSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
    uiint& getNodeID(const uiint& index);
    uiint& getElementID(const uiint& index);
    uiint& getNodeIndex(const uiint& id);   
    uiint& getElementIndex(const uiint& id);
    void constructNodeConnectFEM(const uiint& node_id);
    void getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL);
    void getNodeConnectFEM_Item(uiint itemU[], uiint itemL[]);
    void getNodeConnectFEM_Item_F(iint itemU[], iint itemL[]);
    uiint getNumOfAggregateElement(const uiint& node_id);
    uiint& getAggregateElementID(const uiint& node_id, const uiint& ielem);
    void FinalizeRefine();
    uiint nodetype_s();
    uiint nodetype_v();
    uiint nodetype_sv();
    uiint elemtype_hexa();
    uiint elemtype_hexa2();
    uiint elemtype_tetra();
    uiint elemtype_tetra2();
    uiint elemtype_prism();
    uiint elemtype_prism2();
    uiint elemtype_quad();
    uiint elemtype_quad2();
    uiint elemtype_triangle();
    uiint elemtype_triangle2();
    uiint elemtype_line();
    uiint elemtype_line2();
    uiint fistr_elemtype_hexa();
    uiint fistr_elemtype_hexa2();
    uiint fistr_elemtype_tetra();
    uiint fistr_elemtype_tetra2();
    uiint fistr_elemtype_prism();
    uiint fistr_elemtype_prism2();
    uiint fistr_elemtype_quad();
    uiint fistr_elemtype_quad2();
    uiint fistr_elemtype_triangle();
    uiint fistr_elemtype_triangle2();
    uiint fistr_elemtype_line();
    uiint fistr_elemtype_line2();
    uiint fistr_elemtype_to_mw3_elemtype(const uiint& fistr_elemtype);
    uiint mw3_elemtype_to_fistr_elemtype(const uiint& mw3_elemtype);
    uiint& NumOfIntegPoint(const uiint& shapeType);
    void ShapeFunc_on_pt(const uiint& shapeType, const uiint& igauss, vdouble& N);
    void ShapeFunc_on_pt(uiint shapeType, uiint igauss, double N[]);
    void ShapeFunc(const uiint& shapeType, vvdouble& N);
    double& ShapeFunc_Hexa81(uiint igauss, uiint ishape);
    double& ShapeFunc_Hexa82(uiint igauss, uiint ishape);
    double& ShapeFunc_Hexa201(uiint igauss, uiint ishape);
    double& ShapeFunc_Hexa202(uiint igauss, uiint ishape);
    double& ShapeFunc_Hexa203(uiint igauss, uiint ishape);
    double& ShapeFunc_Tetra41(uiint igauss, uiint ishape);
    double& ShapeFunc_Tetra101(uiint igauss, uiint ishape);
    double& ShapeFunc_Tetra104(uiint igauss, uiint ishape);
    double& ShapeFunc_Tetra1015(uiint igauss, uiint ishape);
    double& ShapeFunc_Prism62(uiint igauss, uiint ishape);
    double& ShapeFunc_Prism156(uiint igauss, uiint ishape);
    double& ShapeFunc_Prism159(uiint igauss, uiint ishape);
    double& ShapeFunc_Prism1518(uiint igauss, uiint ishape);
    double& ShapeFunc_Quad41(uiint igauss, uiint ishape);
    double& ShapeFunc_Quad84(uiint igauss, uiint ishape);
    double& ShapeFunc_Quad89(uiint igauss, uiint ishape);
    double& ShapeFunc_Triangle31(uiint igauss, uiint ishape);
    double& ShapeFunc_Triangle63(uiint igauss, uiint ishape);
    double& ShapeFunc_Line21(uiint igauss, uiint ishape);
    double& ShapeFunc_Line32(uiint igauss, uiint ishape);
    void dNdr_on_pt(const uiint& shapeType, const uiint& igauss, vvdouble& dNdr);
    void dNdr(const uiint& shapeType, vvvdouble& dNdr);
    void dNdr(const uiint& shapeType, double dNdr[]);
    double& dNdr_Hexa81_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Hexa82_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Hexa201_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Hexa202_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Hexa203_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tetra41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tetra101_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tetra104_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tetra1015_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Prism62_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Prism156_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Prism159_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Prism1518_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Quad41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Quad84_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Quad89_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tri31_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Tri63_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    double& dNdr_Line21_on_pt_on_shape(uiint igauss, uiint ishape);
    double& dNdr_Line32_on_pt_on_shape(uiint igauss, uiint ishape);
    void Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index);
    void dNdx_on_pt(const uiint& igauss, vvdouble& dNdx);
    void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdx);
    void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& ielem, double dNdx[]);
    void detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ);
    void Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w);
    uiint shapetype_hexa81();
    uiint shapetype_hexa82();
    uiint shapetype_hexa201();
    uiint shapetype_hexa202();
    uiint shapetype_hexa203();
    uiint shapetype_tetra41();
    uiint shapetype_tetra101();
    uiint shapetype_tetra104();
    uiint shapetype_tetra1015();
    uiint shapetype_prism62();
    uiint shapetype_prism156();
    uiint shapetype_prism159();
    uiint shapetype_prism1518();
    uiint shapetype_quad41();
    uiint shapetype_quad84();
    uiint shapetype_quad89();
    uiint shapetype_tri31();
    uiint shapetype_tri63();
    uiint shapetype_line21();
    uiint shapetype_line32();
    uiint GetNumOfBoundaryNodeMesh();
    uiint GetNumOfBoundaryFaceMesh();
    uiint GetNumOfBoundaryEdgeMesh();
    uiint GetNumOfBoundaryVolumeMesh();
    uiint GetBNDType_BNodeMesh(const uiint& ibmesh);
    uiint GetBNDType_BFaceMesh(const uiint& ibmesh);
    uiint GetBNDType_BEdgeMesh(const uiint& ibmesh);
    uiint GetBNDType_BVolumeMesh(const uiint& ibmesh);
    uiint getNeumannType();
    uiint getDirichletType();
    uiint GetNumOfBNode_BNodeMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BFaceMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BEdgeMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BVolumeMesh(const uiint& ibmesh);
    uiint GetNumOfDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint GetNumOfDOF_BFaceMesh(const uiint& ibmesh);
    uiint GetNumOfDOF_BEdgeMesh(const uiint& ibmesh);
    uiint GetNumOfDOF_BVolumeMesh(const uiint& ibmesh);
    uiint GetDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& idof);
    uiint GetDOF_BFaceMesh(const uiint& ibmesh, const uiint& idof);
    uiint GetDOF_BEdgeMesh(const uiint& ibmesh, const uiint& idof);
    uiint GetDOF_BVolumeMesh(const uiint& ibmesh, const uiint& idof);
    double& GetBNodeValue_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof);
    double& GetBNodeValue_BFaceMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    double& GetBNodeValue_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    double& GetBNodeValue_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    uiint& GetNodeID_BNode_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint GetNumOfBFace(const uiint& ibmesh);
    double& GetBFaceValue(const uiint& ibmesh, const uiint& ibface, const uiint& dof);
    uiint GetNumOfBEdge(const uiint& ibmesh);
    double& GetBEdgeValue(const uiint& ibmesh, const uiint& ibedge, const uiint& dof);
    uiint GetNumOfBVolume(const uiint& ibmesh);
    double& GetBVolumeValue(const uiint& ibmesh, const uiint& ibvol, const uiint& dof);
    uiint GetNumOfNode_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint& GetNodeID_BFace(const uiint& ibmesh, const uiint& ibface, const uiint& ibnode);
    uiint GetNumOfNode_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint& GetNodeID_BEdge(const uiint& ibmesh, const uiint& ibedge, const uiint& ibnode);
    uiint GetNumOfNode_BVolume(const uiint& ibmesh, const uiint& ibvol);
    uiint& GetNodeID_BVolume(const uiint& ibmesh, const uiint& ibvol, const uiint& ibnode);
    uiint GetBNodeMesh_NameLength(const uiint& ibmesh);
    string& GetBNodeMesh_Name(const uiint& ibmesh);
    uiint GetBFaceMesh_NameLength(const uiint& ibmesh);
    string& GetBFaceMesh_Name(const uiint& ibmesh);
    uiint GetBVolumeMesh_NameLength(const uiint& ibmesh);
    string& GetBVolumeMesh_Name(const uiint& ibmesh);
    uiint GetBEdgeMesh_NameLength(const uiint& ibmesh);
    string& GetBEdgeMesh_Name(const uiint& ibmesh);
    uiint GetEdgeID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint GetElemID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint GetFaceID_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint GetElemID_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint GetElemID_BVolume(const uiint& ibmesh, const uiint& ibvol);
    int& GetRank();
    int& GetNumOfProcess();
    int AllReduce(void* sendbuf, void* recvbuf, int buf_size, int datatype, int op, int commworld);
    int Barrier(int commworld);
    int Abort(int commworld, int error);
    int AllGather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    int Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
    int Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    int Bcast(void* buf, int cnt, MPI_Datatype type, int root, MPI_Comm comm);
    uiint GetNumOfNeibPE(const uiint& imesh);
    uiint& GetTransRank(const uiint& imesh, const uiint& ipe);
    void Send_Recv_R(double* buf, const int& num_of_node, const int& dof_size, const int& trans_rank);
    void Send_Recv_I(int* buf, const int& num_of_node, const int& dof_size, const int& trans_rank );
    uiint GetNumOfCommMesh();
    uiint GetNumOfCommNode(const uiint& icmesh);
    uiint& GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode);
    uiint GetNumOfElementGroup();
    uiint GetNumOfElementID(const uiint& iGrp);
    uiint& GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index);
    uiint GetElementGroupName_Length(const uiint& iGrp);
    string& GetElementGroupName(const uiint& iGrp);
    void LoggerMode(const uiint& mode);
    void LoggerDevice(const uiint& mode, const uiint& device);
    void LoggerInfo(const uiint& mode, char* message);
    void LoggerInfo(const uiint& mode, const char* message);
    uiint getErrorMode();
    uiint getWarnMode();
    uiint getInfoMode();
    uiint getDebugMode();
    uiint getDiskDevice();
    uiint getDisplayDevice();
};
#endif
}
