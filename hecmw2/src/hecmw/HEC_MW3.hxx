/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/HEC_MW3.hxx
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
    __declspec(dllexport) static CMW*  Instance(){
	static CMW moMW;
	return &moMW;
    }
private:
    __declspec(dllexport) CMW();
public:
    __declspec(dllexport) virtual ~CMW();
protected:
    FileIO::CFileIO   *mpFileIO;
    Utility::CLogger  *mpLogger;    
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
	uiint mnSolutionType;
public:
    __declspec(dllexport) uiint Initialize( int argc, char** argv); 
    __declspec(dllexport) uiint Initialize_fstr(int argc, char** argv, string& ctrlname);
    __declspec(dllexport) uiint Finalize();  
    __declspec(dllexport) void  Banner_Display();
    __declspec(dllexport) uiint FileRead(string& basename, bool bBinary);
    __declspec(dllexport) uiint FileRead_fstr(bool bBinary);             
    __declspec(dllexport) uiint FileDebugWrite();          
    __declspec(dllexport) uiint FileWriteRes(const uiint& nStep, bool bBinary);
    __declspec(dllexport) uiint SetRestart(const uiint& nStep, bool bBinary);
    __declspec(dllexport) void PrintRlt_Start(const uiint& nStep, bool bBinary);
    __declspec(dllexport) void PrintRlt_P(const uiint& width, const char* format, ... );
    __declspec(dllexport) void PrintRlt_R(const uiint& width, const char* format, ... );
    __declspec(dllexport) void PrintRlt_End();
    __declspec(dllexport) void PrintMicroAVS_Basis();
    __declspec(dllexport) void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    __declspec(dllexport) void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    __declspec(dllexport) void PrintMicroAVS_FEM();  
    __declspec(dllexport) string getFstr_FileName_Mesh();
    __declspec(dllexport) string getFstr_FileName_Control();
    __declspec(dllexport) string getFstr_FileName_Result();
    __declspec(dllexport) string getFstr_FileName_Restart();
    __declspec(dllexport) string getFstr_FileName_PartIn();
    __declspec(dllexport) string getFstr_FileName_PartOut();
    __declspec(dllexport) string getFstr_FileName_VisMesh();
    __declspec(dllexport) string getFstr_FileName_VisIn();
    __declspec(dllexport) string getFstr_FileName_VisOut();
    __declspec(dllexport) void GeneLinearAlgebra(const uiint& nNumOfAlgebra, uiint* vDOF);
    __declspec(dllexport) void SelectAlgebra(const uiint& iequ);
    __declspec(dllexport) uiint Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    __declspec(dllexport) uiint Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double *NodalMatrix);
    __declspec(dllexport) void Matrix_Clear(const uiint& iMesh);
    __declspec(dllexport) void Vector_Clear(const uiint& iMesh);
    __declspec(dllexport) uiint Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& sol_value);
    __declspec(dllexport) uiint Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& rhs_value); 
    __declspec(dllexport) uiint Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                       
    __declspec(dllexport) uiint Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);  
	__declspec(dllexport) void  Sample_Set_BC(uiint iMesh);
    __declspec(dllexport) uiint Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition);
    __declspec(dllexport) void GetSolution_Vector(double* buf, const uiint& imesh);
    __declspec(dllexport) void GetSolution_AssyVector(double* buf);
    __declspec(dllexport) void GetRHS_Vector(double* buf, const uiint& imesh);
    __declspec(dllexport) void GetRHS_AssyVector(double* buf);
    __declspec(dllexport) double& GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    __declspec(dllexport) double& GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    __declspec(dllexport) uiint& GetSolutionVector_DOF();
    __declspec(dllexport) uiint& GetRHSVector_DOF();
    __declspec(dllexport) void multVector(double* vX, double* vB);
    __declspec(dllexport) uiint Refine(const uiint& nNumOfRefine);
    __declspec(dllexport) uiint GetNumOfAssembleModel();
    __declspec(dllexport) void SelectAssembleModel(const uiint& mgLevel);
    __declspec(dllexport) uiint GetNumOfMeshPart();
    __declspec(dllexport) void SelectMeshPart_ID(const uiint& mesh_id);
    __declspec(dllexport) void SelectMeshPart_IX(const uiint& index);
    __declspec(dllexport) void SelectElement_ID(const uiint& elem_id);
    __declspec(dllexport) void SelectElement_IX(const uiint& index);
    __declspec(dllexport) uiint GetElementType();
    __declspec(dllexport) uiint GetNumOfElementVert();
    __declspec(dllexport) void GetElementVertNodeID(iint* vNodeID);
    __declspec(dllexport) uiint GetNumOfElementEdge();
    __declspec(dllexport) void GetElementEdgeNodeID(iint* vNodeID);
    __declspec(dllexport) void GetNodeCoord(const uiint& node_id, double& x, double& y, double& z);
    __declspec(dllexport) uiint GetNumOfDOF(const uiint& node_id);
    __declspec(dllexport) uiint GetNumOfScalar(const uiint& node_id);
    __declspec(dllexport) uiint GetNumOfVector(const uiint& node_id);
    __declspec(dllexport) uiint& GetNodeType(const uiint& node_id);
    __declspec(dllexport) uiint  getNodeSize(){ return mpMesh->getNodeSize();}
    __declspec(dllexport) uiint  getElementSize(){ return mpMesh->getElementSize();}
    __declspec(dllexport) uiint  getNodeSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    __declspec(dllexport) uiint  getElementSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
    __declspec(dllexport) uiint& getNodeID(const uiint& index);
    __declspec(dllexport) uiint& getElementID(const uiint& index);
    __declspec(dllexport) uiint& getNodeIndex(const uiint& id);   
    __declspec(dllexport) uiint& getElementIndex(const uiint& id);
    __declspec(dllexport) void constructNodeConnectFEM(const uiint& node_id);
    __declspec(dllexport) void getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL);
    __declspec(dllexport) void getNodeConnectFEM_Item(uiint itemU[], uiint itemL[]);
    __declspec(dllexport) void getNodeConnectFEM_Item_F(iint itemU[], iint itemL[]);
    __declspec(dllexport) uiint getNumOfAggregateElement(const uiint& node_id);
    __declspec(dllexport) uiint& getAggregateElementID(const uiint& node_id, const uiint& ielem);
    __declspec(dllexport) void FinalizeRefine();
    __declspec(dllexport) uiint nodetype_s();
    __declspec(dllexport) uiint nodetype_v();
    __declspec(dllexport) uiint nodetype_sv();
    __declspec(dllexport) uiint elemtype_hexa();
    __declspec(dllexport) uiint elemtype_hexa2();
    __declspec(dllexport) uiint elemtype_tetra();
    __declspec(dllexport) uiint elemtype_tetra2();
    __declspec(dllexport) uiint elemtype_prism();
    __declspec(dllexport) uiint elemtype_prism2();
    __declspec(dllexport) uiint elemtype_quad();
    __declspec(dllexport) uiint elemtype_quad2();
    __declspec(dllexport) uiint elemtype_triangle();
    __declspec(dllexport) uiint elemtype_triangle2();
    __declspec(dllexport) uiint elemtype_line();
    __declspec(dllexport) uiint elemtype_line2();
    __declspec(dllexport) uiint fistr_elemtype_hexa();
    __declspec(dllexport) uiint fistr_elemtype_hexa2();
    __declspec(dllexport) uiint fistr_elemtype_tetra();
    __declspec(dllexport) uiint fistr_elemtype_tetra2();
    __declspec(dllexport) uiint fistr_elemtype_prism();
    __declspec(dllexport) uiint fistr_elemtype_prism2();
    __declspec(dllexport) uiint fistr_elemtype_quad();
    __declspec(dllexport) uiint fistr_elemtype_quad2();
    __declspec(dllexport) uiint fistr_elemtype_triangle();
    __declspec(dllexport) uiint fistr_elemtype_triangle2();
    __declspec(dllexport) uiint fistr_elemtype_line();
    __declspec(dllexport) uiint fistr_elemtype_line2();
    __declspec(dllexport) uiint fistr_elemtype_to_mw3_elemtype(const uiint& fistr_elemtype);
    __declspec(dllexport) uiint mw3_elemtype_to_fistr_elemtype(const uiint& mw3_elemtype);
    __declspec(dllexport) uiint& NumOfIntegPoint(const uiint& shapeType);
    __declspec(dllexport) void ShapeFunc_on_pt(const uiint& shapeType, const uiint& igauss, vdouble& N);
    __declspec(dllexport) void ShapeFunc_on_pt(uiint shapeType, uiint igauss, double N[]);
    __declspec(dllexport) void ShapeFunc(const uiint& shapeType, vvdouble& N);
    __declspec(dllexport) double& ShapeFunc_Hexa81(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa82(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa201(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa202(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa203(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra41(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra101(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra104(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra1015(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Prism62(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Prism156(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Prism159(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Prism1518(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Quad41(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Quad84(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Quad89(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Triangle31(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Triangle63(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Line21(uiint igauss, uiint ishape);
    __declspec(dllexport) double& ShapeFunc_Line32(uiint igauss, uiint ishape);
    __declspec(dllexport) void dNdr_on_pt(const uiint& shapeType, const uiint& igauss, vvdouble& dNdr);
    __declspec(dllexport) void dNdr(const uiint& shapeType, vvvdouble& dNdr);
    __declspec(dllexport) void dNdr(const uiint& shapeType, double dNdr[]);
    __declspec(dllexport) double& dNdr_Hexa81_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Hexa82_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Hexa201_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Hexa202_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Hexa203_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tetra41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tetra101_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tetra104_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tetra1015_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Prism62_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Prism156_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Prism159_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Prism1518_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Quad41_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Quad84_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Quad89_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tri31_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Tri63_on_pt_on_shape(uiint igauss, uiint ishape, uiint iaxis);
    __declspec(dllexport) double& dNdr_Line21_on_pt_on_shape(uiint igauss, uiint ishape);
    __declspec(dllexport) double& dNdr_Line32_on_pt_on_shape(uiint igauss, uiint ishape);
    __declspec(dllexport) void Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index);
    __declspec(dllexport) void dNdx_on_pt(const uiint& igauss, vvdouble& dNdx);
    __declspec(dllexport) void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdx);
    __declspec(dllexport) void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& ielem, double dNdx[]);
    __declspec(dllexport) void detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ);
    __declspec(dllexport) void Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w);
    __declspec(dllexport) uiint shapetype_hexa81();
    __declspec(dllexport) uiint shapetype_hexa82();
    __declspec(dllexport) uiint shapetype_hexa201();
    __declspec(dllexport) uiint shapetype_hexa202();
    __declspec(dllexport) uiint shapetype_hexa203();
    __declspec(dllexport) uiint shapetype_tetra41();
    __declspec(dllexport) uiint shapetype_tetra101();
    __declspec(dllexport) uiint shapetype_tetra104();
    __declspec(dllexport) uiint shapetype_tetra1015();
    __declspec(dllexport) uiint shapetype_prism62();
    __declspec(dllexport) uiint shapetype_prism156();
    __declspec(dllexport) uiint shapetype_prism159();
    __declspec(dllexport) uiint shapetype_prism1518();
    __declspec(dllexport) uiint shapetype_quad41();
    __declspec(dllexport) uiint shapetype_quad84();
    __declspec(dllexport) uiint shapetype_quad89();
    __declspec(dllexport) uiint shapetype_tri31();
    __declspec(dllexport) uiint shapetype_tri63();
    __declspec(dllexport) uiint shapetype_line21();
    __declspec(dllexport) uiint shapetype_line32();
    __declspec(dllexport) uiint GetNumOfBoundaryNodeMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryFaceMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryEdgeMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryVolumeMesh();
    __declspec(dllexport) uiint GetBNDType_BNodeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BVolumeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint getNeumannType();
    __declspec(dllexport) uiint getDirichletType();
    __declspec(dllexport) uiint GetNumOfBNode_BNodeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BVolumeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfDOF_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfDOF_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfDOF_BVolumeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BFaceMesh(const uiint& ibmesh, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BEdgeMesh(const uiint& ibmesh, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BVolumeMesh(const uiint& ibmesh, const uiint& idof);
    __declspec(dllexport) double& GetBNodeValue_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof);
    __declspec(dllexport) double& GetBNodeValue_BFaceMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) uiint& GetNodeID_BNode_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfBFace(const uiint& ibmesh);
    __declspec(dllexport) double& GetBFaceValue(const uiint& ibmesh, const uiint& ibface, const uiint& dof);
    __declspec(dllexport) uiint GetNumOfBEdge(const uiint& ibmesh);
    __declspec(dllexport) double& GetBEdgeValue(const uiint& ibmesh, const uiint& ibedge, const uiint& dof);
    __declspec(dllexport) uiint GetNumOfBVolume(const uiint& ibmesh);
    __declspec(dllexport) double& GetBVolumeValue(const uiint& ibmesh, const uiint& ibvol, const uiint& dof);
    __declspec(dllexport) uiint GetNumOfNode_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint& GetNodeID_BFace(const uiint& ibmesh, const uiint& ibface, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfNode_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint& GetNodeID_BEdge(const uiint& ibmesh, const uiint& ibedge, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfNode_BVolume(const uiint& ibmesh, const uiint& ibvol);
    __declspec(dllexport) uiint& GetNodeID_BVolume(const uiint& ibmesh, const uiint& ibvol, const uiint& ibnode);
    __declspec(dllexport) uiint GetBNodeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBNodeMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBFaceMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBFaceMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBVolumeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBVolumeMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBEdgeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBEdgeMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetEdgeID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint GetElemID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint GetFaceID_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint GetElemID_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint GetElemID_BVolume(const uiint& ibmesh, const uiint& ibvol);
    __declspec(dllexport) int& GetRank(); 
    __declspec(dllexport) int& GetNumOfProcess();    
    __declspec(dllexport) int AllReduce(void* sendbuf, void* recvbuf, int buf_size, int datatype, int op, int commworld);
    __declspec(dllexport) int Barrier(int commworld);
    __declspec(dllexport) int Abort(int commworld, int error);
    __declspec(dllexport) int AllGather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    __declspec(dllexport) int Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    __declspec(dllexport) int Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    __declspec(dllexport) int Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
    __declspec(dllexport) int Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    __declspec(dllexport) int Bcast(void* buf, int cnt, MPI_Datatype type, int root, MPI_Comm comm);
    __declspec(dllexport) uiint GetNumOfNeibPE(const uiint& imesh);
    __declspec(dllexport) uiint& GetTransRank(const uiint& imesh, const uiint& ipe);
    __declspec(dllexport) void Send_Recv_R(double* buf, const int& num_of_node, const int& dof_size, const int& trans_rank);
    __declspec(dllexport) void Send_Recv_I(int* buf, const int& num_of_node, const int& dof_size, const int& trans_rank );
    __declspec(dllexport) uiint GetNumOfCommMesh();
    __declspec(dllexport) uiint GetNumOfCommNode(const uiint& icmesh);
    __declspec(dllexport) uiint& GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode);
    __declspec(dllexport) uiint GetNumOfElementGroup();
    __declspec(dllexport) uiint GetNumOfElementID(const uiint& iGrp);
    __declspec(dllexport) uiint& GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index);
    __declspec(dllexport) uiint GetElementGroupName_Length(const uiint& iGrp);
    __declspec(dllexport) string& GetElementGroupName(const uiint& iGrp);
    __declspec(dllexport) void LoggerMode(const uiint& mode);
    __declspec(dllexport) void LoggerDevice(const uiint& mode, const uiint& device);
    __declspec(dllexport) void LoggerInfo(const uiint& mode, char* message);
    __declspec(dllexport) void LoggerInfo(const uiint& mode, const char* message);
    __declspec(dllexport) uiint getErrorMode();
    __declspec(dllexport) uiint getWarnMode();
    __declspec(dllexport) uiint getInfoMode();
    __declspec(dllexport) uiint getDebugMode();
    __declspec(dllexport) uiint getDiskDevice();
    __declspec(dllexport) uiint getDisplayDevice();
};
#endif
}
