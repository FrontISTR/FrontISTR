/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/HEC_MW3.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
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
#include "CodeType.h"
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


#ifdef REVOCAP_REFINE
#include "rcapRefiner.h"
#include "rcapRefinerMacros.h"
#endif

#include <cstdio>  // printf
#include <cstdlib> // calloc

#include <cstdarg>
#include <iomanip>
#include <cstring>
typedef pair<uiint,uiint> integPair;

namespace pmw
{
#ifndef PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#define PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
class CMW
{
public:
    static CMW*  Instance() {
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

    uiint mnSolutionType;//----
    uiint mnCodeType;//--------

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
    CContactMesh *mpConMesh;

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

    vuint mv_ItemL;
    vuint mv_ItemU;

    uiint BaseName_BCast(uiint& nLength, string& sName, uiint nType);

    bool mbRefine_Use;//--------- リファイン使用判定

    // API_Fortran 定数
    MPI_Datatype mnMPI_UIINT;
    MPI_Datatype mnMPI_IINT;
    int  mnMyRank;
    int  mnNumOfProcess;

public:
    uiint Initialize( int argc, char** argv);
    uiint Initialize_fstr( int argc, char** argv, string& ctrlname);
    uiint Finalize();
    void  Banner_Display();


#ifdef REVOCAP_REFINE //==========================================REVOCAP_REFINE
    //----
    // REVOCAP_Refiner
    //----
    uiint RevocapRefine(string& filename, const uiint& nRefine);
private:
    void Rcap_BndParamSetup(uiint nElemType, uiint& nNumOfVol,vuint& vBndID, vstring& vBndName, vuint& vBndBType,
                            vuint& vBndDOF, vvuint& vvBndDOF, vuint& vBndType, vuint& vBndCount, vuint& vBndTypeSeq, uiint& nSeq,
                            CBoundaryVolumeMesh *pVolMesh);//境界条件管理変数セットアップ:BndVolume

    void Rcap_CrsElemNodes(size_t nNNode, vector<CElement*> vElem, CIndexBucket *pBucket, int32_t *elemNodes);

    void Rcap_NodeGene(uiint iLevel, uiint nMaxLevel, size_t crsNumOfNode, size_t refineNodeCount, float64_t* resultCoord, CMesh *pCrsMesh, CMesh *pProgMesh);
    void Rcap_ElemGene(int8_t nType, uiint& nIDCount, size_t refineCount, int32_t *refineNodes, CMesh *pProgMesh);

    void Rcap_SetBFaceMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh);
    void Rcap_SetBEdgeMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh);
    void Rcap_SetBVolMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, vuint vBndBType, vstring vBndName, vuint vNumDOF, vvuint vvDOF, CMesh *pProgMesh, CMesh *pMesh);

    vvuint Rcap_NodeIX_SortMerge(size_t nNumOfID, size_t *nEntityCount, size_t nNumOfEntityNode, int32_t **refineNodes);
    void Rcap_BNodeGene_VolMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex);
    void Rcap_BNodeGene_FaceMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex);
    void Rcap_BNodeGene_EdgeMesh(uiint iLevel, uiint nMaxLevel, vuint vBndID, CMesh *pMesh, CMesh *pProgMesh, vvuint& vvNodeIndex);
    void Rcap_NodeNum2Index(vector<map<uiint, uiint> >& vmNodeNum2Index, vvuint& NodeIndex);//逆引き生成

    void Rcap_BVolGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t *nVolCount, uiint nShape, size_t& nNumOfVolNode,
                       int32_t **refineNodes, CMesh *pProgMesh, CMesh *pCrsMesh);
    void Rcap_BFaceGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t *nFaceCount, uiint nShape, size_t& nNumOfFaceNode,
                        int32_t **refineNodes, CMesh *pProgMesh, CMesh *pCrsMesh);
    void Rcap_BEdgeGene(uiint iLevel, vuint& vBndID, vector<map<uiint, uiint> >& vmID2Index, size_t *nEdgeCount, uiint nShape, size_t& nNumOfEdgeNode,
                        int32_t **refineNodes, CMesh *pProgMesh, CMesh *pCrsMesh);

    void Rcap_ElemSearch(uiint& nElementID, vector<CNode*>& vNode, CMesh *pProgMesh);
    void Rcap_ElemFaceSearch(uiint& nElementID, uiint& nFaceID, vector<CNode*>& vNode, CMesh *pProgMesh);
    void Rcap_ElemEdgeSearch(uiint& nElementID, uiint& nEdgeID, vector<CNode*>& vNode, CMesh *pProgMesh);
    void Rcap_ElemPointSearch(uiint& nElementID, uiint& nLocalID, CNode* pNode, CMesh *pProgMesh);

    void Rcap_BNodeValueDist(uiint iLevel, CBoundaryVolumeMesh *pBVolMesh, vdouble& vVal, vuint& vBNodeIndex, map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfEntityNode, int32_t *crsNodes);// Dirichlet値 分配
    void Rcap_BNodeValueDist(uiint iLevel, CBoundaryFaceMesh *pBFaceMesh, vdouble& vVal, vuint& vBNodeIndex, map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfEntityNode, int32_t *crsNodes);// Dirichlet値 分配
    void Rcap_BNodeValueDist(uiint iLevel, CBoundaryEdgeMesh *pBEdgeMesh, vdouble& vVal, vuint& vBNodeIndex, map<uiint,uiint>& NodeNum2BNodeNum, uiint nNumOfEntityNode, int32_t *crsNodes);// Dirichlet値 分配
    void Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryFaceMesh *pBFaceMesh, CBoundaryFace *pBFace);
    void Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryEdgeMesh *pBEdgeMesh, CBoundaryEdge *pBEdge);
    void Rcap_EquivalentNodalForce(uiint iLevel, uiint dof, CBoundaryVolumeMesh *pBVolMesh, CBoundaryVolume *pBVol);

    void Rcap_BFaceMeshDebug(vuint& vBndID, CMesh *pProgMesh);

    //
    // Comm
    //
    void Rcap_SetCommMesh2(uiint iLevel, map<uiint,uiint> mComID, vvuint& vvNum4EachType, CMesh *pCrsMesh, CMesh *pProgMesh);
    void Rcap_NodesPoint2vuint(map<uiint,uiint>& mComID, uiint nNumOfFaceNode, size_t *nFaceCount, int32_t **refineNodes, vvuint& vvRefineNodes);

    // Comm:形状別リファイン・ノード配列を一つにまとめる
    void Rcap_SumRefineNodes( vvuint& vvRefineNodes, vvuint& vvQuad, vvuint& vvQuad2, vvuint& vvTri, vvuint& vvTri2, vvuint& vvBeam, vvuint& vvBeam2, vvuint& vvPoint);
    void Rcap_CommFaceCount( vvuint& vvNum4EachType, map<uiint,uiint> mComID, size_t* commQuadCount, size_t* commQuad2Count,
                             size_t* commTriCount, size_t* commTri2Count, size_t* commBeamCount, size_t* commBeam2Count,
                             size_t* commPointCount);
    void Rcap_CommNodeGene( map<uiint,uiint> mComID, vvuint& vNum4EachType, vvvuint& vvvElemNum, vector<map<uiint, vuint> >& vmaElemNIndex,
                            CMesh *pProgMesh, vvuint& vvRefineNodes, vector<map<uiint, uiint> >& maNNum2CommNNum);
    void Rcap_CommFaceGene( uiint iLevel, map<uiint,uiint>& mComID, vvvuint& vvvElemNum, vector<map<uiint, vuint> >& vmaElemNIndex,
                            vvuint& vvRefineNodes, CMesh *pProgMesh, vector<map<uiint, uiint> >& maNNum2CommNNum );

#endif //=========================================================REVOCAP_REFINE
private:
    void SortMerge(vuint& vec);

public:
    void clearMW3(const uiint& mgLevel);


    uiint FileRead(string& basename, bool bBinary);
    uiint FileRead_fstr(bool bBinary);
    uiint FileDebugWrite();

    //--
    // リスタート
    //--
    uiint FileWriteRes(const uiint& nStep, bool bBinary);
    uiint SetRestart(const uiint& nStep, uiint& nAppNumOfLevel, uiint& nAppNumEquation, bool bBinary);

    void PrintRlt_Start(const uiint& nStep, bool bBinary);
    void PrintRlt(const uiint& width, const char* format, ... );
    void PrintRlt_End();

    void PrintMicroAVS_Basis(const uiint& ieq);
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void PrintMicroAVS_FEM();

    void recVTK_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recVTK_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void PrintVTK_FEM();

    void recUNS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recUNS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void PrintUNS_FEM();

    string getFstr_FileName_Mesh();
    string getFstr_FileName_Control();
    string getFstr_FileName_Result();
    string getFstr_FileName_Restart();
    string getFstr_FileName_PartIn();
    string getFstr_FileName_PartOut();
    string getFstr_FileName_VisMesh();
    string getFstr_FileName_VisIn();
    string getFstr_FileName_VisOut();
    string getFstr_RefineCADFitName();
    uiint  getFstr_RefineNum();
    uiint  getFstr_RefineType();

    //--
    // 線形方程式の生成
    //--
    void GeneLinearAlgebra(const uiint& nNumOfAlgebra, const uiint& nGlobalNumOfMesh, uiint** vvDOF, double** vvTransCoeff);
    void GeneLinearAlgebra(const uiint& nNumOfAlgebra, const uiint& nGlobalNumOfMesh, uiint** vvDOF);//--- 方程式生成Base
    void SelectAlgebra(const uiint& iequ);

    //--
    // 行列・ベクトル(線形方程式)
    //--
    uiint Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    uiint Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double *NodalMatrix);
    void Matrix_Clear(const uiint& iMesh);
    void Vector_Clear(const uiint& iMesh);
    void AssyMatrix_Clear();//----------------------------------11.11.01
    void AssyVector_Clear();//----------------------------------11.11.01
    uiint Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& sol_value);
    uiint Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& rhs_value);

    uiint Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);
    uiint Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);
    uiint Set_BC_NL_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);//---- 節点集中荷重として扱う:Rank大には値は入らない.
    uiint Add_BC_NL_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);//---- 節点集中荷重として扱う:Rank大には値は入らない.

    //--
    // 線形ソルバー
    //--
    uiint Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition);

    void GetSolution_Vector(double* buf, const uiint& imesh);
    void GetSolution_AssyVector(double* buf);
    void GetRHS_Vector(double* buf, const uiint& imesh);
    void GetRHS_AssyVector(double* buf);
    void GetRHS_Load(double* buf, const uiint& imesh);//--- 非線形構造解析の残差力:右辺ベクトルをsumupしたベクトル
    void GetRHS_AssyLoad(double* buf);//------------------- 非線形構造解析の残差力:右辺ベクトルをsumupしたベクトル
    double& GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    double& GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    uiint& GetSolutionVector_DOF(const uiint& imesh);
    uiint& GetRHSVector_DOF(const uiint& imesh);
    uiint& GetSolutionVector_DOF(const uiint& iLevel, const uiint& imesh, const uiint& ieq);
    uiint& GetRHSVector_DOF(const uiint& iLevel, const uiint& imesh, const uiint& ieq);

    void dump_AssyMatrix();//debug用 行列ダンプ(CRT表示)
    void dump_RHSAssyVector();//debug用 行列ダンプ(CRT表示)

////    //--
////    // 分散ベクトル 分散行列 : アセンブルモデル間の通信 => 未実装
////    // --------------------------------
////    // # ユーザー定義ベクトル の通信界面操作:Update,SumupはMatVecで利用
////    // --------------------------------
////    //--
////    void Update_Add(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF);//------通信界面のベクトル値 加算 : アセンブルモデル間の通信 => 未実装
    void Update(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF);//----------通信界面のベクトル値 (rank小から大へ値をコピー) : アセンブルモデル間の通信 => 未実装
    void Sumup(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF);//-----------通信界面のベクトル値の加算(rank大から小に値を加算)
////    void Update_Average(const uiint& iLevel, const uiint& imesh, double* vec, const uiint& nDOF);//--通信界面のベクトル値 平均 : アセンブルモデル間の通信 => 未実装
////
////    void Update_Assy_Add(const uiint& iLevel, double* assy_vec, uiint* vDOF);
////    void Update_Assy(const uiint& iLevel, double* assy_vec, uiint* vDOF);
////    void Sumup_Assy(const uiint& iLevel, double* assy_vec, uiint* vDOF);
////    void Update_Assy_Average(const uiint& iLevel, double* assy_vec, uiint* vDOF);

    // void MatVec(const uiint& iLevel, const uiint& imesh, double* mat, double* x, double* y, const uiint& nDOF);//汎用版 Ax=y
    void MatVec_Assy(const uiint& iLevel, const uiint& ieq, double* assy_x, double* assy_y);     //AssyMatrix版 Ax=y
    void MatVec(const uiint& iLevel, const uiint& imesh, const uiint& ieq, double* x, double* y);//MatrixBCRS版 Ax=y


    uiint  Refine(const uiint& nNumOfRefine);


    uiint GetNumOfAssembleModel();
    void SelectAssembleModel(const uiint& mgLevel);
    uiint GetNumOfMeshPart();
    void SelectMeshPart_ID(const uiint& mesh_id);
    void SelectMeshPart_IX(const uiint& index);
    uiint GetMeshID_Num(const uiint& index);//---------- Mesh 通し番号(Index)から、ID
    uiint GetMeshIndex_Num(const uiint& id);//---------- Mesh ID番号から、Index

    void SelectElement_ID(const uiint& elem_id);
    void SelectElement_IX(const uiint& index);
    uiint GetElementType();
    uiint GetNumOfElementVert();
    void GetElementVertNodeID(iint* vNodeID);
    void GetElementVertNodeIndex(iint* vNodeIndex);
    uiint GetNumOfElementEdge();
    uiint GetNumOfElementFace();
    void GetElementEdgeNodeID(iint* vNodeID);
    uiint GetElementFaceElementID(uiint faceIndex);
    uiint GetNumOfElementEdgeElement(uiint edgeIndex);
    uiint GetElementEdgeElementID(uiint edgeIndex, uiint i);
    void GetNodeCoord(const uiint& node_id, double& x, double& y, double& z);
    uiint GetNumOfDOF(const uiint& node_id);
    uiint GetNumOfScalar(const uiint& node_id);
    uiint GetNumOfVector(const uiint& node_id);
    uiint& GetNodeType(const uiint& node_id);
    uiint  getNodeSize();
    uiint  getElementSize();
    uiint  getNodeSize(uiint iMesh);
    uiint  getElementSize(uiint iMesh);
    uiint& getNodeID(const uiint& index);
    uiint& getElementID(const uiint& index);
    uiint& getNodeIndex(const uiint& id);
    uiint& getElementIndex(const uiint& id);

    uiint  getNumOfParentNode(const uiint& id, const uiint& nLevel);
    uiint  getParentNodeID(const uiint& id, const uiint& nLevel, const uiint& index);

    void constructNodeConnectFEM(const uiint& node_id);
    void getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL);
    void getNodeConnectFEM_Item(uiint itemU[], uiint itemL[]);
    void getNodeConnectFEM_Item_F(iint itemU[], iint itemL[]);

    uiint getNumOfAggregateElement(const uiint& node_id);
    uiint& getAggregateElementID(const uiint& node_id, const uiint& ielem);

    void setupNeighbors();

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

    //--
    // 境界条件
    //--
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

    //--
    // Test
    //--
    double& GetBNode_X_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);//Test
    double& GetBNode_Y_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);//Test
    double& GetBNode_Z_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);//Test

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

    //--
    // MPI
    //--
    int GetRank();
    int GetNumOfProcess();
    iint AllReduce(void* sendbuf, void* recvbuf, iint buf_size, MPI_Datatype datatype, MPI_Op op, MPI_Comm commworld);
    iint Barrier(MPI_Comm commworld);
    iint Abort(MPI_Comm commworld, int error);
    iint AllGather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    iint Gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    iint Scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    iint Recv(void* buf, iint count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
    iint Send(void* buf, iint count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    iint Bcast(void* buf, iint cnt, MPI_Datatype type, int root, MPI_Comm comm);
    iint GetNumOfNeibPE(const uiint& imesh);
    iint GetTransRank(const uiint& imesh, const uiint& ipe);
    void Send_Recv_R(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
    void Send_Recv_I(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank );
    void Sendrecv_R(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
    void Sendrecv_I(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
    MPI_Datatype MPI_UIINT();
    MPI_Datatype MPI_IINT();


    //--
    // 通信テーブル:CommMesh2
    //--
    uiint GetNumOfCommMesh();
    uiint GetNumOfCommNode(const uiint& icmesh);
    uiint& GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode);



    //--
    // 接合面:ContactMesh
    //--
    uiint GetNumOfContactMesh();
    uiint GetContactMeshID(const uiint& icont);//---- ContactMesh ID番号




    uiint GetNumOfElementGroup();
    uiint GetNumOfElementID(const uiint& iGrp);
    uiint& GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index);
    uiint GetElementGroupName_Length(const uiint& iGrp);
    string& GetElementGroupName(const uiint& iGrp);

    void LoggerMode(const uiint& mode);
    void LoggerDevice(const uiint& mode, const uiint& device);
    void LoggerInfoMssg(const uiint& mode, const char* message);
    void LoggerInfo(const uiint& mode, const char* format, ... );
    uiint getErrorMode();
    uiint getWarnMode();
    uiint getInfoMode();
    uiint getDebugMode();
    uiint getDiskDevice();
    uiint getDisplayDevice();
};
#endif
}// namespace pmw;

////#ifdef __cplusplus
////extern "C" {
////#endif
////    MPI_Datatype mpi_uiint();
////    MPI_Datatype mpi_iint();
////    int get_rank();
////    int get_num_of_process();
////
////    iint allreduce(void* sendbuf, void* recvbuf, iint buf_size, MPI_Datatype datatype, MPI_Op op, MPI_Comm commworld);
////    iint barrier(MPI_Comm commworld);
////    iint mpi_abort(MPI_Comm commworld, int error);
////    iint allgather(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
////    iint gather(void* sendbuf , iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
////    iint scatter(void* sendbuf, iint sendcnt, MPI_Datatype sendtype, void* recvbuf, iint recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
////    iint recv(void* buf, iint count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
////    iint send(void* buf, iint count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
////    iint bcast(void* buf, iint cnt, MPI_Datatype type, int root, MPI_Comm comm);
////    iint get_num_of_neib_pe(const uiint& imesh);
////    iint get_trans_rank(const uiint& imesh, const uiint& ipe);
////    void send_recv_r(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
////    void send_recv_i(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank );
////    void sendrecv_r(double* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
////    void sendrecv_i(iint* buf, const iint& num_of_node, const iint& dof_size, const int& trans_rank);
////#ifdef __cplusplus
////}
////#endif


