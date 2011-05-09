//
// HEC_MW3.h
//
// pMW
//  => Main & Interface
// 
// (MeshFactory)
// (GMGModel)
// 
//			2009.05.08
//                      2008.11.19
//			k.Takeda
#include "TypeDef.h"//typedef文
#include "HEC_MPI.h"//MPIラッパー

//Utility
#include "FileIO.h"
#include "Logger.h"

//SolutionType
#include "SolutionType.h"

//メッシュモデル
#include "MeshFactory.h"
#include "GMGModel.h"
#include "AssyMatrix.h"

//形状関数
#include "ShapeHexa.h"
#include "ShapeHexaNic.h"
#include "ShapeTetra.h"
#include "ShapePrism.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeLine.h"
//形状関数カタログ
#include "ShapeFunctionCatalog.h"
//形状関数 番号 <=変換=> MW3の要素形状の辺番号
#include "ISTR2Edge.h"
#include "Edge2ISTR.h"

//Jacobian
#include "Jacobian.h"

//boost文字列処理
#include <boost/lexical_cast.hpp>

typedef pair<uint,uint> integPair;

namespace pmw{
#ifndef PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#define PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
class CMW
{
// Singleton Constructor
//
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

    uint mnSolutionType;

    bool mb_file;//file読み込み後？

    string msInputFileName;
    string msOutputFileName;    
    
    pmw::CGMGModel    *mpGMGModel;
    pmw::CMeshFactory  *mpFactory;
    pmw::CHecMPI      *mpMPI;


    //APIで使用する作業用変数
    pmw::CAssyModel  *mpAssy;
    pmw::CAssyMatrix *mpAssyMatrix;
    pmw::CAssyVector *mpRHSAssyVector;//右辺ベクトル
    pmw::CAssyVector *mpSolAssyVector;//解ベクトル
    pmw::CMesh       *mpMesh;
    pmw::CElement    *mpElement;

    //API 形状関数・導関数
    pmw::CShapeHexa     *mpShapeHexa;
    pmw::CShapeHexaNic  *mpShapeHexaNic;
    pmw::CShapeTetra    *mpShapeTetra;
    pmw::CShapePrism    *mpShapePrism;
    pmw::CShapeQuad     *mpShapeQuad;
    pmw::CShapeTriangle *mpShapeTriangle;
    pmw::CShapeLine     *mpShapeLine;

    //API 形状関数_種類カタログ
    pmw::CShapeFunctionCatalog *mpShapeCatalog;
    
    //形状関数 番号 => 辺番号 変換
    pmw::CISTR2Edge *mpISTR2Edge;
    pmw::CEdge2ISTR *mpEdge2ISTR;


    //dNdxのテンポラリー変数(calc_dNdxコール時に一旦,保持する)
    vvvdouble mvdNdx;

public:
    //----
    // 初期化・終了処理
    //----
    int Initialize(int argc, char** argv, const char* path);//引数:MPI_Init用, *.cntファイルパス
    int Finalize();  // Loggerのログファイルclose
    
    //----
    // File関連
    //----
    int FileRead(); // ファイル入力(MW3書式)
    int FileWrite();// ファイル出力

    //----
    // Solver
    //----
    void GeneLinearAlgebra(const uint& nNumOfAlgebra, uint* vDOF);//全Levelに方程式を生成 vDOF:方程式ごとのDOF
    void SelectAlgebra(const uint& iequ);// mpAssyMatrix等にiequ番目の方程式をロード{ Levelは選択されているとする }

    int Matrix_Add_Elem(const uint& iMesh, const uint& iElem, double *ElemMatrix);
    int Matrix_Add_Node(const uint& iMesh, const uint& iNodeID, const uint& jNodeID, double *NodalMatrix);

    void Matrix_Clear(const uint& iMesh);// Matrix all 0 clear
    void Vector_Clear(const uint& iMesh);// Matrix all 0 clear

    //int Set_BC_Mat_SolVec(uint& iMesh, uint& iNode, uint& nDOF, double& value1, double& value2);//行列対角項, 解ベクトル
    //int SetZero_NonDiag(uint& iMesh, uint& iNode, uint& nDOF);
    int Set_BC_Mat_RHS2(uint& iMesh, uint& iNode, uint& nDOF, double& diag_value, double& rhs_value);//対角項=Val,非対角項=0、右辺ベクトル
    int Set_BC_Mat_RHS(uint& iMesh, uint& iNode, uint& nDOF, double& diag_value, double& rhs_value); //行列対角項, 右辺ベクトル
    int Set_BC_RHS(uint& iMesh, uint& iNode, uint& nDOF, double& value);                        //右辺ベクトル
    int Add_BC_RHS(uint& iMesh, uint& iNode, uint& nDOF, double& value);                        //右辺ベクトルへの加算
    
    void  Sample_Set_BC(uint iMesh);
    
    int Solve(uint& iter_max, double& tolerance, uint& method, uint& precondition);

    //--
    // 解ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    void GetSolution_Vector(double* buf, const uint& imesh);
    void GetSolution_AssyVector(double* buf);
    //--
    // 右辺ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    void GetRHS_Vector(double* buf, const uint& imesh);
    void GetRHS_AssyVector(double* buf);

    //--
    // AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
    //--
    void multVector(double* vX, double* vB);

    
    //----
    // MG constructor (Refiner)
    //----
    int  Refine();//ファイル指定のマルチグリッドMeshデータ生成


    
    // ----
    // Fortran_API用 : 基本void(Fortran_API に合わせる為) : 引数のconst以外=>出力用
    // ----
    // *Fortran用にオブジェクトを選択する.(C++使用者には不要)
    //  選択順序;
    //   1. Assemble Model(階層構造から必要な階層(Level)のアセンブル・モデルを取得)
    //   2. Mesh (Assemble Model から,必要なMeshパーツを取得)
    //   3. Element (Meshから,必要な要素を取得)
    //   4. Node (要素の構成NodeIDからNodeを取得)
    // ----
    // Assemble Modelの選択
    uint GetNumOfAssembleModel();//Assemble Modelの個数==階層数(mMGLevel+1)
    void SelectAssembleModel(const uint& mgLevel);
    // Meshの選択
    uint GetNumOfMeshPart();//Meshパーツの個数
    void SelectMeshPart_ID(const uint& mesh_id);
    void SelectMeshPart_IX(const uint& index);
    // 要素の選択
    void SelectElement_ID(const uint& elem_id);
    void SelectElement_IX(const uint& index);
    // 選択されたMeshの要素の情報
    uint GetElementType();
    uint GetNumOfElementVert();//要素の頂点数
    void GetElementVertNodeID(int* vNodeID);//要素の頂点のノードID
    uint GetNumOfElementEdge();//要素の辺数
    void GetElementEdgeNodeID(int* vNodeID);//要素の辺のノードID
    // 選択されたMeshのノード(Scalar or Vector)
    void GetNodeCoord(const uint& node_id, double& x, double& y, double& z);//ノードの座標
    uint GetNumOfDOF(const uint& node_id);
    uint GetNumOfScalar(const uint& node_id);
    uint GetNumOfVector(const uint& node_id);
    uint& GetNodeType(const uint& node_id);
    void SetNodeValue(const uint& node_id, double value[]);
    void SetNodeValue(const uint& node_id, const uint& idof, const double& value);
    void GetNodeValue(const uint& node_id, double value[]);
    double& GetNodeValue(const uint& node_id, const uint& idof);
    // 選択されたMeshのノード(ScalarVector)
    void SetSVNodeValue(const uint& node_id, double v_value[], double s_value[]);
    void SetSVNodeValue(const uint& node_id, const uint& v_dof, const double& v_value, const uint& s_dof, const double& s_value);
    void GetSVNodeValue(const uint& node_id, double v_value[], double s_value[]);
    void GetSVNodeValue(const uint& node_id, const uint& v_dof, double& v_value, const uint& s_dof, double& s_value);
    // 選択されたMeshのノード数、要素数
    uint  getNodeSize(){ return mpMesh->getNodeSize();}
    uint  getElementSize(){ return mpMesh->getElementSize();}
    uint  getNodeSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    uint  getElementSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
    // ID
    uint& getNodeID(const uint& index);
    uint& getElementID(const uint& index);


    //----
    // メモリーの解放(Node::AggElemIDの解放):Mesh操作によるデータ構築の最後に実行
    //----
    void FinalizeRefine();


    //----
    // 節点タイプ for Fortran
    //----
    uint nodetype_s();
    uint nodetype_v();
    uint nodetype_sv();
    //----
    // 要素タイプ for Fortran
    //----
    uint elemtype_hexa();
    uint elemtype_hexa2();
    uint elemtype_tetra();
    uint elemtype_tetra2();
    uint elemtype_prism();
    uint elemtype_prism2();
    uint elemtype_quad();
    uint elemtype_quad2();
    uint elemtype_triangle();
    uint elemtype_triangle2();
    uint elemtype_line();
    uint elemtype_line2();
    //----
    // FrontISTR 要素タイプ
    //----
    uint fistr_elemtype_hexa();
    uint fistr_elemtype_hexa2();
    uint fistr_elemtype_tetra();
    uint fistr_elemtype_tetra2();
    uint fistr_elemtype_prism();
    uint fistr_elemtype_prism2();
    uint fistr_elemtype_quad();
    uint fistr_elemtype_quad2();
    uint fistr_elemtype_triangle();
    uint fistr_elemtype_triangle2();
    uint fistr_elemtype_line();
    uint fistr_elemtype_line2();
    //----
    // FrontISTR 要素タイプ　=> MW3 要素タイプ 変換
    //----
    uint fistr_elemtype_to_mw3_elemtype(const uint& fistr_elemtype);
    //----
    // MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
    //----
    uint mw3_elemtype_to_fistr_elemtype(const uint& mw3_elemtype);


    // --
    // 積分点数を返す
    // --
    uint& NumOfIntegPoint(const uint& shapeType);
    
    //--
    // 形状関数
    //--
    // N(積分点ごと)
    void ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N);//積分点での形状関数:N[igauss] を返す.
    void ShapeFunc_on_pt(uint shapeType, uint igauss, double N[]);// Fortran API
    // N(丸ごと)
    void ShapeFunc(const uint& shapeType, vvdouble& N);//全積分点の形状関数:N を丸ごと返す.
    // ----
    // Fortran API用 形状関数
    // ----
    // N[igauss][ishape]
    double& ShapeFunc_Hexa81(int igauss, int ishape);
    double& ShapeFunc_Hexa82(int igauss, int ishape);
    double& ShapeFunc_Hexa201(int igauss, int ishape);
    double& ShapeFunc_Hexa202(int igauss, int ishape);
    double& ShapeFunc_Hexa203(int igauss, int ishape);
    double& ShapeFunc_Tetra41(int igauss, int ishape);
    double& ShapeFunc_Tetra101(int igauss, int ishape);
    double& ShapeFunc_Tetra104(int igauss, int ishape);
    double& ShapeFunc_Tetra1015(int igauss, int ishape);
    double& ShapeFunc_Prism62(int igauss, int ishape);
    double& ShapeFunc_Prism156(int igauss, int ishape);
    double& ShapeFunc_Prism159(int igauss, int ishape);
    double& ShapeFunc_Prism1518(int igauss, int ishape);
    double& ShapeFunc_Quad41(int igauss, int ishape);
    double& ShapeFunc_Quad84(int igauss, int ishape);
    double& ShapeFunc_Quad89(int igauss, int ishape);
    double& ShapeFunc_Triangle31(int igauss, int ishape);
    double& ShapeFunc_Triangle63(int igauss, int ishape);
    double& ShapeFunc_Line21(int igauss, int ishape);
    double& ShapeFunc_Line32(int igauss, int ishape);

    //--
    // 自然座標での導関数
    //--
    // dN/dr(積分点ごと)
    void dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr);
    // dN/dr(まるごと)
    void dNdr(const uint& shapeType, vvvdouble& dNdr);
    // ----
    // Fortran API用 導関数
    // ----
    // dNdr[igauss][ishape][iaxis]
    void dNdr(const uint& shapeType, double dNdr[]);
    double& dNdr_Hexa81_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Hexa82_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Hexa201_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Hexa202_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Hexa203_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tetra41_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tetra101_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tetra104_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tetra1015_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Prism62_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Prism156_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Prism159_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Prism1518_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Quad41_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Quad84_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Quad89_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tri31_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Tri63_on_pt_on_shape(int igauss, int ishape, int iaxis);
    double& dNdr_Line21_on_pt_on_shape(int igauss, int ishape);
    double& dNdr_Line32_on_pt_on_shape(int igauss, int ishape);


    //--
    // 空間座標での導関数
    //--
    // dN/dx 計算のみ(クラス・メンバー:mvdNdxに代入)
    void Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index);
    // dN/dx(積分点ごと): mvdNdx[igaus] ,Calculate_dNdx(...)実行後に使用
    void dNdx_on_pt(const uint& igauss, vvdouble& dNdx);
    // dN/dx(まるごと)
    void dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdx);
    // dN/dx(1次元配列  Fortran用途
    void dNdx(const uint& elemType, const uint& numOfInteg, const uint& ielem, double dNdx[]);

    //--
    // J行列式 |J|
    //--
    //これを呼び出す直前に使用した,dNdx計算の detJ の値
    //--
    void detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ);
    
    //--
    // Gauss積分点の重み:Weight
    //--
    void Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w);

    //----
    // 形状関数タイプ for Fortran
    //----
    uint shapetype_hexa81();
    uint shapetype_hexa82();
    uint shapetype_hexa201();
    uint shapetype_hexa202();
    uint shapetype_hexa203();
    uint shapetype_tetra41();
    uint shapetype_tetra101();
    uint shapetype_tetra104();
    uint shapetype_tetra1015();
    uint shapetype_prism62();
    uint shapetype_prism156();
    uint shapetype_prism159();
    uint shapetype_prism1518();
    uint shapetype_quad41();
    uint shapetype_quad84();
    uint shapetype_quad89();
    uint shapetype_tri31();
    uint shapetype_tri63();
    uint shapetype_line21();
    uint shapetype_line32();


    //--
    // Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
    //--
    // BoundaryMesh数
    uint GetNumOfBoundaryNodeMesh();
    uint GetNumOfBoundaryFaceMesh();
    uint GetNumOfBoundaryEdgeMesh();
    uint GetNumOfBoundaryVolumeMesh();
    // BoundaryType { Neumann || Dirichlet }
    uint GetBNDType_BNodeMesh(const uint& ibmesh);
    uint GetBNDType_BFaceMesh(const uint& ibmesh);
    uint GetBNDType_BEdgeMesh(const uint& ibmesh);
    uint GetBNDType_BVolumeMesh(const uint& ibmesh);
    // BoundaryTypeを表す型(定数)
    uint getNeumannType();
    uint getDirichletType();
    // 境界節点数
    uint GetNumOfBNode_BNodeMesh(const uint& ibmesh);
    uint GetNumOfBNode_BFaceMesh(const uint& ibmesh);
    uint GetNumOfBNode_BEdgeMesh(const uint& ibmesh);
    uint GetNumOfBNode_BVolumeMesh(const uint& ibmesh);
    // BoundaryのDOF数
    uint GetNumOfDOF_BNodeMesh(const uint& ibmesh, const uint& ibnode);
    uint GetNumOfDOF_BFaceMesh(const uint& ibmesh);
    uint GetNumOfDOF_BEdgeMesh(const uint& ibmesh);
    uint GetNumOfDOF_BVolumeMesh(const uint& ibmesh);
    // Boundary DOF番号("DOFインデックス"に対応するDOF番号)
    uint GetDOF_BNodeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof);
    uint GetDOF_BFaceMesh(const uint& ibmesh, const uint& idof);
    uint GetDOF_BEdgeMesh(const uint& ibmesh, const uint& idof);
    uint GetDOF_BVolumeMesh(const uint& ibmesh, const uint& idof);


    //--
    // BoundaryNode 境界値
    //--
    double& GetBNodeValue_BNodeMesh(const uint& ibmesh, const uint& ibnode, const uint& dof);
    double& GetBNodeValue_BFaceMesh(const uint& ibmesh, const uint& ibnode, const uint& dof, const uint& mgLevel);
    double& GetBNodeValue_BEdgeMesh(const uint& ibmesh, const uint& ibnode, const uint& dof, const uint& mgLevel);
    double& GetBNodeValue_BVolumeMesh(const uint& ibmesh, const uint& ibnode, const uint& dof, const uint& mgLevel);
    uint& GetNodeID_BNode_BNodeMesh(const uint& ibmesh, const uint& ibnode);
    uint& GetNodeID_BNode_BFaceMesh(const uint& ibmesh, const uint& ibnode);
    uint& GetNodeID_BNode_BEdgeMesh(const uint& ibmesh, const uint& ibnode);
    uint& GetNodeID_BNode_BVolumeMesh(const uint& ibmesh, const uint& ibnode);
    //--
    // Face, Edge, Volume の境界値
    //--
    uint GetNumOfBFace(const uint& ibmesh);
    double& GetBFaceValue(const uint& ibmesh, const uint& ibface, const uint& dof);
    uint GetNumOfBEdge(const uint& ibmesh);
    double& GetBEdgeValue(const uint& ibmesh, const uint& ibedge, const uint& dof);
    uint GetNumOfBVolume(const uint& ibmesh);
    double& GetBVolumeValue(const uint& ibmesh, const uint& ibvol, const uint& dof);
    //--
    // Boundaryの名称
    //--
    uint GetBNodeMesh_NameLength(const uint& ibmesh);
    string& GetBNodeMesh_Name(const uint& ibmesh);
    uint GetBFaceMesh_NameLength(const uint& ibmesh);
    string& GetBFaceMesh_Name(const uint& ibmesh);
    uint GetBVolumeMesh_NameLength(const uint& ibmesh);
    string& GetBVolumeMesh_Name(const uint& ibmesh);
    uint GetBEdgeMesh_NameLength(const uint& ibmesh);
    string& GetBEdgeMesh_Name(const uint& ibmesh);
    //--
    // Entity_ID of BoundaryMesh (for FrontISTR)
    //--
    uint GetEdgeID_BEdge(const uint& ibmesh, const uint& ibedge);
    uint GetElemID_BEdge(const uint& ibmesh, const uint& ibedge);
    uint GetFaceID_BFace(const uint& ibmesh, const uint& ibface);
    uint GetElemID_BFace(const uint& ibmesh, const uint& ibface);
    uint GetElemID_BVolume(const uint& ibmesh, const uint& ibvol);


    //--
    // MPI (直接呼び出したい人向け)
    //--
    int AllReduce(void* sendbuf, void* recvbuf, int buf_size, int datatype, int op, int commworld);
    int Barrier(int commworld);
    int Abort(int commworld, int error);
    int AllGather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    int Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int& GetRank();//自分のプロセス-ランクを取得
    int& GetNumOfProcess();
    // 以下の３メソッドは、ペア
    uint GetNumOfNeibPE(const uint& imesh);//Meshパーツが通信する相手の数
    uint& GetTransRank(const uint& imesh, const uint& ipe);//通信Mesh毎のランク番号
    void Send_Recv_R(double* buf, const uint& num_of_node, const uint& dof_size, const uint& trans_rank);//bufの値を送信-受信

    void Send_Recv_R(double* buf, const uint& nDOF);// bufの値を送信, 受信値をNodeとbufに代入.   bufのサイズ == NumOfCommNode * dof_size
    void Send_Recv_R();                          // 通信Nodeの値を入れ替えて更新





    //--
    // グループ { select された AssyModel,Meshを対象 }
    //--
    uint GetNumOfElementGroup();
    uint GetNumOfElementID(const uint& iGrp);
    uint& GetElementID_with_ElementGroup(const uint& iGrp, const uint& index);
    uint GetElementGroupName_Length(const uint& iGrp);
    string& GetElementGroupName(const uint& iGrp);


    //--
    // Logger
    //--
    void LoggerMode(const uint& mode);
    void LoggerDevice(const uint& mode, const uint& device);
    void LoggerInfo(const uint& mode, char* message);
    //--
    // Logger_Parameter for Fortran
    //--
    uint getErrorMode();
    uint getWarnMode();
    uint getInfoMode();
    uint getDebugMode();
    uint getDiskDevice();
    uint getDisplayDevice();
};
#endif
}


