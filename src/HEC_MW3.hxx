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

    bool mb_file;//file読み込み後？

    string msInputFileName;
    string msOutputFileName;    
    
    CGMGModel    *mpGMGModel;
    CMeshFactory  *mpFactory;
    CHecMPI      *mpMPI;


    //APIで使用する作業用変数
    CAssyModel  *mpAssy;
    CAssyMatrix *mpAssyMatrix;
    CAssyVector *mpRHSAssyVector;//右辺ベクトル
    CAssyVector *mpSolAssyVector;//解ベクトル
    CMesh       *mpMesh;
    CElement    *mpElement;

    //API 形状関数・導関数
    CShapeHexa     *mpShapeHexa;
    CShapeHexaNic  *mpShapeHexaNic;
    CShapeTetra    *mpShapeTetra;
    CShapePrism    *mpShapePrism;
    CShapeQuad     *mpShapeQuad;
    CShapeTriangle *mpShapeTriangle;
    CShapeLine     *mpShapeLine;

    //API 形状間数種類カタログ
    CShapeFunctionCatalog *mpShapeCatalog;
    
    //形状関数 番号 => 辺番号 変換
    CISTR2Edge *mpISTR2Edge;
    CEdge2ISTR *mpEdge2ISTR;


    //dNdxのテンポラリー変数(calc_dNdxコール時に一旦,保持する)
    vvvdouble mvdNdx;

public:
    //----
    // 初期化・終了処理
    //----
    __declspec(dllexport) int Initialize(int argc, char** argv, const char* path);//引数:MPI_Init用, *.cntファイルパス
    __declspec(dllexport) int Finalize();  // Loggerのログファイルclose
    
    //----
    // File関連
    //----
    __declspec(dllexport) int FileRead(); // ファイル入力(MW3書式)
    __declspec(dllexport) int FileWrite();// ファイル出力

    //----
    // Solver
    //----
    __declspec(dllexport) void GeneLinearAlgebra(const uint& nNumOfAlgebra, uint* vDOF);
    __declspec(dllexport) void SelectAlgebra(const uint& iequ);// mpAssyMatrix等にロード

    __declspec(dllexport) int Matrix_Add_Elem(const uint& iMesh, const uint& iElem, double *ElemMatrix);
    __declspec(dllexport) int Matrix_Add_Node(const uint& iMesh, const uint& iNodeID, const uint& jNodeID, double *NodalMatrix);

    __declspec(dllexport) void Matrix_Clear(const uint& iMesh);// Matrix 0 clear
    __declspec(dllexport) void Vector_Clear(const uint& iMesh);// Matrix 0 clear

    __declspec(dllexport) int Set_BC_Mat_SolVec(int iMesh, int iNode, int iDof, double value1, double value2);//行列対角項, 解ベクトル
    __declspec(dllexport) int Set_BC_Mat_RHS(int iMesh, int iNode, int iDof, double value1, double value2);   //行列対角項, 右辺ベクトル
    __declspec(dllexport) int Set_BC_RHS(int iMesh, int iNode, int iDof, double value);                       //右辺ベクトル
    
    __declspec(dllexport) void Sample_Set_BC(int iMesh);
    
    __declspec(dllexport) int Solve(uint iter_max, double tolerance, uint method, uint precondition);
    
    
    //----
    // Refiner
    //----
    __declspec(dllexport) int Refine();//ファイル指定のマルチグリッドMeshデータ生成


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
    __declspec(dllexport) uint GetNumOfAssembleModel();//Assemble Modelの個数==階層数(mMGLevel+1)
    __declspec(dllexport) void SelectAssembleModel(const uint& mgLevel);
    // Meshの選択
    __declspec(dllexport) uint GetNumOfMeshPart();//Meshパーツの個数
    __declspec(dllexport) void SelectMeshPart_ID(const uint& mesh_id);
    __declspec(dllexport) void SelectMeshPart_IX(const uint& index);
    // 要素の選択
    __declspec(dllexport) void SelectElement_ID(const uint& elem_id);
    __declspec(dllexport) void SelectElement_IX(const uint& index);
    // 選択されたMeshの要素の情報
    __declspec(dllexport) uint GetElementType();
    __declspec(dllexport) uint GetNumOfElementVert();//要素の頂点数
    __declspec(dllexport) void GetElementVertNodeID(int* vNodeID);//要素の頂点のノードID
    __declspec(dllexport) uint GetNumOfElementEdge();//要素の辺数
    __declspec(dllexport) void GetElementEdgeNodeID(int* vNodeID);//要素の辺のノードID
    // 選択されたMeshのノード(Scalar or Vector)
    __declspec(dllexport) void GetNodeCoord(const uint& node_id, double& x, double& y, double& z);//ノードの座標
    __declspec(dllexport) uint GetNumOfDOF(const uint& node_id);
    __declspec(dllexport) uint GetNumOfScalar(const uint& node_id);
    __declspec(dllexport) uint GetNumOfVector(const uint& node_id);
    __declspec(dllexport) uint& GetNodeType(const uint& node_id);
    __declspec(dllexport) void SetNodeValue(const uint& node_id, double value[]);
    __declspec(dllexport) void SetNodeValue(const uint& node_id, const uint& idof, const double& value);
    __declspec(dllexport) void GetNodeValue(const uint& node_id, double value[]);
    __declspec(dllexport) double& GetNodeValue(const uint& node_id, const uint& idof);
    // 選択されたMeshのノード(ScalarVector)
    __declspec(dllexport) void SetSVNodeValue(const uint& node_id, double v_value[], double s_value[]);
    __declspec(dllexport) void SetSVNodeValue(const uint& node_id, const uint& v_dof, const double& v_value, const uint& s_dof, const double& s_value);
    __declspec(dllexport) void GetSVNodeValue(const uint& node_id, double v_value[], double s_value[]);
    __declspec(dllexport) void GetSVNodeValue(const uint& node_id, const uint& v_dof, double& v_value, const uint& s_dof, double& s_value);
    // 選択されたMeshのノード数、要素数
    __declspec(dllexport) uint  getNodeSize(){ return mpMesh->getNodeSize();}
    __declspec(dllexport) uint  getElementSize(){ return mpMesh->getElementSize();}
    __declspec(dllexport) uint  getNodeSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    __declspec(dllexport) uint  getElementSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}


    //----
    // メモリーの解放(Node::AggElemIDの解放):Mesh操作によるデータ構築の最後に実行
    //----
    __declspec(dllexport) void FinalizeRefine();


    //----
    // 節点タイプ for Fortran
    //----
    __declspec(dllexport) uint nodetype_s();
    __declspec(dllexport) uint nodetype_v();
    __declspec(dllexport) uint nodetype_sv();
    //----
    // 要素タイプ for Fortran
    //----
    __declspec(dllexport) uint elemtype_hexa();
    __declspec(dllexport) uint elemtype_hexa2();
    __declspec(dllexport) uint elemtype_tetra();
    __declspec(dllexport) uint elemtype_tetra2();
    __declspec(dllexport) uint elemtype_prism();
    __declspec(dllexport) uint elemtype_prism2();
    __declspec(dllexport) uint elemtype_quad();
    __declspec(dllexport) uint elemtype_quad2();
    __declspec(dllexport) uint elemtype_triangle();
    __declspec(dllexport) uint elemtype_triangle2();
    __declspec(dllexport) uint elemtype_line();
    __declspec(dllexport) uint elemtype_line2();



    // --
    // 積分点数を返す
    // --
    __declspec(dllexport) uint& NumOfIntegPoint(const uint& shapeType);

    //--
    // 形状関数
    //--
    // N(積分点ごと)
    __declspec(dllexport) void ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N);//積分点での形状関数:N[igauss] を返す.
    __declspec(dllexport) void ShapeFunc_on_pt(uint shapeType, uint igauss, double N[]);// Fortran API
    // N(丸ごと)
    __declspec(dllexport) void ShapeFunc(const uint& shapeType, vvdouble& N);//全積分点の形状関数:N を丸ごと返す.
    // ----
    // Fortran API用 形状関数
    // ----
    // N[igauss][ishape]
    __declspec(dllexport) double& ShapeFunc_Hexa81(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa82(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa201(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa202(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Hexa203(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra41(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra101(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra104(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Tetra1015(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Prism62(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Prism156(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Prism159(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Prism1518(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Quad41(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Quad84(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Quad89(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Triangle31(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Triangle63(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Line21(int igauss, int ishape);
    __declspec(dllexport) double& ShapeFunc_Line32(int igauss, int ishape);

    //--
    // 自然座標での導関数
    //--
    // dN/dr(積分点ごと)
    __declspec(dllexport) void dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr);
    // dN/dr(まるごと)
    __declspec(dllexport) void dNdr(const uint& shapeType, vvvdouble& dNdr);
    // ----
    // Fortran API用 導関数
    // ----
    // dNdr[igauss][ishape][iaxis]
    __declspec(dllexport) void dNdr(const uint& shapeType, double dNdr[]);
    __declspec(dllexport) double& dNdr_Hexa81_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Hexa82_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Hexa201_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Hexa202_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Hexa203_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tetra41_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tetra101_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tetra104_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tetra1015_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Prism62_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Prism156_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Prism159_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Prism1518_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Quad41_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Quad84_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Quad89_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tri31_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Tri63_on_pt_on_shape(int igauss, int ishape, int iaxis);
    __declspec(dllexport) double& dNdr_Line21_on_pt_on_shape(int igauss, int ishape);
    __declspec(dllexport) double& dNdr_Line32_on_pt_on_shape(int igauss, int ishape);


    //--
    // 空間座標での導関数
    //--
    // dN/dx 計算のみ(クラス・メンバー:mvdNdxに代入)
    __declspec(dllexport) void Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index);
    // dN/dx(積分点ごと): mvdNdx[igaus] ,Calculate_dNdx(...)実行後に使用
    __declspec(dllexport) void dNdx_on_pt(const uint& igauss, vvdouble& dNdx);
    // dN/dx(まるごと)
    __declspec(dllexport) void dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdx);
    // dN/dx(1次元配列  Fortran用途
    __declspec(dllexport) void dNdx(const uint& elemType, const uint& numOfInteg, const uint& ielem, double dNdx[]);

    //--
    // J行列式 |J|
    //--
    //これを呼び出す直前に使用した,dNdx計算の detJ の値
    //--
    __declspec(dllexport) void detJacobian(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& detJ);

    //--
    // Gauss積分点の重み:Weight
    //--
    __declspec(dllexport) void Weight(const uint& elemType, const uint& numOfInteg, const uint& igauss, double& w);

    //----
    // 形状関数タイプ for Fortran
    //----
    __declspec(dllexport) uint shapetype_hexa81();
    __declspec(dllexport) uint shapetype_hexa82();
    __declspec(dllexport) uint shapetype_hexa201();
    __declspec(dllexport) uint shapetype_hexa202();
    __declspec(dllexport) uint shapetype_hexa203();
    __declspec(dllexport) uint shapetype_tetra41();
    __declspec(dllexport) uint shapetype_tetra101();
    __declspec(dllexport) uint shapetype_tetra104();
    __declspec(dllexport) uint shapetype_tetra1015();
    __declspec(dllexport) uint shapetype_prism62();
    __declspec(dllexport) uint shapetype_prism156();
    __declspec(dllexport) uint shapetype_prism159();
    __declspec(dllexport) uint shapetype_prism1518();
    __declspec(dllexport) uint shapetype_quad41();
    __declspec(dllexport) uint shapetype_quad84();
    __declspec(dllexport) uint shapetype_quad89();
    __declspec(dllexport) uint shapetype_tri31();
    __declspec(dllexport) uint shapetype_tri63();
    __declspec(dllexport) uint shapetype_line21();
    __declspec(dllexport) uint shapetype_line32();


    //--
    // Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
    //--
    __declspec(dllexport) uint GetNumOfBoundaryNodeMesh();
    __declspec(dllexport) uint GetNumOfBoundaryFaceMesh();
    __declspec(dllexport) uint GetNumOfBoundaryEdgeMesh();
    __declspec(dllexport) uint GetNumOfBoundaryVolumeMesh();
    __declspec(dllexport) uint GetNumOfBNode_BNodeMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfBNode_BFaceMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfBNode_BEdgeMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfBNode_BVolumeMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfDOF_BNodeMesh(const uint& ibmesh, const uint& ibnode);
    __declspec(dllexport) uint GetNumOfDOF_BFaceMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfDOF_BEdgeMesh(const uint& ibmesh);
    __declspec(dllexport) uint GetNumOfDOF_BVolumeMesh(const uint& ibmesh);
    //--
    // BoundaryNode 境界値
    //--
    __declspec(dllexport) double& GetBNodeValue_BNodeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof);
    __declspec(dllexport) double& GetBNodeValue_BFaceMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BEdgeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BVolumeMesh(const uint& ibmesh, const uint& ibnode, const uint& idof, const uint& mgLevel);
    __declspec(dllexport) uint& GetNodeID_BNode_BNodeMesh(const uint& ibmesh, const uint& ibnode);
    __declspec(dllexport) uint& GetNodeID_BNode_BFaceMesh(const uint& ibmesh, const uint& ibnode);
    __declspec(dllexport) uint& GetNodeID_BNode_BEdgeMesh(const uint& ibmesh, const uint& ibnode);
    __declspec(dllexport) uint& GetNodeID_BNode_BVolumeMesh(const uint& ibmesh, const uint& ibnode);
    //--
    // Face, Edge, Volume の境界値
    //--
    __declspec(dllexport) uint GetNumOfBFace(const uint& ibmesh);
    __declspec(dllexport) double& GetBFaceValue(const uint& ibmesh, const uint& ibface, const uint& idof);
    __declspec(dllexport) uint GetNumOfBEdge(const uint& ibmesh);
    __declspec(dllexport) double& GetBEdgeValue(const uint& ibmesh, const uint& ibedge, const uint& idof);
    __declspec(dllexport) uint GetNumOfBVolume(const uint& ibmesh);
    __declspec(dllexport) double& GetBVolumeValue(const uint& ibmesh, const uint& ibvol, const uint& idof);


    //--
    // MPI (直接呼び出したい人向け)
    //--
    __declspec(dllexport) int AllReduce(void* sendbuf, void* recvbuf, int buf_size, int datatype, int op, int commworld);
    __declspec(dllexport) int Barrier(int commworld);
    __declspec(dllexport) int Abort(int commworld, int error);
    __declspec(dllexport) int AllGather(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, MPI_Comm comm);
    __declspec(dllexport) int Gather(void* sendbuf , int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
    __declspec(dllexport) int Scatter(void* sendbuf, int sendcnt, MPI_Datatype sendtype, void* recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    __declspec(dllexport) int& GetRank(); //自分のプロセス-ランクを取得
    __declspec(dllexport) int& GetNumOfProcess();
    __declspec(dllexport) void Send_Recv_R(double* buf, int dof_size);// bufの値を送信, 受信値をNodeとbufに代入.   bufのサイズ == NumOfCommNode * dof_size
    __declspec(dllexport) void Send_Recv_R();// 通信Nodeの値を入れ替えて更新


    //--
    // グループ { select された AssyModel,Meshを対象 }
    //--
    __declspec(dllexport) uint GetNumOfElementGroup();
    __declspec(dllexport) uint GetNumOfElementID(const uint& iGrp);
    __declspec(dllexport) uint& GetElementID_with_ElementGroup(const uint& iGrp, const uint& index);
    __declspec(dllexport) uint GetElementGroupName_Length(const uint& iGrp);
    __declspec(dllexport) string& GetElementGroupName(const uint& iGrp);

    
    //--
    // Logger
    //--
    __declspec(dllexport) void LoggerMode(const uint& mode);
    __declspec(dllexport) void LoggerDevice(const uint& mode, const uint& device);
    __declspec(dllexport) void LoggerInfo(const uint& mode, char* message);
    //--
    // Logger_Parameter for Fortran
    //--
    __declspec(dllexport) uint getErrorMode();
    __declspec(dllexport) uint getWarnMode();
    __declspec(dllexport) uint getInfoMode();
    __declspec(dllexport) uint getDebugMode();
    __declspec(dllexport) uint getDiskDevice();
    __declspec(dllexport) uint getDisplayDevice();
};
#endif
}


