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

////boost文字列処理
//#include <boost/lexical_cast.hpp>

//REVOCAP_Refiner
#include "rcapRefiner.h"

#include <cstdarg>
#include <iomanip>
#include <cstring>//strlen

typedef pair<uiint,uiint> integPair;

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

    // --
    // rank番号を付加するファイル : HEC_MW3の管理
    // --
    string msInputFileName;  //メッシュファイル名
    string msOutputFileName; //Data_Check ファイル名(debug)
    //-----------------------//解析制御ファイル名
    string msResFileName;    //リスタートファイル名
    string msRltFileName;    //リザルトファイル名
    //-----------------------//パーティショナー入力(単一)メッシュファイル名
    string msPartOutFileName;//パーティショナー出力(分散)メッシュファイル名
    //-----------------------//VisMesh:メッシュファイル名
    //-----------------------//VisIn  :リザルトファイル名
    //-----------------------//VisOut :可視化ファイル名
    
    CGMGModel     *mpGMGModel;
    CMeshFactory  *mpFactory;
    CHecMPI       *mpMPI;


    //APIで使用する作業用変数
    CAssyModel  *mpAssy;
    CAssyMatrix *mpAssyMatrix;
    CAssyVector *mpRHSAssyVector;//右辺ベクトル
    CAssyVector *mpSolAssyVector;//解ベクトル
    CMesh       *mpMesh;
    CElement    *mpElement;
    CCommMesh2  *mpComMesh;

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

    //線形方程式数, 各方程式のDOF(リスタートで利用)
    vuint mvAlgebraEquation;

    //NodeConnectFEMのテンポラリー変数
    vuint mv_ItemL;
    vuint mv_ItemU;
    
    //プロセス０から、分散プロセスへファイル名をセットアップ
    uiint BaseName_BCast(int& nLength, string& sName, int nType);
public:
    //----
    // 初期化・終了処理
    //----
    __declspec(dllexport) uiint Initialize( int argc, char** argv); // 標準 [全体制御ファイルは使用しない]
    __declspec(dllexport) uiint Initialize_fstr(int argc, char** argv, string& ctrlname);//Init(標準), ctrlファイル名
    __declspec(dllexport) uiint Finalize();  // Loggerのログファイルclose
    //----
    // 横断幕
    //----
    __declspec(dllexport) void  Banner_Display();


    //----
    // REVOCAP_Refiner
    //----
    __declspec(dllexport) void set_RevocapMesh(const uiint& nRefine);//REVOCAP_Refinerへのメッシュの登録
    
    

    //----
    // File関連
    //----
    __declspec(dllexport) uiint FileRead(string& basename, bool bBinary);//メッシュ(*.msh) 読み込み:標準
    __declspec(dllexport) uiint FileRead_fstr(bool bBinary);             //メッシュ(*.msh) 読み込み:ctrlファイル使用
    __declspec(dllexport) uiint FileDebugWrite();          // データ・チェック(*.out)
    __declspec(dllexport) uiint FileWriteRes(const uiint& nStep, bool bBinary);// Resファイル出力(*.res)
    // --
    // Restart data のセット(ファイル読み込み => 解ベクトルにセット)
    // --
    __declspec(dllexport) uiint SetRestart(const uiint& nStep, bool bBinary);// Resファイル => Algebra生成
    // --
    // リザルト・ファイル出力
    // --
    __declspec(dllexport) void PrintRlt_Start(const uiint& nStep, bool bBinary);
    __declspec(dllexport) void PrintRlt_P(const uiint& width, const char* format, ... );//リザルト出力(可変長引数:ポインター)
    __declspec(dllexport) void PrintRlt_R(const uiint& width, const char* format, ... );//リザルト出力(可変長引数:値)
    __declspec(dllexport) void PrintRlt_End();
    // *.inp 出力
    __declspec(dllexport) void PrintMicroAVS_Basis();// 基礎変数の出力
    __declspec(dllexport) void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);//出力ラベルの登録
    __declspec(dllexport) void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);//出力変数の登録
    __declspec(dllexport) void PrintMicroAVS_FEM();  // 登録変数の出力
    
    // --
    // ファイル名(fstr)
    // --
    __declspec(dllexport) string getFstr_FileName_Mesh();
    __declspec(dllexport) string getFstr_FileName_Control();
    __declspec(dllexport) string getFstr_FileName_Result();
    __declspec(dllexport) string getFstr_FileName_Restart();
    __declspec(dllexport) string getFstr_FileName_PartIn();
    __declspec(dllexport) string getFstr_FileName_PartOut();
    __declspec(dllexport) string getFstr_FileName_VisMesh();
    __declspec(dllexport) string getFstr_FileName_VisIn();
    __declspec(dllexport) string getFstr_FileName_VisOut();

    
    //----
    // Solver
    //----
    __declspec(dllexport) void GeneLinearAlgebra(const uiint& nNumOfAlgebra, uiint* vDOF);
    __declspec(dllexport) void SelectAlgebra(const uiint& iequ);// mpAssyMatrix等にロード

    __declspec(dllexport) uiint Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    __declspec(dllexport) uiint Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double *NodalMatrix);

    __declspec(dllexport) void Matrix_Clear(const uiint& iMesh);// Matrix 0 clear
    __declspec(dllexport) void Vector_Clear(const uiint& iMesh);// Matrix 0 clear


    __declspec(dllexport) uiint Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& rhs_value);//対角項=Val,非対角項=0、右辺ベクトル
    __declspec(dllexport) uiint Set_BC_Mat_RHS(uiint iMesh, uiint iNodeID, uiint nDOF, double& diag_value, double& rhs_value);//行列対角項, 右辺ベクトル
    __declspec(dllexport) uiint Set_BC_RHS(uiint iMesh, uiint iNodeID, uiint nDOF, double value);                       //右辺ベクトル
    __declspec(dllexport) uiint Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                //右辺ベクトルへの加算

    
    __declspec(dllexport) uiint Solve(iint iter_max, double tolerance, iint method, iint precondition);

    //--
    // 解ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    __declspec(dllexport) void GetSolution_Vector(double* buf, const uiint& imesh);
    __declspec(dllexport) void GetSolution_AssyVector(double* buf);
    //--
    // 右辺ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    __declspec(dllexport) void GetRHS_Vector(double* buf, const uiint& imesh);
    __declspec(dllexport) void GetRHS_AssyVector(double* buf);
    //--
    // ベクトル値を取得
    //--
    __declspec(dllexport) double& GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    __declspec(dllexport) double& GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    //--
    // ベクトルDOFを取得
    //--
    __declspec(dllexport) uiint& GetSolutionVector_DOF();
    __declspec(dllexport) uiint& GetRHSVector_DOF();
    //--
    // AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
    //--
    __declspec(dllexport) void multVector(double* vX, double* vB);

    
    
    //----
    // Refiner
    //----
    __declspec(dllexport) uiint Refine(const uiint& nNumOfRefine);//マルチグリッドMeshデータ生成


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
    __declspec(dllexport) uiint GetNumOfAssembleModel();//Assemble Modelの個数==階層数(mMGLevel+1)
    __declspec(dllexport) void SelectAssembleModel(const uiint& mgLevel);
    // Meshの選択
    __declspec(dllexport) uiint GetNumOfMeshPart();//Meshパーツの個数
    __declspec(dllexport) void SelectMeshPart_ID(const uiint& mesh_id);
    __declspec(dllexport) void SelectMeshPart_IX(const uiint& index);
    // 要素の選択
    __declspec(dllexport) void SelectElement_ID(const uiint& elem_id);
    __declspec(dllexport) void SelectElement_IX(const uiint& index);
    // 選択されたMeshの要素の情報
    __declspec(dllexport) uiint GetElementType();
    __declspec(dllexport) uiint GetNumOfElementVert();//要素の頂点数
    __declspec(dllexport) void GetElementVertNodeID(iint* vNodeID);//要素の頂点のノードID
    __declspec(dllexport) uiint GetNumOfElementEdge();//要素の辺数
    __declspec(dllexport) void GetElementEdgeNodeID(iint* vNodeID);//要素の辺のノードID
    // 選択されたMeshのノード(Scalar or Vector)
    __declspec(dllexport) void GetNodeCoord(const uiint& node_id, double& x, double& y, double& z);//ノードの座標
    __declspec(dllexport) uiint GetNumOfDOF(const uiint& node_id);
    __declspec(dllexport) uiint GetNumOfScalar(const uiint& node_id);
    __declspec(dllexport) uiint GetNumOfVector(const uiint& node_id);
    __declspec(dllexport) uiint& GetNodeType(const uiint& node_id);
////    __declspec(dllexport) void SetNodeValue(const uiint& node_id, double value[]);
////    __declspec(dllexport) void SetNodeValue(const uiint& node_id, const uiint& idof, const double& value);
////    __declspec(dllexport) void GetNodeValue(const uiint& node_id, double value[]);
////    __declspec(dllexport) double& GetNodeValue(const uiint& node_id, const uiint& idof);
    // 選択されたMeshのノード(ScalarVector)
////    __declspec(dllexport) void SetSVNodeValue(const uiint& node_id, double v_value[], double s_value[]);
////    __declspec(dllexport) void SetSVNodeValue(const uiint& node_id, const uiint& v_dof, const double& v_value, const uiint& s_dof, const double& s_value);
////    __declspec(dllexport) void GetSVNodeValue(const uiint& node_id, double v_value[], double s_value[]);
////    __declspec(dllexport) void GetSVNodeValue(const uiint& node_id, const uiint& v_dof, double& v_value, const uiint& s_dof, double& s_value);
    // 選択されたMeshのノード数、要素数
    __declspec(dllexport) uiint  getNodeSize(){ return mpMesh->getNodeSize();}
    __declspec(dllexport) uiint  getElementSize(){ return mpMesh->getElementSize();}
    __declspec(dllexport) uiint  getNodeSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    __declspec(dllexport) uiint  getElementSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
    // ID  &&  index
    __declspec(dllexport) uiint& getNodeID(const uiint& index);
    __declspec(dllexport) uiint& getElementID(const uiint& index);
    __declspec(dllexport) uiint& getNodeIndex(const uiint& id);   //ID -> index
    __declspec(dllexport) uiint& getElementIndex(const uiint& id);//ID -> index
    //--
    // 節点コネクティビティ itemU:上三角のNodeID, itemL:下三角のNodeID
    //--
    __declspec(dllexport) void constructNodeConnectFEM(const uiint& node_id);
    __declspec(dllexport) void getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL);//node_id : 上三角のNodeID数, 下三角のNodeID数
    __declspec(dllexport) void getNodeConnectFEM_Item(uiint itemU[], uiint itemL[]);//node_id : 上三角のNodeID配列, 下三角のNodeID配列
    __declspec(dllexport) void getNodeConnectFEM_Item_F(iint itemU[], iint itemL[]);//FortranAPI用 (引数タイプ : signed int )
    //--
    // 節点周囲 要素
    //--
    __declspec(dllexport) uiint getNumOfAggregateElement(const uiint& node_id);
    __declspec(dllexport) uiint& getAggregateElementID(const uiint& node_id, const uiint& ielem);




    //----
    // メモリーの解放(Node::AggElemIDの解放):Mesh操作によるデータ構築の最後に実行
    //----
    __declspec(dllexport) void FinalizeRefine();


    //----
    // 節点タイプ for Fortran
    //----
    __declspec(dllexport) uiint nodetype_s();
    __declspec(dllexport) uiint nodetype_v();
    __declspec(dllexport) uiint nodetype_sv();
    //----
    // 要素タイプ for Fortran
    //----
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
    //----
    // FrontISTR 要素タイプ
    //----
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
    //----
    // FrontISTR 要素タイプ　=> MW3 要素タイプ 変換
    //----
    __declspec(dllexport) uiint fistr_elemtype_to_mw3_elemtype(const uiint& fistr_elemtype);
    //----
    // MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
    //----
    __declspec(dllexport) uiint mw3_elemtype_to_fistr_elemtype(const uiint& mw3_elemtype);


    // --
    // 積分点数を返す
    // --
    __declspec(dllexport) uiint& NumOfIntegPoint(const uiint& shapeType);

    //--
    // 形状関数
    //--
    // N(積分点ごと)
    __declspec(dllexport) void ShapeFunc_on_pt(const uiint& shapeType, const uiint& igauss, vdouble& N);//積分点での形状関数:N[igauss] を返す.
    __declspec(dllexport) void ShapeFunc_on_pt(uiint shapeType, uiint igauss, double N[]);// Fortran API
    // N(丸ごと)
    __declspec(dllexport) void ShapeFunc(const uiint& shapeType, vvdouble& N);//全積分点の形状関数:N を丸ごと返す.
    // ----
    // Fortran API用 形状関数
    // ----
    // N[igauss][ishape]
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

    //--
    // 自然座標での導関数
    //--
    // dN/dr(積分点ごと)
    __declspec(dllexport) void dNdr_on_pt(const uiint& shapeType, const uiint& igauss, vvdouble& dNdr);
    // dN/dr(まるごと)
    __declspec(dllexport) void dNdr(const uiint& shapeType, vvvdouble& dNdr);
    // ----
    // Fortran API用 導関数
    // ----
    // dNdr[igauss][ishape][iaxis]
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


    //--
    // 空間座標での導関数
    //--
    // dN/dx 計算のみ(クラス・メンバー:mvdNdxに代入)
    __declspec(dllexport) void Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index);
    // dN/dx(積分点ごと): mvdNdx[igaus] ,Calculate_dNdx(...)実行後に使用
    __declspec(dllexport) void dNdx_on_pt(const uiint& igauss, vvdouble& dNdx);
    // dN/dx(まるごと)
    __declspec(dllexport) void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdx);
    // dN/dx(1次元配列  Fortran用途
    __declspec(dllexport) void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& ielem, double dNdx[]);

    //--
    // J行列式 |J|
    //--
    //これを呼び出す直前に使用した,dNdx計算の detJ の値
    //--
    __declspec(dllexport) void detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ);

    //--
    // Gauss積分点の重み:Weight
    //--
    __declspec(dllexport) void Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w);

    //----
    // 形状関数タイプ for Fortran
    //----
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


    //--
    // Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
    //--
    // BoundaryMesh数
    __declspec(dllexport) uiint GetNumOfBoundaryNodeMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryFaceMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryEdgeMesh();
    __declspec(dllexport) uiint GetNumOfBoundaryVolumeMesh();
    // BoundaryType { Neumann || Dirichlet }
    __declspec(dllexport) uiint GetBNDType_BNodeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBNDType_BVolumeMesh(const uiint& ibmesh);
    // BoundaryTypeを表す型(定数)
    __declspec(dllexport) uiint getNeumannType();
    __declspec(dllexport) uiint getDirichletType();
    // 境界節点数
    __declspec(dllexport) uiint GetNumOfBNode_BNodeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfBNode_BVolumeMesh(const uiint& ibmesh);
    // BoundaryのDOF数
    __declspec(dllexport) uiint GetNumOfDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfDOF_BFaceMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfDOF_BEdgeMesh(const uiint& ibmesh);
    __declspec(dllexport) uiint GetNumOfDOF_BVolumeMesh(const uiint& ibmesh);
    // Boundary DOF番号("DOFインデックス"に対応するDOF番号)
    __declspec(dllexport) uiint GetDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BFaceMesh(const uiint& ibmesh, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BEdgeMesh(const uiint& ibmesh, const uiint& idof);
    __declspec(dllexport) uiint GetDOF_BVolumeMesh(const uiint& ibmesh, const uiint& idof);
    //--
    // BoundaryNode 境界値
    //--
    __declspec(dllexport) double& GetBNodeValue_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof);
    __declspec(dllexport) double& GetBNodeValue_BFaceMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) double& GetBNodeValue_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    __declspec(dllexport) uiint& GetNodeID_BNode_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode);
    __declspec(dllexport) uiint& GetNodeID_BNode_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode);
    //--
    // Face, Edge, Volume の境界値
    //--
    __declspec(dllexport) uiint GetNumOfBFace(const uiint& ibmesh);
    __declspec(dllexport) double& GetBFaceValue(const uiint& ibmesh, const uiint& ibface, const uiint& dof);
    __declspec(dllexport) uiint GetNumOfBEdge(const uiint& ibmesh);
    __declspec(dllexport) double& GetBEdgeValue(const uiint& ibmesh, const uiint& ibedge, const uiint& dof);
    __declspec(dllexport) uiint GetNumOfBVolume(const uiint& ibmesh);
    __declspec(dllexport) double& GetBVolumeValue(const uiint& ibmesh, const uiint& ibvol, const uiint& dof);
    //--
    // Face, Edge, Volume のNodeID
    //--
    __declspec(dllexport) uiint GetNumOfNode_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint& GetNodeID_BFace(const uiint& ibmesh, const uiint& ibface, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfNode_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint& GetNodeID_BEdge(const uiint& ibmesh, const uiint& ibedge, const uiint& ibnode);
    __declspec(dllexport) uiint GetNumOfNode_BVolume(const uiint& ibmesh, const uiint& ibvol);
    __declspec(dllexport) uiint& GetNodeID_BVolume(const uiint& ibmesh, const uiint& ibvol, const uiint& ibnode);
    //--
    // Boundaryの名称
    //--
    __declspec(dllexport) uiint GetBNodeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBNodeMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBFaceMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBFaceMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBVolumeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBVolumeMesh_Name(const uiint& ibmesh);
    __declspec(dllexport) uiint GetBEdgeMesh_NameLength(const uiint& ibmesh);
    __declspec(dllexport) string& GetBEdgeMesh_Name(const uiint& ibmesh);
    //--
    // Entity_ID of BoundaryMesh (for FrontISTR)
    //--
    __declspec(dllexport) uiint GetEdgeID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint GetElemID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    __declspec(dllexport) uiint GetFaceID_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint GetElemID_BFace(const uiint& ibmesh, const uiint& ibface);
    __declspec(dllexport) uiint GetElemID_BVolume(const uiint& ibmesh, const uiint& ibvol);



    //--
    // MPI (直接呼び出したい人向け)
    //--
    __declspec(dllexport) int& GetRank(); //自分のプロセス-ランクを取得
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
    // 以下の３メソッドは、ペア
    __declspec(dllexport) uiint GetNumOfNeibPE(const uiint& imesh);//Meshパーツが通信する相手の数
    __declspec(dllexport) uiint& GetTransRank(const uiint& imesh, const uiint& ipe);//通信Mesh毎のランク番号
    __declspec(dllexport) void Send_Recv_R(double* buf, const int& num_of_node, const int& dof_size, const int& trans_rank);//bufの値を送信-受信
    __declspec(dllexport) void Send_Recv_I(int* buf, const int& num_of_node, const int& dof_size, const int& trans_rank );
    ////    __declspec(dllexport) void Send_Recv_R(double* buf, const uiint& nDOF);// bufの値を送信, 受信値をNodeとbufに代入.   bufのサイズ == NumOfCommNode * dof_size
    ////    __declspec(dllexport) void Send_Recv_R();// 通信Nodeの値を入れ替えて更新


    //--
    // CommMesh2 (通信Mesh) { select された AssyModel,Meshを対象 }
    // CommNode (通信Node) :  Visualizer 用途
    //--
    __declspec(dllexport) uiint GetNumOfCommMesh();
    __declspec(dllexport) uiint GetNumOfCommNode(const uiint& icmesh);
    __declspec(dllexport) uiint& GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode);//MeshのNodeID


    //--
    // グループ { select された AssyModel,Meshを対象 }
    //--
    __declspec(dllexport) uiint GetNumOfElementGroup();
    __declspec(dllexport) uiint GetNumOfElementID(const uiint& iGrp);
    __declspec(dllexport) uiint& GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index);
    __declspec(dllexport) uiint GetElementGroupName_Length(const uiint& iGrp);
    __declspec(dllexport) string& GetElementGroupName(const uiint& iGrp);

    
    //--
    // Logger
    //--
    __declspec(dllexport) void LoggerMode(const uiint& mode);
    __declspec(dllexport) void LoggerDevice(const uiint& mode, const uiint& device);
    __declspec(dllexport) void LoggerInfo(const uiint& mode, char* message);
    __declspec(dllexport) void LoggerInfo(const uiint& mode, const char* message);
    //--
    // Logger_Parameter for Fortran
    //--
    __declspec(dllexport) uiint getErrorMode();
    __declspec(dllexport) uiint getWarnMode();
    __declspec(dllexport) uiint getInfoMode();
    __declspec(dllexport) uiint getDebugMode();
    __declspec(dllexport) uiint getDiskDevice();
    __declspec(dllexport) uiint getDisplayDevice();
};
#endif
}


