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

    //API 形状関数_種類カタログ
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
    uiint Initialize( int argc, char** argv);// 標準 [全体制御ファイルは使用しない]
    uiint Initialize_fstr( int argc, char** argv, string& ctrlname);//Init(標準), ctrlファイル名
    uiint Finalize();  // MPI_Finalize, Loggerのログファイルclose
    //----
    // 横断幕
    //----
    void  Banner_Display();

    //----
    // REVOCAP_Refiner
    //----
    void set_RevocapMesh(const uiint& nRefine);//REVOCAP_Refinerへのメッシュの登録

    
    //----
    // File関連
    //----
    uiint FileRead(string& basename, bool bBinary);// メッシュ(*.msh) 読み込み:標準
    uiint FileRead_fstr(bool bBinary);             // メッシュ(*.msh) 読み込み:ctrlファイル使用
    uiint FileDebugWrite();          // データ・チェック(*.out)
    uiint FileWriteRes(const uiint& nStep, bool bBinary);// Resファイル出力(*.res)
    // --
    // Restart data のセット(ファイル読み込み => 解ベクトルにセット)
    // --
    uiint SetRestart(const uiint& nStep, bool bBinary);// Resファイル => Algebra生成
    // --
    // リザルト・ファイル出力
    // --
    void PrintRlt_Start(const uiint& nStep, bool bBinary);
    void PrintRlt_P(const uiint& width, const char* format, ... );//リザルト出力(可変長引数:ポインター)
    void PrintRlt_R(const uiint& width, const char* format, ... );//リザルト出力(可変長引数:値)
    void PrintRlt_End();
    // *.inp 出力 
    void PrintMicroAVS_Basis();// 基礎変数の出力
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);//出力ラベルの登録
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);//出力変数の登録
    void PrintMicroAVS_FEM();  // 登録変数の出力
    
    // --
    // ファイル名(fstr)
    // --
    string getFstr_FileName_Mesh();
    string getFstr_FileName_Control();
    string getFstr_FileName_Result();
    string getFstr_FileName_Restart();
    string getFstr_FileName_PartIn();
    string getFstr_FileName_PartOut();
    string getFstr_FileName_VisMesh();
    string getFstr_FileName_VisIn();
    string getFstr_FileName_VisOut();


    //----
    // Solver
    //----
    void GeneLinearAlgebra(const uiint& nNumOfAlgebra, uiint* vDOF);//全Levelに方程式を生成 vDOF:方程式ごとのDOF
    void SelectAlgebra(const uiint& iequ);// mpAssyMatrix等にiequ番目の方程式をロード{ Levelは選択されているとする }

    uiint Matrix_Add_Elem(const uiint& iMesh, const uiint& iElem, double *ElemMatrix);
    uiint Matrix_Add_Node(const uiint& iMesh, const uiint& inode, const uiint& jnode, double *NodalMatrix);// inode,jnode:節点インデックス

    void Matrix_Clear(const uiint& iMesh);// Matrix all 0 clear
    void Vector_Clear(const uiint& iMesh);// Matrix all 0 clear
    
    uiint Set_BC_Mat_RHS2(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& sol_value);//対角項=Val,非対角項=0、右辺ベクトル
    uiint Set_BC_Mat_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& diag_value, double& rhs_value); //行列対角項, 右辺ベクトル
    uiint Set_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                        //右辺ベクトル
    uiint Add_BC_RHS(uiint& iMesh, uiint& iNodeID, uiint& nDOF, double& value);                        //右辺ベクトルへの加算
    
    void  Sample_Set_BC(uiint iMesh);
    
    uiint Solve(iint& iter_max, double& tolerance, iint& method, iint& precondition);

    //--
    // 解ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    void GetSolution_Vector(double* buf, const uiint& imesh);
    void GetSolution_AssyVector(double* buf);
    //--
    // 右辺ベクトルをbufへコピー (SelectされているLevel && Selectされている方程式番号)
    //--
    void GetRHS_Vector(double* buf, const uiint& imesh);
    void GetRHS_AssyVector(double* buf);
    //--
    // ベクトル値を取得
    //--
    double& GetSolutionVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    double& GetRHSVector_Val(const uiint& imesh, const uiint& inode, const uiint& idof);
    //--
    // ベクトルDOFを取得
    //--
    uiint& GetSolutionVector_DOF();
    uiint& GetRHSVector_DOF();
    //--
    // AssyMatrix * vX = vB , vector_size == NumOfMesh * NumOfNode * DOF
    //--
    void multVector(double* vX, double* vB);

    
    //----
    // MG constructor (Refiner)
    //----
    uiint  Refine(const uiint& nNumOfRefine);//ファイル指定のマルチグリッドMeshデータ生成


    
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
    uiint GetNumOfAssembleModel();//Assemble Modelの個数==階層数(mMGLevel+1)
    void SelectAssembleModel(const uiint& mgLevel);
    // Meshの選択
    uiint GetNumOfMeshPart();//Meshパーツの個数
    void SelectMeshPart_ID(const uiint& mesh_id);
    void SelectMeshPart_IX(const uiint& index);
    // 要素の選択
    void SelectElement_ID(const uiint& elem_id);
    void SelectElement_IX(const uiint& index);
    // 選択されたMeshの要素の情報
    uiint GetElementType();
    uiint GetNumOfElementVert();//要素の頂点数
    void GetElementVertNodeID(iint* vNodeID);//要素の頂点のノードID
    uiint GetNumOfElementEdge();//要素の辺数
    void GetElementEdgeNodeID(iint* vNodeID);//要素の辺のノードID
    // 選択されたMeshのノード(Scalar or Vector)
    void GetNodeCoord(const uiint& node_id, double& x, double& y, double& z);//ノードの座標
    uiint GetNumOfDOF(const uiint& node_id);
    uiint GetNumOfScalar(const uiint& node_id);
    uiint GetNumOfVector(const uiint& node_id);
    uiint& GetNodeType(const uiint& node_id);
////    void SetNodeValue(const uiint& node_id, double value[]);
////    void SetNodeValue(const uiint& node_id, const uiint& idof, const double& value);
////    void GetNodeValue(const uiint& node_id, double value[]);
////    double& GetNodeValue(const uiint& node_id, const uiint& idof);
    // 選択されたMeshのノード(ScalarVector)
////    void SetSVNodeValue(const uiint& node_id, double v_value[], double s_value[]);
////    void SetSVNodeValue(const uiint& node_id, const uiint& v_dof, const double& v_value, const uiint& s_dof, const double& s_value);
////    void GetSVNodeValue(const uiint& node_id, double v_value[], double s_value[]);
////    void GetSVNodeValue(const uiint& node_id, const uiint& v_dof, double& v_value, const uiint& s_dof, double& s_value);
    // 選択されたMeshのノード数、要素数
    uiint  getNodeSize(){ return mpMesh->getNodeSize();}
    uiint  getElementSize(){ return mpMesh->getElementSize();}
    uiint  getNodeSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    uiint  getElementSize(uiint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}
    // ID && index
    uiint& getNodeID(const uiint& index);
    uiint& getElementID(const uiint& index);
    uiint& getNodeIndex(const uiint& id);   //ID -> index
    uiint& getElementIndex(const uiint& id);//ID -> index
    //--
    // 節点コネクティビティ itemU:上三角のNodeID, itemL:下三角のNodeID
    //--
    void constructNodeConnectFEM(const uiint& node_id);
    void getNodeConnectFEM_Size(uiint& nNumOfItemU, uiint& nNumOfItemL);//node_id : 上三角のNodeID数, 下三角のNodeID数
    void getNodeConnectFEM_Item(uiint itemU[], uiint itemL[]);//node_id : 上三角のNodeID配列, 下三角のNodeID配列
    void getNodeConnectFEM_Item_F(iint itemU[], iint itemL[]);//FortranAPI用 (引数タイプ : signed int )
    //--
    // 節点周囲 要素
    //--
    uiint getNumOfAggregateElement(const uiint& node_id);
    uiint& getAggregateElementID(const uiint& node_id, const uiint& ielem);

    

    //----
    // メモリーの解放(Node::AggElemIDの解放):Mesh操作によるデータ構築の最後に実行
    //----
    void FinalizeRefine();


    //----
    // 節点タイプ for Fortran
    //----
    uiint nodetype_s();
    uiint nodetype_v();
    uiint nodetype_sv();
    //----
    // 要素タイプ for Fortran
    //----
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
    //----
    // FrontISTR 要素タイプ
    //----
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
    //----
    // FrontISTR 要素タイプ　=> MW3 要素タイプ 変換
    //----
    uiint fistr_elemtype_to_mw3_elemtype(const uiint& fistr_elemtype);
    //----
    // MW3 要素タイプ　=> FrontISTR 要素タイプ 変換
    //----
    uiint mw3_elemtype_to_fistr_elemtype(const uiint& mw3_elemtype);


    // --
    // 積分点数を返す
    // --
    uiint& NumOfIntegPoint(const uiint& shapeType);
    
    //--
    // 形状関数
    //--
    // N(積分点ごと)
    void ShapeFunc_on_pt(const uiint& shapeType, const uiint& igauss, vdouble& N);//積分点での形状関数:N[igauss] を返す.
    void ShapeFunc_on_pt(uiint shapeType, uiint igauss, double N[]);// Fortran API
    // N(丸ごと)
    void ShapeFunc(const uiint& shapeType, vvdouble& N);//全積分点の形状関数:N を丸ごと返す.
    // ----
    // Fortran API用 形状関数
    // ----
    // N[igauss][ishape]
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

    //--
    // 自然座標での導関数
    //--
    // dN/dr(積分点ごと)
    void dNdr_on_pt(const uiint& shapeType, const uiint& igauss, vvdouble& dNdr);
    // dN/dr(まるごと)
    void dNdr(const uiint& shapeType, vvvdouble& dNdr);
    // ----
    // Fortran API用 導関数
    // ----
    // dNdr[igauss][ishape][iaxis]
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


    //--
    // 空間座標での導関数
    //--
    // dN/dx 計算のみ(クラス・メンバー:mvdNdxに代入)
    void Calculate_dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index);
    // dN/dx(積分点ごと): mvdNdx[igaus] ,Calculate_dNdx(...)実行後に使用
    void dNdx_on_pt(const uiint& igauss, vvdouble& dNdx);
    // dN/dx(まるごと)
    void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& elem_index, vvvdouble& dNdx);
    // dN/dx(1次元配列  Fortran用途
    void dNdx(const uiint& elemType, const uiint& numOfInteg, const uiint& ielem, double dNdx[]);

    //--
    // J行列式 |J|
    //--
    //これを呼び出す直前に使用した,dNdx計算の detJ の値
    //--
    void detJacobian(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& detJ);
    
    //--
    // Gauss積分点の重み:Weight
    //--
    void Weight(const uiint& elemType, const uiint& numOfInteg, const uiint& igauss, double& w);

    //----
    // 形状関数タイプ for Fortran
    //----
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
    // Boundary :: 各Meshが所有するBoundaryMeshから境界値を取得
    //--
    // BoundaryMesh数
    uiint GetNumOfBoundaryNodeMesh();
    uiint GetNumOfBoundaryFaceMesh();
    uiint GetNumOfBoundaryEdgeMesh();
    uiint GetNumOfBoundaryVolumeMesh();
    // BoundaryType { Neumann || Dirichlet }
    uiint GetBNDType_BNodeMesh(const uiint& ibmesh);
    uiint GetBNDType_BFaceMesh(const uiint& ibmesh);
    uiint GetBNDType_BEdgeMesh(const uiint& ibmesh);
    uiint GetBNDType_BVolumeMesh(const uiint& ibmesh);
    // BoundaryTypeを表す型(定数)
    uiint getNeumannType();
    uiint getDirichletType();
    // 境界節点数
    uiint GetNumOfBNode_BNodeMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BFaceMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BEdgeMesh(const uiint& ibmesh);
    uiint GetNumOfBNode_BVolumeMesh(const uiint& ibmesh);
    // BoundaryのDOF数
    uiint GetNumOfDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint GetNumOfDOF_BFaceMesh(const uiint& ibmesh);
    uiint GetNumOfDOF_BEdgeMesh(const uiint& ibmesh);
    uiint GetNumOfDOF_BVolumeMesh(const uiint& ibmesh);
    // Boundary DOF番号("DOFインデックス"に対応するDOF番号)
    uiint GetDOF_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& idof);
    uiint GetDOF_BFaceMesh(const uiint& ibmesh, const uiint& idof);
    uiint GetDOF_BEdgeMesh(const uiint& ibmesh, const uiint& idof);
    uiint GetDOF_BVolumeMesh(const uiint& ibmesh, const uiint& idof);


    //--
    // BoundaryNode 境界値
    //--
    double& GetBNodeValue_BNodeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof);
    double& GetBNodeValue_BFaceMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    double& GetBNodeValue_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    double& GetBNodeValue_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode, const uiint& dof, const uiint& mgLevel);
    uiint& GetNodeID_BNode_BNodeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BFaceMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BEdgeMesh(const uiint& ibmesh, const uiint& ibnode);
    uiint& GetNodeID_BNode_BVolumeMesh(const uiint& ibmesh, const uiint& ibnode);
    //--
    // Face, Edge, Volume の境界値
    //--
    uiint GetNumOfBFace(const uiint& ibmesh);
    double& GetBFaceValue(const uiint& ibmesh, const uiint& ibface, const uiint& dof);
    uiint GetNumOfBEdge(const uiint& ibmesh);
    double& GetBEdgeValue(const uiint& ibmesh, const uiint& ibedge, const uiint& dof);
    uiint GetNumOfBVolume(const uiint& ibmesh);
    double& GetBVolumeValue(const uiint& ibmesh, const uiint& ibvol, const uiint& dof);
    //--
    // Face, Edge, Volume のNodeID
    //--
    uiint GetNumOfNode_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint& GetNodeID_BFace(const uiint& ibmesh, const uiint& ibface, const uiint& ibnode);
    uiint GetNumOfNode_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint& GetNodeID_BEdge(const uiint& ibmesh, const uiint& ibedge, const uiint& ibnode);
    uiint GetNumOfNode_BVolume(const uiint& ibmesh, const uiint& ibvol);
    uiint& GetNodeID_BVolume(const uiint& ibmesh, const uiint& ibvol, const uiint& ibnode);
    //--
    // Boundaryの名称
    //--
    uiint GetBNodeMesh_NameLength(const uiint& ibmesh);
    string& GetBNodeMesh_Name(const uiint& ibmesh);
    uiint GetBFaceMesh_NameLength(const uiint& ibmesh);
    string& GetBFaceMesh_Name(const uiint& ibmesh);
    uiint GetBVolumeMesh_NameLength(const uiint& ibmesh);
    string& GetBVolumeMesh_Name(const uiint& ibmesh);
    uiint GetBEdgeMesh_NameLength(const uiint& ibmesh);
    string& GetBEdgeMesh_Name(const uiint& ibmesh);
    //--
    // Entity_ID of BoundaryMesh (for FrontISTR)
    //--
    uiint GetEdgeID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint GetElemID_BEdge(const uiint& ibmesh, const uiint& ibedge);
    uiint GetFaceID_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint GetElemID_BFace(const uiint& ibmesh, const uiint& ibface);
    uiint GetElemID_BVolume(const uiint& ibmesh, const uiint& ibvol);


    //--
    // MPI (直接呼び出したい人向け)
    //--
    int& GetRank();//自分のプロセス-ランクを取得
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
    // 以下の３メソッドは、ペア
    //
    uiint GetNumOfNeibPE(const uiint& imesh);//Meshパーツが通信する相手の数
    uiint& GetTransRank(const uiint& imesh, const uiint& ipe);//通信Mesh毎のランク番号
    void Send_Recv_R(double* buf, const int& num_of_node, const int& dof_size, const int& trans_rank);//bufの値を送信-受信
    void Send_Recv_I(int* buf, const int& num_of_node, const int& dof_size, const int& trans_rank );

    ////    void Send_Recv_R(double* buf, const uiint& nDOF);// bufの値を送信, 受信値をNodeとbufに代入.   bufのサイズ == NumOfCommNode * dof_size
    ////    void Send_Recv_R();// 通信Nodeの値を入れ替えて更新

    //--
    // CommMesh2 (通信Mesh) { select された AssyModel,Meshを対象 }
    // CommNode (通信Node) :  Visualizer 用途
    //--
    uiint GetNumOfCommMesh();
    uiint GetNumOfCommNode(const uiint& icmesh);
    uiint& GetNodeID_CommNode(const uiint& icmesh, const uiint& icnode);//MeshのNodeID


    //--
    // グループ { select された AssyModel,Meshを対象 }
    //--
    uiint GetNumOfElementGroup();
    uiint GetNumOfElementID(const uiint& iGrp);
    uiint& GetElementID_with_ElementGroup(const uiint& iGrp, const uiint& index);
    uiint GetElementGroupName_Length(const uiint& iGrp);
    string& GetElementGroupName(const uiint& iGrp);


    //--
    // Logger
    //--
    void LoggerMode(const uiint& mode);
    void LoggerDevice(const uiint& mode, const uiint& device);
    void LoggerInfo(const uiint& mode, char* message);
    void LoggerInfo(const uiint& mode, const char* message);
    //--
    // Logger_Parameter for Fortran
    //--
    uiint getErrorMode();
    uiint getWarnMode();
    uiint getInfoMode();
    uiint getDebugMode();
    uiint getDiskDevice();
    uiint getDisplayDevice();
};
#endif
}


