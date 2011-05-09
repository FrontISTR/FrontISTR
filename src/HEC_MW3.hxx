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
#ifndef PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513
#define PMW_MAIN_HH_C8955A70_0EE9_4f3e_82D1_FB32C9626513

//typedef文
#include "TypeDef.h"

//MPIラッパー
#include "HEC_MPI.h"

//Utility
#include "FileIO.h"
#include "Logger.h"

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
    __declspec(dllexport) CMW(void);
public:
    __declspec(dllexport) virtual ~CMW(void);

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
    CAssyVector *mpAssyVector;
    CAssyVector *mpAssyVector2;
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
    __declspec(dllexport) int FileRead(); // ファイル入力(MW3書式):AssyModel階層の構築
    __declspec(dllexport) int FileWrite();// ファイル出力

    //----
    // Solver
    //----
    __declspec(dllexport) int Initialize_Matrix(); // 行列の初期化
    __declspec(dllexport) int Initialize_Vector(); // 行列の初期化
    __declspec(dllexport) int Matrix_Add_Elem(int iMesh, int iElem, double *ElemMatrix);
    __declspec(dllexport) int Set_BC(int iMesh, int iNode, int iDof, double value1, double value2);
    __declspec(dllexport) int Set_BC(int iMesh, int iNode, int iDof, double value);
	__declspec(dllexport) void Sample_Set_BC(int iMesh);
    __declspec(dllexport) int Solve(uint iter_max, double tolerance, uint method, uint precondition); 

    //----
    // Refiner
    //----
    __declspec(dllexport) int Refine();//ファイル指定のマルチグリッドMeshデータ生成

    //----
    // FEM : 基本void(Fortran_API に合わせる為) : 引数のconst以外=>出力用
    //----
    // 作業用変数へ指定オブジェクトを選択:形状関数,導関数のループ処理用
    __declspec(dllexport) void SelectAssyModel(const uint& mgLevel);
    __declspec(dllexport) void SelectMesh_ID(const uint& mesh_id);
    __declspec(dllexport) void SelectMesh_IX(const uint& index);
//    void SelectElement_ID(const uint& elem_id);
//    void SelectElement_IX(const uint& index);
    __declspec(dllexport) uint  getNodeSize(){ return mpMesh->getNodeSize();}
    __declspec(dllexport) uint  getElementSize(){ return mpMesh->getElementSize();}
    __declspec(dllexport) uint  getNodeSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getNodeSize();}
    __declspec(dllexport) uint  getElementSize(uint iMesh){ return mpAssy->getMesh(iMesh)->getElementSize();}

	__declspec(dllexport) void StoreMatrix(){
		mpAssy->setAssyMatrix(mpAssyMatrix);
		mpAssy->setAssyVector(mpAssyVector);
		mpAssy->setAssyVector2(mpAssyVector2);
	}
	
	__declspec(dllexport) void LoadMatrix(){
		mpAssyMatrix = mpAssy->getAssyMatrix();
		mpAssyVector = mpAssy->getAssyVector();
		mpAssyVector2 = mpAssy->getAssyVector2();
	}

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
    __declspec(dllexport) void GetNumOfAssembleModel(uint& numOfAssembleModel);//Assemble Modelの個数==階層数(mMGLevel+1)
    __declspec(dllexport) void SelectAssembleModel(const uint& mgLevel);
    // Meshの選択
    __declspec(dllexport) void GetNumOfMeshPart(uint& numOfMeshPart);//Meshパーツの個数
    __declspec(dllexport) void SelectMeshPart_ID(const uint& mesh_id);
    __declspec(dllexport) void SelectMeshPart_IX(const uint& index);
    // 要素の選択
    __declspec(dllexport) void SelectElement_ID(const uint& elem_id);
    __declspec(dllexport) void SelectElement_IX(const uint& index);
    // 選択されたMeshの要素の情報
    __declspec(dllexport) void GetElementType(uint& elemType);
    __declspec(dllexport) void GetNumOfElementVert(uint& numOfVert);//要素の頂点数
    __declspec(dllexport) void GetElementVertNodeID(vuint& vNodeID);//要素の頂点のノードID
    __declspec(dllexport) void GetNumOfElementEdge(uint& numOfEdge);//要素の辺数
    __declspec(dllexport) void GetElementEdgeNodeID(vuint& vNodeID);//要素の辺のノードID
    // 選択されたMeshのノードの座標
    __declspec(dllexport) void GetNodeCoord(const uint& node_id, double& x, double& y, double& z);//ノードの座標

    // --
    // 積分点数を返す
    // --
    __declspec(dllexport) void NumOfIntegPoint(const uint& shapeType, uint& numOfInteg);
    
    //--
    // 形状関数
    //--
    // N(積分点ごと)
    __declspec(dllexport) void ShapeFunc_on_pt(const uint& shapeType, const uint& igauss, vdouble& N);//積分点での形状関数:N[igauss] を返す.
    // N(まるごと)
    __declspec(dllexport) void ShapeFunc(const uint& shapeType, vvdouble& N);//全積分点の形状関数:N を丸ごと返す.
    
    //--
    // 自然座標での導関数
    //--
    // dN/dr(積分点ごと)
    __declspec(dllexport) void dNdr_on_pt(const uint& shapeType, const uint& igauss, vvdouble& dNdr);
    // dN/dr(まるごと)
    __declspec(dllexport) void dNdr(const uint& shapeType, vvvdouble& dNdr);



    //--
    // 空間座標での導関数
    //--
    // dN/dx 計算のみ(クラス・メンバー:mvdNdxに代入)
    __declspec(dllexport) void Calculate_dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index);
    // dN/dx(積分点ごと): mvdNdx[igaus] ,Calculate_dNdx(...)実行後に使用
    __declspec(dllexport) void dNdx_on_pt(const uint& igauss, vvdouble& dNdx);
    // dN/dx(まるごと)
    __declspec(dllexport) void dNdx(const uint& elemType, const uint& numOfInteg, const uint& elem_index, vvvdouble& dNdx);

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
};
}

#endif
