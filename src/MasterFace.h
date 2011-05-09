//
//  MasterFace.h
//
//  MPCマスター面群の面:親クラス==CSkinFace
//
//				2009.10.15
//				2009.01.08
//				k.Takeda
#ifndef MASTER_FACE_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC
#define MASTER_FACE_HH_CC30A703_CF33_4eaf_98E1_91A1310095BC


#include "SkinFace.h"
#include "MPCValueType.h"//CalcSlaveで補間するConNodeの値の種類指定

namespace pmw{
class CMasterFace:public CSkinFace{
public:
    CMasterFace();
    virtual ~CMasterFace();

private:
    double mDotValue;//dot product 参照変数
    vdouble mvVector;//cross product 参照変数 ::内外判定に用いるスレーブ点と頂点2点で構成される面ベクトル
    vdouble& CrossProduct(CContactNode* pConNode,const uint& ivert, const uint& jvert);//外積(ConNodeと頂点2点の外積)
    double& DotProduct();//内外判定用の内積(mvNormalVectorとmvVectorの内積):面ベクトルは,ContactMesh::setupSPointOnMFace()で計算がコールされている.
    double VectorLength(const vdouble& vec);
    double VectorLength(CContactNode *pConNode0, CContactNode *pConNode1);
    double VectorLength(CContactNode *pConNode, const vdouble& vec);
    double VectorLength(const vdouble& vBasePoint, const vdouble& vPoint);
    double VectorLength(const vdouble& vec, CContactNode *pConNode);

    // 線形補間 :BaseNode-ConNodeのidof番目の値の線形補間
    double LinearInter(const double& coeff, const uint& idof, CContactNode *pConBaseNode, CContactNode *pConNode ,const uint& valType);
    double LinearInter(const double& coeff, const uint& idof, const vdouble& vBaseVal, const vdouble& vVal);

    // 2直線ベクトル L1,L2　の最近接点
    // -----
    // L1とL2ベクトルの最近接点を表すパラメータ計算, 最近接点を表すスカラー値s,t
    // 
    double mS,mT, mParam_A, mParam_B, mParam_C, mParam_D, mParam_E, mParam_F;
    vdouble mvParam_R;

    double& s_NearCrossP();// L1ベクトルの基点からの距離 s
    double& t_NearCrossP();// L2ベクトルの基点からの距離 t
    // 線形代数方程式の各パラメータ
    bool NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, const vdouble& L2P0);//L1とL2の最近接点を表す行列要素a,b,c,d,e,fの計算
    bool NearCrossP(const vdouble& vL1, const vdouble& vL2, CContactNode* pL1P0, CContactNode* pL2P0);//L1とL2の最近接点を表す行列要素a,b,c,d,e,fの計算
    bool NearCrossP(const vdouble& vL1, const vdouble& vL2, const vdouble& L1P0, CContactNode* pL2P0);//L1とL2の最近接点を表す行列要素a,b,c,d,e,fの計算
    // 各パラメータの個別計算関数
    double& a_NearCrossP(const vdouble& vL1);// パラメータa
    double& b_NearCrossP(const vdouble& vL1, const vdouble& vL2);// パラメータb
    double& c_NearCrossP(const vdouble& vL1);// パラメータc
    double& d_NearCrossP();// パラメータd
    double& e_NearCrossP(const vdouble& vL2);// パラメータe
    double& f_NearCrossP(const vdouble& vL2);// パラメータf
    vdouble& r_NearCrossP(const vdouble& L1P0, const vdouble& L2P0);// パラメータr
    vdouble& r_NearCrossP(CContactNode *pL1P0, CContactNode *pL2P0);// パラメータr
    vdouble& r_NearCrossP(const vdouble& L1P0, CContactNode *pL2P0);// パラメータr


    // 点-点とベクトルvVの内積(投影した距離) => 現状,使っていない.
    double& s_ProjecP(const vdouble& vV, CContactNode *pP,CContactNode *pP3);//点P-P3ベクトルからベクトルVへの投影(内積),最近接点パラメータs
    double& t_ProjecP(const vdouble& vV, CContactNode *pP,CContactNode *pP2);//点P-P2ベクトル

    // ベクトルの正規化
    void Normalized(vdouble& vV);
    //    vdouble mNzVector;//正規化ベクトル計算時の一時的な入れ物(計算用テンポラリー変数)
    //    vdouble& Normalized(CContactNode* pConBaseNode, CContactNode* pConNode);

protected:
    // Slaveノード
    vector<CContactNode*> mvSlaveConNode;

    virtual CSkinFace* generateFace();//自身の生成,SkinFace::refine()で使用
public:
    virtual const char* getName(){ return "MasterFace";}

    // Slaveノード : ContactMeshで,Slaveノードをセット
    // --
    virtual void addSlaveNode(CContactNode* pConNode);//内外判定を行い,内ならば追加
    virtual CContactNode* getSlaveNode(const uint& index){ return mvSlaveConNode[index];}
    virtual uint getNumOfSlaveNode(){ return mvSlaveConNode.size();}
    // EQUATION :: valType => MPCValueType
    virtual void CalcSlave(const uint& islave, const uint& valType);//マスター面の変位 => スレーブ点の変位/スカラー値, 変位/スカラーはContactNodeが所有.
    // ----
    // virtual void CalcSlaveDisplacement(const uint& islave);//マスター面の変位 => スレーブ点の変位, 変位はContactNodeが所有.
    // virtual void CalcSlaveScalar(const uint& islave);//EQUATION::マスター面の温度 => スレーブ点の温度, 温度はContactNodeが所有.
};
}
#endif

