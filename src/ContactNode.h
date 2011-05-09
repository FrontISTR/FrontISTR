/* 
 * File:   ContactNode.h
 * Author: ktakeda
 *
 * MPC接合メッシュ(ContactMesh)のノード
 *
 * 親クラス:CVertex
 *
 * Created on 2009/10/26, 14:25
 */
#include "Vertex.h"

namespace pmw{
#ifndef _CONTACTNODE_H
#define	_CONTACTNODE_H
class CContactNode:public CVertex{
public:
    CContactNode();
    virtual ~CContactNode();

    // --
    // 自身のID:class CVertex
    // SkinFaceIDの節点集合関連:class CVertex
    // --

protected:
    // Nodeを特定するパラメータ
    // --
    uint myRank; //ランク
    uint mMeshID;//メッシュ・ノードが所属するMesh <-- 計算領域に接合面が無い場合は,不定
    uint mNodeID;//メッシュ・ノードのID <-- 計算領域に接合面が無い場合は,不定
    uint mLevel;//MultiGridレベル

    // 接合面のノードが計算領域に存在するか否か？
    // ・ ランクだけでは判別できないのでbool値を持たす
    //--
    bool mbMesh;//Meshが自身の計算領域データか.(true == 自分の計算量域内)
    bool mbNode;//Nodeが自身の計算領域データか.(true == 自分の計算量域内)

    // 変位:EQUATION計算の為
    vdouble mvDisplacement;//接合面Nodeの変位:Arbitrary DOF (任意自由度)
    vdouble mvScalar;//接合面の温度等のスカラー量
public:
    void markingSelfMesh(){ mbMesh=true;}
    void markingSelfNode(){ mbNode=true;}

    // Nodeを特定 <= Meshの節点
    // --
    void  setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
    void  setRank(const uint& rank){ myRank= rank;}
    uint& getRank(){ return myRank;}
    void  setMeshID(const uint& meshID){ mMeshID= meshID;}
    uint& getMeshID(){ return mMeshID;}
    void  setNodeID(const uint& nodeID){ mNodeID= nodeID;}
    uint& getNodeID(){ return mNodeID;}


    // EQUATION関連
    // ---
    // 変位
    void resizeDisp(const uint& dof);//変位の自由度配列の確保
    void initDisp();//変位の初期化
    void setDisp(const uint& idof, const double& disp){ mvDisplacement[idof]= disp;}
    double& getDisp(const uint& idof){ return mvDisplacement[idof];}
    uint getNumOfDisp(){return mvDisplacement.size();}
    // スカラー
    void resizeScalar(const uint& numOfScalar);//スカラー量のパラメータ数分の配列確保
    void initScalar();//スカラーの初期化
    void setScalar(const uint& i, const double& val){ mvScalar[i]=val;}
    double& getScalar(const uint& i){ return mvScalar[i];}
    uint getNumOfScalar(){return mvScalar.size();}
};
#endif	/* _CONTACTNODE_H */
}


