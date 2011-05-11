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
#include "TypeDef.h"
#include "Vertex.h"

#include <vector>
#include <iostream>

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
    uiint myRank; //ランク
    uiint mMeshID;//メッシュ・ノードが所属するMesh <-- 計算領域に接合面が無い場合は,不定
    uiint mNodeID;//メッシュ・ノードのID <-- 計算領域に接合面が無い場合は,不定
    uiint mLevel;//MultiGridレベル

    // 接合面のノードが計算領域に存在するか否か？
    // ・ ランクだけでは判別できないのでbool値を持たす
    //--
    bool mbMesh;//Meshが自身の計算領域データか.(true == 自分の計算量域内)
    bool mbNode;//Nodeが自身の計算領域データか.(true == 自分の計算量域内)

    // 変位:EQUATION計算でついでに計算 -> Valueのカップリング
    vdouble mvDisplacement;//接合面Nodeの変位:Arbitrary DOF (任意自由度)
    vdouble mvScalar;//接合面の温度等のスカラー量

    // MPCソルバーの従属行列からのアクセスの為
    uiint mMFaceID;//自身がSlave点だった場合の所属Master面ID => mMFaceID 没
    bool mbSlave; //自身がSlave点か.(true == スレーブ点)
    bool mbMarkingMFace;//Master面IDをセットしたことをマーキング



    // 各Levelでのマスター面ID '10.03.05
    // ---------------------
    map<uiint, uiint, less<uiint> > mmMasterFaceID;// [iLevel]= MasterFaceID
    vector<bool>                mvbMasterFaceID;// [iLevel]= true .or. false

    // 八分木
    // -----
    vector<uiint> mvKnotID;

public:
    void markingSelfMesh();
    void markingSelfNode();

    // Nodeを特定 <= Meshの節点
    // --
    void  setLevel(const uiint& level){ mLevel= level;}
    uiint& getLevel(){ return mLevel;}
    void pushLevelMarking(){  mvbMasterFaceID.push_back(false);}

    void  setRank(const uiint& rank){ myRank= rank;}
    uiint& getRank(){ return myRank;}
    void  setMeshID(const uiint& meshID){ mMeshID= meshID;}
    uiint& getMeshID(){ return mMeshID;}
    void  setNodeID(const uiint& nodeID){ mNodeID= nodeID;}
    uiint& getNodeID(){ return mNodeID;}


    // EQUATION関連
    // ---
    // 変位
    void resizeDisp(const uiint& dof);//変位の自由度配列の確保
    void initDisp();//変位の初期化
    void setDisp(const uiint& idof, const double& disp){ mvDisplacement[idof]= disp;}
    double& getDisp(const uiint& idof){ return mvDisplacement[idof];}
    uiint getNumOfDisp(){return mvDisplacement.size();}
    // スカラー
    void resizeScalar(const uiint& numOfScalar);//スカラー量のパラメータ数分の配列確保
    void initScalar();//スカラーの初期化
    void setScalar(const uiint& i, const double& val){ mvScalar[i]=val;}
    double& getScalar(const uiint& i){ return mvScalar[i];}
    uiint getNumOfScalar(){return mvScalar.size();}
    
    // ----
    // Slave点の場合の所属マスター面ID管理
    // ----
    //  コメントアウト '10.03.05
    //  ---------------------
    //  void setMasterFaceID(const uint& faceID){ mMFaceID=faceID; mbMarkingMFace=true;}
    //  uint& getMasterFaceID(){ return mMFaceID;}
    //  bool have_MasterFaceID(){ return mbMarkingMFace;}//MasterFaceIDを所有しているか？
    void markingSlave();
    bool isSlave(){ return mbSlave;}


    
    // 新バージョン '10.03.05
    // ----
    // 各Levelでのマスター面ID
    // ----
    void setMasterFaceID(const uiint& faceID, const uiint& level);
    uiint& getMasterFaceID(const uiint& level);
    bool have_MasterFaceID(const uiint& level);



    // 八分木
    // -----
    void  resizeOctreeID(const uiint& res_size);
    void  setOctreeID(const uiint& layer, const uiint& knot_id);
    uiint& getOctreeID(const uiint& layer){ return mvKnotID[layer];}

};
#endif	/* _CONTACTNODE_H */
}













