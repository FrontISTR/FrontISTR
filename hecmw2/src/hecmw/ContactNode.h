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
    uint myRank; //ランク
    uint mMeshID;//メッシュ・ノードが所属するMesh <-- 計算領域に接合面が無い場合は,不定
    uint mNodeID;//メッシュ・ノードのID <-- 計算領域に接合面が無い場合は,不定
    uint mLevel;//MultiGridレベル

    // 接合面のノードが計算領域に存在するか否か？
    // ・ ランクだけでは判別できないのでbool値を持たす
    //--
    bool mbMesh;//Meshが自身の計算領域データか.(true == 自分の計算量域内)
    bool mbNode;//Nodeが自身の計算領域データか.(true == 自分の計算量域内)

    // 変位:EQUATION計算でついでに計算 -> Valueのカップリング
    vdouble mvDisplacement;//接合面Nodeの変位:Arbitrary DOF (任意自由度)
    vdouble mvScalar;//接合面の温度等のスカラー量

    // MPCソルバーの従属行列からのアクセスの為
    uint mMFaceID;//自身がSlave点だった場合の所属Master面ID => mMFaceID 没
    bool mbSlave; //自身がSlave点か.(true == スレーブ点)
    bool mbMarkingMFace;//Master面IDをセットしたことをマーキング



    // 各Levelでのマスター面ID '10.03.05
    // ---------------------
    map<uint, uint, less<uint> > mmMasterFaceID;// [iLevel]= MasterFaceID
    vector<bool>                mvbMasterFaceID;// [iLevel]= true .or. false

    // 八分木
    // -----
    vector<uint> mvKnotID;

public:
    void markingSelfMesh();
    void markingSelfNode();

    // Nodeを特定 <= Meshの節点
    // --
    void  setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
    void pushLevelMarking(){  mvbMasterFaceID.push_back(false);}

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
    void setMasterFaceID(const uint& faceID, const uint& level);
    uint& getMasterFaceID(const uint& level);
    bool have_MasterFaceID(const uint& level);



    // 八分木
    // -----
    void  resizeOctreeID(const uint& res_size);
    void  setOctreeID(const uint& layer, const uint& knot_id);
    uint& getOctreeID(const uint& layer){ return mvKnotID[layer];}

};
#endif	/* _CONTACTNODE_H */
}













