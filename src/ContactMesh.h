//
//  ContactMesh.h
//
//  MPC接続面メッシュ
//				2009.10.15
//				2009.01.08
//				k.Takeda
#ifndef CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F
#define CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F


#include "SkinFace.h"
#include "MasterFace.h"
#include "BoundingBox.h"//スレーブ点の絞り込み検索Box
#include "MPCValueType.h"//MPC面で補間する値の種類

#include <map>
#include <utility>//pair
typedef std::pair<pmw::CContactNode*, pmw::CContactNode*> PairConNode;

namespace pmw{
class CContactMesh{
public:
    CContactMesh();
    virtual ~CContactMesh();

protected:
    uint mContactID;//MPC接続面ID
    uint mLevel;//MultiGridレベル
    //--
    //マスター スレーブのMeshIDは,Refine時に新Faceを取得する際に必要になる.
    //--
    vuint mvMasterMeshID;//複数のMeshの要素面を集めて,マスター面を構成しているのでMeshIDは複数存在
    vuint mvSlaveMeshID; //上に同じ
    map<uint, uint, less<uint> >  mmMasterMeshID2Index;
    map<uint, uint, less<uint> >  mmSlaveMeshID2Index;

    // 接合メッシュ全体のノード
    vector<CContactNode*> mvConNode;//接合メッシュ全体のContactNode配列
    map<uint, uint, less<uint> >  mmConNodeID2Index;

    // Master
    vector<CContactNode*> mvMasterConNode;// マスターNode
    vector<CSkinFace*> mvFace;// マスター面(インスタンスは,MasterFace)
    map<uint, uint, less<uint> > mmMasterConNodeID2Index;

    // Slave
    vector<CContactNode*> mvSlaveConNode; // スレーブNode
    vector<CSkinFace*> mvSlaveFace;// スレーブ面(Refine時に必要)
    map<uint, uint, less<uint> > mmSlaveConNodeID2Index;

    // BoundingBox(検索Box)
    CBoundingBox mBoundingBox;

public:
    void setID(const uint& id){ mContactID = id;}
    uint& getID(){ return mContactID;}

    void setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}

    // 接合メッシュ全体のContactNode
    // --
    void reserveConNode(const uint& res_size){ mvConNode.reserve(res_size);}
    void resizeConNode(const uint& res_size){ mvConNode.resize(res_size);}
    void addConNode(CContactNode* pConNode, const uint& id);
    CContactNode* getContactNode(const uint& index){ return mvConNode[index];}
    CContactNode* getContactNode_ID(const uint& id){ uint i=mmConNodeID2Index[id]; return mvConNode[i];}
    uint getNumOfConNode(){ return mvConNode.size();}

    // 接合しているMesh
    // --
    // マスターメッシュIDs
    void addMasterMeshID(const uint& id);
    uint& getMasterMeshID(const uint& index){return mvMasterMeshID[index];}
    uint getNumOfMasterMesh(){ return mvMasterMeshID.size();}
    // スレーブメッシュIDs
    void addSlaveMeshID(const uint& id);
    uint& getSlaveMeshID(const uint& index){ return mvSlaveMeshID[index];}
    uint getNumOfSlaveMesh(){ return mvSlaveMeshID.size();}


    // マスターFace
    // --
    void reserveMasterFace(const uint& res_size){ mvFace.reserve(res_size);}
    void resizeMasterFace(const uint& res_size){ mvFace.resize(res_size);}
    void addMasterFace(CSkinFace* pFace){ mvFace.push_back(pFace);}
    void addMasterFace(vector<CSkinFace*>& vface);
    void setMasterFace(CSkinFace* pFace, const uint& index){ mvFace[index] = pFace;}
    uint getNumOfMasterFace(){ return mvFace.size();}
    CSkinFace* getMasterFace(const uint& index){ return mvFace[index];}
    // マスターConNode
    void addMasterConNode(CContactNode* pConNode, const uint& id);
    void reserveMasterConNode(const uint& res_size){ mvMasterConNode.reserve(res_size);}
    void resizeMasterConNode(const uint& res_size){ mvMasterConNode.resize(res_size);}
    CContactNode* getMasterConNode(const uint& index){ return mvMasterConNode[index];}
    CContactNode* getMasterConNode_ID(const uint& id){ uint i=mmMasterConNodeID2Index[id]; return mvMasterConNode[i];}
    uint getNumOfMasterConNode(){ return mvMasterConNode.size();}

    // スレーブFace
    // --
    void reserveSlaveFace(const uint& res_size){ mvSlaveFace.reserve(res_size);}
    void resizeSlaveFace(const uint& res_size){ mvSlaveFace.resize(res_size);}
    void addSlaveFace(CSkinFace* pSlaveFace){ mvSlaveFace.push_back(pSlaveFace);}
    void addSlaveFace(vector<CSkinFace*>& vface);
    void setSlaveFace(CSkinFace* pSlaveFace, const uint& index){ mvSlaveFace[index]= pSlaveFace;}
    uint getNumOfSlaveFace(){ mvSlaveFace.size();}
    CSkinFace* getSlaveFace(const uint& index){ return mvSlaveFace[index];}
    // スレーブConNode
    void addSlaveConNode(CContactNode* pConNode, const uint& id);
    void reserveSlaveConNode(const uint& res_size){ mvSlaveConNode.reserve(res_size);}
    void resizeSlaveConNode(const uint& res_size){ mvSlaveConNode.resize(res_size);}
    CContactNode* getSlaveConNode(const uint& index){ return mvSlaveConNode[index];}
    CContactNode* getSlaveConNode_ID(const uint& id){ uint i=mmSlaveConNodeID2Index[id]; return mvSlaveConNode[i];}
    uint getNumOfSlaveConNode(){ return mvSlaveConNode.size();}


    // Refine関連
    // --
    void setupAggSkinFace();// ConNodeへ接続しているSkinFaceIDの集合をConNodeへセット.
    void setupEdgeConNode(CContactMesh *pProgConMesh, uint& countID);// 辺に接続するSkinFaceのセットと辺ノードの生成,設置
    void setupFaceConNode(CContactMesh *pProgConMesh, uint& countID);// 面中心に配置する面ノードの生成,設置


    // マスター&スレーブ関連  :: Equationは,MasterFaceにつける予定
    // --
    void setupSPointOnMFace();// マスター面にスレーブ点を設置(BBoxのみ利用バージョン)
    
};
}
#endif
