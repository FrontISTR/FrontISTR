/* 
 * File:   CommMesh2.h
 * Author: ktakeda
 *
 *  共有節点による通信Mesh
 *   => 通信界面{Solid(面),Shell(辺),Beam(点)}メッシュモデル
 *
 * Created on 2010/03/02, 17:18
 */
#include "TypeDef.h"
#include <map>

#include "CommNode.h"
#include "CommFace.h"

namespace pmw{
#ifndef _COMM_MESH2_H
#define	_COMM_MESH2_H
class CCommMesh2{
public:
    CCommMesh2();
    virtual ~CCommMesh2();

protected:
    uint mID;
    uint myRank;
    uint mTransmitRank;

    uint mMGLevel;

    vector<CCommNode*> mvCommNode;    //頂点ノード集合(1次)、頂点と辺ノード集合(2次)
    //vector<CCommNode*> mvEdgeCommNode;//辺ノード集合
    //vector<CCommNode*> mvFaceCommNode;//面ノード集合
    uint mEdgeCount;//辺に生成するノード数のカウンター

    vector<CCommFace*> mvCommFace;    //界面(Meshの幾何学的界面)

    map<uint, uint, less<uint> > mmCommNodeID2Index;
    map<uint, uint, less<uint> > mmCommFaceID2Index;

public:
    void  setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}

    void setLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getLevel(){ return mMGLevel;}

    // Mesh通信界面のランク
    void  setRank(const uint& rank){ myRank= rank;}
    void  setTransmitRank(const uint& trans_rank){ mTransmitRank= trans_rank;}
    uint& getRank(){ return myRank;}
    uint& getTrasmitRank(){ return mTransmitRank;}

public:
    // 通信界面 節点・面の追加
    void reserveCommNode(const uint& res_size){ mvCommNode.reserve(res_size);}
    void addCommNode(CCommNode* pCommNode);

    void reserveCommFace(const uint& res_size){ mvCommFace.reserve(res_size);}
    void addCommFace(CCommFace* pCommFace);

    // 通信界面 節点 - 面の取得
    //CCommNode* getCommVertNode(const uint& id);
    CCommNode* getCommNode(const uint& id);
    //CCommNode* getCommVertNodeIX(const uint& index){ return mvCommNode[index];}
    CCommNode* getCommNodeIX(const uint& index){ return mvCommNode[index];}
    //CCommNode* getCommEdgeNodeIX(const uint& index){ return mvEdgeCommNode[index];}
    //CCommNode* getCommFaceNodeIX(const uint& index){ return mvFaceCommNode[index];}

    CCommFace* getCommFace(const uint& id);
    CCommFace* getCommFaceIX(const uint& index){ return mvCommFace[index];}

    // サイズ
    //uint getCommVertNodeSize(){ return mvCommNode.size();}
    uint getCommNodeSize(){ return mvCommNode.size();}
    uint getCommFaceSize(){ return mvCommFace.size();}

public:
    void setupAggFace();
    void setupEdgeCommNode(CCommMesh2* pProgCommMesh, const uint& nLevel);
    void setupFaceCommNode(CCommMesh2* pProgCommMesh);

    
    //void setupVertCommNode(CCommMesh2* pProgCommMesh);
    void setupCommNode(CCommMesh2* pProgCommMesh);// 下位のノードを,そのまま上位にセット

    //Refine後のvector解放(辺ノード)
    void deleteProgData();
};
#endif	/* _COMM_MESH2_H */
}



