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
    uiint mID;
    uiint myRank;
    uiint mTransmitRank;

    uiint mMGLevel;

    vector<CCommNode*> mvCommNode;    //頂点ノード集合(1次)、頂点と辺ノード集合(2次)
    //vector<CCommNode*> mvEdgeCommNode;//辺ノード集合
    //vector<CCommNode*> mvFaceCommNode;//面ノード集合
    uiint mEdgeCount;//辺に生成するノード数のカウンター

    vector<CCommFace*> mvCommFace;    //界面(Meshの幾何学的界面)

    map<uiint, uiint, less<uiint> > mmCommNodeID2Index;
    map<uiint, uiint, less<uiint> > mmCommFaceID2Index;

public:
    void  setID(const uiint& id){ mID= id;}
    uiint& getID(){ return mID;}

    void setLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getLevel(){ return mMGLevel;}

    // Mesh通信界面のランク
    void  setRank(const uiint& rank){ myRank= rank;}
    void  setTransmitRank(const uiint& trans_rank){ mTransmitRank= trans_rank;}
    uiint& getRank(){ return myRank;}
    uiint& getTrasmitRank(){ return mTransmitRank;}

public:
    // 通信界面 節点・面の追加
    void reserveCommNode(const uiint& res_size){ mvCommNode.reserve(res_size);}
    void addCommNode(CCommNode* pCommNode);

    void reserveCommFace(const uiint& res_size){ mvCommFace.reserve(res_size);}
    void addCommFace(CCommFace* pCommFace);

    // 通信界面 節点 - 面の取得
    //CCommNode* getCommVertNode(const uint& id);
    CCommNode* getCommNode(const uiint& id);
    //CCommNode* getCommVertNodeIX(const uint& index){ return mvCommNode[index];}
    CCommNode* getCommNodeIX(const uiint& index){ return mvCommNode[index];}
    //CCommNode* getCommEdgeNodeIX(const uint& index){ return mvEdgeCommNode[index];}
    //CCommNode* getCommFaceNodeIX(const uint& index){ return mvFaceCommNode[index];}

    CCommFace* getCommFace(const uiint& id);
    CCommFace* getCommFaceIX(const uiint& index){ return mvCommFace[index];}

    // サイズ
    //uint getCommVertNodeSize(){ return mvCommNode.size();}
    uiint getCommNodeSize(){ return mvCommNode.size();}
    uiint getCommFaceSize(){ return mvCommFace.size();}

public:
    void setupAggFace();
    void setupEdgeCommNode(CCommMesh2* pProgCommMesh, const uiint& nLevel);
    void setupFaceCommNode(CCommMesh2* pProgCommMesh);

    
    //void setupVertCommNode(CCommMesh2* pProgCommMesh);
    void setupCommNode(CCommMesh2* pProgCommMesh);// 下位のノードを,そのまま上位にセット

    //Refine後のvector解放(辺ノード)
    void deleteProgData();
};
#endif	/* _COMM_MESH2_H */
}



