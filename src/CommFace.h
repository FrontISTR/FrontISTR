/* 
 * File:   CommFace.h
 * Author: ktakeda
 *
 * Created on 2010/03/02, 17:28
 */
#include "TypeDef.h"
#include <utility>//pair
#include "EdgeTree.h"

#include "CommNode.h"
#include "Element.h"

#include "Logger.h"

namespace pmw{
typedef pair<CCommNode*,CCommNode*> PairCommNode;
#ifndef _COMMFACE_H
#define	_COMMFACE_H
class CCommFace{
public:
    CCommFace();
    virtual ~CCommFace();
private:
    uint mID;  //Face固有のID

    uint mMGLevel;//MultiGrid Level

    // 親Elementの表面
    uint mMeshID;       //Faceが所属するMeshID
    uint mElementID;    //Faceが所属する要素ID
    uint mElementFaceID;//要素のFace番号(局所番号)

    // 接合面の形
    uint mNumOfEdge;//接合面形状は,1辺〜多辺(Polygon)
    uint mFaceType; //接合面の形をElementTypeで表現


    // 面の構成データ
    vector<CCommNode*>  mvCommNode;    // 面の頂点Node
    vector<CCommNode*>  mvEdgeCommNode;// 辺に存在するNode (Refine用途)
    CCommNode*          mpFaceCommNode;// 面に存在するNode (Refine用途)
    vector<CCommFace*>  mvEdgeCommFace;// 辺に隣接するFace,<<<<<-自分自身は除外.
    vector<bool>  mvbEdgeMarking;//辺にノードを生成済みか.(辺ノード生成マーキング)

    // 辺の両端ペア
    PairCommNode        mPairCommNode;//辺の両端ノード
    uint& getEdgeIndex(PairCommNode& pairCommNode);//両端のノードから辺番号を特定


    // Refin時の分割された新CommFace
    vector<CCommFace*> mvProgCommFace;

public:
    // 初期化(Type, vector_resize)
    void initialize(const uint& numOfVert, const uint& numOfEdge);
    
    // ID
    void  setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}

    
    // MultiGrid Level
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}


    // 形状
    uint& getType(){ return mFaceType;}
    uint& getNumOfEdge(){ return mNumOfEdge;}


    // ラップしているMesh要素
    // ----
    void setMeshID(const uint& meshID){ mMeshID= meshID;}
    void setElementID(const uint& elemID){ mElementID= elemID;}
    void setElementFaceID(const uint& iface){ mElementFaceID= iface;}
    uint& getMeshID(){ return mMeshID;}
    uint& getElementID(){ return mElementID;}
    uint& getElementFaceID(){ return mElementFaceID;}




    // 頂点のノード
    // ----
    void setVertCommNode(const uint& ivert, CCommNode* pCommNode){ mvCommNode[ivert]= pCommNode;}
    uint getVertCommNodeSize(){ return mvCommNode.size();}
    CCommNode* getVertCommNode(const uint& ivert){ return mvCommNode[ivert];}
    vector<CCommNode*>& getVertCommNode(){ return mvCommNode;}


    // 辺の両端のCommNode
    // ----
    PairCommNode& getEdgePairCommNode(const uint& iedge);

    // 辺ノード
    // ----
    CCommNode* getEdgeCommNode(const uint& iedge){  return  mvEdgeCommNode[iedge];}

    
    // Edge集合処理のFaceへの処理
    // ----
    void setEdgeCommFace(CCommFace* pNeibFace, const uint& iedge);
    void setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, const uint& iedge);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode);
    void markingEdgeNode(const uint& iedge);
    void markingEdgeNode(PairCommNode& pairCommNode);
    bool isEdgeNodeMarking(const uint& iedge){ return mvbEdgeMarking[iedge];}



    // 面ノード
    // ----
    void setFaceCommNode(CCommNode* pFaceCommNode){ mpFaceCommNode= pFaceCommNode;}
    CCommNode* getFaceCommNode(){ return mpFaceCommNode;}


    // 面の分割
    // ----
    vector<CCommFace*>& refine(CElement* pElement);

};
#endif	/* _COMMFACE_H */
}



