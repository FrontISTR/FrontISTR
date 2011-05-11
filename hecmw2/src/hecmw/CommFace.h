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
    uiint mID;  //Face固有のID

    uiint mMGLevel;//MultiGrid Level

    // 親Elementの表面
    uiint mMeshID;       //Faceが所属するMeshID
    uiint mElementID;    //Faceが所属する要素ID
    uiint mElementFaceID;//要素のFace番号(局所番号)

    // 接合面の形
    uiint mNumOfEdge;//接合面形状は,1辺〜多辺(Polygon)
    uiint mFaceType; //接合面の形をElementTypeで表現



    // 面の構成データ
    vector<CCommNode*>  mvCommNode;    // 面の頂点Node
    vector<CCommNode*>  mvEdgeCommNode;// 辺に存在するNode (Refine用途)
    CCommNode*          mpFaceCommNode;// 面に存在するNode (Refine用途)
    vector<CCommFace*>  mvEdgeCommFace;// 辺に隣接するFace,<<<<<-自分自身は除外.
    bool*  mvbEdgeMarking;//辺にノードを生成済みか.(辺ノード生成マーキング)

    // 辺の両端ペア
    //PairCommNode        mPairCommNode;//辺の両端ノード
    uiint& getEdgeIndex(PairCommNode& pairCommNode);//両端のノードから辺番号を特定


    // Refin時の分割された新CommFace
    vector<CCommFace*> mvProgCommFace;

public:
    // 初期化(Type, vector_resize)
    void initialize(const uiint& numOfVert, const uiint& numOfEdge, const uiint& nOrder);
    
    // ID
    void  setID(const uiint& id){ mID= id;}
    uiint& getID(){ return mID;}

    
    // MultiGrid Level
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}


    // 形状
    uiint& getType(){ return mFaceType;}
    uiint& getNumOfEdge(){ return mNumOfEdge;}
    uiint  getNumOfVert();

    // 1次 || 2次
    uiint getOrder();

    // ラップしているMesh要素
    // ----
    void setMeshID(const uiint& meshID){ mMeshID= meshID;}
    void setElementID(const uiint& elemID){ mElementID= elemID;}
    void setElementFaceID(const uiint& iface){ mElementFaceID= iface;}
    uiint& getMeshID(){ return mMeshID;}
    uiint& getElementID(){ return mElementID;}
    uiint& getElementFaceID(){ return mElementFaceID;}




    // 通信ノード全体
    // ----
    //void setVertCommNode(const uint& ivert, CCommNode* pCommNode){ mvCommNode[ivert]= pCommNode;}
    void setCommNode(const uiint& index, CCommNode* pCommNode){ mvCommNode[index]= pCommNode;}
    //uint getVertCommNodeSize(){ return mvCommNode.size();}
    uiint getCommNodeSize(){ return mvCommNode.size();}
    //CCommNode* getVertCommNode(const uint& ivert){ return mvCommNode[ivert];}
    //vector<CCommNode*>& getVertCommNode(){ return mvCommNode;}
    CCommNode* getCommNode(const uiint& index){ return mvCommNode[index];}
    vector<CCommNode*>& getCommNode(){ return mvCommNode;}


    // 辺の両端のCommNode
    // ----
    PairCommNode getEdgePairCommNode(const uiint& iedge);

    // 辺ノード
    // ----
    CCommNode* getEdgeCommNode(const uiint& iedge){  return  mvEdgeCommNode[iedge];}

    
    // Edge集合処理のFaceへの処理
    // ----
    void setEdgeCommFace(CCommFace* pNeibFace, const uiint& iedge);
    void setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, const uiint& iedge);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode);
    void markingEdgeNode(const uiint& iedge);
    void markingEdgeNode(PairCommNode& pairCommNode);
    bool isEdgeNodeMarking(const uiint& iedge){ return mvbEdgeMarking[iedge];}

    // 2次要素面のときの辺ノード => mvCommNodeへ移設
    //
    void replaceEdgeCommNode();


    // 面ノード
    // ----
    void setFaceCommNode(CCommNode* pFaceCommNode){ mpFaceCommNode= pFaceCommNode;}
    CCommNode* getFaceCommNode(){ return mpFaceCommNode;}


    // 面の分割
    // ----
    vector<CCommFace*>& refine(CElement* pElement);


    // Refine後のvector解放
    // ----
    void deleteProgData();

};
#endif	/* _COMMFACE_H */
}



