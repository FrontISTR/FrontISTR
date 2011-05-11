/* 
 * File:   SkinFace.h
 * Author: ktakeda
 *
 * ・MPCに使用するモデル表面の面ベースクラス
 * ・SlaveFaceとしても利用
 *
 * Created on 2009/10/15, 14:57
 */
#include "TypeDef.h"
#include "ContactNode.h"

#include "ElementType.h"//要素形状定数
#include "ElementProperty.h"

#include "Node.h"
#include "Element.h" //refineで利用

#include "EdgeTree.h"

#include <utility>//pair
typedef std::pair<pmw::CContactNode*, pmw::CContactNode*> PairConNode;

#include "Logger.h"


namespace pmw{
#ifndef _SKINFACE_H
#define	_SKINFACE_H
class CSkinFace{
public:
    CSkinFace();
    virtual ~CSkinFace();

protected:
    uiint mID;//SkinFace固有のID
    uiint mRank;//SkinFaceのrank
    uiint mLevel;//MG Level

    bool mbSelfDom;//Mesh自身の計算領域内のSkinFaceか.
    uiint mMeshID;//Faceが所属するMeshID
    uiint mElementID;    //Faceが所属する要素ID
    uiint mElementFaceID;//要素のFace番号(局所番号)

    // 面の構成ContactNode
    vector<CContactNode*> mvConNode;

    // 辺に生成するContactNode (Refine用途)
    //
    vector<CContactNode*> mvEdgeNode;//辺に生成するContactNode
    bool* mvbEdgeMarking;            //辺にノードを生成済みか.

    vector<CSkinFace*>    mvEdgeFace;//辺に隣接するSkinFace,<<<<<<<<<<<<<-自分自身は除外.

    uiint& getEdgeIndex(PairConNode& pairConNode);//両端のノードから辺番号を特定

    // 面中心のContactNode (Refine用途)
    CContactNode *mpFaceNode;

    // Refineされた上位レベルのSkinFace
    vector<CSkinFace*> mvProgFace;

    // Shape(形状タイプ)
    uiint mNumOfEdge;
    uiint mShapeType;//ElementTypeによる定数

    uiint mnOrder;//1次｜2次 次数

    // refine()からコールされる関数
    // --
    virtual CSkinFace* generateFace();//refine()でMasterFace,SkinFaceの型違い生成に使用
    CSkinFace* mpOtherFace;//自分以外のSkinFace(MasterFace),返り値-用途


    // 面ベクトル
    vdouble mvNormalVector;
    
public:
    virtual const char* getName(){return "SkinFace";}

    //Level
    void setLevel(const uiint& level){ mLevel= level;}
    uiint& getLevel(){ return mLevel;}

    //SkinFaceのrank
    void setRank(const uiint& rank){ mRank= rank;}
    uiint& getRank(){ return mRank;}

    //SkinFace固有のID
    void setID(const uiint& id){ mID= id;}
    uiint& getID(){ return mID;}

    // 所属Mesh-Element-ElemFace関連
    // 
    void markingSelf(){ mbSelfDom=true;}//計算領域内に存在するSkinFaceであることをマーキング
    bool isSelf(){ return mbSelfDom;}   //計算領域内に存在するSkinFaceか？
    
    // 所属MeshID
    void setMeshID(const uiint& meshID){ mMeshID= meshID;}
    uiint& getMeshID(){ return mMeshID;}
    // 所属要素のID
    void  setElementID(const uiint& id){ mElementID = id;}
    uiint& getElementID(){ return mElementID;}
    // 面番号
    void  setFaceID(const uiint& faceID){ mElementFaceID= faceID;}
    uiint& getFaceID(){ return mElementFaceID;}


    // Face構成のノード
    // --
    void resizeNode(const uiint& size){  mvConNode.resize(size);}
    void addNode(CContactNode* pConNode){ mvConNode.push_back(pConNode);}
    void setNode(CContactNode* pConNode, const uiint& index){ mvConNode[index] = pConNode;}
    void setNode(vector<CContactNode*>& vConNode){ mvConNode = vConNode;} // "="代入と同じ
    void operator=(vector<CContactNode*>& vConNode){ mvConNode = vConNode;}// ":=" setNode(vNode)

    uiint   getNumOfNode() { return mvConNode.size();}
    vector<CContactNode*>&  getNodes(){ return mvConNode;}
    CContactNode*  getNode(const uiint& index){ return mvConNode[index];}


    // 形状タイプ
    // --
    void setShapeType(const uiint& shapeType);
    uiint& getShapeType(){ return mShapeType;}
    uiint& getNumOfEdge(){ return mNumOfEdge;}
    uiint  getNumOfVert();
    uiint& getOrder(){ return mnOrder;}
    


    // --
    // MasterFace用途: CSkinFaceでは無効
    //
    // <> Edge集合のFaceを集める為に"型"を同一のものとして扱う必要があった.
    // <> MasterFace,SkinFaceともにSkinFace*で扱う.
    // --
    virtual void addSlaveNode(CContactNode* pConNode);// スレーブ点の追加, 警告メッセージあり
    virtual CContactNode* getSlaveNode(const uiint& index){ return NULL;}
    virtual uiint getNumOfSlaveNode(){ return 0;}
    virtual void CalcSlave(const uiint& islave, const uiint& valType);// EQUATION::マスター面の値 => スレーブ点の値, SkinFaceなので警告メッセージあり
    virtual double& getCoef(const uiint& slaveID, const uiint& ivert);
    
    
    // Refine関連(Refineのための辺ノード設置)
    // --
    PairConNode getEdgePairNode(const uiint& iedge);
    void setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode);//隣接するSkinFaceのセット
    void setEdgeFace(CSkinFace* pFace, const uiint& iedge);       //隣接するSkinFaceのセット
    void setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode);//mvEdgeNode:辺ノード(ContactNode)のセット
    void setEdgeConNode(CContactNode* pEdgeConNode, const uiint& iedge);       //mvEdgeNode:辺ノード(ContactNode)のセット
    bool isEdgeNodeMarking(const uiint& iedge){ return mvbEdgeMarking[iedge];}
    bool isEdgeNodeMarking(PairConNode& pairConNode);
    void markingEdgeNode(const uiint& iedge){ mvbEdgeMarking[iedge]=true;}
    void markingEdgeNode(PairConNode& pairConNode);
    void setFaceConNode(CContactNode* pFaceConNode){ mpFaceNode= pFaceConNode;}//mpFaceNode:面中心ノード(ContactNode)のセット

    CContactNode* getEdgeConNode(const uiint& iedge){ return mvEdgeNode[iedge];}
    CContactNode* getEdgeConNode(PairConNode& pairConNode);
    CContactNode* getFaceConNode(){ return mpFaceNode;}
    
    vector<CSkinFace*>& getProgFace(){ return mvProgFace;}
    CSkinFace* getProgFace(const uiint& ivert){ return mvProgFace[ivert];}
    void refine(CElement* pElem, uiint& faceID);////MasterFace, SlaveFaceの"再分割"
protected:
    //辺,面に生成されたConNodeへ,メッシュのNodeIDをセットする(対象：progFace)
    void setupNodeID_progFace(CElement* pElem, const uiint& numOfVert);
    //辺に生成されたConNodeへ,メッシュのNodeIDをセット
    void setupEdgeNodeID(CElement* pElem, const uiint& numOfVert);
public:
    //辺に生成されたConNodeへ、メッシュのNodeIDをセットする(2次要素+最終Level専用) (対象：自身)
    void setupNodeID_2nd_LastLevel(CElement* pElem);


public:
    // 面ベクトル
    // --
    vdouble& CalcNzNormalVector();//面ベクトル(正規化)計算
    vdouble& getNzNormalVector();//面ベクトル(正規化)取得のみ
    
    // 辺ノードを mvConNode へ移動 (2次要素面の場合)
    //
    void replaceEdgeNode();

    // 辺ノードvectorのメモリー解放
    //
    void deleteProgData();
};
#endif	/* _SKINFACE_H */
}


