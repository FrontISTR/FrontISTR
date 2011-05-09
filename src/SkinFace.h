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
    uint mID;//SkinFace固有のID

    bool mbSelfDom;//Mesh自身の計算領域内のSkinFaceか.
    uint mMeshID;//Faceが所属するMeshID
    uint mElementID;    //Faceが所属する要素ID
    uint mElementFaceID;//要素のFace番号(局所番号)

    // 面の構成ContactNode
    vector<CContactNode*> mvConNode;

    // 辺に生成するContactNode (Refine用途)
    vector<CContactNode*> mvEdgeNode;//辺に生成するContactNode
    vector<bool>          mvbEdgeMarking;//辺にノードを生成済みか.
    vector<CSkinFace*>    mvEdgeFace;//辺に隣接するSkinFace,<<<<<<<<<<<<<-自分自身は除外.
    PairConNode           mPairConNode;//辺の両端ノードのreturn value用
    uint& getEdgeIndex(PairConNode& pairConNode);//両端のノードから辺番号を特定

    // 面中心のContactNode (Refine用途)
    CContactNode *mpFaceNode;

    // Refineされた上位レベルのSkinFace
    vector<CSkinFace*> mvProgFace;

    // Shape(形状タイプ)
    uint mNumOfEdge;
    uint mShapeType;//ElementTypeによる定数

    // refine()からコールされる関数
    // --
    virtual CSkinFace* generateFace();//refine()でMasterFace,SkinFaceの型違い生成に使用
    CSkinFace* mpOtherFace;//自分以外のSkinFace(MasterFace),返り値-用途
    void setupNodeID_progFace(CElement* pElem, const uint& numOfConNode);//辺,面に生成されたConNodeへ,メッシュのNodeIDをセットする

    // 面ベクトル
    vdouble mvNormalVector;
    
    // Logger
    Utility::CLogger* mpLogger;
public:
    virtual const char* getName(){return "SkinFace";}

    //SkinFace固有のID
    void setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}

    // 所属Mesh-Element-ElemFace関連
    // 
    void markingSelf(){ mbSelfDom=true;}//計算領域内に存在するSkinFaceであることをマーキング
    bool isSelf(){ return mbSelfDom;}   //計算領域内に存在するSkinFaceか？
    
    // 所属MeshID
    void setMeshID(const uint& meshID){ mMeshID= meshID;}
    uint& getMeshID(){ return mMeshID;}
    // 所属要素のID
    void  setElementID(const uint& id){ mElementID = id;}
    uint& getElementID(){ return mElementID;}
    // 面番号
    void  setFaceID(const uint& faceID){ mElementFaceID= faceID;}
    uint& getFaceID(){ return mElementFaceID;}


    // Face構成のノード
    // --
    void resizeNode(const uint& size){  mvConNode.resize(size);}
    void addNode(CContactNode* pConNode){ mvConNode.push_back(pConNode);}
    void setNode(CContactNode* pConNode, const uint& index){ mvConNode[index] = pConNode;}
    void setNode(vector<CContactNode*>& vConNode){ mvConNode = vConNode;} // "="代入と同じ
    void operator=(vector<CContactNode*>& vConNode){ mvConNode = vConNode;}// ":=" setNode(vNode)

    uint   getNumOfNode() { return mvConNode.size();}
    vector<CContactNode*>&  getNodes(){ return mvConNode;}
    CContactNode*  getNode(const uint& index){ return mvConNode[index];}

    // 形状タイプ
    // --
    void setShapeType(const uint& shapeType);
    uint& getShapeType(){ return mShapeType;}
    uint& getNumOfEdge(){ return mNumOfEdge;}
    


    // --
    // MasterFace用途: CSkinFaceでは無効
    //
    // <> Edge集合のFaceを集める為に"型"を同一のものとして扱う必要があった.
    // <> MasterFace,SkinFaceともにSkinFace*で扱う.
    // --
    virtual void addSlaveNode(CContactNode* pConNode);// スレーブ点の追加, 警告メッセージあり
    virtual CContactNode* getSlaveNode(const uint& index){ return NULL;}
    virtual uint getNumOfSlaveNode(){ return 0;}
    virtual void CalcSlave(const uint& islave, const uint& valType);// EQUATION::マスター面の値 => スレーブ点の値, SkinFaceなので警告メッセージあり

    
    
    // Refine関連(Refineのための辺ノード設置)
    // --
    PairConNode& getEdgePairNode(const uint& iedge);
    void setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode);//隣接するSkinFaceのセット
    void setEdgeFace(CSkinFace* pFace, const uint& iedge);       //隣接するSkinFaceのセット
    void setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode);//mvEdgeNode:辺ノード(ContactNode)のセット
    void setEdgeConNode(CContactNode* pEdgeConNode, const uint& iedge);       //mvEdgeNode:辺ノード(ContactNode)のセット
    bool isEdgeNodeMarking(const uint& iedge){ return mvbEdgeMarking[iedge];}
    bool isEdgeNodeMarking(PairConNode& pairConNode);
    void markingEdgeNode(const uint& iedge){ mvbEdgeMarking[iedge]=true;}
    void markingEdgeNode(PairConNode& pairConNode);
    void setFaceConNode(CContactNode* pFaceConNode){ mpFaceNode= pFaceConNode;}//mpFaceNode:面中心ノード(ContactNode)のセット

    CContactNode* getEdgeConNode(const uint& iedge){ return mvEdgeNode[iedge];}
    CContactNode* getEdgeConNode(PairConNode& pairConNode);
    CContactNode* getFaceConNode(){ return mpFaceNode;}
    
    vector<CSkinFace*>& getProgFace(){ return mvProgFace;}
    CSkinFace* getProgFace(const uint& ivert){ return mvProgFace[ivert];}
    void refine(CElement* pElem, uint& faceID);////MasterFace, SlaveFaceの"再分割"


    // 面ベクトル
    // --
    vdouble& CalcNzNormalVector();//面ベクトル(正規化)計算
    vdouble& getNzNormalVector();//面ベクトル(正規化)取得のみ
    
};
#endif	/* _SKINFACE_H */
}


