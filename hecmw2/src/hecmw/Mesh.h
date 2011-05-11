//
//  Mesh.h
//			2010.02.28
//			2008.11.10
//			k.Takeda
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"

#include "Node.h"
#include "Element.h"


// Boundary
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "BoundaryEdgeMesh.h"
#include "BoundaryVolumeMesh.h"

#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryEdge.h"
#include "BoundaryVolume.h"

// BoundaryGroup
// FaceMesh,EdgeMesh, VolumeMesh
// それぞれの"境界種類"の管理テンプレートとして使用
#include "BoundaryGroup.h"
// BoundaryNodeMesh
#include "BNodeMeshGrp.h"

// Node around Elements
#include "AggregateElement.h"

// Node around Nodes
#include "AggregateNode.h"

// 通信領域(要素共有型)
#include <map>
#include "CommunicationMesh.h"
// ID->Index (頂点Node,Element)
#include "IndexBucket.h"
// 通信要素の並べ替え
#include "QuickSort.h"


// 通信節点界面(節点共有)
#include "CommMesh2.h"

// FEM, FVM
#include "SolutionType.h"


// グループ
#include "ElementGroup.h"

// MPI
#include "HEC_MPI.h"

// Edge
typedef std::pair<pmw::CNode*,pmw::CNode*> PairNode;

namespace pmw{
#ifndef _MESH_HH_C77868BF_2116_41b3_
#define _MESH_HH_C77868BF_2116_41b3_
class CMesh
{
private:
    Utility::CLogger *mpLogger;
public:
    CMesh(void);
    CMesh(const uiint& numofnode, const uiint& numofelem);
    virtual ~CMesh(void);

protected:
    // Mesh
    // --
    uiint            mnMeshID;   // Mesh ID
    vector<CNode*>    mvNode;   // Node
    vector<CElement*> mvElement;// Element

    // SolutionType
    uiint   mnSolutionType;// FEM, FVM

    // Property : 構造:0 || 流体:1 || ...未定
    uiint mnProp;


    // 計算領域の終端を表す番号 '09.09.09
    // --
    uiint mNodeEndIndex;//mvNode内の計算領域の終端Index(サイズ)  <= 計算に使用する配列数
    uiint mElemEndIndex;//mvElement内の計算領域の終端Index(サイズ) <= 計算に使用する配列数


    //    // Node, Elementのハッシュ "ID => Index"  <- IndexBucketがあるので不要
    //    // --
    //    map<uint, uint, less<uint> > mmNodeIndex;
    //    map<uint, uint, less<uint> > mmElementIndex;


    // 節点集合Node, 節点集合Element
    // --
    vector<CAggregateElement*> mvAggElement;// Aggregate-Elemens(ノード周囲の要素集合を,ノード数分確保)
    vector<CAggregateNode*>    mvAggNode;   // Aggregate-Nodes(ノード周囲のノード集合を,ノード数ぶん確保)


    // restriction(下位_Mesh)のノード管理
    // 自分自身のノードではない -> 同じ型Meshの下位にくるrestriction_Meshのノードのdeleteを回避
    //    -> デストラクタで使用
    uiint mMGLevel;//Mesh自身のMultiGrid_Level
    uiint mMaxMGLevel;//GMGModel全体での最大Level


    // 通信領域(CommMesh)
    // --
    vector<CCommMesh*>     mvCommMesh;//Meshが通信する,通信領域を管理(通信領域Mesh=>CommMesh)
    map<uiint, uiint, less<uiint> > mmCommIndex;//CommID(通信番号:グローバル番号) => mvCommMesh Index
    
    
    // 通信界面(CommMesh2):節点共有型
    // --
    vector<CCommMesh2*> mvCommMesh2;
    map<uiint, uiint, less<uiint> > mmComm2Index;//CommID(通信番号:グローバル番号) => mvCommMesh2 Index
    
    
    
    
    // IndexBucket for "Node & Element"
    // --
    CIndexBucket     moBucket;

    uiint mNumOfNode, mNumOfElement;// Number of Node, Element(ファイル入力時 => いらないかも)
    uiint maxNodeID,  maxElementID;// max ID(ファイル入力時 => いらないかも)
    uiint minNodeID,  minElementID;// min ID(ファイル入力時 => いらないかも)

    // BoundaryGroup_Template エンティティ種類別に管理(節点,面,辺,体積)
    // --
    CBNodeMeshGrp                       *mpBNodeMeshGrp;
    BoundaryGroup<CBoundaryFaceMesh*>   mGrpBndFaceMesh;
    BoundaryGroup<CBoundaryEdgeMesh*>   mGrpBndEdgeMesh;
    BoundaryGroup<CBoundaryVolumeMesh*> mGrpBndVolumeMesh;


    //prolongater用-親Nodeセット
    void setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode);//辺から生成される子Nodeに,親Nodeにセット
    void setupParentNode(vector<CNode*>& vNode, CNode* inNode);       //面,体から生成される子Nodeに,1親Nodeをセット
    
    
    // グループ{ Element }
    // --
    vector<CElementGroup*> mvElementGroup;
    map<uiint, uiint, less<uiint> > mmElemGrpID2IX;//GrpID => Index


public:
    // Mesh ID
    //
    void setMeshID(const uiint& id){ mnMeshID = id;}
    uiint& getMeshID(){ return mnMeshID;}

    // MultiGrid Level
    //
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){return mMGLevel;}
    void setMaxMGLevel(const uiint& maxLevel){ mMaxMGLevel= maxLevel;}

    // SolutionType:: FEM .or. FVM
    //
    void setSolutionType(const uiint& nSolutionType){ mnSolutionType = nSolutionType;}

    // Property : 構造:0 || 流体:1 || ...
    //
    void  setProp(const uiint& nProp){ mnProp = nProp;}
    uiint& getProp(){ return mnProp;}


    // Bucket
    //
    void initBucketNode(const uiint& max_id, const uiint& min_id);
    void setupBucketNodeIndex(const uiint& id, const uiint& index);
    void setupBucketNode();// 一括処理
    void initBucketElement(const uiint& max_id, const uiint& min_id);

    void setupBucketElementIndex(const uiint& id, const uiint& index);
    void setupBucketElement();// 一括処理

    virtual CIndexBucket* getBucket(){ return &moBucket;}

    uiint& getMaxNodeID(){ return moBucket.getMaxNodeID();}
    uiint& getMaxElementID(){ return moBucket.getMaxElementID();}
    uiint& getMinNodeID(){ return moBucket.getMinNodeID();}
    uiint& getMinElementID(){ return moBucket.getMinElementID();}


    // Node
    //
    void reserveNode(const uiint& num_of_node);
    void setNode(CNode* pNode);
    CNode* getNode(const uiint& id);     //IDを指定してNode*を取得
    CNode* getNodeIX(const uiint& index);//直接Indexを指定してNode*を取得

    uiint  getNodeSize(){ return mvNode.size();}
    uiint& getNumOfNode(){return mNumOfNode;}
    //void setNumOfNode(const uint& num){ numOfNode= num;}// <= こちらを使用する場合は,setupNumOfNode()は不使用
    void setupNumOfNode(){ mNumOfNode= mvNode.size();}   // <= こちらを使用する場合は,setNumOfNode()は不使用

    // <<< 計算に使用するNode配列数 >>>
    // --
    uiint& getEfficientNodeSize(){ return mNodeEndIndex;}//計算に有効なNode配列数


    // Element
    //
    void reserveElement(const uiint& num_of_elem);
    void setElement(CElement* pElement);
    CElement* getElement(const uiint& id);     //IDを指定してElement*を取得
    CElement* getElementIX(const uiint& index);//直接Indexを指定してElement*を取得
    vector<CElement*>& getElements(){ return mvElement;}

    uiint  getElementSize(){ return mvElement.size();}
    uiint& getNumOfElement(){ return mNumOfElement;}
    //void setNumOfElement(const uint& num){ numOfElement= num;}// <= こちらを使用する場合は,setupNumOfElement()は不使用
    void setupNumOfElement(){ mNumOfElement= mvElement.size();}// <= こちらを使用する場合は,setNumOfElement()は不使用
    
    // <<< 計算に使用するElement配列数 >>>
    // --
    uiint& getEfficientElementSize(){ return mElemEndIndex;}//計算に有効なElement配列数





    // Aggregate
    // ---
    // Aggregate-Element & Aggregate-Node (Node数と同じだけ存在)
    // ---
    void resizeAggregate(const uiint& res_size);// mvAggElement,mvAggNodeのリザーブ
    void setupAggregate(const uiint& nLevel);   // mvAggElement,mvAggNodeのセットアップ:nLevelは、2次要素の初期(Level=0)処理のために使用

    void setAggNode(CAggregateNode* pAggNode, const uiint& inode);
    void setAggElement(CAggregateElement* pAggElem, const uiint& inode);

    vector<CAggregateElement*>& getAggElems(){ return mvAggElement;}
    vector<CAggregateNode*>&    getAggNodes(){ return mvAggNode;}
    CAggregateElement* getAggElem(const uiint& node_id);
    CAggregateElement* getAggElemIX(const uiint& inode);
    CAggregateNode*    getAggNode(const uiint& node_id);
    CAggregateNode*    getAggNodeIX(const uiint& inode);

    
    // prolongationの準備
    // --    
    // 1.EdgeElementのセットアップ,辺中間ノードのセット
    // 2.FaceElementのセットアップ,面中間ノードのセット
    //   引数：pProgMeshは、prolongationされた上位Mesh
    //
    void presetProgMesh(CMesh* pProgMesh);//prolongationノード&要素をセットする前に,土台Meshのデータをプリセット.
    void setupEdgeElement(CMesh* pProgMesh, const uiint& nLevel);//新たに生成したNodeは,pProgMeshにセット. EdgeElementは,MeshのElementにセット.
    void setupFaceElement(CMesh* pProgMesh);//新たに生成したNodeは,pProgMeshにセット. FaceElementは,MeshのElementにセット.
    void setupVolumeNode(CMesh* pProgMesh); //新たに生成したNodeは,pProgMeshにセット
    void replaceEdgeNode();//要素が二次要素の場合、要素クラス内で辺ノードをmvNodeに移し替え
protected:
    CNode* GeneInterNode(CNode *pNode);
    void avgCoord(vector<CNode*> vCnvNode, CNode* pNode);//複数ノード座標の平均をpNodeにセット
    
public:
    // 境界条件Mesh
    // --
    // Boundary Node
    void setBNodeMeshGrp(CBNodeMeshGrp *pBNodeMeshGrp){ mpBNodeMeshGrp= pBNodeMeshGrp;}
    CBNodeMeshGrp* getBNodeMeshGrp(){ return mpBNodeMeshGrp;}

    void reserveBndNodeMesh(const uiint& res_size){ mpBNodeMeshGrp->reserveBndNodeMesh(res_size);}
    void setBndNodeMesh(CBoundaryNodeMesh *pBNodeMesh){
        //cout << "Mesh::setBndNodeMesh, id " << pBNodeMesh->getID() << endl;
        mpBNodeMeshGrp->setBndNodeMesh(pBNodeMesh);
    }
    CBoundaryNodeMesh* getBndNodeMeshIX(const uiint& index){ return mpBNodeMeshGrp->getBndNodeMeshIX(index);}
    CBoundaryNodeMesh* getBndNodeMeshID(const uiint& id){ return mpBNodeMeshGrp->getBndNodeMeshID(id);}
    uiint getNumOfBoundaryNodeMesh(){ return mpBNodeMeshGrp->getNumOfBoundaryNodeMesh();}


    // Boundary Face
    void reserveBndFaceMesh(const uiint& res_size){ mGrpBndFaceMesh.reserve(res_size);}
    void setBndFaceMesh(CBoundaryFaceMesh *pBFaceMesh){ mGrpBndFaceMesh.push(pBFaceMesh);}
    CBoundaryFaceMesh* getBndFaceMeshIX(const uiint& index){ return mGrpBndFaceMesh.get_withIndex(index);}
    CBoundaryFaceMesh* getBndFaceMeshID(const uiint& id){ return mGrpBndFaceMesh.get_withID(id);}
    uiint getNumOfBoundaryFaceMesh(){ return mGrpBndFaceMesh.NumOfBoundary();}

    // Boundary Volume
    void reserveBndVolumeMesh(const uiint& res_size){ mGrpBndVolumeMesh.reserve(res_size);}
    void setBndVolumeMesh(CBoundaryVolumeMesh *pBVolMesh){ mGrpBndVolumeMesh.push(pBVolMesh);}
    CBoundaryVolumeMesh* getBndVolumeMeshIX(const uiint& index){ return mGrpBndVolumeMesh.get_withIndex(index);}
    CBoundaryVolumeMesh* getBndVolumeMeshID(const uiint& id){ return mGrpBndVolumeMesh.get_withID(id);}
    uiint getNumOfBoundaryVolumeMesh(){ return mGrpBndVolumeMesh.NumOfBoundary();}

    // Boundary Edge
    void reserveBndEdgeMesh(const uiint& res_size){ mGrpBndEdgeMesh.reserve(res_size);}
    void setBndEdgeMesh(CBoundaryEdgeMesh *pBEdgeMesh){ mGrpBndEdgeMesh.push(pBEdgeMesh);}
    CBoundaryEdgeMesh* getBndEdgeMeshIX(const uiint& index){ return mGrpBndEdgeMesh.get_withIndex(index);}
    CBoundaryEdgeMesh* getBndEdgeMeshID(const uiint& id){ return mGrpBndEdgeMesh.get_withID(id);}
    uiint getNumOfBoundaryEdgeMesh(){ return mGrpBndEdgeMesh.NumOfBoundary();}




    // 通信(CCommMesh):CommID => 通信番号
    // --
    void reserveCommMesh(const uiint& res_size){ mvCommMesh.reserve(res_size);}//通信対象領域数に応じた,通信テーブル領域確保
    void setCommMesh(CCommMesh* pCommMesh);// Hashデータ:mmCommIndex[commID] => vector Index もセット(CommMeshにはプロパティを全てセットしておくこと)
    uiint getNumOfCommMesh(){return mvCommMesh.size();}
    CCommMesh* getCommMesh(const uiint& comID);
    CCommMesh* getCommMeshIX(const uiint& index){ return mvCommMesh[index];}

    // 通信で不要になったNode,Elementの並び替え
    // --
    // CommMesh内での名称：不要Node=>mvDNode, 不要Element=>mvDElement
    //  => prolongation時に発生 -> Factory::refineMesh()からコール
    // --
    void sortMesh();
    
    
    
    // 通信(CommMesh2):(節点共有型)
    // --
    void setCommMesh2(CCommMesh2 *pCommMesh2);
    CCommMesh2* getCommMesh2(const uiint& comID);
    CCommMesh2* getCommMesh2IX(const uiint& index){ return mvCommMesh2[index];}
    uiint getCommMesh2Size(){ return mvCommMesh2.size();}


    // Refine後の処理
    // --
    // 辺-面 Node*の削除
    // 辺-面 Element*の削除
    // --
    void deleteProgData();

    //  CMW::FinalizeRefine()からコール
    //
    void deleteAggregate_on_Node();//vertexのAggElemを削除 <= 全てのMesh処理が終わった後で呼び出す


    // グループ { Element }
    // --
    void addElemGrp(CElementGroup* pElemGrp);
    uiint getNumOfElemGrp();
    CElementGroup* getElemGrpIX(const uiint& index);
    CElementGroup* getElemGrpID(const uiint& nGrpID);


////    // ID 決定の為のメソッド { 下位グリッドからの増加節点 }
////    //
////    uint increaseNode(CMesh *pRestMesh);
};
#endif
}

















