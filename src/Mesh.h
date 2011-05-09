//
//
//			2008.09.28
//			2008.11.10
//			k.Takeda
#ifndef MESH_HH_C77868BF_2116_41b3_B667_59BA392FCC5C
#define MESH_HH_C77868BF_2116_41b3_B667_59BA392FCC5C

#include "CommonStd.h"
#include "TypeDef.h"

#include "Logger.h"

#include "Node.h"
#include "Element.h"


// Boundary(Node,Face,Volume)
#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryVolume.h"
// BoundaryGroup(BoundaryNode,Face,Volumeの管理テンプレート)
#include "BoundaryGroup.h"

// Node around Elements
#include "AggregateElement.h"

// Node around Nodes
#include "AggregateNode.h"

// 通信領域
#include <map>
#include "CommunicationMesh.h"

#include "IndexBucket.h"

#include "QuickSort.h"

// Edge
typedef std::pair<pmw::CNode*,pmw::CNode*> PairNode;

namespace pmw{
class CMesh
{
private:
    Utility::CLogger *mpLogger;
public:
    CMesh(void);
    CMesh(const uint& numofnode, const uint& numofelem);
    virtual ~CMesh(void);

protected:
    // Mesh
    // --
    uint            mnMeshID;   // Mesh ID
    vector<CNode*>    mvNode;   // Node
    vector<CElement*> mvElement;// Element

    // 計算領域の終端を表す番号 '09.09.09
    // --
    uint mNodeEndIndex;//mvNode内の計算領域の終端Index(サイズ)  <= 計算に使用する配列数
    uint mElemEndIndex;//mvElement内の計算領域の終端Index(サイズ) <= 計算に使用する配列数

    // Node, Elementのハッシュ "ID => Index" '09.09.28
    // --
    map<uint, uint, less<uint> > mmNodeIndex;
    map<uint, uint, less<uint> > mmElementIndex;


    // 節点集合Node, 節点集合Element
    // --
    vector<CAggregateElement*> mvAggElement;// Aggregate-Elemens(ノード周囲の要素集合を,ノード数分確保)
    vector<CAggregateNode*>    mvAggNode;   // Aggregate-Nodes(ノード周囲のノード集合を,ノード数ぶん確保)


    // restriction(下位_Mesh)のノード管理
    // 自分自身のノードではない -> 同じ型Meshの下位にくるrestriction_Meshのノードのdeleteを回避
    //    -> デストラクタで使用
    uint mMGLevel;//Mesh自身のMultiGrid_Level


    // 通信領域(CommMesh)
    // --
    vector<CCommMesh*>     mvCommMesh;//Meshが通信する,通信領域を管理(通信領域Mesh=>CommMesh)
    map<uint, uint, less<uint> > mmCommIndex;//CommID(通信番号:グローバル番号) => mvCommMesh Index
    
    
    // IndexBucket for "Node & Element"
    // --
    CIndexBucket     moBucket;

    uint mNumOfNode, mNumOfElement;// Number of Node, Element(ファイル入力時 => いらないかも)
    uint maxNodeID, maxElementID;// max ID(ファイル入力時 => いらないかも)
    uint minNodeID, minElementID;// min ID(ファイル入力時 => いらないかも)

    // BoundaryGroup_Template 
    // --
    BoundaryGroup<CBoundaryNode*>   mBoundaryNodes;
    BoundaryGroup<CBoundaryFace*>   mBoundaryFaces;
    BoundaryGroup<CBoundaryVolume*> mBoundaryVolumes;



    //debug 用途 dummy counter
    uint mnDummyCount;

public:
    // Mesh ID
    //
    void setMeshID(const uint& id){ mnMeshID = id;}
    uint& getMeshID(){ return mnMeshID;}

    // MultiGrid Level
    //
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){return mMGLevel;}


    // Bucket
    //
    void initBucketNode(const uint& max_id, const uint& min_id);
    void setupBucketNodeIndex(const uint& id, const uint& index);
    void setupBucketNode();// 一括処理
    void initBucketElement(const uint& max_id, const uint& min_id);

    void setupBucketElementIndex(const uint& id, const uint& index);
    void setupBucketElement();// 一括処理

    virtual CIndexBucket* getBucket(){ return &moBucket;}


    // Node
    //
    void reserveNode(const uint& num_of_node);
    void setNode(CNode* pNode);
    CNode* getNode(const uint& id);     //IDを指定してNode*を取得
    CNode* getNodeIX(const uint& index);//直接Indexを指定してNode*を取得

    uint  getNodeSize(){ return mvNode.size();}
    uint& getNumOfNode(){return mNumOfNode;}
    //void setNumOfNode(const uint& num){ numOfNode= num;}// <= こちらを使用する場合は,setupNumOfNode()は不使用
    void setupNumOfNode(){ mNumOfNode= mvNode.size();}   // <= こちらを使用する場合は,setNumOfNode()は不使用

    // <<< 計算に使用するNode配列数 >>>
    // --
    uint& getEfficientNodeSize(){ return mNodeEndIndex;}//計算に有効なNode配列数


    // Element
    //
    void reserveElement(const uint& num_of_elem);
    void setElement(CElement* pElement);
    CElement* getElement(const uint& id);     //IDを指定してElement*を取得
    CElement* getElementIX(const uint& index);//直接Indexを指定してElement*を取得
    vector<CElement*>& getElements(){ return mvElement;}

    uint  getElementSize(){ return mvElement.size();}
    uint& getNumOfElement(){ return mNumOfElement;}
    //void setNumOfElement(const uint& num){ numOfElement= num;}// <= こちらを使用する場合は,setupNumOfElement()は不使用
    void setupNumOfElement(){ mNumOfElement= mvElement.size();}// <= こちらを使用する場合は,setNumOfElement()は不使用
    
    // <<< 計算に使用するElement配列数 >>>
    // --
    uint& getEfficientElementSize(){ return mElemEndIndex;}//計算に有効なElement配列数





    // Aggregate
    // ---
    // Aggregate-Element & Aggregate-Node (Node数と同じだけ存在)
    // ---
    void reserveAggregate(const uint& res_size);// AggElement,AggNodeのリザーブ
    void setAggNode(CAggregateNode* pAggNode);
    void setAggElement(CAggregateElement* pAggElem);
    void setupAggregate();// mvAggElement,mvAggNodeのセットアップ
    vector<CAggregateElement*>& getAggElems(){ return mvAggElement;}
    vector<CAggregateNode*>& getAggNodes(){ return mvAggNode;}

    
    // prolongationの準備
    // --    
    // 1.EdgeElementのセットアップ,辺中間ノードのセット
    // 2.FaceElementのセットアップ,面中間ノードのセット
    //   引数：pProgMeshは、prolongationされた上位Mesh
    //
    void presetProgMesh(CMesh* pProgMesh);//prolongationノード&要素をセットする前に,土台Meshのデータをプリセット.
    void setupEdgeElement(CMesh* pProgMesh);//新たに生成したNodeは,pProgMeshにセット. EdgeElementは,MeshのElementにセット.
    void setupFaceElement(CMesh* pProgMesh);//新たに生成したNodeは,pProgMeshにセット. FaceElementは,MeshのElementにセット.
    void setupVolumeNode(CMesh* pProgMesh); //新たに生成したNodeは,pProgMeshにセット
protected:
    CNode* GeneInterNode(CNode *pNode);
    void avgCoord(vector<CNode*> vCnvNode, CNode* pNode);//複数ノード座標の平均をpNodeにセット

public:
    // 境界条件
    // --
    // Boundary Nodes
    void reserveBoundaryNode(const uint& res_size){ mBoundaryNodes.reserve(res_size);}
    void setBoundaryNode(CBoundaryNode *pBNode){ mBoundaryNodes.push(pBNode);}
    CBoundaryNode* getBoundaryNode_withID(const uint& id){ return mBoundaryNodes.get_withID(id);}
    CBoundaryNode* getBoundaryNode_withIndex(const uint& index){ return mBoundaryNodes.get_withIndex(index);}
    uint getNumOfBoundaryNode(){ return mBoundaryNodes.NumOfBoundary();}

    // Boundary Faces
    void reserveBoundaryFace(const uint& res_size){ mBoundaryFaces.reserve(res_size);}
    void setBoundaryFace(CBoundaryFace *pBFace){ mBoundaryFaces.push(pBFace);}
    CBoundaryFace* getBoundaryFace_withID(const uint& id){ return mBoundaryFaces.get_withID(id);}
    CBoundaryFace* getBoundaryFace_withIndex(const uint& index){ return mBoundaryFaces.get_withIndex(index);}
    uint getNumOfBoundaryFace(){ return mBoundaryFaces.NumOfBoundary();}

    // Boundary Volumes
    void reserveBoundaryVolume(const uint& res_size){ mBoundaryVolumes.reserve(res_size);}
    void setBoundaryVolume(CBoundaryVolume *pBVolume){ mBoundaryVolumes.push(pBVolume);}
    CBoundaryVolume* getBoundaryVolume_withID(const uint& id){ return mBoundaryVolumes.get_withID(id);}
    CBoundaryVolume* getBoundaryVolume_withIndex(const uint& index){ return mBoundaryVolumes.get_withIndex(index);}
    uint getNumOfBoundaryVolume(){ return mBoundaryVolumes.NumOfBoundary();}


    // 通信(CCommMesh):CommID => 通信番号
    // --
    void reserveCommMesh(const uint& res_size){ mvCommMesh.reserve(res_size);}//通信対象領域数に応じた,通信テーブル領域確保
    void setCommMesh(CCommMesh* pCommMesh);// Hashデータ:mmCommIndex[commID] => vector Index もセット(CommMeshにはプロパティを全てセットしておくこと)
    uint getNumOfCommMesh(){return mvCommMesh.size();}
    CCommMesh* getCommMesh(const uint& comID);
    CCommMesh* getCommMeshIX(const uint& index){ return mvCommMesh[index];}

    // 通信で不要になったNode,Elementの並び替え
    // --
    // CommMesh内での名称：不要Node=>mvDNode, 不要Element=>mvDElement
    //  => prolongation時に発生 -> Factory::refineMesh()からコール
    // --
    void sortMesh();
};
}
#endif
