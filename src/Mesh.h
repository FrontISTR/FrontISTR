//
//
//			2008.05.15
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


#include "IndexBucket.h"

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
    //
    uint            mnMeshID;   // Mesh ID
    vector<CNode*>    mvNode;   // Node
    vector<CElement*> mvElement;// Element
    
    vector<CAggregateElement*> mvAggElement;// Aggregate-Elemens(ノード周囲の要素集合を,ノード数分確保)
    vector<CAggregateNode*>    mvAggNode;   // Aggregate-Nodes(ノード周囲のノード集合を,ノード数ぶん確保)

    // restriction(下位_Mesh)のノード管理
    // 自分自身のノードではない -> 同じ型Meshの下位にくるrestriction_Meshのノードのdeleteを回避
    //    -> デストラクタで使用
    uint mMGLevel;//Mesh自身のMultiGrid_Level


    // IndexBucket for "Node & Element"
    //
    CIndexBucket     moBucket;

    uint numOfNode, numOfElement;// Number of Node, Element
    uint maxNodeID, maxElementID;// max ID
    uint minNodeID, minElementID;// min ID

    // BoundaryGroup_Template 
    //
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
    virtual void reserveNode(const uint& num_of_node);
    virtual void setNode(CNode* pNode);
    virtual CNode* getNode_inID(const uint& id);
    virtual CNode* getNode(const uint& index);

    uint  getNodeSize(){ return mvNode.size();}
    uint& getNumOfNode(){return numOfNode;}
    //void setNumOfNode(const uint& num){ numOfNode= num;}// <= こちらを使用する場合は,setupNumOfNode()は不使用
    void setupNumOfNode(){ numOfNode= mvNode.size();}   // <= こちらを使用する場合は,setNumOfNode()は不使用


    // Element
    //
    virtual void reserveElement(const uint& num_of_elem);
    virtual void setElement(CElement* pElement);
    virtual CElement* getElement_inID(const uint& id);
    virtual CElement* getElement(const uint& index);

    uint  getElementSize(){ return mvElement.size();}
    uint& getNumOfElement(){ return numOfElement;}
    //void setNumOfElement(const uint& num){ numOfElement= num;}// <= こちらを使用する場合は,setupNumOfElement()は不使用
    void setupNumOfElement(){ numOfElement= mvElement.size();}// <= こちらを使用する場合は,setNumOfElement()は不使用


    // Aggregate
    // ---
    // Aggregate Elements (Node数と同じだけ存在)
    //
    //void reserveAggElement(const uint& res_size);
    //void setAggElement(CAggregateElement* pAggElem);
    //void setupAggElement(); // mvAggElementのセットアップ
    //
    // Aggregate-Element & Aggregate-Node
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
};
}
#endif
