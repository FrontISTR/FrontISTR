/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Mesh.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef MESH_HH_C77868BF_2116_41b3_B667_59BA392FCC5C
#define MESH_HH_C77868BF_2116_41b3_B667_59BA392FCC5C
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
#include "Node.h"
#include "Element.h"
#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryVolume.h"
#include "BoundaryGroup.h"
#include "AggregateElement.h"
#include "AggregateNode.h"
#include <map>
#include "CommunicationMesh.h"
#include "IndexBucket.h"
#include "QuickSort.h"
#include "CommMesh2.h"
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
    uint            mnMeshID;   
    vector<CNode*>    mvNode;   
    vector<CElement*> mvElement;
    uint mNodeEndIndex;
    uint mElemEndIndex;
    vector<CNode*> mvEdgeNode; 
    map<uint, uint, less<uint> > mmEdgeNodeID2IX;
    vector<CAggregateElement*> mvAggElement;
    vector<CAggregateNode*>    mvAggNode;   
    uint mMGLevel;
    vector<CCommMesh*>     mvCommMesh;
    map<uint, uint, less<uint> > mmCommIndex;
    vector<CCommMesh2*> mvCommMesh2;
    map<uint, uint, less<uint> > mmComm2Index;
    CIndexBucket     moBucket;
    uint mNumOfNode, mNumOfElement;
    uint maxNodeID, maxElementID;
    uint minNodeID, minElementID;
    BoundaryGroup<CBoundaryNode*>   mBoundaryNodes;
    BoundaryGroup<CBoundaryFace*>   mBoundaryFaces;
    BoundaryGroup<CBoundaryVolume*> mBoundaryVolumes;
    void setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode);
    void setupParentNode(vector<CNode*>& vNode, CNode* inNode);       
    void setupChildNode(CNode* pNode0, CNode* pNode1, CNode* inNode);
    void setupChildNode(vector<CNode*>& vNode, CNode* inNode);       
    uint mnDummyCount;
public:
    void setMeshID(const uint& id){ mnMeshID = id;}
    uint& getMeshID(){ return mnMeshID;}
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){return mMGLevel;}
    void initBucketNode(const uint& max_id, const uint& min_id);
    void setupBucketNodeIndex(const uint& id, const uint& index);
    void setupBucketNode();
    void initBucketElement(const uint& max_id, const uint& min_id);
    void setupBucketElementIndex(const uint& id, const uint& index);
    void setupBucketElement();
    virtual CIndexBucket* getBucket(){ return &moBucket;}
    void reserveNode(const uint& num_of_node);
    void setNode(CNode* pNode);
    CNode* getNode(const uint& id);     
    CNode* getNodeIX(const uint& index);
    uint  getNodeSize(){ return mvNode.size();}
    uint& getNumOfNode(){return mNumOfNode;}
    void setupNumOfNode(){ mNumOfNode= mvNode.size();}   
    uint& getEfficientNodeSize(){ return mNodeEndIndex;}
    void reserveElement(const uint& num_of_elem);
    void setElement(CElement* pElement);
    CElement* getElement(const uint& id);     
    CElement* getElementIX(const uint& index);
    vector<CElement*>& getElements(){ return mvElement;}
    uint  getElementSize(){ return mvElement.size();}
    uint& getNumOfElement(){ return mNumOfElement;}
    void setupNumOfElement(){ mNumOfElement= mvElement.size();}
    uint& getEfficientElementSize(){ return mElemEndIndex;}
    void reserveAggregate(const uint& res_size);
    void setAggNode(CAggregateNode* pAggNode);
    void setAggElement(CAggregateElement* pAggElem);
    void setupAggregate();
    vector<CAggregateElement*>& getAggElems(){ return mvAggElement;}
    vector<CAggregateNode*>& getAggNodes(){ return mvAggNode;}
    CAggregateElement* getAggElem(const uint& node_id);
    CAggregateNode*    getAggNode(const uint& node_id);
    void presetProgMesh(CMesh* pProgMesh);
    void setupEdgeElement(CMesh* pProgMesh);
    void setupFaceElement(CMesh* pProgMesh);
    void setupVolumeNode(CMesh* pProgMesh); 
protected:
    CNode* GeneInterNode(CNode *pNode);
    void avgCoord(vector<CNode*> vCnvNode, CNode* pNode);
public:
    void reserveBoundaryNode(const uint& res_size){ mBoundaryNodes.reserve(res_size);}
    void setBoundaryNode(CBoundaryNode *pBNode){ mBoundaryNodes.push(pBNode);}
    CBoundaryNode* getBoundaryNode_withID(const uint& id){ return mBoundaryNodes.get_withID(id);}
    CBoundaryNode* getBoundaryNode_withIndex(const uint& index){ return mBoundaryNodes.get_withIndex(index);}
    uint getNumOfBoundaryNode(){ return mBoundaryNodes.NumOfBoundary();}
    void reserveBoundaryFace(const uint& res_size){ mBoundaryFaces.reserve(res_size);}
    void setBoundaryFace(CBoundaryFace *pBFace){ mBoundaryFaces.push(pBFace);}
    CBoundaryFace* getBoundaryFace_withID(const uint& id){ return mBoundaryFaces.get_withID(id);}
    CBoundaryFace* getBoundaryFace_withIndex(const uint& index){ return mBoundaryFaces.get_withIndex(index);}
    uint getNumOfBoundaryFace(){ return mBoundaryFaces.NumOfBoundary();}
    void reserveBoundaryVolume(const uint& res_size){ mBoundaryVolumes.reserve(res_size);}
    void setBoundaryVolume(CBoundaryVolume *pBVolume){ mBoundaryVolumes.push(pBVolume);}
    CBoundaryVolume* getBoundaryVolume_withID(const uint& id){ return mBoundaryVolumes.get_withID(id);}
    CBoundaryVolume* getBoundaryVolume_withIndex(const uint& index){ return mBoundaryVolumes.get_withIndex(index);}
    uint getNumOfBoundaryVolume(){ return mBoundaryVolumes.NumOfBoundary();}
    void reserveCommMesh(const uint& res_size){ mvCommMesh.reserve(res_size);}
    void setCommMesh(CCommMesh* pCommMesh);
    uint getNumOfCommMesh(){return mvCommMesh.size();}
    CCommMesh* getCommMesh(const uint& comID);
    CCommMesh* getCommMeshIX(const uint& index){ return mvCommMesh[index];}
    void sortMesh();
    void setCommMesh2(CCommMesh2 *pCommMesh2);
    CCommMesh2* getCommMesh2(const uint& comID);
    CCommMesh2* getCommMesh2IX(const uint& index){ return mvCommMesh2[index];}
    uint getCommMesh2Size(){ return mvCommMesh2.size();}
    void setEdgeNode(CNode* pNode);
    CNode* getEdgeNode(const uint& id);
    uint getEdgeNodeSize(){ return mvEdgeNode.size();}
};
}
#endif
