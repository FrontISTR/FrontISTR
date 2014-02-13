/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Mesh.h
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "CommonStd.h"
#include "TypeDef.h"
#include "Logger.h"
#include "Node.h"
#include "Element.h"
#include "BoundaryNodeMesh.h"
#include "BoundaryFaceMesh.h"
#include "BoundaryEdgeMesh.h"
#include "BoundaryVolumeMesh.h"
#include "BoundaryNode.h"
#include "BoundaryFace.h"
#include "BoundaryEdge.h"
#include "BoundaryVolume.h"
#include "BoundaryGroup.h"
#include "BNodeMeshGrp.h"
#include "AggregateElement.h"
#include "AggregateNode.h"
#include <map>
#include "CommunicationMesh.h"
#include "IndexBucket.h"
#include "QuickSort.h"
#include "CommMesh2.h"

#include "SolutionType.h"
#include "CodeType.h"

#include "ElementGroup.h"
#include "HEC_MPI.h"

typedef std::pair<pmw::CNode*,pmw::CNode*> PairNode;
namespace pmw
{
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
    uiint            mnMeshID;
    vector<CNode*>    mvNode;
    vector<CElement*> mvElement;

    uiint mnSolutionType;// FEM, FVM
    uiint mnProp;// CodeType

    uiint mNodeEndIndex;
    uiint mElemEndIndex;

    vector<CAggregateElement*> mvAggElement;
    vector<CAggregateNode*>    mvAggNode;

    uiint mMGLevel;
    uiint mMaxMGLevel;

    vector<CCommMesh*>     mvCommMesh;
    map<uiint, uiint, less<uiint> > mmCommIndex;
    vector<CCommMesh2*> mvCommMesh2;
    map<uiint, uiint, less<uiint> > mmComm2Index;

    CIndexBucket     moBucket;
    uiint mNumOfNode, mNumOfElement;
    uiint maxNodeID,  maxElementID;
    uiint minNodeID,  minElementID;

    CBNodeMeshGrp                       *mpBNodeMeshGrp;
    BoundaryGroup<CBoundaryFaceMesh*>   mGrpBndFaceMesh;
    BoundaryGroup<CBoundaryEdgeMesh*>   mGrpBndEdgeMesh;
    BoundaryGroup<CBoundaryVolumeMesh*> mGrpBndVolumeMesh;

    void setupParentNode(CNode* pNode0, CNode* pNode1, CNode* inNode);   //辺ノード
    void setupParentNode2nd(CNode* pNode0, CNode* pNode1, CNode* inNode);//辺ノード(最終Level 2次要素)
    void setupParentNode(vector<CNode*>& vNode, CNode* inNode);          //面ノード、体積ノード

    bool* mvMarkingBNode;    //-- Node全体での境界マーキング
    bool* mvMarkingDirichlet;//-- Node全体でのディレクレ・マーキング
    bool* mvMarkingNeumann;  //-- Node全体でのノイマン・マーキング


    bool* mvMarkingLargeRankNode; //-- Node全体での"Rank大"の通信ノード・マーキング


    vector<CElementGroup*> mvElementGroup;
    map<uiint, uiint, less<uiint> > mmElemGrpID2IX;
public:
    void setMeshID(const uiint& id) {
        mnMeshID = id;
    }
    uiint& getMeshID() {
        return mnMeshID;
    }

    void setMGLevel(const uiint& mgLevel) {
        mMGLevel= mgLevel;
    }
    uiint& getMGLevel() {
        return mMGLevel;
    }
    void setMaxMGLevel(const uiint& maxLevel) {
        mMaxMGLevel= maxLevel;
    }

    void setSolutionType(const uiint& nSolutionType) {
        mnSolutionType = nSolutionType;
    }

    void  setProp(const uiint& nProp) {
        mnProp = nProp;
    }
    uiint& getProp() {
        return mnProp;
    }

    void initBucketNode(const uiint& max_id, const uiint& min_id);
    void setupBucketNodeIndex(const uiint& id, const uiint& index);
    void setupBucketNode();
    void initBucketElement(const uiint& max_id, const uiint& min_id);
    void setupBucketElementIndex(const uiint& id, const uiint& index);
    void setupBucketElement();
    virtual CIndexBucket* getBucket() {
        return &moBucket;
    }

    uiint& getMaxNodeID() {
        return moBucket.getMaxNodeID();
    }
    uiint& getMaxElementID() {
        return moBucket.getMaxElementID();
    }
    uiint& getMinNodeID() {
        return moBucket.getMinNodeID();
    }
    uiint& getMinElementID() {
        return moBucket.getMinElementID();
    }

    void reserveNode(const uiint& num_of_node);
    void setNode(CNode* pNode);
    CNode* getNode(const uiint& id);
    CNode* getNodeIX(const uiint& index);
    uiint  getNodeSize() {
        return mvNode.size();
    }
    uiint& getNumOfNode() {
        return mNumOfNode;
    }
    void setupNumOfNode() {
        mNumOfNode= mvNode.size();
    }
    uiint& getEfficientNodeSize() {
        return mNodeEndIndex;
    }
    void reserveElement(const uiint& num_of_elem);
    void setElement(CElement* pElement);
    CElement* getElement(const uiint& id);
    CElement* getElementIX(const uiint& index);
    vector<CElement*>& getElements() {
        return mvElement;
    }
    uiint  getElementSize() {
        return mvElement.size();
    }
    uiint& getNumOfElement() {
        return mNumOfElement;
    }
    void setupNumOfElement() {
        mNumOfElement= mvElement.size();
    }
    uiint& getEfficientElementSize() {
        return mElemEndIndex;
    }

    void resizeAggregate(const uiint& res_size);
    void setupAggregate(const uiint& nLevel);
    void setAggNode(CAggregateNode* pAggNode, const uiint& inode);
    void setAggElement(CAggregateElement* pAggElem, const uiint& inode);
    vector<CAggregateElement*>& getAggElems() {
        return mvAggElement;
    }
    vector<CAggregateNode*>&    getAggNodes() {
        return mvAggNode;
    }
    CAggregateElement* getAggElem(const uiint& node_id);
    CAggregateElement* getAggElemIX(const uiint& inode);
    CAggregateNode*    getAggNode(const uiint& node_id);
    CAggregateNode*    getAggNodeIX(const uiint& inode);

    void presetProgMesh(CMesh* pProgMesh);
    void setupEdgeElement(CMesh* pProgMesh, const uiint& nLevel);
    void setupFaceElement(CMesh* pProgMesh);
    void setupFaceElement2(CMesh* pProgMesh);
    void setupVolumeNode(CMesh* pProgMesh);

    void replaceEdgeNode();
protected:
    CNode* GeneInterNode(CNode *pNode);
    void avgCoord(vector<CNode*> vCnvNode, CNode* pNode);
public:
    void setBNodeMeshGrp(CBNodeMeshGrp *pBNodeMeshGrp) {
        mpBNodeMeshGrp= pBNodeMeshGrp;
    }
    CBNodeMeshGrp* getBNodeMeshGrp() {
        return mpBNodeMeshGrp;
    }
    void reserveBndNodeMesh(const uiint& res_size) {
        mpBNodeMeshGrp->reserveBndNodeMesh(res_size);
    }
    void setBndNodeMesh(CBoundaryNodeMesh *pBNodeMesh) {
        mpBNodeMeshGrp->setBndNodeMesh(pBNodeMesh);
    }
    CBoundaryNodeMesh* getBndNodeMeshIX(const uiint& index) {
        return mpBNodeMeshGrp->getBndNodeMeshIX(index);
    }
    CBoundaryNodeMesh* getBndNodeMeshID(const uiint& id) {
        return mpBNodeMeshGrp->getBndNodeMeshID(id);
    }
    uiint getNumOfBoundaryNodeMesh() {
        return mpBNodeMeshGrp->getNumOfBoundaryNodeMesh();
    }

    void reserveBndFaceMesh(const uiint& res_size) {
        mGrpBndFaceMesh.reserve(res_size);
    }
    void setBndFaceMesh(CBoundaryFaceMesh *pBFaceMesh) {
        mGrpBndFaceMesh.push(pBFaceMesh);
    }
    CBoundaryFaceMesh* getBndFaceMeshIX(const uiint& index) {
        return mGrpBndFaceMesh.get_withIndex(index);
    }
    CBoundaryFaceMesh* getBndFaceMeshID(const uiint& id) {
        return mGrpBndFaceMesh.get_withID(id);
    }
    uiint getNumOfBoundaryFaceMesh() {
        return mGrpBndFaceMesh.NumOfBoundary();
    }

    void reserveBndVolumeMesh(const uiint& res_size) {
        mGrpBndVolumeMesh.reserve(res_size);
    }
    void setBndVolumeMesh(CBoundaryVolumeMesh *pBVolMesh) {
        mGrpBndVolumeMesh.push(pBVolMesh);
    }
    CBoundaryVolumeMesh* getBndVolumeMeshIX(const uiint& index) {
        return mGrpBndVolumeMesh.get_withIndex(index);
    }
    CBoundaryVolumeMesh* getBndVolumeMeshID(const uiint& id) {
        return mGrpBndVolumeMesh.get_withID(id);
    }
    uiint getNumOfBoundaryVolumeMesh() {
        return mGrpBndVolumeMesh.NumOfBoundary();
    }

    void reserveBndEdgeMesh(const uiint& res_size) {
        mGrpBndEdgeMesh.reserve(res_size);
    }
    void setBndEdgeMesh(CBoundaryEdgeMesh *pBEdgeMesh) {
        mGrpBndEdgeMesh.push(pBEdgeMesh);
    }
    CBoundaryEdgeMesh* getBndEdgeMeshIX(const uiint& index) {
        return mGrpBndEdgeMesh.get_withIndex(index);
    }
    CBoundaryEdgeMesh* getBndEdgeMeshID(const uiint& id) {
        return mGrpBndEdgeMesh.get_withID(id);
    }
    uiint getNumOfBoundaryEdgeMesh() {
        return mGrpBndEdgeMesh.NumOfBoundary();
    }

    //--
    // 境界Node判定 関数 : MG境界処理
    //--
    bool isBNode(uiint& inode);  //--全Node中のBNode判定
    bool* getBNodeMarkingArray();//--全Node中のBNode判定のbool配列
    bool isDirichletBNode(uiint& inode);  //--全Node中のディレクレBNode判定
    bool* getDirichletBNodeMarkingArray();//--全Node中のディレクレBNode判定のbool配列
    bool isNeumannBNode(uiint& inode);//--全Node中のノイマンBNode判定
    bool* getNeumannBNodeArray();     //--全Node中のノイマンBNode判定のbool配列

    void setupBNodeMarking();//-----------境界条件Node(All,Dirichlet,Neumann) マーキング配列の生成


    void reserveCommMesh(const uiint& res_size) {
        mvCommMesh.reserve(res_size);
    }
    void setCommMesh(CCommMesh* pCommMesh);
    uiint getNumOfCommMesh() {
        return mvCommMesh.size();
    }
    CCommMesh* getCommMesh(const uiint& comID);
    CCommMesh* getCommMeshIX(const uiint& index) {
        return mvCommMesh[index];
    }
    void sortMesh();

    void setCommMesh2(CCommMesh2 *pCommMesh2);
    CCommMesh2* getCommMesh2(const uiint& comID);
    CCommMesh2* getCommMesh2IX(const uiint& index) {
        return mvCommMesh2[index];
    }
    uiint getCommMesh2Size() {
        return mvCommMesh2.size();
    }

    //--
    // Nodeの内、Rankが大の通信ノード(荷重,境界値を与えてはいけないNode)
    // # "通常のNode"と"最小Rankの通信ノード" 以外ということ.
    //--
    void setupLargeRank_CommNodeMarking();//--------- Rank大の通信Nodeマーキング

    bool isLargeRankCommNode(const uiint& inode);//--- 引数:Node全体の通し番号


    void deleteProgData();
    void deleteAggregate_on_Node();

    void addElemGrp(CElementGroup* pElemGrp);
    uiint getNumOfElemGrp();
    CElementGroup* getElemGrpIX(const uiint& index);
    CElementGroup* getElemGrpID(const uiint& nGrpID);

    void clear();

};
#endif
}
