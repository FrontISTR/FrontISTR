/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   SkinFace.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "TypeDef.h"
#include "ContactNode.h"
#include "ElementType.h"
#include "ElementProperty.h"
#include "Node.h"
#include "Element.h" 
#include "EdgeTree.h"
#include <utility>
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
    uint mID;
    uint mRank;
    uint mLevel;
    bool mbSelfDom;
    uint mMeshID;
    uint mElementID;    
    uint mElementFaceID;
    vector<CContactNode*> mvConNode;
    vector<CContactNode*> mvEdgeNode;
    vector<bool>          mvbEdgeMarking;
    vector<CSkinFace*>    mvEdgeFace;
    PairConNode           mPairConNode;
    uint& getEdgeIndex(PairConNode& pairConNode);
    CContactNode *mpFaceNode;
    vector<CSkinFace*> mvProgFace;
    uint mNumOfEdge;
    uint mShapeType;
    virtual CSkinFace* generateFace();
    CSkinFace* mpOtherFace;
    void setupNodeID_progFace(CElement* pElem, const uint& numOfConNode);
    vdouble mvNormalVector;
    Utility::CLogger* mpLogger;
public:
    virtual const char* getName(){return "SkinFace";}
    void setLevel(const uint& level){ mLevel= level;}
    uint& getLevel(){ return mLevel;}
    void setRank(const uint& rank){ mRank= rank;}
    uint& getRank(){ return mRank;}
    void setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}
    void markingSelf(){ mbSelfDom=true;}
    bool isSelf(){ return mbSelfDom;}   
    void setMeshID(const uint& meshID){ mMeshID= meshID;}
    uint& getMeshID(){ return mMeshID;}
    void  setElementID(const uint& id){ mElementID = id;}
    uint& getElementID(){ return mElementID;}
    void  setFaceID(const uint& faceID){ mElementFaceID= faceID;}
    uint& getFaceID(){ return mElementFaceID;}
    void resizeNode(const uint& size){  mvConNode.resize(size);}
    void addNode(CContactNode* pConNode){ mvConNode.push_back(pConNode);}
    void setNode(CContactNode* pConNode, const uint& index){ mvConNode[index] = pConNode;}
    void setNode(vector<CContactNode*>& vConNode){ mvConNode = vConNode;} 
    void operator=(vector<CContactNode*>& vConNode){ mvConNode = vConNode;}
    uint   getNumOfNode() { return mvConNode.size();}
    vector<CContactNode*>&  getNodes(){ return mvConNode;}
    CContactNode*  getNode(const uint& index){ return mvConNode[index];}
    void setShapeType(const uint& shapeType);
    uint& getShapeType(){ return mShapeType;}
    uint& getNumOfEdge(){ return mNumOfEdge;}
    virtual void addSlaveNode(CContactNode* pConNode);
    virtual CContactNode* getSlaveNode(const uint& index){ return NULL;}
    virtual uint getNumOfSlaveNode(){ return 0;}
    virtual void CalcSlave(const uint& islave, const uint& valType);
    virtual double& getCoef(const uint& slaveID, const uint& ivert);
    PairConNode& getEdgePairNode(const uint& iedge);
    void setEdgeFace(CSkinFace* pFace, PairConNode& pairConNode);
    void setEdgeFace(CSkinFace* pFace, const uint& iedge);       
    void setEdgeConNode(CContactNode* pEdgeConNode, PairConNode& pairConNode);
    void setEdgeConNode(CContactNode* pEdgeConNode, const uint& iedge);       
    bool isEdgeNodeMarking(const uint& iedge){ return mvbEdgeMarking[iedge];}
    bool isEdgeNodeMarking(PairConNode& pairConNode);
    void markingEdgeNode(const uint& iedge){ mvbEdgeMarking[iedge]=true;}
    void markingEdgeNode(PairConNode& pairConNode);
    void setFaceConNode(CContactNode* pFaceConNode){ mpFaceNode= pFaceConNode;}
    CContactNode* getEdgeConNode(const uint& iedge){ return mvEdgeNode[iedge];}
    CContactNode* getEdgeConNode(PairConNode& pairConNode);
    CContactNode* getFaceConNode(){ return mpFaceNode;}
    vector<CSkinFace*>& getProgFace(){ return mvProgFace;}
    CSkinFace* getProgFace(const uint& ivert){ return mvProgFace[ivert];}
    void refine(CElement* pElem, uint& faceID);
    vdouble& CalcNzNormalVector();
    vdouble& getNzNormalVector();
};
#endif	/* _SKINFACE_H */
}
