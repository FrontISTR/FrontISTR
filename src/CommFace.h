/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommFace.h
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
#include <utility>
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
    uint mID;  
    uint mMGLevel;
    uint mMeshID;       
    uint mElementID;    
    uint mElementFaceID;
    uint mNumOfEdge;
    uint mFaceType; 
    vector<CCommNode*>  mvCommNode;    
    vector<CCommNode*>  mvEdgeCommNode;
    CCommNode*          mpFaceCommNode;
    vector<CCommFace*>  mvEdgeCommFace;
    vector<bool>  mvbEdgeMarking;
    PairCommNode        mPairCommNode;
    uint& getEdgeIndex(PairCommNode& pairCommNode);
    vector<CCommFace*> mvProgCommFace;
public:
    void initialize(const uint& numOfVert, const uint& numOfEdge);
    void  setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}
    uint& getType(){ return mFaceType;}
    uint& getNumOfEdge(){ return mNumOfEdge;}
    void setMeshID(const uint& meshID){ mMeshID= meshID;}
    void setElementID(const uint& elemID){ mElementID= elemID;}
    void setElementFaceID(const uint& iface){ mElementFaceID= iface;}
    uint& getMeshID(){ return mMeshID;}
    uint& getElementID(){ return mElementID;}
    uint& getElementFaceID(){ return mElementFaceID;}
    void setVertCommNode(const uint& ivert, CCommNode* pCommNode){ mvCommNode[ivert]= pCommNode;}
    uint getVertCommNodeSize(){ return mvCommNode.size();}
    CCommNode* getVertCommNode(const uint& ivert){ return mvCommNode[ivert];}
    vector<CCommNode*>& getVertCommNode(){ return mvCommNode;}
    PairCommNode& getEdgePairCommNode(const uint& iedge);
    CCommNode* getEdgeCommNode(const uint& iedge){  return  mvEdgeCommNode[iedge];}
    void setEdgeCommFace(CCommFace* pNeibFace, const uint& iedge);
    void setEdgeCommFace(CCommFace* pNeibFace, PairCommNode& pairCommNode);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, const uint& iedge);
    void setEdgeCommNode(CCommNode* pEdgeCommNode, PairCommNode& pairCommNode);
    void markingEdgeNode(const uint& iedge);
    void markingEdgeNode(PairCommNode& pairCommNode);
    bool isEdgeNodeMarking(const uint& iedge){ return mvbEdgeMarking[iedge];}
    void setFaceCommNode(CCommNode* pFaceCommNode){ mpFaceCommNode= pFaceCommNode;}
    CCommNode* getFaceCommNode(){ return mpFaceCommNode;}
    vector<CCommFace*>& refine(CElement* pElement);
};
#endif	/* _COMMFACE_H */
}
