/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Element.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef ELEMENT_HH_05B90007_6581_437f_AD9D_7C9227A8D7B6
#define ELEMENT_HH_05B90007_6581_437f_AD9D_7C9227A8D7B6
#include "CommonStd.h"
#include "TypeDef.h"
#include <map>
#include <utility> 
#include "Node.h"
#include "ElementType.h"
typedef std::pair<pmw::CNode*,pmw::CNode*> PairNode;
#include "EdgeTree.h"    
#include "FaceTree.h"    
#include "EdgeFaceTree.h"
#include "NodeConnectNodeTree.h"
#include "Logger.h"
namespace pmw{
class CElement 
{
public:
    CElement(void);
    virtual ~CElement(void);
protected:
    uint mnID;
    vector<CNode*> mvNode;
    uint mMGLevel;
    vector<vector<CNode*> > mvvFaceCnvNode;
    vuint mvPairNodeLocalNum;
    PairNode mEdgeNode;
    map<uint, uint, less<uint> > mmIDLocal;
    vector<CNode*> mvEdgeInterNode;
    vector<vector<CElement*> > mvvEdgeElement;
    vector<bool>   mvb_edge;
    vector<CElement*> mvFaceElement;
    vector<CNode*> mvFaceNode;
    vector<bool>   mvb_face;  
    CNode* mpCenterNode;
    uint parentID;
    virtual bool IndexCheck(const uint& propType, const uint& index, string& method_name)=0;
    vuint& getLocalNodeNum(CNode* pNode0, CNode* pNode1);
    virtual uint& getLocalFaceNum(const vuint& vLocalNodeNum)=0;
    vector<CElement*> mvProgElement;
    bool mbComm; 
    bool mbDComm;
    bool mbRComm;
    int mCommID; 
    bool mbMaster;
    bool mbSlave;
    vector<bool> mvbMPCFace;
    bool mbCommMesh2;          
    vector<bool> mvbCommEntity;
    uint mnCommEntity;
public:
    void setID(const uint& id){ mnID = id;}
    uint& getID(){ return mnID;}
    void setParentID(const uint& id){ parentID= id;}
    uint& getParentID(){ return parentID;}
    void setProgElem(CElement* pProgElem, const uint& ivert){ mvProgElement[ivert]= pProgElem;}
    CElement* getProgElem(const uint& ivert){ return mvProgElement[ivert];}
    CElement* getProgElem_NodeID(const uint& nodeID);
    uint& getLocalVertNum(const uint& nodeID){ return mmIDLocal[nodeID]; }
    void interCommElem(){ mbComm= true;}    
    bool& isInterCommElem(){ return mbComm;}
    void setCommID(const int& comID){ mCommID= comID;}
    int& getCommID(){ return mCommID;}
    void interDCommElem(){ mbDComm= true;}
    bool& isInterDCommElem(){ return mbDComm;}
    void interRCommElem(){ mbRComm= true;}
    bool& isInterRCommElem(){ return mbRComm;}
    void setMGLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getMGLevel(){ return mMGLevel;}
    void  setNode(CNode* pNode,const uint& local_id);
    vector<CNode*>& getNode() { return mvNode;}
    CNode* getNode(const uint& local_id){ return mvNode[local_id];}
    virtual const uint& getType()=0;
    virtual const uint& getNumOfFace()=0;
    virtual const uint& getNumOfEdge()=0;
    virtual const uint& getNumOfNode()=0;
    virtual const uint& getEntityType()=0;
    virtual PairNode& getPairNode(const uint& edgeIndex)=0;
    virtual void     getPairNode(vint& pairNodeIndex, const uint& edgeIndex)=0;
    void reserveEdgeElement(const uint& edgeIndex, const uint& numOfElem);
    void setEdgeElement(const uint& edgeIndex, CElement* pElem);
    void setEdgeAggElement(const uint& edgeIndex, vector<CElement*> vElement);
    vector<CElement*>& getEdgeElement(const uint& edgeIndex){ return mvvEdgeElement[edgeIndex];}
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1)=0;
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1)=0;
    virtual uint& getEdgeIndex(CNode* pNode0, CNode* pNode1)=0;
    virtual uint& getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1)=0;
    void setEdgeInterNode(CNode* pNode, const uint& edgeIndex);
    CNode* getEdgeInterNode(const uint& edgeIndex){ return mvEdgeInterNode[edgeIndex];}
    void setFaceElement(CElement* pElem, const uint& faceIndex){ mvFaceElement[faceIndex]= pElem;}
    vector<CNode*>& getFaceCnvNodes(const uint& faceIndex){ return  mvvFaceCnvNode[faceIndex];}
    virtual void setupFaceCnvNodes()=0;
    void setFaceNode(CNode* pNode, const uint& faceIndex){ mvFaceNode[faceIndex]= pNode;}
    CNode* getFaceNode(const uint& faceIndex){ return mvFaceNode[faceIndex];}
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);   
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    virtual uint& getFaceIndex(const uint& edge0, const uint& edge1)=0;
    virtual uint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)=0;
    void setVolumeNode(CNode* pNode){ mpCenterNode= pNode;}
    CNode* getVolumeNode(){return mpCenterNode;}
    virtual vector<CNode*> getConnectNode(CNode* pNode)=0;
    void markingMPCMaster(){ mbMaster= true;}
    bool& isMPCMaster(){ return mbMaster;}
    void markingMPCSlave(){ mbSlave= true;}
    bool& isMPCSlave(){ return mbSlave;}
    void markingMPCFace(const uint& iface){ mvbMPCFace[iface]= true;}
    bool isMPCFace(const uint& iface){ return mvbMPCFace[iface];}
    void markingCommMesh2(){ mbCommMesh2= true;}
    bool& isCommMesh2(){ return mbCommMesh2;}
    void markingCommEntity(const uint& ient){ mnCommEntity= ient;  mvbCommEntity[ient]= true;}
    bool isCommEntity(const uint& ient){ return mvbCommEntity[ient];}
    uint& getCommEntityID(){ return  mnCommEntity;}
};
}
#endif
