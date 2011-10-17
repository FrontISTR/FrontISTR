/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/Element.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
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
#include "EdgeTree.h"    
#include "FaceTree.h"    
#include "EdgeFaceTree.h"
#include "NodeConnectNodeTree.h"
#include "Logger.h"
typedef std::pair<pmw::CNode*,pmw::CNode*> PairNode;
namespace pmw{
class CElement 
{
public:
    CElement(void);
    virtual ~CElement(void);
protected:
    uiint mnID;
    vector<CNode*> mvNode;
    map<uiint, uiint, less<uiint> > mmIDLocal;
    uiint mMGLevel;
    vector<CNode*> mvEdgeInterNode;
    vector<vector<CElement*> > mvvEdgeElement;
    bool*   mvb_edge;
    vector<CElement*> mvFaceElement;
    vector<CNode*> mvFaceNode;
    bool*   mvb_face;  
    CNode* mpCenterNode;
    virtual bool IndexCheck(const uiint& propType, const uiint& index, string& method_name)=0;
    vuint getLocalNodeNum(CNode* pNode0, CNode* pNode1);
    virtual uiint& getLocalFaceNum(const vuint& vLocalNodeNum)=0;
    vector<CElement*> mvProgElement;
    bool mbComm; 
    bool mbDComm;
    bool mbRComm;
    int mCommID; 
    bool mbMaster;
    bool mbSlave;
    bool* mvbMPCFace;
    bool mbCommMesh2;   
    bool* mvbCommEntity;
    uiint mnCommEntity;  
public:
    virtual void initialize()=0;
    void setID(const uiint& id){ mnID = id;}
    uiint& getID(){ return mnID;}
    void setProgElem(CElement* pProgElem, const uiint& ivert){ mvProgElement[ivert]= pProgElem;}
    CElement* getProgElem(const uiint& ivert){ return mvProgElement[ivert];}
    CElement* getProgElem_NodeID(const uiint& nodeID);
    uiint& getLocalVertNum(const uiint& nodeID){ return mmIDLocal[nodeID]; }
    void interCommElem(){ mbComm= true;}    
    bool& isInterCommElem(){ return mbComm;}
    void setCommID(const int& comID){ mCommID= comID;}
    int& getCommID(){ return mCommID;}
    void interDCommElem(){ mbDComm= true;}
    bool& isInterDCommElem(){ return mbDComm;}
    void interRCommElem(){ mbRComm= true;}
    bool& isInterRCommElem(){ return mbRComm;}
    void setMGLevel(const uiint& mgLevel){ mMGLevel= mgLevel;}
    uiint& getMGLevel(){ return mMGLevel;}
    void  setNode(CNode* pNode,const uiint& local_id);
    vector<CNode*>& getNode() { return mvNode;}
    CNode* getNode(const uiint& local_id){ return mvNode[local_id];}
    void ReInit_IDLocal();
    virtual const uiint& getType()=0;
    virtual const uiint& getOrder()=0;
    virtual const uiint& getNumOfFace()=0;
    virtual const uiint& getNumOfEdge()=0;
    virtual const uiint& getNumOfNode()=0;
    virtual const uiint& getNumOfVert()=0;
    virtual const uiint& getEntityType()=0;
    virtual PairNode getPairNode(const uiint& iedge)=0;
    virtual void     getPairNode(vint& pairNodeIndex, const uiint& iedge)=0;
    void reserveEdgeElement(const uiint& edgeIndex, const uiint& numOfElem);
    void setEdgeElement(const uiint& edgeIndex, CElement* pElem);
    void setEdgeAggElement(const uiint& edgeIndex, vector<CElement*> vElement);
    vector<CElement*>& getEdgeElement(const uiint& edgeIndex){ return mvvEdgeElement[edgeIndex];}
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1)=0;
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1)=0;
    virtual uiint& getEdgeIndex(CNode* pNode0, CNode* pNode1)=0;
    virtual uiint& getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1)=0;
    void setEdgeInterNode(CNode* pNode, const uiint& edgeIndex);
    CNode* getEdgeInterNode(const uiint& edgeIndex){ return mvEdgeInterNode[edgeIndex];}
    uiint getEdgeInterNodeSize(){ return mvEdgeInterNode.size();}
    virtual void replaseEdgeNode(){ ;}
    void setFaceElement(CElement* pElem, const uiint& faceIndex){ mvFaceElement[faceIndex]= pElem;}
    virtual vector<CNode*> getFaceCnvNodes(const uiint& iface)=0;
    void setFaceNode(CNode* pNode, const uiint& faceIndex){ mvFaceNode[faceIndex]= pNode;}
    CNode* getFaceNode(const uiint& faceIndex){ return mvFaceNode[faceIndex];}
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);     
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    uiint getFaceNodeSize(){ return mvFaceNode.size();}
    virtual uiint& getFaceIndex(const uiint& edge0, const uiint& edge1)=0;       
    virtual uiint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2)=0;
    void setVolumeNode(CNode* pNode){ mpCenterNode= pNode;}
    CNode* getVolumeNode(){return mpCenterNode;}
    virtual vector<CNode*> getConnectNode(CNode* pNode)=0;
    void markingMPCMaster(){ mbMaster= true;}
    bool& isMPCMaster(){ return mbMaster;}
    void markingMPCSlave(){ mbSlave= true;}
    bool& isMPCSlave(){ return mbSlave;}
    void markingMPCFace(const uiint& iface){ mvbMPCFace[iface]= true;}
    bool isMPCFace(const uiint& iface){ return mvbMPCFace[iface];}
    void markingCommMesh2(){ mbCommMesh2= true;}
    bool& isCommMesh2(){ return mbCommMesh2;}
    void markingCommEntity(const uiint& ient){ mnCommEntity= ient;  mvbCommEntity[ient]= true;}
    bool isCommEntity(const uiint& ient){ return mvbCommEntity[ient];}
    uiint& getCommEntityID(){ return  mnCommEntity;}
    virtual void deleteProgData()=0;
};
}
#endif
