/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Quad.h
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
#ifndef QUAD_HH_6BF7327A_6740_43e8_ADA1_B90414CCA216
#define QUAD_HH_6BF7327A_6740_43e8_ADA1_B90414CCA216
#include "ShellBase.h"
namespace pmw
{
class CQuad:public CShellBase
{
public:
    CQuad(void);
    virtual ~CQuad(void);
private:
    static uiint mnElemType;
    static uiint mnElemOrder;
    static uiint mNumOfFace;
    static uiint mNumOfEdge;
    static uiint mNumOfNode;
    static uiint mNumOfVert;
    uiint mTempo;
protected:
    virtual bool IndexCheck(const uiint& propType, const uiint& index, string& method_name);
    virtual uiint& getLocalFaceNum(const vuint& vLocalNodeNum);
public:
    virtual void initialize();
public:
    virtual const uiint& getType();
    virtual const uiint& getOrder();
    virtual const uiint& getNumOfFace() {
        return mNumOfFace;
    }
    virtual const uiint& getNumOfEdge() {
        return mNumOfEdge;
    }
    virtual const uiint& getNumOfNode() {
        return mNumOfNode;
    }
    virtual const uiint& getNumOfVert() {
        return mNumOfVert;
    }
    virtual PairNode getPairNode(const uiint& iedge);
    virtual void getPairNode(vint& pairNodeIndex, const uiint& iedge);
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual uiint& getEdgeIndex(CNode* pNode0, CNode* pNode1);
    virtual uiint& getEdgeIndex(const uiint& nodeID_0, const uiint& nodeID_1);
    virtual uiint& getFaceIndex(const uiint& edge0, const uiint& edge1);
    virtual uiint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    virtual vector<CNode*> getFaceCnvNodes(const uiint& iface);
    virtual vector<CNode*> getConnectNode(CNode* pNode);
    virtual void deleteProgData();
};
}
#endif
