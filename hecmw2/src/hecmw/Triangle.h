/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Triangle.h
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
#ifndef TRIANGLE_HH_A3EF00A3_5876_4b0b_A196_DB3BBE2795A5
#define TRIANGLE_HH_A3EF00A3_5876_4b0b_A196_DB3BBE2795A5
#include "ShellBase.h"
namespace pmw
{
class CTriangle:public CShellBase
{
public:
    CTriangle(void);
    virtual ~CTriangle(void);
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
