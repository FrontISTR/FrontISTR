/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Prism.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _PRISM_H_4a9dd9b8_e22d_43b7_b22e_67c0612b96ba
#define	_PRISM_H_4a9dd9b8_e22d_43b7_b22e_67c0612b96ba
#include "SolidBase.h"
namespace pmw{
class CPrism:public CSolidBase
{
public:
    CPrism();
    virtual ~CPrism();
private:
    static uint mnElemType;
    static uint mnElemType2;
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;
protected:
    virtual bool IndexCheck(const uint& propType, const uint& index, string& method_name);
    virtual uint& getLocalFaceNum(const vuint& vLocalNodeNum);
public:
    virtual const uint& getType();
    virtual const uint& getNumOfFace(){ return mNumOfFace;}
    virtual const uint& getNumOfEdge(){ return mNumOfEdge;}
    virtual const uint& getNumOfNode(){ return mNumOfNode;}
    virtual PairNode& getPairNode(const uint& edgeIndex);
    virtual void getPairNode(vint& pairNodeIndex, const uint& edgeIndex);
    virtual bool isEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual void setBoolEdgeElem(CNode* pNode0, CNode* pNode1);
    virtual uint& getEdgeIndex(CNode* pNode0, CNode* pNode1);
    virtual uint& getEdgeIndex(const uint& nodeID_0, const uint& nodeID_1);
    virtual uint& getFaceIndex(const uint& edge0, const uint& edge1);
    virtual uint& getFaceIndex(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    virtual void setupFaceCnvNodes();
    virtual vector<CNode*> getConnectNode(CNode* pNode);
};
}
#endif	/* _PRISM_H */
