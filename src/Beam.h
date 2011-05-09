/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Beam.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef BEAM_HH_785972AE_015C_439f_AA7A_0464A768D304
#define BEAM_HH_785972AE_015C_439f_AA7A_0464A768D304
#include "BeamBase.h"
namespace pmw{
class CBeam:public CBeamBase
{
public:
    CBeam(void);
    virtual ~CBeam(void);
private:
    static uint mnElemType;
    static uint mnElemType2;
    static uint mNumOfFace;
    static uint mNumOfEdge;
    static uint mNumOfNode;
    uint mDummy;
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
    virtual bool isFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);   
    virtual void setBoolFaceElem(CNode* pNode0, CNode* pNode1, CNode* pNode2);
    virtual void setupFaceCnvNodes();
    virtual vector<CNode*> getConnectNode(CNode* pNode);
};
}
#endif
