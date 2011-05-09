/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   CommMesh2.h
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
#include <map>
#include "CommNode.h"
#include "CommFace.h"
namespace pmw{
#ifndef _COMM_MESH2_H
#define	_COMM_MESH2_H
class CCommMesh2{
public:
    CCommMesh2();
    virtual ~CCommMesh2();
protected:
    uint mID;
    uint myRank;
    uint mTransmitRank;
    uint mMGLevel;
    vector<CCommNode*> mvCommNode;    
    vector<CCommNode*> mvEdgeCommNode;
    vector<CCommNode*> mvFaceCommNode;
    vector<CCommFace*> mvCommFace;    
    map<uint, uint, less<uint> > mmCommNodeID2Index;
    map<uint, uint, less<uint> > mmCommFaceID2Index;
public:
    void  setID(const uint& id){ mID= id;}
    uint& getID(){ return mID;}
    void setLevel(const uint& mgLevel){ mMGLevel= mgLevel;}
    uint& getLevel(){ return mMGLevel;}
    void  setRank(const uint& rank){ myRank= rank;}
    void  setTransmitRank(const uint& trans_rank){ mTransmitRank= trans_rank;}
    uint& getRank(){ return myRank;}
    uint& getTrasmitRank(){ return mTransmitRank;}
public:
    void reserveCommNode(const uint& res_size){ mvCommNode.reserve(res_size);}
    void addCommNode(CCommNode* pCommNode);
    void reserveCommFace(const uint& res_size){ mvCommFace.reserve(res_size);}
    void addCommFace(CCommFace* pCommFace);
    CCommNode* getCommVertNode(const uint& id);
    CCommNode* getCommVertNodeIX(const uint& index){ return mvCommNode[index];}
    CCommNode* getCommEdgeNodeIX(const uint& index){ return mvEdgeCommNode[index];}
    CCommNode* getCommFaceNodeIX(const uint& index){ return mvFaceCommNode[index];}
    CCommFace* getCommFace(const uint& id);
    CCommFace* getCommFaceIX(const uint& index){ return mvCommFace[index];}
    uint getCommVertNodeSize(){ return mvCommNode.size();}
    uint getCommFaceSize(){ return mvCommFace.size();}
public:
    void setupAggFace();
    void setupEdgeCommNode(CCommMesh2* pProgCommMesh);
    void setupFaceCommNode(CCommMesh2* pProgCommMesh);
    void setupVertCommNode(CCommMesh2* pProgCommMesh);
};
#endif	/* _COMMSURFACE_H */
}
