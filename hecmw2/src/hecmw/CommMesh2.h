/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/CommMesh2.h
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
#include "TypeDef.h"
#include <map>
#include "CommNode.h"
#include "CommFace.h"
namespace pmw
{
#ifndef _COMM_MESH2_H
#define	_COMM_MESH2_H
class CCommMesh2
{
public:
    CCommMesh2();
    virtual ~CCommMesh2();
protected:
    uiint mID;
    uiint myRank;
    uiint mTransmitRank;
    uiint mMGLevel;
    vector<CCommNode*> mvCommNode;
    uiint mEdgeCount;
    vector<CCommFace*> mvCommFace;
    map<uiint, uiint, less<uiint> > mmCommNodeID2Index;
    map<uiint, uiint, less<uiint> > mmCommFaceID2Index;
public:
    void  setID(const uiint& id) {
        mID= id;
    }
    uiint& getID() {
        return mID;
    }
    void setLevel(const uiint& mgLevel) {
        mMGLevel= mgLevel;
    }
    uiint& getLevel() {
        return mMGLevel;
    }
    void  setRank(const uiint& rank) {
        myRank= rank;
    }
    void  setTransmitRank(const uiint& trans_rank) {
        mTransmitRank= trans_rank;
    }
    uiint& getRank() {
        return myRank;
    }
    uiint& getTrasmitRank() {
        return mTransmitRank;
    }
public:
    void reserveCommNode(const uiint& res_size) {
        mvCommNode.reserve(res_size);
    }
    void addCommNode(CCommNode* pCommNode);
    void reserveCommFace(const uiint& res_size) {
        mvCommFace.reserve(res_size);
    }
    void addCommFace(CCommFace* pCommFace);
    CCommNode* getCommNode(const uiint& id);
    CCommNode* getCommNodeIX(const uiint& index) {
        return mvCommNode[index];
    }
    CCommFace* getCommFace(const uiint& id);
    CCommFace* getCommFaceIX(const uiint& index) {
        return mvCommFace[index];
    }
    uiint getCommNodeSize() {
        return mvCommNode.size();
    }
    uiint getCommFaceSize() {
        return mvCommFace.size();
    }
public:
    void setupAggFace();
    void setupEdgeCommNode(CCommMesh2* pProgCommMesh, const uiint& nLevel);
    void setupFaceCommNode(CCommMesh2* pProgCommMesh);
    void setupCommNode(CCommMesh2* pProgCommMesh);
    void deleteProgData();
};
#endif	/* _COMM_MESH2_H */
}
