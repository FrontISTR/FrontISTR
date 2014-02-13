/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Node.h
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
#include "CommonStd.h"
#include "TypeDef.h"
#include "NodeType.h"
#include "Logger.h"
#include "Vertex.h"
namespace pmw
{
#ifndef NODE_HH_
#define NODE_HH_
class CNode:public CVertex
{
public:
    CNode();
    virtual ~CNode();
protected:
    vector<vector<CNode*> > mvParentNode;//---- Level別 コースグリッド・ノード(親ノード)

////    bool mbSComm;
////    vector<pair<uiint,uiint> > mvPairRank;//--- 通信テーブル毎のペアランク:1st=myRank,2nd=transRank

public:
    virtual uiint& getType()=0;

////    void markingSCommNode();
////    bool isSCommNode();
////
////    void addRank(uiint myRank, uiint transRank);
////    uiint getNumOfCommMesh();// Nodeが所属(Affiliation)するCommMesh数
////    pair<uiint,uiint> getPairRank(const uiint& index);

    virtual void setScalarDOF(const uiint& nNDOF)=0;
    virtual void setVectorDOF(const uiint& nNDOF)=0;
    virtual uiint& getScalarDOF()=0;
    virtual uiint& getVectorDOF()=0;
    virtual uiint getTotalDOF()=0;

    vector<CNode*>& getParentNode(const uiint& iLevel) {
        return mvParentNode[iLevel];
    }
    CNode* getParentNode(const uiint& iLevel, const uiint& index) {
        return mvParentNode[iLevel][index];
    }

    void resizeGridLevel(const uiint& nNumOfLevel);//---- Levelサイズ分のParentNodeを確保
    void reserveParentNode(const uiint& iLevel, const uiint& res_size) {
        mvParentNode[iLevel].reserve(res_size);
    }
    void addParentNode(const uiint& iLevel, CNode* pNode) {
        mvParentNode[iLevel].push_back(pNode);
    }
    uiint getNumOfParentNode(const uiint& iLevel) {
        return mvParentNode[iLevel].size();
    }

};
#endif
}
