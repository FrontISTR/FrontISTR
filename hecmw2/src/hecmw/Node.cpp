/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/Node.cpp
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
#include <vector>
#include "Node.h"
using namespace pmw;
CNode::CNode(void)
{
    ////mbSComm = false;
}
CNode::~CNode(void)
{
    ;
}
//--
// Levelサイズ分のParentNodeを確保
//--
void CNode::resizeGridLevel(const uiint& nNumOfLevel)
{
    mvParentNode.resize(nNumOfLevel);
}

////void CNode::markingSCommNode()
////{
////    mbSComm = true;
////}
////bool CNode::isSCommNode()
////{
////    return mbSComm;
////}
////void CNode::addRank(uiint myRank, uiint transRank)
////{
////    pair<uiint,uiint> pairRank;
////
////    pairRank.first = myRank;
////    pairRank.second= transRank;
////
////    mvPairRank.push_back(pairRank);
////}
//////
////// Nodeが所属(Affiliation)するCommMesh数
//////
////uiint CNode::getNumOfCommMesh()
////{
////    return mvPairRank.size();
////}
//////
////// Nodeが所属するCommMesh2のRank
//////
////pair<uiint,uiint> CNode::getPairRank(const uiint& index)
////{
////    return mvPairRank[index];
////}




