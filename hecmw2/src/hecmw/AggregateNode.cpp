/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/AggregateNode.cpp
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
#include "HEC_MPI.h"
#include <vector>
#include "AggregateNode.h"
using namespace pmw;
CAggregateNode::CAggregateNode()
{
    ;
}
CAggregateNode::~CAggregateNode()
{
    ;
}
void CAggregateNode::setNode(CNode* pNode)
{
    if(mvNode.size()==0) mvNode.push_back(pNode);
    CNode* mvpNode;
    uiint i;
    bool bCheck(false);
    for(i=0; i< mvNode.size(); i++){
        mvpNode= mvNode[i];
        if(mvpNode->getID()==pNode->getID()) bCheck=true;
    };
    if(!bCheck) mvNode.push_back(pNode);
}
