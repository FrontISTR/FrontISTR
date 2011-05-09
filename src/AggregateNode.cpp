/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   AggregateNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
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
    uint i;
    bool bCheck(false);
    for(i=0; i< mvNode.size(); i++){
        mvpNode= mvNode[i];
        if(mvpNode->getID()==pNode->getID()) bCheck=true;
    };
    if(!bCheck) mvNode.push_back(pNode);
}
