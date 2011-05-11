
#include <vector>

//
//  AggregateNode.cpp
//
//
//
//                  2009.08.11
//                  2009.08.11
//                  k.Takeda
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

// コア・ノード周囲のノード => Aggregate Node をセットする.
// --
// 既に取得済みかどうかの判定の上でpush_back.
// --
void CAggregateNode::setNode(CNode* pNode)
{
    //まったくセットされていない場合.
    if(mvNode.size()==0) mvNode.push_back(pNode);


    CNode* mvpNode;
    uiint i;
    bool bCheck(false);
    // 一致するIDが無ければ,AggregateNodeとしてpush_backする.
    //
    for(i=0; i< mvNode.size(); i++){
        mvpNode= mvNode[i];
        if(mvpNode->getID()==pNode->getID()) bCheck=true;
    };
    
    if(!bCheck) mvNode.push_back(pNode);
}






