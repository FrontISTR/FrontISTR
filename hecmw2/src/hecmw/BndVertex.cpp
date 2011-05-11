//
//  BndVertex.cpp
//
//  座標を持たないVertex
//    Vertexの親
//
//          2010.04.28
//          k.Takeda
#include "BndVertex.h"
using namespace pmw;

CBndVertex::CBndVertex()
{
    ;
}
CBndVertex::~CBndVertex()
{
    ;
}


// AggregateElementの解放
//
// |全てのMesh処理が終わった後で呼び出す|
//
void CBndVertex::deleteAggregate()
{
    vuint vTemp;
    vector<uiint>(vTemp).swap(mvAggElementID);

    map<uiint, uiint, less<uiint> > mTemp;
    map<uiint, uiint, less<uiint> >(mTemp).swap(mmNeibElemVertNum);
}







