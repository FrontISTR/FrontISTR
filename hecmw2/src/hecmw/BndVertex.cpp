/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/BndVertex.cpp
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
#include "BndVertex.h"
using namespace pmw;
CBndVertex::CBndVertex()
{
    mID= 0;
    mvAggElementID.clear();
}
CBndVertex::~CBndVertex()
{
    ;
}
void CBndVertex::deleteAggregate()
{
    vuint vTemp;
    vector<uiint>(vTemp).swap(mvAggElementID);
    map<uiint, uiint, less<uiint> > mTemp;
    map<uiint, uiint, less<uiint> >(mTemp).swap(mmNeibElemVertNum);
}
