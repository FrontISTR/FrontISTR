/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/BndVertex.cpp
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
void CBndVertex::deleteAggregate()
{
    vuint vTemp;
    vector<uiint>(vTemp).swap(mvAggElementID);
    map<uiint, uiint, less<uiint> > mTemp;
    map<uiint, uiint, less<uiint> >(mTemp).swap(mmNeibElemVertNum);
}
