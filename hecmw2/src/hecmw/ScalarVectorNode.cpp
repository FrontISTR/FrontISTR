/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ScalarVectorNode.cpp
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
#include "ScalarVectorNode.h"
using namespace pmw;
#include "Logger.h"
uiint CScalarVectorNode::mnType= NodeType::ScalarVector;
CScalarVectorNode::CScalarVectorNode()
{
    ;
}
CScalarVectorNode::~CScalarVectorNode()
{
}
void CScalarVectorNode::setScalarDOF(const uiint& nNDOF)
{
    mnNumOfSDOF= nNDOF;
}
void CScalarVectorNode::setVectorDOF(const uiint& nNDOF)
{
    mnNumOfVDOF= nNDOF;
}
uiint& CScalarVectorNode::getScalarDOF()
{
    return mnNumOfSDOF;
}
uiint& CScalarVectorNode::getVectorDOF()
{
    return mnNumOfVDOF;
}
uiint CScalarVectorNode::getTotalDOF()
{
    return mnNumOfSDOF + mnNumOfVDOF;
}
