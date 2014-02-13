/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ScalarNode.cpp
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
#include "ScalarNode.h"
using namespace pmw;
#include "Logger.h"
uiint CScalarNode::mnType= NodeType::Scalar;
CScalarNode::CScalarNode()
{
    ;
}
CScalarNode::~CScalarNode()
{
}
void CScalarNode::setScalarDOF(const uiint& nNDOF)
{
    mnNumOfDOF= nNDOF;
}
void CScalarNode::setVectorDOF(const uiint& nNDOF)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method  ScalarNode::setVectorDOF");
}
uiint& CScalarNode::getScalarDOF()
{
    return mnNumOfDOF;
}
uiint& CScalarNode::getVectorDOF()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->setUDummyValue(0);
    return pLogger->getUDummyValue();
}
uiint CScalarNode::getTotalDOF()
{
    return mnNumOfDOF;
}
