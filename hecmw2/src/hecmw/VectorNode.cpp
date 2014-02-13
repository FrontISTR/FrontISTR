/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/VectorNode.cpp
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
#include "VectorNode.h"
using namespace pmw;
#include "Logger.h"
uiint CVectorNode::mnType = NodeType::Vector;
CVectorNode::CVectorNode()
{
    ;
}
CVectorNode::~CVectorNode()
{
}
void CVectorNode::setScalarDOF(const uiint& nNDOF)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method  VectorNode::setScalarDOF");
}
void CVectorNode::setVectorDOF(const uiint& nNDOF)
{
    mnNumOfDOF= nNDOF;
}
uiint& CVectorNode::getScalarDOF()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->setUDummyValue(0);
    return pLogger->getUDummyValue();
}
uiint& CVectorNode::getVectorDOF()
{
    return mnNumOfDOF;
}
uiint CVectorNode::getTotalDOF()
{
    return mnNumOfDOF;
}
