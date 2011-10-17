/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/VectorNode.cpp
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
