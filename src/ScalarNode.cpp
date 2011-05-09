/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ScalarNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
uint CScalarNode::mnType= NodeType::Scalar;
CScalarNode::CScalarNode()
{
    ;
}
CScalarNode::~CScalarNode()
{
}
void CScalarNode::Message_getScalar()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::getScalar");
}
void CScalarNode::resizeScalar(const uint& res_size)
{
    mvParam.resize(res_size);
}
void CScalarNode::setScalar(const double& val, const uint& index)
{
    mvParam[index]= val;
}
uint CScalarNode::numOfScalarParam()
{
    return mvParam.size();
}
uint CScalarNode::numOfTotalParam()
{
    return mvParam.size();
}
void CScalarNode::resizeVector(const uint& res_size)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::resizeVector");
}
void CScalarNode::setVector(const double& val, const uint& index)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::setVector");
}
vdouble& CScalarNode::getVector()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::getVector");
    return mvParam;
}
uint CScalarNode::numOfVectorParam()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::numOfVectorParam");
    return 0;
}
