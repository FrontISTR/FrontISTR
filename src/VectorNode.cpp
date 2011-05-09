/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   VectorNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
uint CVectorNode::mnType = NodeType::Vector;
CVectorNode::CVectorNode()
{
    ;
}
CVectorNode::~CVectorNode()
{
}
void CVectorNode::resizeVector(const uint& res_size)
{
    mvParam.resize(res_size);
}
void CVectorNode::setVector(const double& val, const uint& index)
{
    mvParam[index] = val;
}
uint CVectorNode::numOfVectorParam()
{
    return mvParam.size();
}
uint CVectorNode::numOfTotalParam()
{
    return mvParam.size();
}
void CVectorNode::resizeScalar(const uint& res_size)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::resizeScalar");
}
void CVectorNode::setScalar(const double& val, const uint& index)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::setScalar");
}
double& CVectorNode::getScalar(const uint& i)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::getScalar");
    return mvParam[0];
}
uint CVectorNode::numOfScalarParam()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::numOfScalarParam");
    return 0;
}
