/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ScalarVectorNode.cxx
|
|                     Written by T.Takeda,    2010/06/01
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
uint CScalarVectorNode::mnType= NodeType::ScalarVector;
CScalarVectorNode::CScalarVectorNode()
{
    ;
}
CScalarVectorNode::~CScalarVectorNode()
{
}
void CScalarVectorNode::Message_getScalar()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarVectorNode::getScalar");
}
void CScalarVectorNode::resizeScalar(const uint& res_size)
{
    mvScalarParam.resize(res_size);
}
void CScalarVectorNode::resizeVector(const uint& res_size)
{
    mvVectorParam.resize(res_size);
}
void CScalarVectorNode::setScalar(const double& val, const uint& index)
{
    mvScalarParam[index]= val;
}
void CScalarVectorNode::setVector(const double& val, const uint& index)
{
    mvVectorParam[index]= val;
}
uint CScalarVectorNode::numOfScalarParam()
{
    return mvScalarParam.size();
}
uint CScalarVectorNode::numOfVectorParam()
{
    return mvVectorParam.size();
}
uint CScalarVectorNode::numOfTotalParam()
{
    return mvScalarParam.size() + mvVectorParam.size();
}
