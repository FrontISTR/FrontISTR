//
//  ScalarNode.cpp
//
//
//
//                      2009.05.26
//                      2009.05.26
//                      k.Takeda
#include "ScalarNode.h"
using namespace pmw;

#include "Logger.h"

uint CScalarNode::mnType= NodeType::Scalar;

// construct & destruct
//
CScalarNode::CScalarNode()
{
    ;
}
CScalarNode::~CScalarNode()
{
//    //debug
//    cout << "~CScalarNode" << endl;
}


// method
// --
//
//
void CScalarNode::Message_getScalar()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::getScalar");
}

// Scalar Parameter
//

//void CScalarNode::reserveScalar(const uint& res_size)
//{
//    mvParam.reserve(res_size);
//}

void CScalarNode::resizeScalar(const uint& res_size)
{
    mvParam.resize(res_size);
}

//void CScalarNode::setScalar(const double& val)
//{
//    mvParam.push_back(val);
//}


void CScalarNode::setScalar(const double& val, const uint& index)
{
    mvParam[index]= val;
}

uint CScalarNode::numOfScalarParam()
{
    return mvParam.size();
}


// スカラー,ベクトル 総自由度
//
uint CScalarNode::numOfTotalParam()
{
    return mvParam.size();
}


// Vector Parameter
//

//void CScalarNode::reserveVector(const uint& res_size)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::reserveVector");
//}

void CScalarNode::resizeVector(const uint& res_size)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::resizeVector");
}

//void CScalarNode::setVector(const vdouble& vVal)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::setVector");
//}

//void CScalarNode::setVector(const double& val)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::setVector");
//}

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






