//
//  VectorNode.cpp
//
//
//
//                      2009.05.26
//                      2009.05.26
//                      k.Takeda
#include "VectorNode.h"
using namespace pmw;

#include "Logger.h"

uint CVectorNode::mnType = NodeType::Vector;

// コンストラクター&デストラクター
//
CVectorNode::CVectorNode()
{
    ;
}
CVectorNode::~CVectorNode()
{
//    //debug
//    cout << "~CVectorNode" << endl;
}





// Vector Parameter
//
//void CVectorNode::reserveVector(const uint& res_size)
//{
//    mvParam.reserve(res_size);
//}

void CVectorNode::resizeVector(const uint& res_size)
{
    mvParam.resize(res_size);
}


//void CVectorNode::setVector(const vdouble& vVal)
//{
//    mvParam = vVal;
//}

//void CVectorNode::setVector(const double& val)
//{
//    mvParam.push_back(val);
//}

void CVectorNode::setVector(const double& val, const uint& index)
{
    mvParam[index] = val;
}

uint CVectorNode::numOfVectorParam()
{
    return mvParam.size();
}

// スカラー,ベクトル 総自由度
//
uint CVectorNode::numOfTotalParam()
{
    return mvParam.size();
}

// Scalar Parameter
//

//void CVectorNode::reserveScalar(const uint& res_size)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::reserveScalar");
//}

void CVectorNode::resizeScalar(const uint& res_size)
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::resizeScalar");
}

//void CVectorNode::setScalar(const double& val)
//{
//    Utility::CLogger *pLogger = Utility::CLogger::Instance();
//    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::setScalar");
//}

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












