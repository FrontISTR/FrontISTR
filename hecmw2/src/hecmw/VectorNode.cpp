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

uiint CVectorNode::mnType = NodeType::Vector;

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



// Vector Parameter
//
////void CVectorNode::resizeVector(const uiint& res_size)
////{
////    mvParam.resize(res_size);
////}
////
////void CVectorNode::setVector(const double& val, const uiint& index)
////{
////    mvParam[index] = val;
////}
////
////uiint CVectorNode::numOfVectorParam()
////{
////    return mvParam.size();
////}
////
////// スカラー,ベクトル 総自由度
//////
////uiint CVectorNode::numOfTotalParam()
////{
////    return mvParam.size();
////}
////
////// Scalar Parameter
//////
////void CVectorNode::resizeScalar(const uiint& res_size)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::resizeScalar");
////}
////
////void CVectorNode::setScalar(const double& val, const uiint& index)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::setScalar");
////}
////
////double& CVectorNode::getScalar(const uiint& i)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::getScalar");
////
////    return mvParam[0];
////}
////
////uiint CVectorNode::numOfScalarParam()
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at VectorNode::numOfScalarParam");
////
////    return 0;
////}












