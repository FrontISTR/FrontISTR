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

uiint CScalarNode::mnType= NodeType::Scalar;

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

// Scalar Parameter
//
////void CScalarNode::resizeScalar(const uiint& res_size)
////{
////    mvParam.resize(res_size);
////}
////
////void CScalarNode::setScalar(const double& val, const uiint& index)
////{
////    mvParam[index]= val;
////}
////
////uiint CScalarNode::numOfScalarParam()
////{
////    return mvParam.size();
////}
// スカラー,ベクトル 総自由度
//
////uiint CScalarNode::numOfTotalParam()
////{
////    return mvParam.size();
////}
// Vector Parameter
//
////void CScalarNode::resizeVector(const uiint& res_size)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::resizeVector");
////}
////void CScalarNode::setVector(const double& val, const uiint& index)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::setVector");
////}
////
////double& CScalarNode::getVector(const uiint& i)
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::getVector");
////
////    return mvParam[i];
////}
////
////
////uiint CScalarNode::numOfVectorParam()
////{
////    Utility::CLogger *pLogger = Utility::CLogger::Instance();
////    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarNode::numOfVectorParam");
////
////    return 0;
////}






