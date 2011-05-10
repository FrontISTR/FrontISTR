//
//  ScalarVectorNode.cpp
//
//
//
//                      2009.05.26
//                      2009.05.26
//                      k.Takeda
#include "ScalarVectorNode.h"
using namespace pmw;

#include "Logger.h"

uint CScalarVectorNode::mnType= NodeType::ScalarVector;

// construct & destruct
//
CScalarVectorNode::CScalarVectorNode()
{
    ;
}
CScalarVectorNode::~CScalarVectorNode()
{
//    //debug
//    cout << "~CScalarVectorNode" << endl;
}




void CScalarVectorNode::Message_getScalar()
{
    Utility::CLogger *pLogger = Utility::CLogger::Instance();
    pLogger->Info(Utility::LoggerMode::Warn,"invalid method, at ScalarVectorNode::getScalar");
}


// Parameter Reserve
//

//void CScalarVectorNode::reserveScalar(const uint& res_size)
//{
//    mvScalarParam.reserve(res_size);
//}

void CScalarVectorNode::resizeScalar(const uint& res_size)
{
    mvScalarParam.resize(res_size);
}

//void CScalarVectorNode::reserveVector(const uint& res_size)
//{
//    mvVectorParam.reserve(res_size);
//}

void CScalarVectorNode::resizeVector(const uint& res_size)
{
    mvVectorParam.resize(res_size);
}



// Parameter setting
//

//void CScalarVectorNode::setScalar(const double& val)
//{
//    mvScalarParam.push_back(val);
//}

void CScalarVectorNode::setScalar(const double& val, const uint& index)
{
    mvScalarParam[index]= val;
}

//void CScalarVectorNode::setVector(const double& val)
//{
//    mvVectorParam.push_back(val);
//}

//void CScalarVectorNode::setVector(const vdouble& vVal)
//{
//    mvVectorParam = vVal;
//}

void CScalarVectorNode::setVector(const double& val, const uint& index)
{
    mvVectorParam[index]= val;
}

//
// get Method -> *.h
//

// Count for Parameter
//
uint CScalarVectorNode::numOfScalarParam()
{
    return mvScalarParam.size();
}

uint CScalarVectorNode::numOfVectorParam()
{
    return mvVectorParam.size();
}


// スカラー,ベクトル 総自由度
//
uint CScalarVectorNode::numOfTotalParam()
{
    return mvScalarParam.size() + mvVectorParam.size();
}


























