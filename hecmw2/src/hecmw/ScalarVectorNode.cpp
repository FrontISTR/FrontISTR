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

uiint CScalarVectorNode::mnType= NodeType::ScalarVector;

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

void CScalarVectorNode::setScalarDOF(const uiint& nNDOF)
{
    mnNumOfSDOF= nNDOF;
}
void CScalarVectorNode::setVectorDOF(const uiint& nNDOF)
{
    mnNumOfVDOF= nNDOF;
}
uiint& CScalarVectorNode::getScalarDOF()
{
    return mnNumOfSDOF;
}
uiint& CScalarVectorNode::getVectorDOF()
{
    return mnNumOfVDOF;
}
uiint CScalarVectorNode::getTotalDOF()
{
    return mnNumOfSDOF + mnNumOfVDOF;
}

//
// Parameter Reserve
//
////void CScalarVectorNode::resizeScalar(const uiint& res_size)
////{
////    mvScalarParam.resize(res_size);
////}
////
////void CScalarVectorNode::resizeVector(const uiint& res_size)
////{
////    mvVectorParam.resize(res_size);
////}
//
// Parameter setting
//
////void CScalarVectorNode::setScalar(const double& val, const uiint& index)
////{
////    mvScalarParam[index]= val;
////}
////
////void CScalarVectorNode::setVector(const double& val, const uiint& index)
////{
////    mvVectorParam[index]= val;
////}
////
//
// Count for Parameter
//
////uiint CScalarVectorNode::numOfScalarParam()
////{
////    return mvScalarParam.size();
////}
////
////uiint CScalarVectorNode::numOfVectorParam()
////{
////    return mvVectorParam.size();
////}
////
////
////// スカラー,ベクトル 総自由度
//////
////uiint CScalarVectorNode::numOfTotalParam()
////{
////    return mvScalarParam.size() + mvVectorParam.size();
////}


























