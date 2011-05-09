//
//  BoundaryNode.cpp
//
//  境界Node
//
//                          2009.05.13
//                          2009.05.13
//                          k.Takeda
#include "BoundaryNode.h"
using namespace pmw;

CBoundaryNode::CBoundaryNode()
{
    ;
}

CBoundaryNode::~CBoundaryNode()
{
//    //debug
//    cout << "~CBoundaryNode" << endl;
}


// 初期化
// Type別のmvValueの領域確保
//
void CBoundaryNode::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;

    switch(mnType){
        case(BoundaryTypeNode::Disp)://変位
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Accel)://加速度
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Load)://荷重
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Temp)://温度
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Velo)://速度
            mvValue.resize(dof);
            break;
        case(BoundaryTypeNode::Thermal_Flux)://集中熱流束
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        default:
            break;
    }
}



