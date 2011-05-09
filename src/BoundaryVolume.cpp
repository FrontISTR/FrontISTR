//
//  BoundaryVolume.cpp
//
//
//
//                          2009.05.18
//                          2009.05.13
//                          k.Takeda
#include "BoundaryVolume.h"
using namespace pmw;

CBoundaryVolume::CBoundaryVolume()
{
    ;
}

CBoundaryVolume::~CBoundaryVolume()
{
//    //debug
//    cout << "~CBoundaryVolume" << endl;
}

// 初期化
// Type別のmvValueの領域確保
//
void CBoundaryVolume::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;

    switch(mnType){
        case(BoundaryTypeVolume::Heat)://発熱
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Gravity)://重力
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Accel)://加速度
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Centrifugal_Force)://遠心力(点1,点2の座標, 角速度)
            if(dof != 7){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
                pLogger->Info(Utility::LoggerMode::Info,"Cent_Force DOF => point_1(x,y,z),point_2(x,y,z), omega");
            }
            mvValue.resize(7);
            break;
        default:
            break;
    }
}


