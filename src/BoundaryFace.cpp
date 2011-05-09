//
//  BoundaryFace.cpp
//
//
//
//                      2009.05.18
//                      2009.05.18
//                      k.Takeda
#include "BoundaryFace.h"
using namespace pmw;

CBoundaryFace::CBoundaryFace()
{
    ;
}
CBoundaryFace::~CBoundaryFace()
{
//    //debug
//    cout << "~CBoundaryFace" << endl;
}

// 初期化
// Type別のmvValueの領域確保
//
void CBoundaryFace::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;

    switch(mnType){
        case(BoundaryTypeFace::Pressure)://圧力
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::TractionVector)://粘性力(面方向力)のx,y,z方向の値
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::Thermal_Flux)://熱流束(面)
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::Temp)://温度
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


