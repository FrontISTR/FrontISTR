/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryVolume.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "BoundaryVolume.h"
using namespace pmw;
CBoundaryVolume::CBoundaryVolume()
{
    ;
}
CBoundaryVolume::~CBoundaryVolume()
{
}
void CBoundaryVolume::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;
    switch(mnType){
        case(BoundaryTypeVolume::Heat):
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Gravity):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Accel):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeVolume::Centrifugal_Force):
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
