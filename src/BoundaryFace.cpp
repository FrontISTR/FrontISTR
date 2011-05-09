/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   BoundaryFace.cxx
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "BoundaryFace.h"
using namespace pmw;
CBoundaryFace::CBoundaryFace()
{
    ;
}
CBoundaryFace::~CBoundaryFace()
{
}
void CBoundaryFace::initialize(const uint& bndType,const uint& dof)
{
    mnType = bndType;
    switch(mnType){
        case(BoundaryTypeFace::Pressure):
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::TractionVector):
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::Thermal_Flux):
            if(dof != 1){
                Utility::CLogger *pLogger = Utility::CLogger::Instance();
                pLogger->Info(Utility::LoggerMode::Warn,"mismatch DOF, input dof =>",dof);
            }
            mvValue.resize(dof);
            break;
        case(BoundaryTypeFace::Temp):
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
